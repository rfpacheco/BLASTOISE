import os
import shutil
import subprocess
import uuid
from pathlib import Path
from typing import Optional, Dict, List

# noinspection PyPackageRequirements
from fastapi import (
    FastAPI,  # Main class to create the app
    File,  # For handling file uploads in POST request. Loads whole file in bytes.
    UploadFile,  # For handling file uploads in POST request. Streams file as an object (doesn't load in memory)
    Form,  # For reading form-encoded data
    BackgroundTasks,  # To run long-running tasks without blocking the HTTP response
    Request  # To pass request object into templates
)
# noinspection PyPackageRequirements
from fastapi.responses import (
    HTMLResponse,  # Returns raw HTML a response (renders web pages)
    FileResponse,  # Send files for download
    PlainTextResponse,  # Sends plain text
    StreamingResponse  # For streaming large responses to avoid loading them entirely in memory  # TODO: remove?
)
# noinspection PyPackageRequirements
from fastapi.staticfiles import StaticFiles
# noinspection PyPackageRequirements
from fastapi.templating import Jinja2Templates


# =====================================================================================
# Constants
# =====================================================================================
APP_ROOT: Path = Path("/app")  # Base directory of the app (used like this because of Docker container)
UPLOADS_DIR: Path = APP_ROOT / "uploads"  # Where uploaded files are stored
OUTPUT_DIR: Path = APP_ROOT / "output"  # Where results are written


# =====================================================================================
# Setting FastAPI
# =====================================================================================
# Create FastAPI application
app: FastAPI = FastAPI(title="BLASTOISE Web", description="Simple FastAPI front-end for BLASTOISE pipeline")

# Template and static configuration
BASE_DIR: Path = Path(__file__).parent
TEMPLATES_DIR: Path = BASE_DIR / "templates"
STATIC_DIR: Path = BASE_DIR / "static"

templates: Jinja2Templates = Jinja2Templates(directory=str(TEMPLATES_DIR))
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

# Global dict to store running jobs and their log outputs
running_jobs: Dict[str, Dict] = {}


# =====================================================================================
# Functions
# =====================================================================================
def _save_upload(dst_path: Path, up: UploadFile) -> None:
    """
    Saves an uploaded file to the specified destination path.

    The function ensures that the parent directories of the destination path
    exist, creating them if necessary, and then writes the content of the
    uploaded file to the specified path.

    Parameters:
    dst_path : Path
        The destination file path where the uploaded content will be saved.
    up : UploadFile
        The uploaded file object that contains the content to be saved.

    Returns:
    None
    """
    # Create parent directories if they don't exist
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    # Open destination file in write binary mode (wb) and copy uploaded file content
    # Uses shutil.copyfileobj for efficient streaming of file data
    # This avoids loading the entire file into memory at once
    with dst_path.open("wb") as f:
        shutil.copyfileobj(up.file, f)  # type: ignore


def _zip_dir(src_dir: Path, zip_path: Path) -> None:  # TODO: too complex maybe
    """
    Create a ZIP file from the contents of the specified directory.

    This function compresses all files (including those in subdirectories) within
    the provided directory into a single ZIP archive. Files are stored in the
    archive with paths relative to the parent of the given source directory.

    Parameters
    ----------
    src_dir : Path
        The source directory containing files and subdirectories to compress.
    zip_path : Path
        The destination path where the ZIP file will be created.

    """
    from zipfile import ZipFile, ZIP_DEFLATED
    with ZipFile(zip_path, 'w', ZIP_DEFLATED) as zipf:  # Open ZIP file in write mode with compression
        for root, _, files in os.walk(src_dir):  # Walk through all files in source directory
            for file in files:  # For each file found
                full_path = Path(root) / file  # Get full path to current file
                arcname = full_path.relative_to(src_dir.parent)  # Get relative path for ZIP structure
                zipf.write(full_path, arcname)  # Add file to ZIP using relative path


def run_blastoise_background(job_id: str, cmd: List[str], log_path: Path):
    """
    Runs a command in the background while capturing its output.

    This function runs a specified command as a subprocess, logs the output to a file,
    and updates the status and output associated with the given job ID. The subprocess
    output is monitored in real-time, and relevant details such as the return code and
    status are updated based on the execution result. If the command succeeds, the job
    output is archived in a zip file; otherwise, the job is marked as failed or errored.

    Parameters
    ----------
    job_id : str
        The unique identifier for the job being executed.
    cmd : list
        The command to be executed as a subprocess, represented as a list of strings.
    log_path : Path
        The path to the file where the subprocess output will be logged.
    """
    running_jobs[job_id]["status"] = "running"
    running_jobs[job_id]["output"] = []

    try:
        # Start the process
        # noinspection PyTypeChecker
        process = subprocess.Popen(  # TODO: check types
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        # Read output line by line
        with log_path.open("w") as logf:
            for line in process.stdout:
                line = line.rstrip()
                if line:  # Only add non-empty lines
                    running_jobs[job_id]["output"].append(line)
                    logf.write(line + "\n")
                    logf.flush()

        # Wait for process to complete
        process.wait()
        running_jobs[job_id]["return_code"] = process.returncode

        if process.returncode == 0:
            running_jobs[job_id]["status"] = "completed"
            # Create zip of results
            zip_path = OUTPUT_DIR / f"{job_id}.zip"
            job_output_dir = OUTPUT_DIR / job_id
            _zip_dir(job_output_dir, zip_path)
        else:
            running_jobs[job_id]["status"] = "failed"

    except Exception as e:
        running_jobs[job_id]["status"] = "error"
        running_jobs[job_id]["error"] = str(e)
        running_jobs[job_id]["output"].append(f"ERROR: {str(e)}")


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/health")
async def health() -> dict:
    return {"status": "ok"}


@app.post("/run", response_class=HTMLResponse)
async def run_blastoise(
        request: Request,
        background_tasks: BackgroundTasks,
        data_file: UploadFile = File(...),
        genome_file: UploadFile = File(...),
        identity: int = Form(60),
        word_size: int = Form(11),
        min_length: int = Form(100),
        extend: int = Form(100),
        limit_len: int = Form(1000),
        jobs: int = Form(-1),
        job_name: Optional[str] = Form(None),
):
    # Create unique job id and paths
    job_id = (job_name.strip().replace(' ', '_') if job_name else str(uuid.uuid4()))
    job_upload_dir = UPLOADS_DIR / job_id
    job_output_dir = OUTPUT_DIR / job_id

    # Save uploaded files
    data_path = job_upload_dir / (data_file.filename or "data.fasta")
    genome_path = job_upload_dir / (genome_file.filename or "genome.fasta")
    _save_upload(data_path, data_file)
    _save_upload(genome_path, genome_file)

    # Ensure output dir exists
    job_output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare command
    cmd = [
        "blastoise",
        "-d", str(data_path),
        "-g", str(genome_path),
        "-o", str(job_output_dir),
        "-i", str(identity),
        "-ws", str(word_size),
        "-min", str(min_length),
        "-ext", str(extend),
        "-lim", str(limit_len),
        "-j", str(jobs),
    ]

    # Initialize job status
    running_jobs[job_id] = {
        "status": "starting",
        "output": [],
        "return_code": None,
        "error": None
    }

    # Start background task
    log_path = job_output_dir / "web_run.log"
    background_tasks.add_task(run_blastoise_background, job_id, cmd, log_path)

    # Return progress page (rendered via template)
    return templates.TemplateResponse("progress.html", {"request": request, "job_id": job_id})


@app.get("/progress/{job_id}")
async def get_progress(job_id: str):
    """Get the current progress of a job"""
    if job_id not in running_jobs:
        return {"status": "not_found", "output": [], "return_code": None}

    job_data = running_jobs[job_id]
    return {
        "status": job_data["status"],
        "output": job_data["output"],
        "return_code": job_data["return_code"],
        "error": job_data.get("error")
    }


@app.get("/download/{job_id}")
async def download(job_id: str, file: Optional[str] = None):
    if file:
        # Send a single file from the job output directory
        target = OUTPUT_DIR / job_id / file
        if not target.exists():
            return PlainTextResponse("File not found", status_code=404)
        return FileResponse(path=str(target), filename=target.name)

    # Otherwise send the prepared zip
    zip_path = OUTPUT_DIR / f"{job_id}.zip"
    if not zip_path.exists():
        return PlainTextResponse("Archive not found", status_code=404)
    return FileResponse(path=str(zip_path), filename=zip_path.name)