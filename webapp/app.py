import os
import shutil
import subprocess
import uuid
import asyncio
import threading
from pathlib import Path
from typing import Optional, Dict

from fastapi import FastAPI, File, UploadFile, Form, BackgroundTasks
from fastapi.responses import HTMLResponse, FileResponse, PlainTextResponse, StreamingResponse

APP_ROOT = Path("/app")
UPLOADS_DIR = APP_ROOT / "uploads"
OUTPUT_DIR = APP_ROOT / "output"

app = FastAPI(title="BLASTOISE Web", description="Simple FastAPI front-end for BLASTOISE pipeline")

# Global dict to store running jobs and their log outputs
running_jobs: Dict[str, Dict] = {}


def _save_upload(dst_path: Path, up: UploadFile) -> None:
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    with dst_path.open("wb") as f:
        shutil.copyfileobj(up.file, f)


def _zip_dir(src_dir: Path, zip_path: Path) -> None:
    from zipfile import ZipFile, ZIP_DEFLATED
    with ZipFile(zip_path, 'w', ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(src_dir):
            for file in files:
                full_path = Path(root) / file
                arcname = full_path.relative_to(src_dir.parent)
                zipf.write(full_path, arcname)


def run_blastoise_background(job_id: str, cmd: list, log_path: Path):
    """Run BLASTOISE in background and capture output line by line"""
    running_jobs[job_id]["status"] = "running"
    running_jobs[job_id]["output"] = []

    try:
        # Start the process
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
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
async def index() -> str:
    return (
        """
        <html>
          <head>
            <title>BLASTOISE Web</title>
            <style>
              body { font-family: Arial, sans-serif; margin: 40px; }
              .box { max-width: 720px; margin: auto; padding: 20px; border: 1px solid #ddd; border-radius: 8px; }
              label { display: block; margin-top: 12px; font-weight: bold; }
              input[type=number], input[type=text] { width: 200px; }
              .row { margin-bottom: 10px; }
              button { margin-top: 16px; padding: 10px 16px; }

              /* Progress page styles */
              .progress-container { max-width: 900px; margin: auto; padding: 20px; }
              .console-output { 
                background: #1e1e1e; 
                color: #00ff00; 
                font-family: 'Courier New', monospace; 
                font-size: 12px;
                height: 400px; 
                overflow-y: auto; 
                padding: 15px; 
                border: 1px solid #333; 
                border-radius: 5px;
                white-space: pre-wrap;
              }
              .status { 
                padding: 10px; 
                margin: 10px 0; 
                border-radius: 5px; 
                font-weight: bold; 
              }
              .status.running { background: #fff3cd; color: #856404; }
              .status.completed { background: #d4edda; color: #155724; }
              .status.failed { background: #f8d7da; color: #721c24; }
              .status.error { background: #f8d7da; color: #721c24; }
              .download-links { margin-top: 20px; }
              .download-links a { 
                display: inline-block; 
                margin: 5px 10px 5px 0; 
                padding: 8px 15px; 
                background: #007bff; 
                color: white; 
                text-decoration: none; 
                border-radius: 4px; 
              }
            </style>
          </head>
          <body>
            <div class="box">
              <h2>BLASTOISE: Run from your browser</h2>
              <form action="/run" method="post" enctype="multipart/form-data">
                <div class="row">
                  <label>Data file (-d)</label>
                  <input type="file" name="data_file" required />
                </div>
                <div class="row">
                  <label>Genome file (-g)</label>
                  <input type="file" name="genome_file" required />
                </div>
                <div class="row"><label>Identity (-i)</label><input type="number" name="identity" value="60" /></div>
                <div class="row"><label>Word size (-ws)</label><input type="number" name="word_size" value="11" /></div>
                <div class="row"><label>Min length (-min)</label><input type="number" name="min_length" value="100" /></div>
                <div class="row"><label>Extend (-ext)</label><input type="number" name="extend" value="100" /></div>
                <div class="row"><label>Limit (-lim)</label><input type="number" name="limit_len" value="1000" /></div>
                <div class="row"><label>Jobs (-j)</label><input type="number" name="jobs" value="-1" /></div>
                <div class="row"><label>Optional job name</label><input type="text" name="job_name" placeholder="my_run" /></div>
                <button type="submit">Run BLASTOISE</button>
              </form>
            </div>
          </body>
        </html>
        """
    )


@app.get("/health")
async def health() -> dict:
    return {"status": "ok"}


@app.post("/run", response_class=HTMLResponse)
async def run_blastoise(
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

    # Return progress page
    return HTMLResponse(content=f"""
        <html>
          <head>
            <title>BLASTOISE Progress - Job: {job_id}</title>
            <style>
              body {{ font-family: Arial, sans-serif; margin: 40px; }}
              .progress-container {{ max-width: 900px; margin: auto; padding: 20px; }}
              .console-output {{ 
                background: #1e1e1e; 
                color: #00ff00; 
                font-family: 'Courier New', monospace; 
                font-size: 12px;
                height: 400px; 
                overflow-y: auto; 
                padding: 15px; 
                border: 1px solid #333; 
                border-radius: 5px;
                white-space: pre-wrap;
              }}
              .status {{ 
                padding: 10px; 
                margin: 10px 0; 
                border-radius: 5px; 
                font-weight: bold; 
              }}
              .status.running {{ background: #fff3cd; color: #856404; }}
              .status.completed {{ background: #d4edda; color: #155724; }}
              .status.failed {{ background: #f8d7da; color: #721c24; }}
              .status.error {{ background: #f8d7da; color: #721c24; }}
              .download-links {{ margin-top: 20px; }}
              .download-links a {{ 
                display: inline-block; 
                margin: 5px 10px 5px 0; 
                padding: 8px 15px; 
                background: #007bff; 
                color: white; 
                text-decoration: none; 
                border-radius: 4px; 
              }}
            </style>
          </head>
          <body>
            <div class="progress-container">
              <h2>BLASTOISE Progress</h2>
              <h3>Job ID: {job_id}</h3>
              <div id="status" class="status running">Status: Starting...</div>

              <h4>Console Output:</h4>
              <div id="console" class="console-output">Initializing...</div>

              <div id="downloads" class="download-links" style="display: none;">
                <a href="/download/{job_id}">Download Results (.zip)</a>
                <a href="/download/{job_id}?file=web_run.log">Download Log</a>
                <a href="/">Start New Job</a>
              </div>
            </div>

            <script>
              const jobId = '{job_id}';
              const statusDiv = document.getElementById('status');
              const consoleDiv = document.getElementById('console');
              const downloadsDiv = document.getElementById('downloads');

              // Poll for updates every 2 seconds
              function updateProgress() {{
                fetch(`/progress/${{jobId}}`)
                  .then(response => response.json())
                  .then(data => {{
                    // Update status
                    statusDiv.className = `status ${{data.status}}`;
                    statusDiv.textContent = `Status: ${{data.status.toUpperCase()}}`;

                    // Update console output
                    if (data.output && data.output.length > 0) {{
                      consoleDiv.textContent = data.output.join('\\n');
                      consoleDiv.scrollTop = consoleDiv.scrollHeight;
                    }}

                    // Show downloads if completed
                    if (data.status === 'completed') {{
                      downloadsDiv.style.display = 'block';
                      clearInterval(pollInterval);
                    }} else if (data.status === 'failed' || data.status === 'error') {{
                      downloadsDiv.innerHTML = `
                        <a href="/download/${{jobId}}?file=web_run.log">Download Error Log</a>
                        <a href="/">Start New Job</a>
                      `;
                      downloadsDiv.style.display = 'block';
                      clearInterval(pollInterval);
                    }}
                  }})
                  .catch(error => {{
                    console.error('Error:', error);
                    statusDiv.className = 'status error';
                    statusDiv.textContent = 'Status: ERROR - Connection failed';
                  }});
              }}

              // Update immediately and then every 2 seconds
              updateProgress();
              const pollInterval = setInterval(updateProgress, 2000);
            </script>
          </body>
        </html>
    """)


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