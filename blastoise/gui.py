"""
BLASTOISE Tkinter GUI
=====================

A simple graphical interface to run the BLASTOISE pipeline without using the
command line. This GUI mirrors the CLI arguments from blastoise.main and
executes the same steps under the hood.

Features
- File pickers for input data (FASTA/FA), genome (FASTA), and output folder
- Parameters with sensible defaults (same as CLI): identity, word size,
  min length, extend, limit, jobs
- Run button executes the pipeline in a background thread to avoid freezing the UI
- Live log panel that shows progress messages emitted by the pipeline
- Result summary dialog with paths to CSV and GFF on success

How it works
- Reuses the public functions from blastoise.main to avoid code duplication:
  setup_workspace, run_initial_blast, repetitive_sider_searcher, finalize_results
- Adds a file logger into the chosen output folder, just like the CLI does

Requirements
- Tkinter is part of the Python standard library; on conda it is provided by the
  "tk" package (added in meta.yaml). On some Linux distributions outside conda,
  you may need to install your system's python3-tk package.

Usage
- From a terminal (within the blastoise environment):
    blastoise-gui
  or
    python -m blastoise.gui

- Fill in the fields, select files and folder, then click Run.

Note: BLAST+ tools must be available in PATH (provided by conda "blast").
"""
from __future__ import annotations

import os
import sys
import threading
import queue
import traceback
from dataclasses import dataclass
from typing import Optional, Tuple

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import logging

# Reuse the pipeline building blocks from the CLI implementation
from blastoise.main import (
    setup_workspace,
    run_initial_blast,
    repetitive_sider_searcher,
    finalize_results,
    print_message_box,
    logger as module_logger,
)


# -------------------------
# Data structure for params
# -------------------------
@dataclass
class Params:
    data: str
    genome: str
    output: str
    identity: int = 60
    word_size: int = 11
    min_length: int = 100
    extend: int = 100
    limit: int = 1000
    jobs: int = -1


# -------------------------
# Logger/printing utilities
# -------------------------
class TextHandler(logging.Handler):
    """A logging handler that writes log records into a Tk Text widget via a queue.

    This avoids direct cross-thread UI updates. Records are placed on a queue
    and the GUI periodically pulls and inserts them into the Text widget.
    """
    def __init__(self, q: "queue.Queue[str]") -> None:
        super().__init__()
        self.q = q

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
        except Exception:
            msg = record.getMessage()
        self.q.put(msg + "\n")


class StdoutInterceptor:
    """Redirects print() output to a queue (for the GUI log panel)."""
    def __init__(self, q: "queue.Queue[str]", original_stream):
        self.q = q
        self.original_stream = original_stream

    def write(self, s: str):
        if s:
            self.q.put(s)
        # Also forward to original stream to preserve normal stdout behavior
        try:
            self.original_stream.write(s)
        except Exception:
            pass

    def flush(self):
        try:
            self.original_stream.flush()
        except Exception:
            pass


# -------------------------
# GUI Application
# -------------------------
class BlastoiseGUI(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("BLASTOISE GUI")
        self.geometry("820x640")
        self.minsize(740, 560)

        # Queue for cross-thread logging
        self.log_queue: "queue.Queue[str]" = queue.Queue()
        self._stdout_interceptor: Optional[StdoutInterceptor] = None
        self._file_handler_added: bool = False

        self._build_widgets()
        self._poll_log_queue()

    # ----- UI construction -----
    def _build_widgets(self) -> None:
        pad = 8

        frm = ttk.Frame(self)
        frm.pack(fill=tk.BOTH, expand=True, padx=pad, pady=pad)

        # Paths section
        paths_labelframe = ttk.LabelFrame(frm, text="Paths")
        paths_labelframe.pack(fill=tk.X, padx=pad, pady=pad)

        # Data file
        self.data_var = tk.StringVar()
        self._add_file_picker(paths_labelframe, "Input Data (FASTA)", self.data_var)

        # Genome file
        self.genome_var = tk.StringVar()
        self._add_file_picker(paths_labelframe, "Reference Genome (FASTA)", self.genome_var)

        # Output directory
        self.output_var = tk.StringVar()
        self._add_dir_picker(paths_labelframe, "Output Directory", self.output_var)

        # Parameters
        params_frame = ttk.LabelFrame(frm, text="Parameters")
        params_frame.pack(fill=tk.X, padx=pad, pady=pad)

        # Row 1
        row1 = ttk.Frame(params_frame)
        row1.pack(fill=tk.X, pady=4)
        self.identity_var = tk.IntVar(value=60)
        self._add_spinbox(row1, "Identity (%)", self.identity_var, 40, 100)
        self.word_size_var = tk.IntVar(value=11)
        self._add_spinbox(row1, "Word Size", self.word_size_var, 7, 64)
        self.min_length_var = tk.IntVar(value=100)
        self._add_spinbox(row1, "Min Length", self.min_length_var, 1, 1000000)

        # Row 2
        row2 = ttk.Frame(params_frame)
        row2.pack(fill=tk.X, pady=4)
        self.extend_var = tk.IntVar(value=100)
        self._add_spinbox(row2, "Extend (nt)", self.extend_var, 0, 1000000)
        self.limit_var = tk.IntVar(value=1000)
        self._add_spinbox(row2, "Length Limit (nt)", self.limit_var, 1, 100000000)
        self.jobs_var = tk.IntVar(value=-1)
        self._add_spinbox(row2, "Jobs (-1=all)", self.jobs_var, -1, 256)

        # Actions
        actions = ttk.Frame(frm)
        actions.pack(fill=tk.X, padx=pad, pady=(pad, 0))
        self.run_btn = ttk.Button(actions, text="Run", command=self._on_run)
        self.run_btn.pack(side=tk.LEFT)
        self.cancel_btn = ttk.Button(actions, text="Cancel", command=self._on_cancel, state=tk.DISABLED)
        self.cancel_btn.pack(side=tk.LEFT, padx=(8, 0))

        # Log output
        log_frame = ttk.LabelFrame(frm, text="Log")
        log_frame.pack(fill=tk.BOTH, expand=True, padx=pad, pady=pad)
        self.log_text = tk.Text(log_frame, wrap=tk.WORD, height=16)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.configure(state=tk.DISABLED)

        # Footer/help
        footer = ttk.Label(frm, anchor="w", justify=tk.LEFT,
                           text=(
                               "Usage: Select input files and output folder, adjust parameters if needed,\n"
                               "then click Run. The process may take a while depending on genome size.\n"
                               "Results (CSV and GFF) will be written into the chosen output directory."
                           ))
        footer.pack(fill=tk.X, padx=pad, pady=(0, pad))

    def _add_file_picker(self, parent, label, var: tk.StringVar) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=4)
        ttk.Label(row, text=label, width=24).pack(side=tk.LEFT)
        entry = ttk.Entry(row, textvariable=var)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(row, text="Browse...", command=lambda: self._browse_file(var)).pack(side=tk.LEFT, padx=6)

    def _add_dir_picker(self, parent, label, var: tk.StringVar) -> None:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=4)
        ttk.Label(row, text=label, width=24).pack(side=tk.LEFT)
        entry = ttk.Entry(row, textvariable=var)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(row, text="Browse...", command=lambda: self._browse_dir(var)).pack(side=tk.LEFT, padx=6)

    def _add_spinbox(self, parent, label, var: tk.IntVar, frm: int, to: int) -> None:
        block = ttk.Frame(parent)
        block.pack(side=tk.LEFT, padx=8)
        ttk.Label(block, text=label).pack(anchor="w")
        spin = ttk.Spinbox(block, textvariable=var, from_=frm, to=to, width=10)
        spin.pack(anchor="w")

    def _browse_file(self, var: tk.StringVar) -> None:
        path = filedialog.askopenfilename(title="Select file")
        if path:
            var.set(path)

    def _browse_dir(self, var: tk.StringVar) -> None:
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            var.set(path)

    # ----- Run/Cancel logic -----
    def _on_run(self) -> None:
        params = self._collect_params()
        if not params:
            return
        self._prepare_logging(params.output)
        self.run_btn.configure(state=tk.DISABLED)
        self.cancel_btn.configure(state=tk.NORMAL)
        self._append_log("Starting BLASTOISE pipeline...\n")
        self._worker_stop = False
        thread = threading.Thread(target=self._run_pipeline_thread, args=(params,), daemon=True)
        thread.start()

    def _on_cancel(self) -> None:
        # We cannot safely kill the pipeline (external processes may run).
        # We just set a flag and inform the user. For full cancellation, stop BLAST manually.
        self._worker_stop = True
        messagebox.showinfo("BLASTOISE", "Cancellation requested. The current step will need to finish before exit.")

    def _collect_params(self) -> Optional[Params]:
        data = self.data_var.get().strip()
        genome = self.genome_var.get().strip()
        output = self.output_var.get().strip()
        if not data or not os.path.isfile(data):
            messagebox.showerror("Missing input", "Please select a valid input data file.")
            return None
        if not genome or not os.path.isfile(genome):
            messagebox.showerror("Missing genome", "Please select a valid genome FASTA file.")
            return None
        if not output:
            messagebox.showerror("Missing output", "Please select an output directory.")
            return None
        return Params(
            data=data,
            genome=genome,
            output=output,
            identity=int(self.identity_var.get()),
            word_size=int(self.word_size_var.get()),
            min_length=int(self.min_length_var.get()),
            extend=int(self.extend_var.get()),
            limit=int(self.limit_var.get()),
            jobs=int(self.jobs_var.get()),
        )

    def _prepare_logging(self, output_dir: str) -> None:
        # Redirect stdout/stderr into the GUI log as well as original streams
        if not self._stdout_interceptor:
            self._stdout_interceptor = StdoutInterceptor(self.log_queue, sys.__stdout__)
            sys.stdout = self._stdout_interceptor
            sys.stderr = self._stdout_interceptor

        # Attach a file handler to the root logger into the output dir (once)
        if not self._file_handler_added:
            try:
                os.makedirs(output_dir, exist_ok=True)
                log_file_path = os.path.join(output_dir, 'blastoise.log')
                file_handler = logging.FileHandler(log_file_path)
                file_handler.setLevel(logging.INFO)
                file_handler.setFormatter(
                    logging.Formatter('%(asctime)s [%(levelname)s] %(name)s:%(funcName)s:%(lineno)d: %(message)s')
                )
                root_logger = logging.getLogger()
                root_logger.setLevel(logging.INFO)
                # Avoid duplicate handlers for same file
                has_same = any(
                    isinstance(h, logging.FileHandler) and getattr(h, 'baseFilename', None) == file_handler.baseFilename
                    for h in root_logger.handlers
                )
                if not has_same:
                    root_logger.addHandler(file_handler)
                self._file_handler_added = True
                module_logger.info("Logging to file: %s", log_file_path)
            except Exception as e:
                module_logger.warning("Failed to attach file handler for logging: %s", e)

        # Also stream logs into the GUI by attaching a TextHandler to root logger
        root_logger = logging.getLogger()
        if not any(isinstance(h, TextHandler) for h in root_logger.handlers):
            th = TextHandler(self.log_queue)
            th.setLevel(logging.INFO)
            th.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
            root_logger.addHandler(th)

    def _run_pipeline_thread(self, p: Params) -> None:
        try:
            # Mirror the CLI logic using the same functions
            print_message_box("Program started (GUI)")

            # 1. Setup workspace
            output_dir, data_path, genome_path = setup_workspace(
                output_dir=os.path.expanduser(p.output),
                data_file=os.path.expanduser(p.data),
                genome_file=os.path.expanduser(p.genome),
            )

            # 2. Initial BLAST
            initial_data, blast_db_path = run_initial_blast(
                data_path=data_path,
                genome_path=genome_path,
                output_dir=output_dir,
                identity=p.identity,
                word_size=p.word_size,
            )

            # 3. Iterative search (if any)
            if initial_data.empty:
                print_message_box("No data to process after initial BLAST. Writing empty results.")
                csv_path, gff_path = finalize_results(output_dir, initial_data, p.data, p.genome)
            else:
                print_message_box("Running iterative search")
                final_data = repetitive_sider_searcher(
                    data_input=initial_data,
                    genome_path=blast_db_path,
                    extend_number=p.extend,
                    word_size=p.word_size,
                    identity=p.identity,
                    min_length=p.min_length,
                    limit_len=p.limit,
                    n_jobs=p.jobs,
                )
                csv_path, gff_path = finalize_results(output_dir, final_data, p.data, p.genome)

            # 4. Done
            print_message_box("END OF THE PROGRAM (GUI)")
            messagebox.showinfo(
                "BLASTOISE",
                f"Finished.\nCSV: {csv_path}\nGFF: {gff_path}",
            )
        except Exception:
            tb = traceback.format_exc()
            self._append_log("\n" + tb + "\n")
            messagebox.showerror("BLASTOISE - Error", tb)
        finally:
            self.run_btn.configure(state=tk.NORMAL)
            self.cancel_btn.configure(state=tk.DISABLED)

    # ----- Log UI helpers -----
    def _poll_log_queue(self) -> None:
        try:
            while True:
                line = self.log_queue.get_nowait()
                self._append_log(line)
        except queue.Empty:
            pass
        self.after(100, self._poll_log_queue)

    def _append_log(self, text: str) -> None:
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, text)
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)


def main() -> None:
    """Launch the BLASTOISE GUI application."""
    app = BlastoiseGUI()
    app.mainloop()


if __name__ == "__main__":
    main()
