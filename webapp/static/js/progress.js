/**
 * progress.js - Job Progress Monitor
 *
 * Polls the server for job status updates and updates the UI accordingly.
 * Handles status display, console output streaming, and download link management.
 */

document.addEventListener('DOMContentLoaded', () => {
  // DOM element references
  const jobId = document.body.dataset.jobId;
  const statusDiv = document.getElementById('status');
  const consoleDiv = document.getElementById('console');
  const downloadsDiv = document.getElementById('downloads');
  const dlZip = document.getElementById('dl-zip');
  const dlLog = document.getElementById('dl-log');

  /**
   * Updates the status display with the given status value.
   * @param {string} status - Job status (e.g., 'running', 'completed', 'failed')
   */
  function setStatus(status) {
    statusDiv.className = `status ${status}`;
    statusDiv.textContent = `Status: ${status.toUpperCase()}`;
  }

  /**
   * Fetches the latest job progress from the server and updates the UI.
   * Handles status changes, console output updates, and completion states.
   * Stops polling when job reaches a terminal state (completed/failed/error).
   */
  function updateProgress() {
    fetch(`/progress/${jobId}`)
      .then(response => response.json())
      .then(data => {
        setStatus(data.status || 'unknown');

        // Update console output if available
        if (Array.isArray(data.output) && data.output.length > 0) {
          consoleDiv.textContent = data.output.join('\n');
          consoleDiv.scrollTop = consoleDiv.scrollHeight;
        }

        // Handle job completion
        if (data.status === 'completed') {
          downloadsDiv.style.display = 'block';
          clearInterval(pollInterval);
        } else if (data.status === 'failed' || data.status === 'error') {
          if (dlLog) dlLog.style.display = 'inline-block';
          downloadsDiv.style.display = 'block';
          clearInterval(pollInterval);
        }
      })
      .catch(err => {
        console.error('Progress poll failed:', err);
        setStatus('error');
      });
  }

  // Initialize download links
  if (dlZip) dlZip.href = `/download/${jobId}`;
  if (dlLog) dlLog.href = `/download/${jobId}?file=web_run.log`;

  // Start polling (2 second interval)
  updateProgress();
  const pollInterval = setInterval(updateProgress, 2000);
});