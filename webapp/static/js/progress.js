document.addEventListener('DOMContentLoaded', () => {
  const jobId = document.body.dataset.jobId;
  const statusDiv = document.getElementById('status');
  const consoleDiv = document.getElementById('console');
  const downloadsDiv = document.getElementById('downloads');
  const dlZip = document.getElementById('dl-zip');
  const dlLog = document.getElementById('dl-log');

  function setStatus(status) {
    statusDiv.className = `status ${status}`;
    statusDiv.textContent = `Status: ${status.toUpperCase()}`;
  }

  function updateProgress() {
    fetch(`/progress/${jobId}`)
      .then(response => response.json())
      .then(data => {
        setStatus(data.status || 'unknown');

        if (Array.isArray(data.output) && data.output.length > 0) {
          consoleDiv.textContent = data.output.join('\n');
          consoleDiv.scrollTop = consoleDiv.scrollHeight;
        }

        if (data.status === 'completed') {
          downloadsDiv.style.display = 'block';
          clearInterval(pollInterval);
        } else if (data.status === 'failed' || data.status === 'error') {
          // Ensure log link is present and visible
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

  // Initialize download links (in case template changes)
  if (dlZip) dlZip.href = `/download/${jobId}`;
  if (dlLog) dlLog.href = `/download/${jobId}?file=web_run.log`;

  // Start polling
  updateProgress();
  const pollInterval = setInterval(updateProgress, 2000);
});
