/* ============================================================
   PhytClust – ui/file.js
   File handling: load, drag-and-drop.
   ============================================================ */

import { newickEl } from "../dom.js";
import { showStatus } from "./status.js";

export function handleFileSelect(evt) {
  var file = evt.target.files ? evt.target.files[0] : null;
  if (!file) return;
  loadFile(file);
}

export function loadFile(file) {
  var fileNameEl = document.getElementById("file-name");
  var fileBadge = document.getElementById("file-badge");
  if (fileNameEl && fileBadge) {
    fileBadge.textContent = file.name;
    fileNameEl.style.display = "flex";
  }
  showStatus("Loading tree...", "info");
  var reader = new FileReader();
  reader.onload = function (e) {
    newickEl.value = (e.target.result || "").trim();
    showStatus("Tree loaded", "success");
  };
  reader.readAsText(file);
}
