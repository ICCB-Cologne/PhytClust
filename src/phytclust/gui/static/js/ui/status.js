/* ============================================================
   PhytClust – ui/status.js
   Status bar updates.
   ============================================================ */

import { statusEl, statusDot } from "../dom.js";
import { showToast } from "./toast.js";

export function showStatus(message, type) {
  type = type || "info";
  if (statusEl) statusEl.textContent = message;
  if (statusDot) {
    statusDot.className = "status-dot";
    if (type === "info") statusDot.classList.add("running");
    if (type === "success") statusDot.classList.add("ready");
    if (type === "danger") statusDot.classList.add("error");
  }
  if (type === "danger") showToast(message, "danger", 6000);
  else if (type === "success") showToast(message, "success", 3000);
}
