/* ============================================================
   PhytClust – ui/tabs.js
   Tab switching logic.
   ============================================================ */

import { state } from "../state.js";
import { drawOptimalK } from "../views/optimal_k.js";
import { drawComparison } from "../views/compare.js";

export function switchTab(tabName) {
  document
    .querySelectorAll(".tab-btn")
    .forEach((btn) =>
      btn.classList.toggle("active", btn.dataset.tab === tabName),
    );
  document
    .querySelectorAll(".viz-panel")
    .forEach((panel) => panel.classList.toggle("active", panel.id === tabName));

  var treeToolbar = document.getElementById("tree-toolbar");
  var optkToolbar = document.getElementById("optk-toolbar");
  var compareToolbar = document.getElementById("compare-toolbar");
  var compareConfigList = document.getElementById("compare-config-list");
  if (treeToolbar)
    treeToolbar.style.display = tabName === "viewer" ? "" : "none";
  if (optkToolbar)
    optkToolbar.style.display = tabName === "viewer-optimal-k" ? "" : "none";
  if (compareToolbar)
    compareToolbar.style.display = tabName === "compare" ? "" : "none";
  if (compareConfigList)
    compareConfigList.classList.toggle("visible", tabName === "compare");

  if (tabName === "viewer-optimal-k" && state.latestOptimalKData)
    drawOptimalK(state.latestOptimalKData);
  if (tabName === "compare") drawComparison();
}
