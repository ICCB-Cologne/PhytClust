/* ============================================================
   PhytClust – views/compare.js
   Comparison mode: multi-bar tree comparisons.
   ============================================================ */

import { state } from "../state.js";
import { showToast } from "../ui/toast.js";
import { escapeHtmlAttr } from "../utils.js";
import { drawComparisonBarsInto } from "../tree/draw.js";

function listTreeLeaves() {
  const out = [];
  function walk(node) {
    if (!node) return;
    if (node.children && node.children.length) {
      for (const c of node.children) walk(c);
      return;
    }
    if (node.name) out.push(node.name);
  }
  walk(state.NEWICK_RAW_TREE);
  return out;
}

export function normalizeClustersForLeaves(rawMap) {
  if (!rawMap) return null;
  const leaves = listTreeLeaves();
  if (!leaves.length) return null;
  const leafSet = new Set(leaves);
  const out = {};
  let found = 0;
  Object.keys(rawMap).forEach((name) => {
    if (!leafSet.has(name)) return;
    const cid = parseInt(rawMap[name], 10);
    if (isNaN(cid)) return;
    out[name] = cid;
    found += 1;
  });
  if (!found) return null;
  for (const name of leaves) {
    if (out[name] == null) out[name] = 0;
  }
  return out;
}

export function getAvailableCompareKs(data) {
  if (!data) return [];
  if (data.all_ks && data.all_ks.length) return data.all_ks.slice();
  return (data.k_values || data.ks || []).slice();
}

export function getClustersForK(k) {
  if (!state.latestApiData) return null;
  if (state.latestApiData.all_clusters && state.latestApiData.all_ks) {
    const idx = state.latestApiData.all_ks.indexOf(k);
    if (idx >= 0) return state.latestApiData.all_clusters[idx];
  }
  const peakKs = state.latestApiData.k_values || state.latestApiData.ks;
  if (peakKs && state.latestApiData.clusters) {
    const idx = peakKs.indexOf(k);
    if (idx >= 0) return state.latestApiData.clusters[idx];
  }
  return null;
}

export function parseClusterTableText(rawText) {
  const lines = String(rawText || "")
    .replace(/\r\n/g, "\n")
    .split("\n")
    .map((s) => s.trim())
    .filter(Boolean);
  if (!lines.length) return { map: null, invalidRows: 0 };

  const delimiter = lines[0].indexOf("\t") >= 0 ? "\t" : ",";
  const map = {};
  let invalidRows = 0;
  let parsedRows = 0;

  lines.forEach((line, idx) => {
    let cols = line.split(delimiter).map((s) => s.trim());
    if (cols.length < 2 && delimiter === ",") {
      cols = line.split(/\s+/).map((s) => s.trim());
    }
    if (cols.length < 2) {
      invalidRows += 1;
      return;
    }

    let leaf = cols[0];
    let cidRaw = cols[1];

    const firstNum = parseInt(cols[0], 10);
    const secondNum = parseInt(cols[1], 10);
    if (!isNaN(firstNum) && isNaN(secondNum)) {
      leaf = cols[1];
      cidRaw = cols[0];
    }

    const cid = parseInt(cidRaw, 10);
    if (isNaN(cid)) {
      if (idx === 0) return;
      invalidRows += 1;
      return;
    }
    map[leaf] = cid;
    parsedRows += 1;
  });

  return { map: parsedRows ? map : null, invalidRows };
}

export function makeCompareConfig() {
  const ks = getAvailableCompareKs(state.latestApiData);
  const k = ks.length ? Number(ks[0]) : null;
  const id = state.COMPARE_CONFIG_NEXT_ID++;
  return {
    id,
    title: "Bar " + id,
    source: "k",
    k,
    fileName: "",
    fileClusters: null,
  };
}

export function ensureCompareConfigs(data) {
  const ks = getAvailableCompareKs(data);
  if (!state.COMPARE_CONFIGS.length) {
    const cfg = makeCompareConfig();
    cfg.title = ks.length ? "k=" + ks[0] : "Bar 1";
    cfg.k = ks.length ? Number(ks[0]) : null;
    state.COMPARE_CONFIGS = [cfg];
  }

  state.COMPARE_CONFIGS.forEach((cfg) => {
    if (
      cfg.source === "k" &&
      (cfg.k == null || ks.indexOf(Number(cfg.k)) < 0)
    ) {
      cfg.k = ks.length ? Number(ks[0]) : null;
    }
    if (cfg.source === "file" && cfg.fileClusters) {
      cfg.fileClusters = normalizeClustersForLeaves(cfg.fileClusters);
    }
  });
}

function updateCompareHeader(comparisons) {
  const leftLabel = document.getElementById("compare-left-label");
  if (!leftLabel) return;
  if (!comparisons.length) {
    leftLabel.textContent =
      "Add comparison bars from k results or uploaded files.";
    return;
  }
  leftLabel.textContent = comparisons.map((c) => c.title).join(" | ");
}

export function renderCompareConfigList() {
  const host = document.getElementById("compare-config-list");
  if (!host) return;
  if (!state.latestApiData) {
    host.innerHTML = "";
    return;
  }

  const ks = getAvailableCompareKs(state.latestApiData);
  host.innerHTML = "";

  state.COMPARE_CONFIGS.forEach((cfg) => {
    const card = document.createElement("div");
    card.className = "compare-config-card";
    card.dataset.cfgId = String(cfg.id);

    const titleGroup = document.createElement("div");
    titleGroup.className = "form-group";
    titleGroup.innerHTML =
      '<label>Title</label><input class="form-input compare-title" type="text" value="' +
      escapeHtmlAttr(cfg.title || "") +
      '" />';

    const sourceGroup = document.createElement("div");
    sourceGroup.className = "form-group";
    sourceGroup.innerHTML =
      '<label>Source</label><select class="form-input compare-source"><option value="k">Run result (k)</option><option value="file">CSV/TSV file</option></select>';

    const removeWrap = document.createElement("div");
    removeWrap.style.alignSelf = "end";
    const removeBtn = document.createElement("button");
    removeBtn.className = "btn-ghost";
    removeBtn.type = "button";
    removeBtn.textContent = "Remove";
    removeWrap.appendChild(removeBtn);

    const kGroup = document.createElement("div");
    kGroup.className = "form-group full-row compare-k-wrap";
    const kLabel = document.createElement("label");
    kLabel.textContent = "k value";
    const kSelect = document.createElement("select");
    kSelect.className = "form-input compare-k";
    ks.forEach((k) => {
      const opt = document.createElement("option");
      opt.value = String(k);
      opt.textContent = "k = " + k;
      kSelect.appendChild(opt);
    });
    if (cfg.k != null) kSelect.value = String(cfg.k);
    kGroup.appendChild(kLabel);
    kGroup.appendChild(kSelect);

    const fileGroup = document.createElement("div");
    fileGroup.className = "form-group full-row compare-file-wrap";
    const fLabel = document.createElement("label");
    fLabel.textContent = "Cluster file";
    const fInput = document.createElement("input");
    fInput.type = "file";
    fInput.className = "form-input compare-file";
    fInput.accept = ".csv,.tsv,.txt";
    const fHint = document.createElement("div");
    fHint.className = "hint";
    fHint.textContent = cfg.fileName
      ? "Loaded: " + cfg.fileName
      : "Two-column file: leaf_name, cluster_id";
    fileGroup.appendChild(fLabel);
    fileGroup.appendChild(fInput);
    fileGroup.appendChild(fHint);

    card.appendChild(titleGroup);
    card.appendChild(sourceGroup);
    card.appendChild(removeWrap);
    card.appendChild(kGroup);
    card.appendChild(fileGroup);
    host.appendChild(card);

    const titleInput = card.querySelector(".compare-title");
    const sourceSel = card.querySelector(".compare-source");
    sourceSel.value = cfg.source;

    function updateVisibility() {
      const useFile = sourceSel.value === "file";
      kGroup.style.display = useFile ? "none" : "";
      fileGroup.style.display = useFile ? "" : "none";
    }
    updateVisibility();

    titleInput.addEventListener("input", function () {
      cfg.title = titleInput.value.trim() || "Bar " + cfg.id;
      drawComparison();
    });

    sourceSel.addEventListener("change", function () {
      cfg.source = sourceSel.value;
      updateVisibility();
      drawComparison();
    });

    kSelect.addEventListener("change", function () {
      cfg.k = parseInt(kSelect.value, 10);
      if (!titleInput.value.trim() || /^k=\d+$/.test(titleInput.value.trim())) {
        cfg.title = "k=" + cfg.k;
        titleInput.value = cfg.title;
      }
      drawComparison();
    });

    fInput.addEventListener("change", function () {
      const file = fInput.files && fInput.files[0] ? fInput.files[0] : null;
      if (!file) return;
      const reader = new FileReader();
      reader.onload = function (e) {
        const txt = e && e.target ? e.target.result : "";
        const parsed = parseClusterTableText(txt);
        const normalized = normalizeClustersForLeaves(parsed.map);
        if (!normalized) {
          showToast(
            "Could not parse a valid leaf/cluster mapping from file.",
            "danger",
            3200,
          );
          return;
        }
        cfg.fileClusters = normalized;
        cfg.fileName = file.name;
        if (
          !titleInput.value.trim() ||
          /^Bar \d+$/.test(titleInput.value.trim())
        ) {
          cfg.title = file.name.replace(/\.[^.]+$/, "");
          titleInput.value = cfg.title;
        }
        if (parsed.invalidRows > 0) {
          showToast(
            "Loaded with " + parsed.invalidRows + " skipped row(s).",
            "info",
            2500,
          );
        }
        renderCompareConfigList();
        drawComparison();
      };
      reader.readAsText(file);
    });

    removeBtn.addEventListener("click", function () {
      state.COMPARE_CONFIGS = state.COMPARE_CONFIGS.filter((x) => x.id !== cfg.id);
      if (!state.COMPARE_CONFIGS.length) state.COMPARE_CONFIGS = [makeCompareConfig()];
      renderCompareConfigList();
      drawComparison();
    });
  });
}

export function populateCompareSelectors(data) {
  ensureCompareConfigs(data);
  renderCompareConfigList();
}

export function drawComparison() {
  if (!state.latestApiData || !state.NEWICK_RAW_TREE) return;

  const leftHost = document.getElementById("compare-tree-left");
  if (!leftHost) return;

  const comparisons = [];
  state.COMPARE_CONFIGS.forEach((cfg) => {
    if (cfg.source === "k") {
      const raw = getClustersForK(Number(cfg.k));
      const clusters = normalizeClustersForLeaves(raw);
      if (clusters) {
        comparisons.push({
          title: cfg.title || (cfg.k != null ? "k=" + cfg.k : "k"),
          clusters,
        });
      }
      return;
    }
    if (cfg.source === "file" && cfg.fileClusters) {
      comparisons.push({
        title: cfg.title || cfg.fileName || "file",
        clusters: cfg.fileClusters,
      });
    }
  });

  updateCompareHeader(comparisons);
  if (!comparisons.length) {
    leftHost.innerHTML =
      '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);font-size:14px;">No comparison bars configured yet.</div>';
    return;
  }
  drawComparisonBarsInto(leftHost, comparisons);
}
