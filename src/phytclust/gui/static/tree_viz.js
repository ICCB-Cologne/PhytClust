/* ============================================================
   PhytClust – tree_viz.js
   ============================================================ */

// ── Globals ──
let HIER_CART = null;
let HIER_CIRC = null;

const EXAMPLE_NEWICK =
  "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

const newickEl = document.getElementById("newick-input");
const resultEl = document.getElementById("result");
const statusEl = document.getElementById("status-message");
const statusDot = document.getElementById("status-dot");
const treeHost = document.getElementById("tree_display");

const BASE_COLORS = [
  "#b84b4b",
  "#849060",
  "#3d7c74",
  "#6e3f8a",
  "#ceb94b",
  "#3f648a",
  "#3f408a",
  "#da63aa",
  "#c06f2e",
  "#2f6f93",
  "#4f8f4a",
  "#ad5c7a",
  "#7a5d3b",
  "#2f8a85",
  "#8f4b7f",
  "#5b6bb3",
];

let SHOW_LEAF_NAMES = true;
let SHOW_INTERNAL_NAMES = false;
let LEAF_NODE_RADIUS = 3.0;
let INTERNAL_NODE_RADIUS = 1.8;
let LABEL_FONT_SIZE = 9;
let AXIS_FONT_SIZE = 10;
let TREE_WIDTH_SCALE = 1.0;
let TREE_HEIGHT_SCALE = 1.0;
let BRANCH_STROKE_WIDTH = 1.2;
let BRANCH_COLOR_OVERRIDE = null;
let CLUSTER_COLORS = [];
let CURRENT_CLUSTERS = {};
let NEWICK_RAW_TREE = null;
let CURRENT_LAYOUT_MODE = "rectangular";
let COLOR_MODE = "bars";
const LABEL_PAD = 6;

let isRunning = false;
let latestOptimalKData = null;
let latestApiData = null;
let latestRunId = null;
let LAST_TREE_SVG = null;
let LAST_TREE_ZOOM = null;
let LAST_TREE_TRANSFORM = d3.zoomIdentity;
let SEARCH_TERM = "";
let SHOW_BOX_LABELS = true;
let SHOW_OUTLIER_BOXES = true;
let BOX_ALPHA = 0.12;
let BOX_PAD_V = 8;
let BOX_PAD_H = 6;
let BOX_CORNER_RADIUS = 6;
let CLUSTER_LABEL_FONT_SIZE = 8;
const SELECTED_CLUSTER_IDS = new Set();
const BOX_LABEL_MAP = {}; // cid -> custom name, editable via context menu
const BOX_ADJUST_MAP = {}; // cid -> {dx, dy, padX, padY} in px

function hasClusterFocus(cid) {
  if (cid == null) return SELECTED_CLUSTER_IDS.size === 0;
  return (
    SELECTED_CLUSTER_IDS.size === 0 || SELECTED_CLUSTER_IDS.has(Number(cid))
  );
}

function getBoxAdjust(cid) {
  const key = String(cid);
  if (!BOX_ADJUST_MAP[key]) {
    BOX_ADJUST_MAP[key] = { dx: 0, dy: 0, padX: 0, padY: 0 };
  }
  return BOX_ADJUST_MAP[key];
}

// Per-node customizations: keyed by data node reference
const NODE_CUSTOM = new WeakMap();

// Element references for params
const extraKEl = document.getElementById("extra-k");
const extraOutgroupEl = document.getElementById("extra-outgroup");
const extraRootTaxonEl = document.getElementById("extra-root-taxon");
const extraTopNEl = document.getElementById("extra-topn");
const extraResolutionEl = document.getElementById("extra-resolution");
const extraBinsEl = document.getElementById("extra-bins");
const extraMaxKEl = document.getElementById("extra-maxk");
const extraMaxKLimitEl = document.getElementById("extra-maxklimit");
const extraLambdaEl = document.getElementById("extra-lambda");
const extraMinClusterEl = document.getElementById("extra-min-cluster-size");
const extraOutlierEl = document.getElementById("extra-outlier");

// D3 tooltip
const d3Tooltip = d3.select("body").append("div").attr("class", "d3-tooltip");

/* ─────────────────────────────────────
   Dark Mode
   ───────────────────────────────────── */
function initTheme() {
  document.documentElement.setAttribute("data-theme", "light");
  localStorage.removeItem("phytclust-theme");
}
initTheme();

function toggleTheme() {
  showToast("Theme switching is disabled in this GUI.", "info", 1800);
}

/* ─────────────────────────────────────
   Toast Notifications
   ───────────────────────────────────── */
function showToast(message, type, duration) {
  type = type || "info";
  duration = duration || 4000;
  const container = document.getElementById("toast-container");
  if (!container) return;
  const toast = document.createElement("div");
  toast.className = "toast " + type;
  toast.textContent = message;
  container.appendChild(toast);
  setTimeout(function () {
    toast.classList.add("fade-out");
    setTimeout(function () {
      toast.remove();
    }, 300);
  }, duration);
}

/* ─────────────────────────────────────
   Status Bar
   ───────────────────────────────────── */
function showStatus(message, type) {
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

/* ─────────────────────────────────────
   Utility Functions
   ───────────────────────────────────── */
function escapeRegExp(str) {
  return str.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function isOutgroupInNewick(newick, outgroup) {
  if (!outgroup) return true;
  var name = escapeRegExp(outgroup);
  return new RegExp("[\\(,]\\s*'?" + name + "'?\\s*(?:[:,\\)])").test(newick);
}

function estimateLeafCount(newick) {
  var inQuote = null,
    commas = 0;
  for (var i = 0; i < newick.length; i++) {
    var ch = newick[i];
    if (inQuote) {
      if (ch === inQuote) inQuote = null;
      continue;
    }
    if (ch === '"' || ch === "'") {
      inQuote = ch;
      continue;
    }
    if (ch === ",") commas++;
  }
  return Math.max(0, commas + 1);
}

function resetExtraParams() {
  [
    "extra-k",
    "extra-outgroup",
    "extra-root-taxon",
    "extra-topn",
    "extra-bins",
    "extra-maxk",
    "extra-maxklimit",
    "extra-lambda",
    "extra-min-cluster-size",
    "extra-min-prominence",
    "extra-outlier-threshold",
    "extra-min-support",
    "extra-support-weight",
  ].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.value = "";
  });
  ["extra-outlier", "extra-optimize-polytomies"].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.checked = true;
  });
  [
    "extra-resolution",
    "extra-compute-all",
    "extra-outlier-prefer-fewer",
    "extra-no-split-zero",
    "extra-use-support",
    "extra-relative-prom",
  ].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.checked = false;
  });
  var rankMode = document.getElementById("extra-ranking-mode");
  if (rankMode) rankMode.value = "adjusted";
  var ratioMode = document.getElementById("extra-outlier-ratio-mode");
  if (ratioMode) ratioMode.value = "exp";
  showToast("Parameters reset to defaults.", "info", 2000);
}

/* ─────────────────────────────────────
   Newick Parser
   ───────────────────────────────────── */
function parseNewick(newick) {
  var tokens = newick.split(/\s*(;|\(|\)|,|:)\s*/);
  var stack = [];
  var current = {};
  for (var i = 0; i < tokens.length; i++) {
    var token = tokens[i];
    switch (token) {
      case "(":
        var child = {};
        current.children = [child];
        stack.push(current);
        current = child;
        break;
      case ",":
        var sibling = {};
        stack[stack.length - 1].children.push(sibling);
        current = sibling;
        break;
      case ")":
        current = stack.pop();
        break;
      case ":":
        var length = parseFloat(tokens[++i]);
        current.length = isNaN(length) ? 0 : length;
        break;
      default:
        if (token && token !== ";") current.name = token;
    }
  }
  return current;
}

function getVisibleChildren(node) {
  if (node && node._collapsed) return null;
  return node && node.children && node.children.length ? node.children : null;
}

function getAllChildren(node) {
  return node && node.children && node.children.length ? node.children : null;
}

/* ─────────────────────────────────────
   Color Generation
   ───────────────────────────────────── */
function shuffle(arr) {
  let a = arr.slice();
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

function adjustLight(hex, factor) {
  const num = parseInt(hex.slice(1), 16);
  let r = (num >> 16) + factor * 255;
  let g = ((num >> 8) & 0xff) + factor * 255;
  let b = (num & 0xff) + factor * 255;
  return `rgb(${Math.min(255, Math.max(0, r)) | 0}, ${Math.min(255, Math.max(0, g)) | 0}, ${Math.min(255, Math.max(0, b)) | 0})`;
}

function withAlpha(color, alpha) {
  if (color.startsWith("rgb"))
    return color.replace("rgb", "rgba").replace(")", `, ${alpha})`);
  const num = parseInt(color.slice(1), 16);
  return `rgba(${num >> 16}, ${(num >> 8) & 0xff}, ${num & 0xff}, ${alpha})`;
}

function generateClusterColors(nClusters) {
  let palette = shuffle(BASE_COLORS);
  if (nClusters <= palette.length) return palette.slice(0, nClusters);
  let colors = [];
  const minAlpha = 0.58;
  const alphaStep = 0.12;
  const lightStep = 0.14;
  const repeats = Math.ceil(nClusters / palette.length);
  for (let r = 0; r < repeats; r++) {
    const factor = r * lightStep;
    const alpha = Math.max(1 - r * alphaStep, minAlpha);
    palette.forEach((hex) => {
      const adjusted = adjustLight(hex, factor);
      colors.push(alpha < 1 ? withAlpha(adjusted, alpha) : adjusted);
    });
  }
  return colors.slice(0, nClusters);
}

/* ─────────────────────────────────────
   Theme-aware colors for D3
   ───────────────────────────────────── */
function getThemeColors() {
  const style = getComputedStyle(document.documentElement);
  const themeBranch =
    style.getPropertyValue("--pc-tree-branch").trim() || "#000000";
  return {
    branch: BRANCH_COLOR_OVERRIDE || themeBranch,
    label: style.getPropertyValue("--pc-tree-label").trim() || "#334155",
    internal: style.getPropertyValue("--pc-tree-internal").trim() || "#64748b",
    text: style.getPropertyValue("--pc-text").trim() || "#241b3d",
    muted: style.getPropertyValue("--pc-text-muted").trim() || "#9a91b2",
    bg: style.getPropertyValue("--pc-bg").trim() || "#f7f5fa",
    border: style.getPropertyValue("--pc-border").trim() || "#e4ddf0",
    accent: style.getPropertyValue("--pc-accent").trim() || "#6b61ac",
  };
}

/* ─────────────────────────────────────
   Tree Layout Computation
   ───────────────────────────────────── */
function accumulateBranchLength(node, length) {
  length = length || 0;
  node._bl = length;
  if (node.children) {
    for (const c of node.children)
      accumulateBranchLength(c, length + (c.length || 0));
  }
}

function computeLayouts() {
  if (!NEWICK_RAW_TREE) return;
  // Always compute on full topology so collapsing does not move node positions.
  HIER_CART = d3.hierarchy(NEWICK_RAW_TREE, getAllChildren);
  HIER_CIRC = d3.hierarchy(NEWICK_RAW_TREE, getAllChildren);

  var width = treeHost.clientWidth || 800;
  var height = treeHost.clientHeight || 500;
  var margin = { top: 20, right: 80, bottom: 20, left: 80 };
  var innerW = width - margin.left - margin.right;
  var innerH = height - margin.top - margin.bottom;
  var radius = Math.min(innerW, innerH) / 2;

  var maxBl = d3.max(HIER_CART.descendants(), (d) => d.data._bl || 0) || 1;
  var blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);
  var blToR = d3.scaleLinear().domain([0, maxBl]).range([0, radius]);

  d3.cluster().size([innerH, 1])(HIER_CART);
  HIER_CART.each((d) => {
    d._x = d.x;
    d._y = blToX(d.data._bl || 0);
  });

  d3.cluster().size([2 * Math.PI, 1])(HIER_CIRC);
  HIER_CIRC.each((d) => {
    d._angle = d.x;
    d._radius = blToR(d.data._bl || 0);
  });
}

/* ─────────────────────────────────────
   Node Context Menu
   ───────────────────────────────────── */
let CTX_TARGET_DATA = null;

function getNodeCustom(dataNode) {
  if (!NODE_CUSTOM.has(dataNode))
    NODE_CUSTOM.set(dataNode, {
      highlighted: false,
      radiusScale: 1,
      renamedTo: null,
      lockOriginalSize: false,
    });
  return NODE_CUSTOM.get(dataNode);
}

function nodeDisplayName(dataNode) {
  const c = NODE_CUSTOM.has(dataNode) ? NODE_CUSTOM.get(dataNode) : null;
  return c && c.renamedTo != null ? c.renamedTo : dataNode.name || "";
}

function countLeafDescendantsData(node) {
  if (!node) return 0;
  if (!node.children || !node.children.length) return 1;
  var total = 0;
  for (var i = 0; i < node.children.length; i++) {
    total += countLeafDescendantsData(node.children[i]);
  }
  return total;
}

function customNodeRadius(d) {
  let base =
    d.data && d.data.children && d.data.children.length
      ? INTERNAL_NODE_RADIUS
      : LEAF_NODE_RADIUS;
  const c = NODE_CUSTOM.has(d.data) ? NODE_CUSTOM.get(d.data) : null;

  // Per-node override: keep default radius even when collapsed.
  if (c && c.lockOriginalSize) {
    return base;
  }

  // When an internal node is collapsed, scale its marker by subtree leaf count.
  if (
    d.data &&
    d.data._collapsed &&
    d.data.children &&
    d.data.children.length
  ) {
    var leafCount = countLeafDescendantsData(d.data);
    if (leafCount > 1) {
      // Log scaling keeps larger collapsed nodes readable without obscuring tree.
      var autoCollapsed =
        INTERNAL_NODE_RADIUS + LEAF_NODE_RADIUS * 0.55 * Math.log2(leafCount);
      base = Math.max(base, Math.min(14, autoCollapsed));
    }
  }

  return base * (c ? c.radiusScale : 1);
}

function showNodeContextMenu(event, d) {
  event.preventDefault();
  event.stopPropagation();
  CTX_TARGET_DATA = d.data;
  const menu = document.getElementById("node-context-menu");
  if (!menu) return;
  menu.style.display = "block";
  menu.style.left = event.clientX + "px";
  menu.style.top = event.clientY + "px";
  // Sync size slider to current node's scale
  var ctxSlider = document.getElementById("ctx-size-slider");
  var ctxSliderVal = document.getElementById("ctx-size-value");
  var c = NODE_CUSTOM.has(d.data)
    ? NODE_CUSTOM.get(d.data)
    : { radiusScale: 1 };
  if (ctxSlider) ctxSlider.value = c.radiusScale;
  if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);

  // Show/hide items based on context
  var isLeaf = !(d.data.children && d.data.children.length);
  var hasClusters =
    CURRENT_CLUSTERS && Object.keys(CURRENT_CLUSTERS).length > 0;
  var isBoxMode = COLOR_MODE === "boxes";

  var copyBtn = document.getElementById("ctx-copy-subtree");
  var originalSizeBtn = document.getElementById("ctx-original-size");
  var clusterBtn = document.getElementById("ctx-select-cluster");
  var renameClusterBtn = document.getElementById("ctx-rename-cluster");
  var adjustBoxBtn = document.getElementById("ctx-adjust-box");
  var c = getNodeCustom(d.data);

  // "Copy subtree names" only for internal nodes
  if (copyBtn) copyBtn.style.display = isLeaf ? "none" : "";
  // "Keep original size" only useful for internal nodes
  if (originalSizeBtn) {
    originalSizeBtn.style.display = isLeaf ? "none" : "";
    originalSizeBtn.textContent = c.lockOriginalSize
      ? "Use auto collapsed size"
      : "Keep original size";
  }
  // "Focus cluster" only when clusters are available
  if (clusterBtn) clusterBtn.style.display = hasClusters ? "" : "none";
  // "Rename cluster box" only in box mode with clusters
  if (renameClusterBtn)
    renameClusterBtn.style.display = isBoxMode && hasClusters ? "" : "none";
  // "Adjust box" only in box mode with clusters
  if (adjustBoxBtn)
    adjustBoxBtn.style.display = isBoxMode && hasClusters ? "" : "none";
}

function hideNodeContextMenu() {
  const menu = document.getElementById("node-context-menu");
  if (menu) menu.style.display = "none";
  CTX_TARGET_DATA = null;
}

document.addEventListener("click", hideNodeContextMenu);
document.addEventListener("contextmenu", function (e) {
  if (!e.target.closest(".tree-node")) hideNodeContextMenu();
});

/* ─────────────────────────────────────
   Collect leaf names from subtree
   ───────────────────────────────────── */
function collectLeafNamesData(node, out) {
  if (!node) return;
  const kids = node.children && node.children.length ? node.children : null;
  if (!kids) {
    if (node.name) out.push(node.name);
    return;
  }
  for (const c of kids) collectLeafNamesData(c, out);
}

function representativeClusterIdFromData(dataNode) {
  if (!dataNode || !CURRENT_CLUSTERS) return null;
  const leaves = [];
  collectLeafNamesData(dataNode, leaves);
  if (!leaves.length) return null;
  let cid = null;
  for (const nm of leaves) {
    const v = CURRENT_CLUSTERS[nm];
    if (v == null) return null;
    if (cid == null) cid = v;
    else if (cid !== v) return null;
  }
  return cid;
}

/* ─────────────────────────────────────
   Search highlighting helper
   ───────────────────────────────────── */
function isSearchMatch(dataNode) {
  if (!SEARCH_TERM) return false;
  const name = nodeDisplayName(dataNode).toLowerCase();
  return name.includes(SEARCH_TERM.toLowerCase());
}

/* ─────────────────────────────────────
   Tree Drawing
   ───────────────────────────────────── */
function clearTree() {
  treeHost.innerHTML = "";
}

function clearAllCollapsedFlags(node) {
  if (!node) return;
  if (node._collapsed) delete node._collapsed;
  if (node.children && node.children.length) {
    for (const c of node.children) clearAllCollapsedFlags(c);
  }
}

function applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform) {
  if (!svg || !zoom || !zoomLayer) return;
  let bbox = null;
  try {
    bbox = zoomLayer.node().getBBox();
  } catch {
    bbox = null;
  }
  const pad = 80;
  if (
    bbox &&
    isFinite(bbox.x) &&
    isFinite(bbox.y) &&
    isFinite(bbox.width) &&
    isFinite(bbox.height) &&
    bbox.width > 0 &&
    bbox.height > 0
  ) {
    zoom.translateExtent([
      [bbox.x - pad, bbox.y - pad],
      [bbox.x + bbox.width + pad, bbox.y + bbox.height + pad],
    ]);
  }

  svg.call(zoom);
  svg.call(zoom.transform, previousTransform || d3.zoomIdentity);
}

function clearDescendantCollapsedFlags(node) {
  if (!node || !node.children || !node.children.length) return;
  for (const c of node.children) {
    if (c._collapsed) delete c._collapsed;
    clearDescendantCollapsedFlags(c);
  }
}

function toggleNodeCollapsedState(dataNode) {
  if (!dataNode || !dataNode.children || !dataNode.children.length) return;
  if (dataNode._collapsed) {
    delete dataNode._collapsed;
    return;
  }
  dataNode._collapsed = true;
  // Keep only the highest collapsed node active in any subtree.
  clearDescendantCollapsedFlags(dataNode);
}

function drawTree() {
  clearTree();
  if (!NEWICK_RAW_TREE) return;
  computeLayouts();

  const layoutMode = CURRENT_LAYOUT_MODE || "rectangular";
  const colorMode = COLOR_MODE || "bars";
  const tc = getThemeColors();

  const container = d3.select("#tree_display");
  const width = treeHost.clientWidth || 800;
  const height = treeHost.clientHeight || 500;
  const margin = { top: 20, right: 80, bottom: 20, left: 80 };

  const svg = container
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");

  const zoomLayer = svg.append("g");
  const g = zoomLayer
    .append("g")
    .attr(
      "transform",
      layoutMode === "circular"
        ? `translate(${width / 2},${height / 2})`
        : `translate(${margin.left},${margin.top})`,
    );

  const zoom = d3
    .zoom()
    .scaleExtent([0.3, 8])
    .on("zoom", (event) => {
      zoomLayer.attr("transform", event.transform);
      LAST_TREE_TRANSFORM = event.transform;
    });
  const previousTransform = LAST_TREE_TRANSFORM || d3.zoomIdentity;
  LAST_TREE_SVG = svg;
  LAST_TREE_ZOOM = zoom;

  function hasOriginalChildren(d) {
    return !!(
      d &&
      d.data &&
      Array.isArray(d.data.children) &&
      d.data.children.length > 0
    );
  }
  function canCollapse(d) {
    return hasOriginalChildren(d);
  }
  function isNodeVisible(d) {
    if (!d) return false;
    var cur = d.parent;
    while (cur) {
      if (cur.data && cur.data._collapsed) return false;
      cur = cur.parent;
    }
    return true;
  }
  function hasVisibleChildren(d) {
    if (!d || !d.children || !d.children.length) return false;
    if (d.data && d.data._collapsed) return false;
    for (var i = 0; i < d.children.length; i++) {
      if (isNodeVisible(d.children[i])) return true;
    }
    return false;
  }
  function isVisibleLeaf(d) {
    return !hasVisibleChildren(d);
  }
  function nodeRadius(d) {
    return customNodeRadius(d);
  }
  function isHighlighted(d) {
    if (isSearchMatch(d.data)) return true;
    return NODE_CUSTOM.has(d.data) && NODE_CUSTOM.get(d.data).highlighted;
  }

  function highlightColor() {
    return "#f59e0b";
  }

  function clusterColor(cid) {
    if (cid == null || !CLUSTER_COLORS.length) return null;
    if (Number(cid) < 0) return tc.internal; // outlier: neutral color
    return CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
  }

  function nodeFill(d) {
    const name = d.data && d.data.name;
    if (colorMode === "bars") return tc.label;
    if (name && !hasOriginalChildren(d)) {
      const cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[name] : null;
      return clusterColor(cid) || tc.label;
    }
    const rep = representativeClusterIdFromData(d.data);
    return clusterColor(rep) || tc.label;
  }

  function nodeClusterId(d) {
    const name = d.data && d.data.name;
    if (name && !hasOriginalChildren(d)) {
      const cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[name] : null;
      return cid == null ? null : Number(cid);
    }
    const rep = representativeClusterIdFromData(d.data);
    return rep == null ? null : Number(rep);
  }

  function clusterVisible(d) {
    if (SELECTED_CLUSTER_IDS.size === 0) return true;
    const cid = nodeClusterId(d);
    return cid != null && SELECTED_CLUSTER_IDS.has(cid);
  }

  // ── CIRCULAR LAYOUT ──
  if (layoutMode === "circular") {
    const root = HIER_CIRC;
    const visibleNodes = root.descendants().filter((d) => isNodeVisible(d));
    const visibleLinks = root
      .links()
      .filter((l) => isNodeVisible(l.source) && isNodeVisible(l.target));
    const maxRadiusPx = Math.min(width, height) / 2 - 30;
    const maxR = d3.max(root.descendants(), (d) => d._radius || 0) || 1;
    const rScale = d3.scaleLinear().domain([0, maxR]).range([0, maxRadiusPx]);

    function polarToXY(angle, radius) {
      const a = angle - Math.PI / 2;
      const rr = rScale(radius);
      return [
        rr * Math.cos(a) * TREE_WIDTH_SCALE,
        rr * Math.sin(a) * TREE_HEIGHT_SCALE,
      ];
    }
    function nodeXY(d) {
      return polarToXY(d._angle || 0, d._radius || 0);
    }

    function radialElbowPath(link) {
      const s = link.source,
        t = link.target;
      const sa = s._angle || 0,
        sr = s._radius || 0;
      const ta = t._angle || 0,
        tr = t._radius || 0;
      const p0 = polarToXY(sa, sr),
        p1 = polarToXY(ta, sr),
        p2 = polarToXY(ta, tr);
      const delta = ta - sa,
        sweep = delta >= 0 ? 1 : 0;
      const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
      const arcR = rScale(sr);
      return `M${p0[0]},${p0[1]}A${arcR},${arcR} 0 ${largeArc},${sweep} ${p1[0]},${p1[1]}L${p2[0]},${p2[1]}`;
    }

    g.append("g")
      .selectAll(".tree-link")
      .data(visibleLinks)
      .enter()
      .append("path")
      .attr("class", "tree-link")
      .attr("fill", "none")
      .attr("stroke", tc.branch)
      .attr("stroke-opacity", 1)
      .attr("stroke-width", BRANCH_STROKE_WIDTH)
      .attr("d", radialElbowPath);

    const node = g
      .append("g")
      .selectAll(".tree-node")
      .data(visibleNodes)
      .enter()
      .append("g")
      .attr(
        "class",
        (d) =>
          "tree-node" + (isSearchMatch(d.data) ? " node-search-match" : ""),
      )
      .attr("transform", (d) => {
        const p = nodeXY(d);
        return `translate(${p[0]},${p[1]})`;
      });

    node
      .append("circle")
      .attr("r", nodeRadius)
      .attr("fill", nodeFill)
      .attr("stroke", (d) => (isHighlighted(d) ? highlightColor() : "none"))
      .attr("stroke-width", (d) => (isHighlighted(d) ? 2.5 : 0))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.22))
      .style("cursor", "pointer")
      .on("click", function (event, d) {
        if (!canCollapse(d)) return;
        toggleNodeCollapsedState(d.data);
        drawTree();
      })
      .on("contextmenu", showNodeContextMenu);

    node
      .filter((d) => isVisibleLeaf(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", LABEL_FONT_SIZE + "px")
      .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.label))
      .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
      .attr("transform", function (d) {
        const a = d._angle || 0;
        let deg = (a * 180) / Math.PI - 90;
        const left = deg > 90 || deg < -90;
        return `rotate(${deg}) translate(${nodeRadius(d) + LABEL_PAD},0) rotate(${left ? 180 : 0})`;
      })
      .style("text-anchor", (d) => {
        const deg = ((d._angle || 0) * 180) / Math.PI - 90;
        return deg > 90 || deg < -90 ? "end" : "start";
      })
      .text((d) => (SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : ""));

    node
      .filter((d) => !isVisibleLeaf(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", LABEL_FONT_SIZE - 1 + "px")
      .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.internal))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
      .attr("transform", (d) => {
        const a = d._angle || 0;
        let deg = (a * 180) / Math.PI - 90;
        return `rotate(${deg}) translate(0,${-(nodeRadius(d) + 5)})`;
      })
      .style("text-anchor", "middle")
      .text((d) => (SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : ""));

    // Show collapsed subtree leaf count next to collapsed node markers.
    node
      .filter((d) => d.data && d.data._collapsed && hasOriginalChildren(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", Math.max(8, LABEL_FONT_SIZE - 2) + "px")
      .attr("font-weight", "700")
      .attr("fill", tc.internal)
      .attr("transform", function (d) {
        const a = d._angle || 0;
        const deg = (a * 180) / Math.PI - 90;
        const left = deg > 90 || deg < -90;
        return `rotate(${deg}) translate(${nodeRadius(d) + 7},0) rotate(${left ? 180 : 0})`;
      })
      .style("text-anchor", function (d) {
        const deg = ((d._angle || 0) * 180) / Math.PI - 90;
        return deg > 90 || deg < -90 ? "end" : "start";
      })
      .text((d) => String(countLeafDescendantsData(d.data)));

    // Circular cluster visualization
    if (CURRENT_CLUSTERS && CLUSTER_COLORS.length) {
      const visLeaves = root
        .descendants()
        .filter((d) => isNodeVisible(d))
        .filter((d) => isVisibleLeaf(d))
        .sort((a, b) => (a._angle || 0) - (b._angle || 0));

      function visLeafCid(d) {
        const nm = d.data && d.data.name;
        const isOrigLeaf = nm && !(d.data.children && d.data.children.length);
        if (isOrigLeaf) return CURRENT_CLUSTERS[nm];
        return representativeClusterIdFromData(d.data);
      }

      if (colorMode === "bars" && visLeaves.length >= 2) {
        const angles = visLeaves.map((d) => d._angle || 0);
        const boundaries = new Array(visLeaves.length + 1);
        for (let i = 1; i < visLeaves.length; i++)
          boundaries[i] = (angles[i - 1] + angles[i]) / 2;
        boundaries[0] = angles[0] - (boundaries[1] - angles[0]);
        boundaries[visLeaves.length] =
          angles[visLeaves.length - 1] +
          (angles[visLeaves.length - 1] - boundaries[visLeaves.length - 1]);

        const leafLabelChars = visLeaves.reduce(
          (mx, d) => Math.max(mx, (d.data.name || "").length),
          0,
        );
        const estRadialLabelWidth = SHOW_LEAF_NAMES
          ? leafLabelChars * LABEL_FONT_SIZE * 0.55 +
            LABEL_PAD +
            LEAF_NODE_RADIUS
          : 10;
        const ringInner = maxRadiusPx + estRadialLabelWidth + 6;
        const ringOuter = ringInner + 14;
        const arc = d3.arc().innerRadius(ringInner).outerRadius(ringOuter);

        let runStart = 0,
          runCid = visLeafCid(visLeaves[0]);
        for (let i = 1; i <= visLeaves.length; i++) {
          const cid =
            i < visLeaves.length ? visLeafCid(visLeaves[i]) : Symbol("END");
          if (cid !== runCid) {
            if (runCid != null) {
              var arcColor = clusterColor(runCid) || tc.internal;
              var isOutlierArc = Number(runCid) < 0;
              var arcEl = g
                .append("path")
                .attr(
                  "d",
                  arc.startAngle(boundaries[runStart]).endAngle(boundaries[i]),
                )
                .attr("fill", isOutlierArc ? "none" : arcColor)
                .attr("stroke", isOutlierArc ? arcColor : "none")
                .attr("stroke-width", isOutlierArc ? 1.5 : 0)
                .attr("opacity", hasClusterFocus(runCid) ? 0.75 : 0.18);
              if (isOutlierArc) arcEl.attr("stroke-dasharray", "4,2");
            }
            runStart = i;
            runCid = i < visLeaves.length ? cid : null;
          }
        }
        g.selectAll(".tree-link, .tree-node").raise();
      }

      // Circular boxes: wedge-shaped arcs from MRCA radius to leaf tips
      if (colorMode === "boxes" && visLeaves.length >= 2) {
        var circClusterGroups = {};
        visLeaves.forEach(function (d) {
          var cid = visLeafCid(d);
          if (cid == null) return;
          if (!SHOW_OUTLIER_BOXES && Number(cid) < 0) return;
          if (!circClusterGroups[cid]) circClusterGroups[cid] = [];
          circClusterGroups[cid].push(d);
        });
        // Find MRCA for circular nodes
        function findMRCACirc(nodes) {
          if (!nodes.length) return null;
          if (nodes.length === 1) return nodes[0];
          function ancestors(n) {
            var a = [];
            var cur = n;
            while (cur) {
              a.push(cur);
              cur = cur.parent;
            }
            return a;
          }
          var common = ancestors(nodes[0]);
          for (var i = 1; i < nodes.length; i++) {
            var anc = new Set(ancestors(nodes[i]));
            common = common.filter(function (n) {
              return anc.has(n);
            });
          }
          return common[0];
        }
        var boxesG = g.append("g").attr("class", "cluster-boxes");
        for (var cid of Object.keys(circClusterGroups)) {
          var cLeaves = circClusterGroups[cid];
          var mrca = findMRCACirc(cLeaves);
          if (!mrca) continue;
          var leafAngles = cLeaves
            .map(function (d) {
              return d._angle || 0;
            })
            .sort(function (a, b) {
              return a - b;
            });
          var mrcaR = mrca._radius || 0;
          var maxLeafR =
            d3.max(cLeaves, function (d) {
              return d._radius || 0;
            }) || 0;

          // Angular extent with padding
          var angPad = BOX_PAD_V * 0.01; // convert px pad to radians approximately
          var minAngle = leafAngles[0] - angPad;
          var maxAngle = leafAngles[leafAngles.length - 1] + angPad;
          // If only one leaf, give a small angular spread
          if (cLeaves.length === 1) {
            var spread = 0.03;
            minAngle -= spread;
            maxAngle += spread;
          }

          var innerR = rScale(mrcaR) - BOX_PAD_H;
          var outerR =
            rScale(maxLeafR) + BOX_PAD_H + (SHOW_LEAF_NAMES ? 40 : 5);
          if (innerR < 0) innerR = 0;

          var isOutlierWedge = Number(cid) < 0;
          var wedgeColor = isOutlierWedge
            ? tc.internal
            : clusterColor(cid) || tc.internal;
          var wedgeArc = d3
            .arc()
            .innerRadius(innerR * TREE_WIDTH_SCALE)
            .outerRadius(outerR * TREE_WIDTH_SCALE)
            .startAngle(minAngle)
            .endAngle(maxAngle)
            .cornerRadius(BOX_CORNER_RADIUS);

          var wedge = boxesG
            .append("path")
            .attr("d", wedgeArc)
            .attr("fill", isOutlierWedge ? "none" : wedgeColor)
            .attr("opacity", hasClusterFocus(cid) ? BOX_ALPHA : 0.02)
            .attr("stroke", wedgeColor)
            .attr("stroke-width", 1.5)
            .attr(
              "stroke-opacity",
              hasClusterFocus(cid) ? Math.min(BOX_ALPHA * 3.5, 0.8) : 0.08,
            );
          if (isOutlierWedge) wedge.attr("stroke-dasharray", "5,3");

          // Label at the outer edge of the wedge
          if (SHOW_BOX_LABELS) {
            var midAngle = (minAngle + maxAngle) / 2;
            var labelR = outerR * TREE_WIDTH_SCALE + 4;
            var lx = labelR * Math.cos(midAngle - Math.PI / 2);
            var ly = labelR * Math.sin(midAngle - Math.PI / 2);
            var deg = (midAngle * 180) / Math.PI - 90;
            var flip = deg > 90 || deg < -90;
            boxesG
              .append("text")
              .attr("x", lx)
              .attr("y", ly)
              .attr("dy", "0.35em")
              .attr("font-size", CLUSTER_LABEL_FONT_SIZE + "px")
              .attr("font-weight", "600")
              .attr("fill", wedgeColor)
              .attr("font-style", isOutlierWedge ? "italic" : "normal")
              .attr("opacity", hasClusterFocus(cid) ? 1 : 0.25)
              .attr(
                "transform",
                "rotate(" +
                  (flip ? deg + 180 : deg) +
                  "," +
                  lx +
                  "," +
                  ly +
                  ")",
              )
              .style("text-anchor", flip ? "end" : "start")
              .text(
                BOX_LABEL_MAP[cid] || (isOutlierWedge ? "outlier" : "C" + cid),
              );
          }
        }
        boxesG.attr("pointer-events", "none");
        boxesG.lower();
      }
    }
    applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform);
    return;
  }

  // ── CARTESIAN (rectangular / cladogram) ──
  const root = HIER_CART;
  const allNodes = root.descendants().filter((d) => isNodeVisible(d));
  const allLinks = root
    .links()
    .filter((l) => isNodeVisible(l.source) && isNodeVisible(l.target));
  const maxYCart = d3.max(allNodes, (d) => d._y || 0) || 0;
  const leafNames = allNodes
    .filter((d) => !d.children || !d.children.length)
    .map((d) => (d.data.name || "").length);
  const maxLabelChars = d3.max(leafNames) || 0;
  const estLabelWidth = SHOW_LEAF_NAMES
    ? maxLabelChars * LABEL_FONT_SIZE * 0.6 + LABEL_PAD + LEAF_NODE_RADIUS
    : 20;
  const labelColumnX = maxYCart * TREE_WIDTH_SCALE + estLabelWidth + 12;

  g.append("g")
    .selectAll(".tree-link")
    .data(allLinks)
    .enter()
    .append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", tc.branch)
    .attr("stroke-opacity", 1)
    .attr("stroke-width", BRANCH_STROKE_WIDTH)
    .attr("d", function (d) {
      if (layoutMode === "rectangular") {
        return (
          "M" +
          d.source._y * TREE_WIDTH_SCALE +
          "," +
          d.source._x * TREE_HEIGHT_SCALE +
          "V" +
          d.target._x * TREE_HEIGHT_SCALE +
          "H" +
          d.target._y * TREE_WIDTH_SCALE
        );
      }
      return (
        "M" +
        d.source._y * TREE_WIDTH_SCALE +
        "," +
        d.source._x * TREE_HEIGHT_SCALE +
        "L" +
        d.target._y * TREE_WIDTH_SCALE +
        "," +
        d.target._x * TREE_HEIGHT_SCALE
      );
    });

  const node = g
    .append("g")
    .selectAll(".tree-node")
    .data(allNodes)
    .enter()
    .append("g")
    .attr(
      "class",
      (d) => "tree-node" + (isSearchMatch(d.data) ? " node-search-match" : ""),
    )
    .attr(
      "transform",
      (d) =>
        `translate(${(d._y || 0) * TREE_WIDTH_SCALE},${(d._x || 0) * TREE_HEIGHT_SCALE})`,
    );

  node
    .append("circle")
    .attr("r", nodeRadius)
    .attr("fill", nodeFill)
    .attr("stroke", (d) => (isHighlighted(d) ? highlightColor() : "none"))
    .attr("stroke-width", (d) => (isHighlighted(d) ? 2.5 : 0))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.22))
    .style("cursor", "pointer")
    .on("click", function (event, d) {
      if (!canCollapse(d)) return;
      toggleNodeCollapsedState(d.data);
      drawTree();
    })
    .on("contextmenu", showNodeContextMenu);

  // Leaf labels
  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", (d) => nodeRadius(d) + LABEL_PAD)
    .style("text-anchor", "start")
    .style("font-size", LABEL_FONT_SIZE + "px")
    .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.label))
    .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
    .text((d) => (SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : ""));

  // Internal node labels — offset below-right to avoid branch overlap
  node
    .filter((d) => !!(d.children && d.children.length))
    .append("text")
    .attr("dy", (d) => nodeRadius(d) + 12)
    .attr("x", (d) => nodeRadius(d) + 3)
    .style("text-anchor", "start")
    .style("font-size", LABEL_FONT_SIZE - 1 + "px")
    .style("font-style", "italic")
    .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.internal))
    .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
    .text((d) => (SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : ""));

  node
    .filter(
      (d) =>
        !!(
          d.data &&
          d.data._collapsed &&
          d.data.children &&
          d.data.children.length
        ),
    )
    .append("text")
    .attr("dy", "0.32em")
    .attr("x", (d) => nodeRadius(d) + 5)
    .style("text-anchor", "start")
    .style("font-size", Math.max(8, LABEL_FONT_SIZE - 2) + "px")
    .attr("font-weight", "700")
    .attr("fill", tc.internal)
    .attr("opacity", (d) => (clusterVisible(d) ? 0.95 : 0.35))
    .text((d) => String(countLeafDescendantsData(d.data)));

  // Cluster side bars
  if (colorMode === "bars") {
    const leafNodes = allNodes
      .filter((d) => isVisibleLeaf(d))
      .sort((a, b) => (a._x || 0) - (b._x || 0));
    const barX = labelColumnX,
      barW = 20;
    const ys = leafNodes.map((d) => (d._x || 0) * TREE_HEIGHT_SCALE);

    if (leafNodes.length === 1) {
      const cid = nodeClusterId(leafNodes[0]);
      if (cid != null) {
        var barColor = clusterColor(cid) || tc.internal;
        var barRect = g
          .append("rect")
          .attr("x", barX)
          .attr("y", ys[0] - 6)
          .attr("width", barW)
          .attr("height", 12)
          .attr("fill", Number(cid) < 0 ? "none" : barColor)
          .attr("stroke", Number(cid) < 0 ? barColor : "none")
          .attr("stroke-width", Number(cid) < 0 ? 1.5 : 0)
          .attr("opacity", 0.75);
        if (Number(cid) < 0) barRect.attr("stroke-dasharray", "4,2");
      }
    } else if (leafNodes.length > 1) {
      const boundaries = new Array(leafNodes.length + 1);
      for (let i = 1; i < leafNodes.length; i++)
        boundaries[i] = (ys[i - 1] + ys[i]) / 2;
      boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
      boundaries[leafNodes.length] =
        ys[leafNodes.length - 1] +
        (ys[leafNodes.length - 1] - boundaries[leafNodes.length - 1]);

      const barsG = g.append("g").attr("class", "cluster-bars");
      function leafCid(i) {
        return nodeClusterId(leafNodes[i]);
      }

      let runStart = 0,
        runCid = leafCid(0);
      for (let i = 1; i <= leafNodes.length; i++) {
        const cid = i < leafNodes.length ? leafCid(i) : Symbol("END");
        if (cid !== runCid) {
          if (runCid != null) {
            const y0 = boundaries[runStart],
              y1 = boundaries[i];
            var bColor = clusterColor(runCid) || tc.internal;
            var isOutlierBar = Number(runCid) < 0;
            var barEl = barsG
              .append("rect")
              .attr("x", barX)
              .attr("y", y0)
              .attr("width", barW)
              .attr("height", Math.max(1, y1 - y0))
              .attr("fill", isOutlierBar ? "none" : bColor)
              .attr("stroke", isOutlierBar ? bColor : "none")
              .attr("stroke-width", isOutlierBar ? 1.5 : 0)
              .attr("opacity", hasClusterFocus(runCid) ? 0.75 : 0.18);
            if (isOutlierBar) barEl.attr("stroke-dasharray", "4,2");
          }
          runStart = i;
          runCid = i < leafNodes.length ? cid : null;
        }
      }
      barsG.attr("pointer-events", "none");
    }
  }

  // Cluster boxes (MRCA mode) — draw coloured rectangles from MRCA to leaf tips
  if (colorMode === "boxes" && CURRENT_CLUSTERS && CLUSTER_COLORS.length) {
    // Group leaves by cluster id
    const clusterGroups = {};
    allNodes.forEach((d) => {
      if (!isVisibleLeaf(d)) return;
      const cid = nodeClusterId(d);
      if (cid == null) return;
      if (!SHOW_OUTLIER_BOXES && Number(cid) < 0) return;
      if (!clusterGroups[cid]) clusterGroups[cid] = [];
      clusterGroups[cid].push(d);
    });
    // Find MRCA for each cluster group
    function findMRCA(nodes) {
      if (!nodes.length) return null;
      if (nodes.length === 1) return nodes[0];
      // Collect ancestors for each node
      function ancestors(n) {
        const a = [];
        let cur = n;
        while (cur) {
          a.push(cur);
          cur = cur.parent;
        }
        return a;
      }
      let common = ancestors(nodes[0]);
      for (let i = 1; i < nodes.length; i++) {
        const anc = new Set(ancestors(nodes[i]));
        common = common.filter((n) => anc.has(n));
      }
      return common[0]; // first (deepest) common ancestor
    }
    // Estimate rightmost extent including leaf labels
    var estLabelPx = SHOW_LEAF_NAMES
      ? maxLabelChars * LABEL_FONT_SIZE * 0.6 +
        LABEL_PAD +
        LEAF_NODE_RADIUS +
        16
      : 16;
    var boxPadV = BOX_PAD_V; // vertical padding around topmost/bottommost leaf
    var boxPadL = BOX_PAD_H; // left padding before MRCA
    var boxesG = g.append("g").attr("class", "cluster-boxes");
    for (var cid of Object.keys(clusterGroups)) {
      var leaves = clusterGroups[cid];
      var mrca = findMRCA(leaves);
      if (!mrca) continue;
      var ys = leaves.map(function (d) {
        return (d._x || 0) * TREE_HEIGHT_SCALE;
      });
      var minY = Math.min.apply(null, ys) - boxPadV;
      var maxY = Math.max.apply(null, ys) + boxPadV;
      var mrcaX = (mrca._y || 0) * TREE_WIDTH_SCALE - boxPadL;
      var maxLeafX = Math.max.apply(
        null,
        leaves.map(function (d) {
          return (d._y || 0) * TREE_WIDTH_SCALE;
        }),
      );
      var boxRight = maxLeafX + estLabelPx;
      var adj = getBoxAdjust(cid);
      minY += adj.dy - adj.padY;
      maxY += adj.dy + adj.padY;
      mrcaX += adj.dx - adj.padX;
      boxRight += adj.dx + adj.padX;
      var boxW = boxRight - mrcaX;
      var isOutlier = Number(cid) < 0;
      var color = isOutlier
        ? tc.internal
        : CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
      var rect = boxesG
        .append("rect")
        .attr("x", mrcaX)
        .attr("y", minY)
        .attr("width", Math.max(boxW, 8))
        .attr("height", Math.max(maxY - minY, 4))
        .attr("rx", BOX_CORNER_RADIUS)
        .attr("ry", BOX_CORNER_RADIUS)
        .attr("fill", isOutlier ? "none" : color)
        .attr("opacity", hasClusterFocus(cid) ? BOX_ALPHA : 0.02)
        .attr("stroke", color)
        .attr("stroke-width", isOutlier ? 1.5 : 1.5)
        .attr(
          "stroke-opacity",
          hasClusterFocus(cid) ? Math.min(BOX_ALPHA * 3.5, 0.8) : 0.08,
        );
      if (isOutlier) {
        rect.attr("stroke-dasharray", "5,3");
      }
      // Cluster label on right edge of box
      if (SHOW_BOX_LABELS) {
        boxesG
          .append("text")
          .attr("x", boxRight + 4)
          .attr("y", (minY + maxY) / 2)
          .attr("dy", "0.35em")
          .attr("font-size", CLUSTER_LABEL_FONT_SIZE + "px")
          .attr("font-weight", "600")
          .attr("fill", color)
          .attr("font-style", isOutlier ? "italic" : "normal")
          .attr("opacity", hasClusterFocus(cid) ? 1 : 0.25)
          .text(BOX_LABEL_MAP[cid] || (isOutlier ? "outlier" : "C" + cid));
      }
    }
    boxesG.attr("pointer-events", "none");
    // Keep MRCA boxes beneath branches and labels.
    boxesG.lower();
  }

  // Branch-length axis
  var axisMaxX = d3.max(allNodes, (d) => d._x || 0) || 0;
  var axisMaxY = d3.max(allNodes, (d) => d._y || 0) || 0;
  var axisMaxBl = d3.max(allNodes, (d) => d.data._bl || 0) || 1;
  var axisYPos = axisMaxX * TREE_HEIGHT_SCALE + 20;
  var blScale = d3
    .scaleLinear()
    .domain([0, axisMaxBl])
    .range([0, axisMaxY * TREE_WIDTH_SCALE]);

  var axisG = g
    .append("g")
    .attr("class", "branch-length-axis")
    .attr("transform", "translate(0, " + axisYPos + ")")
    .call(d3.axisBottom(blScale).ticks(5));
  axisG
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", AXIS_FONT_SIZE + "px");
  axisG.selectAll("line").attr("stroke", tc.branch);
  axisG.selectAll("path").attr("stroke", tc.branch);
  axisG
    .append("text")
    .attr("x", (axisMaxY * TREE_WIDTH_SCALE) / 2)
    .attr("y", AXIS_FONT_SIZE + 20)
    .attr("text-anchor", "middle")
    .attr("font-size", AXIS_FONT_SIZE)
    .attr("fill", tc.internal)
    .text("Branch length");

  svg.attr("height", Math.max(height, axisYPos + 50 + margin.bottom));
  applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform);
}

/* ─────────────────────────────────────
   Draw tree into a specific container
   (used for comparison mode)
   ───────────────────────────────────── */
function drawTreeInto(hostEl, clusters) {
  hostEl.innerHTML = "";
  if (!NEWICK_RAW_TREE) return;

  const tc = getThemeColors();
  const hier = d3.hierarchy(NEWICK_RAW_TREE, getVisibleChildren);
  const width = hostEl.clientWidth || 400;
  const height = hostEl.clientHeight || 400;
  const margin = { top: 16, right: 60, bottom: 16, left: 60 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const maxBl = d3.max(hier.descendants(), (d) => d.data._bl || 0) || 1;
  const blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);

  d3.cluster().size([innerH, 1])(hier);
  hier.each((d) => {
    d._x = d.x;
    d._y = blToX(d.data._bl || 0);
  });

  const nClusters = clusters ? Math.max(0, ...Object.values(clusters)) + 1 : 0;
  const colors = nClusters > 0 ? generateClusterColors(nClusters) : [];

  const svg = d3
    .select(hostEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  const zoomLayer = svg.append("g");
  const g = zoomLayer
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);
  svg.call(
    d3
      .zoom()
      .scaleExtent([0.3, 8])
      .on("zoom", (event) => {
        zoomLayer.attr("transform", event.transform);
      }),
  );

  g.append("g")
    .selectAll(".tree-link")
    .data(hier.links())
    .enter()
    .append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", tc.branch)
    .attr("stroke-opacity", 1)
    .attr("stroke-width", BRANCH_STROKE_WIDTH)
    .attr(
      "d",
      (d) =>
        "M" +
        d.source._y +
        "," +
        d.source._x +
        "V" +
        d.target._x +
        "H" +
        d.target._y,
    );

  const allNodes = hier.descendants();
  const node = g
    .append("g")
    .selectAll(".tree-node")
    .data(allNodes)
    .enter()
    .append("g")
    .attr("class", "tree-node")
    .attr("transform", (d) => `translate(${d._y},${d._x})`);

  node
    .append("circle")
    .attr("r", (d) => (d.children && d.children.length ? 1.5 : 2.5))
    .attr("fill", function (d) {
      const name = d.data && d.data.name;
      if (name && !(d.children && d.children.length) && clusters) {
        const cid = clusters[name];
        if (cid != null && colors.length) return colors[cid % colors.length];
      }
      return tc.label;
    });

  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", 6)
    .style("text-anchor", "start")
    .style("font-size", "8px")
    .attr("fill", tc.label)
    .text((d) => d.data.name || "");

  // Side bars
  if (clusters && colors.length) {
    const leaves = allNodes
      .filter((d) => !(d.children && d.children.length))
      .sort((a, b) => a._x - b._x);
    if (leaves.length > 1) {
      const estLabel =
        d3.max(leaves, (d) => (d.data.name || "").length) * 5 + 12;
      const barX = (d3.max(allNodes, (d) => d._y) || 0) + estLabel;
      const ys = leaves.map((d) => d._x);
      const boundaries = new Array(leaves.length + 1);
      for (let i = 1; i < leaves.length; i++)
        boundaries[i] = (ys[i - 1] + ys[i]) / 2;
      boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
      boundaries[leaves.length] =
        ys[leaves.length - 1] +
        (ys[leaves.length - 1] - boundaries[leaves.length - 1]);

      let runStart = 0,
        runCid = clusters[leaves[0].data.name];
      for (let i = 1; i <= leaves.length; i++) {
        const cid =
          i < leaves.length ? clusters[leaves[i].data.name] : Symbol("END");
        if (cid !== runCid) {
          if (runCid != null) {
            g.append("rect")
              .attr("x", barX)
              .attr("y", boundaries[runStart])
              .attr("width", 14)
              .attr("height", Math.max(1, boundaries[i] - boundaries[runStart]))
              .attr("fill", colors[runCid % colors.length])
              .attr("opacity", 0.7);
          }
          runStart = i;
          runCid = i < leaves.length ? cid : null;
        }
      }
    }
  }
}

/* ─────────────────────────────────────
   Draw one tree with side-by-side bars
   ───────────────────────────────────── */
function drawComparisonBarsInto(hostEl, comparisons) {
  hostEl.innerHTML = "";
  if (!NEWICK_RAW_TREE) return;

  const tc = getThemeColors();
  const hier = d3.hierarchy(NEWICK_RAW_TREE, getVisibleChildren);
  const width = hostEl.clientWidth || 840;
  const height = hostEl.clientHeight || 500;
  const colW = 16;
  const gap = 6;
  const nCols = Math.max(0, (comparisons || []).length);
  const margin = {
    top: 28,
    right: Math.max(170, nCols * (colW + gap) + 46),
    bottom: 64,
    left: 60,
  };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const maxBl = d3.max(hier.descendants(), (d) => d.data._bl || 0) || 1;
  const blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);

  d3.cluster().size([innerH, 1])(hier);
  hier.each((d) => {
    d._x = d.x;
    d._y = blToX(d.data._bl || 0);
  });

  const leaves = hier
    .descendants()
    .filter((d) => !(d.children && d.children.length))
    .sort((a, b) => a._x - b._x);
  const leafLabels = leaves.map((d) => (d.data.name || "").length);
  const maxLabelChars = d3.max(leafLabels) || 0;
  const estLabelWidth = SHOW_LEAF_NAMES
    ? maxLabelChars * LABEL_FONT_SIZE * 0.6 + LABEL_PAD + LEAF_NODE_RADIUS
    : 20;

  const svg = d3
    .select(hostEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  const zoomLayer = svg.append("g");
  const g = zoomLayer
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);
  svg.call(
    d3
      .zoom()
      .scaleExtent([0.3, 8])
      .on("zoom", (event) => {
        zoomLayer.attr("transform", event.transform);
      }),
  );

  g.append("g")
    .selectAll(".tree-link")
    .data(hier.links())
    .enter()
    .append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", tc.branch)
    .attr("stroke-opacity", 1)
    .attr("stroke-width", BRANCH_STROKE_WIDTH)
    .attr(
      "d",
      (d) =>
        "M" +
        d.source._y +
        "," +
        d.source._x +
        "V" +
        d.target._x +
        "H" +
        d.target._y,
    );

  const node = g
    .append("g")
    .selectAll(".tree-node")
    .data(hier.descendants())
    .enter()
    .append("g")
    .attr("class", "tree-node")
    .attr("transform", (d) => `translate(${d._y},${d._x})`);

  node
    .append("circle")
    .attr("r", (d) =>
      d.children && d.children.length ? INTERNAL_NODE_RADIUS : LEAF_NODE_RADIUS,
    )
    .attr("fill", tc.label)
    .attr("stroke", "none");

  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", LEAF_NODE_RADIUS + LABEL_PAD)
    .style("text-anchor", "start")
    .style("font-size", LABEL_FONT_SIZE + "px")
    .attr("fill", tc.label)
    .text((d) => (SHOW_LEAF_NAMES ? d.data.name || "" : ""));

  const barX =
    (d3.max(hier.descendants(), (d) => d._y) || 0) + estLabelWidth + 8;
  const ys = leaves.map((d) => d._x);
  if (leaves.length > 1) {
    const boundaries = new Array(leaves.length + 1);
    for (let i = 1; i < leaves.length; i++)
      boundaries[i] = (ys[i - 1] + ys[i]) / 2;
    boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
    boundaries[leaves.length] =
      ys[leaves.length - 1] +
      (ys[leaves.length - 1] - boundaries[leaves.length - 1]);

    function drawBarColumn(map, x, colIdx) {
      if (!map) return;
      const maxCid = Math.max(-1, ...Object.values(map).map(Number));
      const colors = maxCid >= 0 ? generateClusterColors(maxCid + 1) : [];
      let runStart = 0;
      let runCid = map[leaves[0].data.name];
      for (let i = 1; i <= leaves.length; i++) {
        const cid =
          i < leaves.length ? map[leaves[i].data.name] : Symbol("END");
        if (cid !== runCid) {
          if (runCid != null) {
            const nCid = Number(runCid);
            const isOutlier = nCid < 0;
            const fillColor = isOutlier
              ? tc.internal
              : colors.length
                ? colors[nCid % colors.length]
                : tc.label;
            g.append("rect")
              .attr("x", x)
              .attr("y", boundaries[runStart])
              .attr("width", colW)
              .attr("height", Math.max(1, boundaries[i] - boundaries[runStart]))
              .attr("fill", fillColor)
              .attr("opacity", 0.76)
              .attr("pointer-events", "none");
          }
          runStart = i;
          runCid = i < leaves.length ? cid : null;
        }
      }

      const title =
        comparisons[colIdx] && comparisons[colIdx].title
          ? comparisons[colIdx].title
          : "Bar " + (colIdx + 1);
      const shortTitle = title.length > 14 ? title.slice(0, 12) + ".." : title;
      g.append("text")
        .attr("x", x + colW / 2)
        .attr("y", -8)
        .attr("text-anchor", "middle")
        .attr("fill", tc.internal)
        .attr("font-size", "10px")
        .text(shortTitle);
    }

    for (let i = 0; i < (comparisons || []).length; i++) {
      drawBarColumn(comparisons[i].clusters, barX + i * (colW + gap), i);
    }
  }

  // Branch-length axis (same semantics as main tree view)
  const axisYPos = innerH + 20;
  const axisG = g
    .append("g")
    .attr("class", "branch-length-axis")
    .attr("transform", "translate(0, " + axisYPos + ")")
    .call(d3.axisBottom(blToX).ticks(5));
  axisG
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", AXIS_FONT_SIZE + "px");
  axisG.selectAll("line").attr("stroke", tc.branch);
  axisG.selectAll("path").attr("stroke", tc.branch);
  axisG
    .append("text")
    .attr("x", innerW / 2)
    .attr("y", AXIS_FONT_SIZE + 20)
    .attr("text-anchor", "middle")
    .attr("font-size", AXIS_FONT_SIZE)
    .attr("fill", tc.internal)
    .text("Branch length");
}

/* ─────────────────────────────────────
   API Call
   ───────────────────────────────────── */
async function apiPostJson(url, payload) {
  const res = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const text = await res.text();
  let data = null;
  try {
    data = text ? JSON.parse(text) : null;
  } catch {
    /* ignore */
  }
  if (!res.ok) {
    throw new Error(
      data && data.detail ? data.detail : `Request failed (${res.status})`,
    );
  }
  return data;
}

/* ─────────────────────────────────────
   Get Current Mode
   ───────────────────────────────────── */
function getCurrentMode() {
  const activeBtn = document.querySelector("#mode-selector .mode-btn.active");
  return activeBtn ? activeBtn.dataset.mode : "global";
}

/* ─────────────────────────────────────
   Read a parameter value
   ───────────────────────────────────── */
function readIntParam(id) {
  var el = document.getElementById(id);
  if (!el) return null;
  var raw = (el.value || "").trim();
  if (!raw) return null;
  var v = parseInt(raw, 10);
  return isNaN(v) ? null : v;
}

function readFloatParam(id) {
  var el = document.getElementById(id);
  if (!el) return null;
  var raw = (el.value || "").trim();
  if (!raw) return null;
  var v = parseFloat(raw);
  return isNaN(v) ? null : v;
}

function readCheckParam(id) {
  var el = document.getElementById(id);
  return el ? el.checked : false;
}

function readSelectParam(id) {
  var el = document.getElementById(id);
  return el ? el.value : null;
}

/* ─────────────────────────────────────
   Run PhytClust
   ───────────────────────────────────── */
async function runPhytClust() {
  if (isRunning) return;
  var newickText = (newickEl.value || "").trim();
  if (!newickText) {
    showToast("Please upload or paste a Newick tree.", "danger", 4000);
    return;
  }

  var numSamples = estimateLeafCount(newickText);
  var mode = getCurrentMode();

  // Read params
  var kVal = readIntParam("extra-k");
  var outgroupVal =
    (extraOutgroupEl ? (extraOutgroupEl.value || "").trim() : "") || null;
  var rootTaxonVal =
    (extraRootTaxonEl ? (extraRootTaxonEl.value || "").trim() : "") || null;
  var topNVal = readIntParam("extra-topn");
  var binsVal = readIntParam("extra-bins");
  var maxKVal = readIntParam("extra-maxk");
  var maxKLimitVal = readFloatParam("extra-maxklimit");
  var lambdaVal = readFloatParam("extra-lambda");
  var minClusterVal = readIntParam("extra-min-cluster-size");

  if (outgroupVal && !isOutgroupInNewick(newickText, outgroupVal)) {
    showStatus("Outgroup not found in Newick.", "danger");
    return;
  }

  if (extraResolutionEl) extraResolutionEl.checked = mode === "resolution";

  // Build payload
  var payload = { newick: newickText, mode: mode };
  if (mode === "k") {
    if (kVal === null) {
      showToast("Please enter a value for k.", "danger");
      return;
    }
    payload.k = kVal;
  }
  if (outgroupVal) payload.outgroup = outgroupVal;
  if (rootTaxonVal) payload.root_taxon = rootTaxonVal;
  if (topNVal !== null) payload.top_n = topNVal;
  if (binsVal !== null) payload.num_bins = binsVal;
  if (maxKVal !== null) payload.max_k = maxKVal;
  if (maxKLimitVal !== null) payload.max_k_limit = maxKLimitVal;
  if (lambdaVal !== null) payload.lambda_weight = lambdaVal;
  if (minClusterVal !== null) payload.min_cluster_size = minClusterVal;
  if (mode === "resolution") payload.by_resolution = true;

  // New params
  if (readCheckParam("extra-compute-all")) payload.compute_all_clusters = true;
  if (readCheckParam("extra-use-support")) payload.use_branch_support = true;
  var minSup = readFloatParam("extra-min-support");
  if (minSup !== null) payload.min_support = minSup;
  var supW = readFloatParam("extra-support-weight");
  if (supW !== null) payload.support_weight = supW;

  // Outlier config
  var outlierThresh = readIntParam("extra-outlier-threshold");
  if (outlierThresh !== null) payload.outlier_size_threshold = outlierThresh;
  if (readCheckParam("extra-outlier-prefer-fewer"))
    payload.outlier_prefer_fewer = true;
  var ratioMode = readSelectParam("extra-outlier-ratio-mode");
  if (ratioMode && ratioMode !== "exp") payload.outlier_ratio_mode = ratioMode;

  // Polytomy config
  payload.optimize_polytomies = readCheckParam("extra-optimize-polytomies");
  if (readCheckParam("extra-no-split-zero"))
    payload.no_split_zero_length = true;

  // Peak config
  var rankMode = readSelectParam("extra-ranking-mode");
  if (rankMode) payload.ranking_mode = rankMode;
  var minProm = readFloatParam("extra-min-prominence");
  if (minProm !== null) payload.min_prominence = minProm;
  if (readCheckParam("extra-relative-prom"))
    payload.use_relative_prominence = true;

  showStatus("Running PhytClust...", "info");
  resultEl.textContent = "Running PhytClust...";
  clearTree();

  var runBtn = document.getElementById("btn-run");
  try {
    isRunning = true;
    if (runBtn) {
      runBtn.disabled = true;
      runBtn.innerHTML = '<span class="spinner"></span> Running...';
    }

    const t0 = performance.now();
    const data = await apiPostJson("/api/run", payload);
    const dt = (performance.now() - t0) / 1000;

    latestApiData = data;
    latestRunId = data.run_id || null;
    showStatus(`Finished in ${dt.toFixed(2)}s`, "success");
    resultEl.textContent = JSON.stringify(data, null, 2);

    // Update leaf count label
    var lcLabel = document.getElementById("leaf-count-label");
    if (lcLabel && data.newick)
      lcLabel.textContent = estimateLeafCount(data.newick) + " leaves";

    if (data.newick) {
      NEWICK_RAW_TREE = parseNewick(data.newick);
      accumulateBranchLength(NEWICK_RAW_TREE);
      computeLayouts();
      populateClusterSelector(data);
      try {
        populateCompareSelectors(data);
      } catch (cmpErr) {
        console.warn("Compare panel init error:", cmpErr);
        showToast("Compare panel failed to initialize.", "danger", 2500);
      }

      const clusterMap =
        data.clusters && data.clusters.length > 0 ? data.clusters[0] : {};
      CURRENT_CLUSTERS = clusterMap;
      if (Object.keys(clusterMap).length > 0) {
        CLUSTER_COLORS = generateClusterColors(
          Math.max(...Object.values(clusterMap)) + 1,
        );
      } else {
        CLUSTER_COLORS = [];
        showStatus("No clusters found.", "danger");
      }
      updateClusterEditorAvailability();
      drawTree();
    } else {
      clearTree();
      showStatus("No Newick tree returned by API.", "danger");
    }

    // Draw mini scores panel on tree page
    try {
      drawMiniScores(data);
    } catch (e) {
      console.warn("Mini scores error:", e);
    }

    // Draw optimal k plot
    try {
      if (data && data.scores) {
        var plotHost = document.getElementById("optimalk_plot");
        if (plotHost) plotHost.innerHTML = "";
        latestOptimalKData = data;
        drawOptimalK(latestOptimalKData);
      } else {
        var plotHost2 = document.getElementById("optimalk_plot");
        if (plotHost2)
          plotHost2.innerHTML =
            '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);font-size:14px;">Score plot unavailable in fixed k mode.</div>';
      }
    } catch (plotErr) {
      console.warn("Plot rendering error:", plotErr);
    }
  } catch (e) {
    console.error(e);
    showStatus("Error: " + e.message, "danger");
    resultEl.textContent = "Error: " + e;
    latestOptimalKData = null;
    latestApiData = null;
  } finally {
    isRunning = false;
    if (runBtn) {
      runBtn.disabled = false;
      runBtn.innerHTML =
        '<svg width="14" height="14" viewBox="0 0 16 16" fill="currentColor"><polygon points="3,1 13,8 3,15"/></svg> Run PhytClust';
    }
  }
}

/* ─────────────────────────────────────
   Cluster Selector (multiple results)
   ───────────────────────────────────── */
let CLUSTER_VIEW_MODE = "peaks";
let COMPARE_CONFIGS = [];
let COMPARE_CONFIG_NEXT_ID = 1;

function getCurrentClusterIds() {
  if (!CURRENT_CLUSTERS) return [];
  const ids = new Set();
  Object.keys(CURRENT_CLUSTERS).forEach((leaf) => {
    const cid = Number(CURRENT_CLUSTERS[leaf]);
    if (isNaN(cid)) return;
    ids.add(cid);
  });
  return Array.from(ids).sort((a, b) => a - b);
}

function defaultClusterLabel(cid) {
  return Number(cid) < 0 ? "outlier" : "C" + cid;
}

function updateClusterEditorAvailability() {
  const btn = document.getElementById("btn-edit-clusters");
  if (!btn) return;
  const hasAnyClusters = getCurrentClusterIds().length > 0;
  btn.style.display = hasAnyClusters ? "inline-flex" : "none";
}

function renderClusterEditorRows() {
  const body = document.getElementById("cluster-editor-body");
  const filterInput = document.getElementById("cluster-filter-text");
  const showOutliersCb = document.getElementById(
    "cluster-filter-show-outliers",
  );
  if (!body) return;
  const ids = getCurrentClusterIds();
  const filter = (filterInput && filterInput.value ? filterInput.value : "")
    .toLowerCase()
    .trim();
  const includeOutliers = showOutliersCb ? showOutliersCb.checked : true;

  const visibleIds = ids.filter((cid) => {
    if (!includeOutliers && Number(cid) < 0) return false;
    if (!filter) return true;
    const label = (BOX_LABEL_MAP[cid] || defaultClusterLabel(cid))
      .toLowerCase()
      .trim();
    return String(cid).includes(filter) || label.includes(filter);
  });

  body.innerHTML = "";

  if (!visibleIds.length) {
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = 2;
    td.className = "cluster-editor-empty";
    td.textContent = ids.length
      ? "No clusters match this filter."
      : "No clusters available. Run clustering first.";
    tr.appendChild(td);
    body.appendChild(tr);
    return;
  }

  visibleIds.forEach((cid) => {
    const tr = document.createElement("tr");

    const tdId = document.createElement("td");
    tdId.textContent = String(cid);

    const tdLabel = document.createElement("td");
    const input = document.createElement("input");
    input.type = "text";
    input.setAttribute("data-cid", String(cid));
    input.value = BOX_LABEL_MAP[cid] || defaultClusterLabel(cid);
    tdLabel.appendChild(input);

    tr.appendChild(tdId);
    tr.appendChild(tdLabel);
    body.appendChild(tr);
  });
}

function openClusterEditor() {
  const modal = document.getElementById("clusterEditorModal");
  const bg = document.getElementById("clusterEditorBackdrop");
  const labelSizeInput = document.getElementById("cluster-label-size-input");
  const filterInput = document.getElementById("cluster-filter-text");
  const showOutliersCb = document.getElementById(
    "cluster-filter-show-outliers",
  );
  if (filterInput) filterInput.value = "";
  if (showOutliersCb) showOutliersCb.checked = true;
  renderClusterEditorRows();
  if (labelSizeInput) {
    labelSizeInput.value = String(Math.round(CLUSTER_LABEL_FONT_SIZE));
  }
  if (modal) modal.classList.add("show");
  if (bg) bg.classList.add("show");
}

function closeClusterEditor() {
  const modal = document.getElementById("clusterEditorModal");
  const bg = document.getElementById("clusterEditorBackdrop");
  if (modal) modal.classList.remove("show");
  if (bg) bg.classList.remove("show");
}

function saveClusterEditor() {
  const body = document.getElementById("cluster-editor-body");
  const labelSizeInput = document.getElementById("cluster-label-size-input");
  if (!body) return;

  body.querySelectorAll("input[data-cid]").forEach((input) => {
    const cid = Number(input.getAttribute("data-cid"));
    const val = (input.value || "").trim();
    const fallback = defaultClusterLabel(cid);
    if (!val || val === fallback) delete BOX_LABEL_MAP[cid];
    else BOX_LABEL_MAP[cid] = val;
  });

  if (labelSizeInput) {
    const v = parseFloat(labelSizeInput.value);
    if (!isNaN(v)) {
      CLUSTER_LABEL_FONT_SIZE = Math.min(36, Math.max(6, v));
    }
  }

  drawTree();
  closeClusterEditor();
  showToast("Cluster labels updated.", "success", 1800);
}

function populateClusterSelector(data) {
  const controls = document.getElementById("cluster-select-controls");
  const selectEl = document.getElementById("cluster-select");
  const labelEl = document.getElementById("cluster-select-label");
  const toggleEl = document.getElementById("cluster-view-toggle");
  if (!controls || !selectEl) return;

  const hasPeaks = (data.clusters || []).length > 1;
  const hasAll = !!(data.all_clusters && data.all_clusters.length > 1);

  if (toggleEl) toggleEl.style.display = hasAll ? "inline-flex" : "none";

  const peakKs = data.k_values || data.ks;

  const hasAny =
    (data.clusters && data.clusters.length > 0) ||
    (data.all_clusters && data.all_clusters.length > 0);
  if (!hasAny) {
    controls.classList.remove("visible");
    if (toggleEl) toggleEl.style.display = "none";
    updateClusterEditorAvailability();
    return;
  }

  if (hasAll && CLUSTER_VIEW_MODE === "all") {
    _fillClusterSelect(
      selectEl,
      labelEl,
      data.all_clusters,
      data.all_ks,
      false,
    );
  } else {
    if (data.clusters && data.clusters.length) {
      CLUSTER_VIEW_MODE = "peaks";
      _fillClusterSelect(selectEl, labelEl, data.clusters, peakKs, true);
    } else {
      CLUSTER_VIEW_MODE = "all";
      _fillClusterSelect(
        selectEl,
        labelEl,
        data.all_clusters || [],
        data.all_ks || [],
        false,
      );
    }
  }
  controls.classList.add("visible");
  updateClusterEditorAvailability();
}

function _fillClusterSelect(selectEl, labelEl, clusters, ks, showRank) {
  selectEl.innerHTML = "";
  for (let i = 0; i < clusters.length; i++) {
    const opt = document.createElement("option");
    opt.value = i;
    const kLabel = ks && ks[i] != null ? ks[i] : i + 1;
    opt.textContent = showRank
      ? `k = ${kLabel} (rank ${i + 1})`
      : `k = ${kLabel}`;
    selectEl.appendChild(opt);
  }
  if (labelEl) labelEl.textContent = `1 / ${clusters.length}`;
}

function switchCluster(idx) {
  if (!latestApiData) return;
  let clusters, ks;
  if (CLUSTER_VIEW_MODE === "all" && latestApiData.all_clusters) {
    clusters = latestApiData.all_clusters;
    ks = latestApiData.all_ks;
  } else {
    clusters = latestApiData.clusters || [];
    ks = latestApiData.k_values || latestApiData.ks || [];
  }
  if (idx < 0 || idx >= clusters.length) return;
  CURRENT_CLUSTERS = clusters[idx];
  if (Object.keys(CURRENT_CLUSTERS).length > 0) {
    CLUSTER_COLORS = generateClusterColors(
      Math.max(...Object.values(CURRENT_CLUSTERS)) + 1,
    );
  } else {
    CLUSTER_COLORS = [];
  }
  updateClusterEditorAvailability();
  drawTree();
  const labelEl = document.getElementById("cluster-select-label");
  if (labelEl) labelEl.textContent = `${idx + 1} / ${clusters.length}`;
}

function cycleCluster(delta) {
  var selectEl = document.getElementById("cluster-select");
  if (!selectEl || !selectEl.options.length) return;
  var cur = selectEl.selectedIndex;
  var next = cur + delta;
  if (next < 0) next = selectEl.options.length - 1;
  if (next >= selectEl.options.length) next = 0;
  selectEl.selectedIndex = next;
  switchCluster(parseInt(selectEl.value, 10));
}

function toggleClusterViewMode() {
  if (!latestApiData) return;
  CLUSTER_VIEW_MODE = CLUSTER_VIEW_MODE === "peaks" ? "all" : "peaks";
  populateClusterSelector(latestApiData);
  switchCluster(0);
  const toggleEl = document.getElementById("cluster-view-toggle");
  if (toggleEl)
    toggleEl.textContent =
      CLUSTER_VIEW_MODE === "all" ? "Show peaks only" : "Show all k";
}

/* ─────────────────────────────────────
   Comparison Mode
   ───────────────────────────────────── */
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
  walk(NEWICK_RAW_TREE);
  return out;
}

function normalizeClustersForLeaves(rawMap) {
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

function getAvailableCompareKs(data) {
  if (!data) return [];
  if (data.all_ks && data.all_ks.length) return data.all_ks.slice();
  return (data.k_values || data.ks || []).slice();
}

function getClustersForK(k) {
  if (!latestApiData) return null;
  if (latestApiData.all_clusters && latestApiData.all_ks) {
    const idx = latestApiData.all_ks.indexOf(k);
    if (idx >= 0) return latestApiData.all_clusters[idx];
  }
  const peakKs = latestApiData.k_values || latestApiData.ks;
  if (peakKs && latestApiData.clusters) {
    const idx = peakKs.indexOf(k);
    if (idx >= 0) return latestApiData.clusters[idx];
  }
  return null;
}

function parseClusterTableText(rawText) {
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
      // Gracefully skip one potential header row.
      if (idx === 0) return;
      invalidRows += 1;
      return;
    }
    map[leaf] = cid;
    parsedRows += 1;
  });

  return { map: parsedRows ? map : null, invalidRows };
}

function makeCompareConfig() {
  const ks = getAvailableCompareKs(latestApiData);
  const k = ks.length ? Number(ks[0]) : null;
  const id = COMPARE_CONFIG_NEXT_ID++;
  return {
    id,
    title: "Bar " + id,
    source: "k",
    k,
    fileName: "",
    fileClusters: null,
  };
}

function ensureCompareConfigs(data) {
  const ks = getAvailableCompareKs(data);
  if (!COMPARE_CONFIGS.length) {
    const cfg = makeCompareConfig();
    cfg.title = ks.length ? "k=" + ks[0] : "Bar 1";
    cfg.k = ks.length ? Number(ks[0]) : null;
    COMPARE_CONFIGS = [cfg];
  }

  COMPARE_CONFIGS.forEach((cfg) => {
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

function escapeHtmlAttr(value) {
  return String(value || "")
    .replace(/&/g, "&amp;")
    .replace(/"/g, "&quot;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

function renderCompareConfigList() {
  const host = document.getElementById("compare-config-list");
  if (!host) return;
  if (!latestApiData) {
    host.innerHTML = "";
    return;
  }

  const ks = getAvailableCompareKs(latestApiData);
  host.innerHTML = "";

  COMPARE_CONFIGS.forEach((cfg) => {
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
      COMPARE_CONFIGS = COMPARE_CONFIGS.filter((x) => x.id !== cfg.id);
      if (!COMPARE_CONFIGS.length) COMPARE_CONFIGS = [makeCompareConfig()];
      renderCompareConfigList();
      drawComparison();
    });
  });
}

function populateCompareSelectors(data) {
  ensureCompareConfigs(data);
  renderCompareConfigList();
}

function drawComparison() {
  if (!latestApiData || !NEWICK_RAW_TREE) return;

  const leftHost = document.getElementById("compare-tree-left");
  if (!leftHost) return;

  const comparisons = [];
  COMPARE_CONFIGS.forEach((cfg) => {
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

/* ─────────────────────────────────────
   Optimal k Plot
   ───────────────────────────────────── */
function drawOptimalK(data) {
  if (!data) data = latestOptimalKData;
  var plotEl = document.getElementById("optimalk_plot");
  if (!plotEl) return;

  var scores = data && data.scores ? data.scores : [];
  var peaks = data && data.peaks ? data.peaks : [];
  const tc = getThemeColors();

  if (!Array.isArray(scores) || scores.length === 0) {
    plotEl.innerHTML =
      '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);">No scores available to plot.</div>';
    window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
    return;
  }

  plotEl.innerHTML = "";
  var dataPoints = [];
  for (let i = 0; i < scores.length - 1; i++)
    dataPoints.push({ k: i + 2, score: scores[i + 1] });

  var peakPoints = peaks
    .map((k) => {
      var idx = k - 2;
      return idx >= 0 && idx < dataPoints.length
        ? { k, score: dataPoints[idx].score }
        : null;
    })
    .filter((d) => d);

  var width = plotEl.clientWidth || 700;
  var height = plotEl.clientHeight || 420;
  var margin = { top: 24, right: 24, bottom: 48, left: 58 };
  var innerWidth = width - margin.left - margin.right;
  var innerHeight = height - margin.top - margin.bottom;

  var svg = d3
    .select(plotEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  var g = svg
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  var axisModeEl = document.getElementById("axis-mode");
  var defaultMode = scores.length > 50 ? "log" : "normal";
  var mode = defaultMode;
  if (axisModeEl) {
    if (!axisModeEl.__initialized) {
      axisModeEl.value = defaultMode;
      axisModeEl.__initialized = true;
    }
    if (axisModeEl.__userOverride) mode = axisModeEl.value || defaultMode;
    else axisModeEl.value = defaultMode;
  }

  var xScale, xAxis;
  if (mode === "log") {
    xScale = d3
      .scaleLog()
      .domain([2, d3.max(dataPoints, (d) => d.k)])
      .range([0, innerWidth]);
    xAxis = d3
      .axisBottom(xScale)
      .tickValues(dataPoints.map((d) => d.k))
      .tickFormat(d3.format("d"));
  } else {
    xScale = d3
      .scaleBand()
      .domain(dataPoints.map((d) => d.k))
      .range([0, innerWidth])
      .padding(0.2);
    var nPts = dataPoints.length;
    var dtick = Math.max(1, Math.ceil(nPts / 10));
    xAxis = d3
      .axisBottom(xScale)
      .tickValues(dataPoints.map((d) => d.k).filter((d, i) => i % dtick === 0))
      .tickFormat(d3.format("d"));
  }

  const yScale = d3
    .scaleLinear()
    .domain([
      d3.min(dataPoints, (d) => d.score),
      d3.max(dataPoints, (d) => d.score),
    ])
    .nice()
    .range([innerHeight, 0]);

  // Grid
  g.append("g")
    .attr("class", "grid")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(d3.axisBottom(xScale).tickSize(-innerHeight).tickFormat(""))
    .selectAll("line")
    .attr("stroke", tc.border)
    .attr("stroke-dasharray", "2,2");
  g.append("g")
    .attr("class", "grid")
    .call(d3.axisLeft(yScale).tickSize(-innerWidth).tickFormat(""))
    .selectAll("line")
    .attr("stroke", tc.border)
    .attr("stroke-dasharray", "2,2");
  g.selectAll(".grid .domain").remove();

  // Axes
  var xAxisG = g
    .append("g")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(xAxis);
  var yAxisG = g.append("g").call(d3.axisLeft(yScale).ticks(6));
  [xAxisG, yAxisG].forEach((ag) => {
    ag.selectAll("text").attr("fill", tc.internal).style("font-size", "11px");
    ag.selectAll("line").attr("stroke", tc.branch);
    ag.selectAll("path").attr("stroke", tc.branch);
  });

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 38)
    .attr("text-anchor", "middle")
    .attr("font-size", 12)
    .attr("fill", tc.internal)
    .text("k (number of clusters)");
  g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("x", -innerHeight / 2)
    .attr("y", -42)
    .attr("text-anchor", "middle")
    .attr("font-size", 12)
    .attr("fill", tc.internal)
    .text("CalBow Score");

  // Line + points
  const lineColor = tc.accent;
  const line = d3
    .line()
    .x((d) =>
      mode === "log" ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2,
    )
    .y((d) => yScale(d.score))
    .curve(d3.curveMonotoneX);
  g.append("path")
    .datum(dataPoints)
    .attr("fill", "none")
    .attr("stroke", lineColor)
    .attr("stroke-width", 2)
    .attr("d", line);

  g.selectAll(".score-point")
    .data(dataPoints)
    .enter()
    .append("circle")
    .attr("class", "score-point")
    .attr("cx", (d) =>
      mode === "log" ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2,
    )
    .attr("cy", (d) => yScale(d.score))
    .attr("r", 3.5)
    .attr("fill", lineColor)
    .attr("stroke", tc.bg)
    .attr("stroke-width", 1.5)
    .on("mouseover", (event, d) =>
      d3Tooltip
        .style("opacity", 1)
        .html(`k = ${d.k}<br/>score = ${d.score.toFixed(4)}`),
    )
    .on("mousemove", (event) =>
      d3Tooltip
        .style("left", event.pageX + 12 + "px")
        .style("top", event.pageY + 12 + "px"),
    )
    .on("mouseout", () => d3Tooltip.style("opacity", 0));

  // Peaks
  g.selectAll(".peak-point")
    .data(peakPoints)
    .enter()
    .append("path")
    .attr("class", "peak-point")
    .attr(
      "transform",
      (d) =>
        `translate(${mode === "log" ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2},${yScale(d.score)})`,
    )
    .attr("d", d3.symbol().type(d3.symbolDiamond).size(100))
    .attr("fill", "#ef4444")
    .attr("stroke", tc.bg)
    .attr("stroke-width", 1.5)
    .on("mouseover", (event, d) =>
      d3Tooltip
        .style("opacity", 1)
        .html(
          `<strong>Optimal k = ${d.k}</strong><br/>score = ${d.score.toFixed(4)}`,
        ),
    )
    .on("mousemove", (event) =>
      d3Tooltip
        .style("left", event.pageX + 12 + "px")
        .style("top", event.pageY + 12 + "px"),
    )
    .on("mouseout", () => d3Tooltip.style("opacity", 0));

  // Legend
  const legend = g
    .append("g")
    .attr("transform", `translate(${innerWidth - 130}, 8)`);
  legend
    .append("rect")
    .attr("x", -8)
    .attr("y", -8)
    .attr("width", 140)
    .attr("height", 44)
    .attr("fill", tc.bg)
    .attr("rx", 6)
    .attr("stroke", tc.border)
    .attr("opacity", 0.9);
  legend
    .append("line")
    .attr("x1", 0)
    .attr("y1", 6)
    .attr("x2", 18)
    .attr("y2", 6)
    .attr("stroke", lineColor)
    .attr("stroke-width", 2);
  legend
    .append("circle")
    .attr("cx", 9)
    .attr("cy", 6)
    .attr("r", 3)
    .attr("fill", lineColor);
  legend
    .append("text")
    .attr("x", 24)
    .attr("y", 10)
    .attr("font-size", 11)
    .attr("fill", tc.internal)
    .text("Score");
  legend
    .append("path")
    .attr("transform", "translate(9,26)")
    .attr("d", d3.symbol().type(d3.symbolDiamond).size(80))
    .attr("fill", "#ef4444");
  legend
    .append("text")
    .attr("x", 24)
    .attr("y", 30)
    .attr("font-size", 11)
    .attr("fill", tc.internal)
    .text("Optimal k");

  if (axisModeEl && !axisModeEl.__wired) {
    axisModeEl.__wired = true;
    axisModeEl.addEventListener("change", function () {
      axisModeEl.__userOverride = true;
      drawOptimalK(data);
    });
  }

  const okSvgNode = document.querySelector("#optimalk_plot svg");
  window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = okSvgNode
    ? new XMLSerializer().serializeToString(okSvgNode)
    : null;
}

/* ─────────────────────────────────────
   Mini Scores Panel (overlay on tree page)
   ───────────────────────────────────── */
function drawMiniScores(data) {
  var panel = document.getElementById("mini-scores-panel");
  var plotEl = document.getElementById("mini-scores-plot");
  if (!panel || !plotEl) return;

  var scores = data && data.scores ? data.scores : [];
  var peaks = data && data.peaks ? data.peaks : [];
  if (!Array.isArray(scores) || scores.length === 0) {
    panel.classList.remove("visible");
    return;
  }

  panel.classList.add("visible");
  plotEl.innerHTML = "";
  var tc = getThemeColors();

  var dataPoints = [];
  for (var i = 0; i < scores.length - 1; i++)
    dataPoints.push({ k: i + 2, score: scores[i + 1] });
  if (!dataPoints.length) return;

  var width = 252,
    height = 124;
  var margin = { top: 8, right: 8, bottom: 20, left: 36 };
  var innerW = width - margin.left - margin.right;
  var innerH = height - margin.top - margin.bottom;

  var svg = d3
    .select(plotEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  var g = svg
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var xScale = d3
    .scaleLinear()
    .domain([d3.min(dataPoints, (d) => d.k), d3.max(dataPoints, (d) => d.k)])
    .range([0, innerW]);
  var yScale = d3
    .scaleLinear()
    .domain([
      d3.min(dataPoints, (d) => d.score),
      d3.max(dataPoints, (d) => d.score),
    ])
    .nice()
    .range([innerH, 0]);

  // Axes (minimal)
  g.append("g")
    .attr("transform", "translate(0," + innerH + ")")
    .call(d3.axisBottom(xScale).ticks(4).tickFormat(d3.format("d")))
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", "9px");
  g.append("g")
    .call(d3.axisLeft(yScale).ticks(3))
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", "9px");
  g.selectAll(".domain").attr("stroke", tc.branch);
  g.selectAll(".tick line").attr("stroke", tc.branch);

  // Line
  var line = d3
    .line()
    .x((d) => xScale(d.k))
    .y((d) => yScale(d.score))
    .curve(d3.curveMonotoneX);
  g.append("path")
    .datum(dataPoints)
    .attr("fill", "none")
    .attr("stroke", tc.accent)
    .attr("stroke-width", 1.5)
    .attr("d", line);

  // Clickable points
  g.selectAll(".mini-point")
    .data(dataPoints)
    .enter()
    .append("circle")
    .attr("cx", (d) => xScale(d.k))
    .attr("cy", (d) => yScale(d.score))
    .attr("r", (d) => (peaks.indexOf(d.k) >= 0 ? 4 : 2.5))
    .attr("fill", (d) => (peaks.indexOf(d.k) >= 0 ? "#ef4444" : tc.accent))
    .attr("stroke", tc.bg)
    .attr("stroke-width", 1)
    .style("cursor", "pointer")
    .on("click", function (event, d) {
      switchToK(d.k);
    })
    .on("mouseover", function (event, d) {
      d3Tooltip
        .style("opacity", 1)
        .html("k=" + d.k + " score=" + d.score.toFixed(3));
    })
    .on("mousemove", function (event) {
      d3Tooltip
        .style("left", event.pageX + 10 + "px")
        .style("top", event.pageY + 10 + "px");
    })
    .on("mouseout", function () {
      d3Tooltip.style("opacity", 0);
    });
}

function switchToK(k) {
  if (!latestApiData) return;
  // Try all_clusters first, then peaks
  if (latestApiData.all_clusters && latestApiData.all_ks) {
    var idx = latestApiData.all_ks.indexOf(k);
    if (idx >= 0) {
      // Ensure we're in "all" view mode
      if (CLUSTER_VIEW_MODE !== "all") {
        CLUSTER_VIEW_MODE = "all";
        populateClusterSelector(latestApiData);
        var toggleEl = document.getElementById("cluster-view-toggle");
        if (toggleEl) toggleEl.textContent = "Show peaks only";
      }
      var selectEl = document.getElementById("cluster-select");
      if (selectEl) selectEl.selectedIndex = idx;
      switchCluster(idx);
      return;
    }
  }
  var peakKs = latestApiData.k_values || latestApiData.ks;
  if (peakKs && latestApiData.clusters) {
    var idx2 = peakKs.indexOf(k);
    if (idx2 >= 0) {
      if (CLUSTER_VIEW_MODE !== "peaks") {
        CLUSTER_VIEW_MODE = "peaks";
        populateClusterSelector(latestApiData);
        var toggleEl2 = document.getElementById("cluster-view-toggle");
        if (toggleEl2) toggleEl2.textContent = "Show all k";
      }
      var selectEl2 = document.getElementById("cluster-select");
      if (selectEl2) selectEl2.selectedIndex = idx2;
      switchCluster(idx2);
      return;
    }
  }
  showToast(
    "k=" + k + " not cached. Enable 'Compute all clusters' to access any k.",
    "info",
    3000,
  );
}

/* ─────────────────────────────────────
   File Handling
   ───────────────────────────────────── */
function handleFileSelect(evt) {
  var file = evt.target.files ? evt.target.files[0] : null;
  if (!file) return;
  loadFile(file);
}

function loadFile(file) {
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

/* ─────────────────────────────────────
   Save / Export
   ───────────────────────────────────── */
function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 100);
}

function getTightSvgSnapshot(svgNode, paddingPx) {
  if (!svgNode) return null;
  const padding = Math.max(0, Number(paddingPx || 10));

  // The first <g> is the zoom/pan layer. We measure *its* bbox (not its
  // first child's) because ``SVGGraphicsElement.getBBox()`` returns a bbox
  // in the element's **local** coordinate system and accounts for descendant
  // transforms. Measuring the zoom layer therefore captures every drawn
  // child — tree, leaf labels, cluster side-bars, box overlays, legends —
  // in the root svg's coordinate frame, which is the same frame the final
  // ``viewBox`` uses. The previous implementation measured the inner
  // content <g> (which carries a ``translate(margin.left, margin.top)``
  // transform), so the bbox was in that <g>'s shifted local space; when
  // it was applied as ``viewBox`` on the root svg, content drawn by g's
  // still-translated children landed outside the viewBox on the right and
  // bottom — cropping the color bars and the right edge of the tree.
  let zoomLayer = null;
  const firstChild = svgNode.firstElementChild;
  if (firstChild && firstChild.tagName.toLowerCase() === "g") {
    zoomLayer = firstChild;
  }

  // If we have a zoom layer, measure it directly. Its own transform
  // (current pan/zoom) is NOT included in the returned bbox because
  // getBBox reports in the element's local coord system — which is
  // exactly what we want: a viewport-invariant measurement.
  let target = zoomLayer || svgNode.querySelector("g") || svgNode;

  let bbox = null;
  try {
    bbox = target.getBBox();
  } catch {
    bbox = null;
  }

  const fallbackW =
    parseFloat(svgNode.getAttribute("width")) || svgNode.clientWidth || 800;
  const fallbackH =
    parseFloat(svgNode.getAttribute("height")) || svgNode.clientHeight || 600;

  const x = bbox && isFinite(bbox.x) ? bbox.x - padding : 0;
  const y = bbox && isFinite(bbox.y) ? bbox.y - padding : 0;
  const w =
    bbox && isFinite(bbox.width) && bbox.width > 0
      ? bbox.width + padding * 2
      : fallbackW;
  const h =
    bbox && isFinite(bbox.height) && bbox.height > 0
      ? bbox.height + padding * 2
      : fallbackH;

  const clone = svgNode.cloneNode(true);
  if (!clone.getAttribute("xmlns")) {
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  }

  // Normalize zoom state in the exported clone so bounds and drawing match.
  const cloneFirstChild = clone.firstElementChild;
  if (cloneFirstChild && cloneFirstChild.tagName.toLowerCase() === "g") {
    cloneFirstChild.removeAttribute("transform");
  }

  clone.setAttribute("viewBox", `${x} ${y} ${w} ${h}`);
  clone.setAttribute("width", String(w));
  clone.setAttribute("height", String(h));

  return {
    svgData: new XMLSerializer().serializeToString(clone),
    width: w,
    height: h,
  };
}

function exportSvgFromEl(selector, filename) {
  const svgNode = document.querySelector(selector + " svg");
  if (!svgNode) {
    showToast("No visualization to export.", "danger");
    return;
  }
  const snap = getTightSvgSnapshot(svgNode, 10);
  if (!snap) {
    showToast("SVG export failed.", "danger");
    return;
  }
  downloadBlob(
    new Blob([snap.svgData], {
      type: "image/svg+xml",
    }),
    filename,
  );
  showToast("Downloaded " + filename, "success", 2000);
}

async function exportTSV() {
  try {
    const res = await fetch("/api/export_tsv", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ top_n: 1, outlier: true, run_id: latestRunId }),
    });
    if (!res.ok) {
      const err = await res.json().catch(() => ({}));
      throw new Error(err.detail || "Export failed");
    }
    downloadBlob(
      new Blob([await res.text()], { type: "text/tab-separated-values" }),
      "phytclust_results.tsv",
    );
    showToast("Downloaded TSV", "success", 2000);
  } catch (e) {
    showToast("TSV export failed: " + e.message, "danger");
  }
}

async function saveToServer() {
  var el = document.getElementById("output-dir");
  const dir = el && el.value ? el.value : "results";
  try {
    await apiPostJson("/api/save", {
      results_dir: dir,
      top_n: 1,
      outlier: true,
    });
    showToast("Saved to " + dir, "success", 3000);
  } catch (e) {
    showToast("Save failed: " + e.message, "danger");
  }
}

function exportPngFromEl(selector, filename, dpi) {
  const svgNode = document.querySelector(selector + " svg");
  if (!svgNode) {
    showToast("No visualization to export.", "danger");
    return;
  }
  const snap = getTightSvgSnapshot(svgNode, 10);
  if (!snap) {
    showToast("PNG export failed.", "danger");
    return;
  }
  dpi = dpi || 300;
  const scale = dpi / 96; // browsers render at 96 DPI
  const svgData = snap.svgData;
  const svgW = snap.width;
  const svgH = snap.height;
  const canvas = document.createElement("canvas");
  canvas.width = Math.ceil(svgW * scale);
  canvas.height = Math.ceil(svgH * scale);
  const ctx = canvas.getContext("2d");
  // Fill with white background for publication
  const bgColor =
    getComputedStyle(document.documentElement)
      .getPropertyValue("--pc-bg")
      .trim() || "#ffffff";
  ctx.fillStyle = bgColor;
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(scale, scale);
  const img = new Image();
  const blob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  img.onload = function () {
    ctx.drawImage(img, 0, 0, svgW, svgH);
    URL.revokeObjectURL(url);
    canvas.toBlob(function (pngBlob) {
      if (!pngBlob) {
        showToast("PNG export failed.", "danger");
        return;
      }
      downloadBlob(pngBlob, filename);
      showToast(
        "Downloaded " + filename + " at " + dpi + " DPI",
        "success",
        3000,
      );
    }, "image/png");
  };
  img.onerror = function () {
    URL.revokeObjectURL(url);
    showToast("PNG rendering failed.", "danger");
  };
  img.src = url;
}

function copyRasterToClipboard(selector, dpi, mimeType) {
  var svgNode = document.querySelector(selector + " svg");
  if (!svgNode) {
    showToast("Nothing to copy.", "danger");
    return;
  }
  var snap = getTightSvgSnapshot(svgNode, 10);
  if (!snap) {
    showToast("Copy failed.", "danger");
    return;
  }
  dpi = dpi || 300;
  mimeType = mimeType || "image/png";
  var fmtLabel = mimeType === "image/jpeg" ? "JPG" : "PNG";
  var scale = dpi / 96;
  var svgData = snap.svgData;
  var svgW = snap.width;
  var svgH = snap.height;
  var canvas = document.createElement("canvas");
  canvas.width = Math.ceil(svgW * scale);
  canvas.height = Math.ceil(svgH * scale);
  var ctx = canvas.getContext("2d");
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(scale, scale);
  var img = new Image();
  var blob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" });
  var url = URL.createObjectURL(blob);
  img.onload = function () {
    ctx.drawImage(img, 0, 0, svgW, svgH);
    URL.revokeObjectURL(url);
    // Clipboard API only supports image/png, so for JPG we download instead
    if (mimeType === "image/jpeg") {
      canvas.toBlob(
        function (jpgBlob) {
          if (!jpgBlob) {
            showToast("JPG render failed.", "danger");
            return;
          }
          downloadBlob(jpgBlob, "phytclust_tree.jpg");
          showToast("Downloaded JPG at " + dpi + " DPI", "success", 2000);
        },
        "image/jpeg",
        0.95,
      );
      return;
    }
    canvas.toBlob(function (pngBlob) {
      if (!pngBlob) {
        showToast(fmtLabel + " copy failed.", "danger");
        return;
      }
      try {
        var item = new ClipboardItem({ "image/png": pngBlob });
        navigator.clipboard.write([item]).then(
          function () {
            showToast(
              fmtLabel + " copied to clipboard (" + dpi + " DPI)",
              "success",
              2000,
            );
          },
          function () {
            showToast(
              "Clipboard write failed — try the Export menu instead.",
              "danger",
            );
          },
        );
      } catch (e) {
        showToast(
          "Clipboard API not supported — try the Export menu instead.",
          "danger",
        );
      }
    }, "image/png");
  };
  img.onerror = function () {
    URL.revokeObjectURL(url);
    showToast(fmtLabel + " rendering failed.", "danger");
  };
  img.src = url;
}

function copySvgToClipboard(selector) {
  const svgNode = document.querySelector(selector + " svg");
  if (!svgNode) {
    showToast("Nothing to copy.", "danger");
    return;
  }
  const snap = getTightSvgSnapshot(svgNode, 10);
  if (!snap) {
    showToast("Copy failed", "danger");
    return;
  }
  navigator.clipboard.writeText(snap.svgData).then(
    () => showToast("SVG copied to clipboard", "success", 2000),
    () => showToast("Copy failed", "danger"),
  );
}

/* ─────────────────────────────────────
   Tab System
   ───────────────────────────────────── */
function switchTab(tabName) {
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

  if (tabName === "viewer-optimal-k" && latestOptimalKData)
    drawOptimalK(latestOptimalKData);
  if (tabName === "compare") drawComparison();
}

/* ─────────────────────────────────────
   Collapsible section helper
   ───────────────────────────────────── */
function wireCollapse(toggleId, bodyId) {
  var toggle = document.getElementById(toggleId);
  var body = document.getElementById(bodyId);
  if (!toggle || !body) return;
  toggle.addEventListener("click", function () {
    var isOpen = body.classList.contains("open");
    body.classList.toggle("open", !isOpen);
    toggle.setAttribute("aria-expanded", !isOpen);
    var chev = toggle.querySelector(".chevron");
    if (chev) chev.innerHTML = isOpen ? "&#x25B6;" : "&#x25BC;";
  });
}

/* ─────────────────────────────────────
   DOM Ready – Wire Everything
   ───────────────────────────────────── */
document.addEventListener("DOMContentLoaded", function () {
  var appLoading = document.getElementById("app-loading");
  if (appLoading) appLoading.style.display = "none";
  showStatus("Ready", "success");

  // ── Branch options collapse ──
  var branchOptsBtn = document.getElementById("branch-options-toggle");
  var branchOptsPanel = document.getElementById("branch-options-panel");
  var branchOptsWrap = document.querySelector(".branch-options-wrap");
  if (branchOptsBtn && branchOptsPanel) {
    var closeBranchOptionsPanel = function () {
      branchOptsPanel.classList.remove("open");
      branchOptsBtn.setAttribute("aria-expanded", "false");
    };

    var openBranchOptionsPanel = function () {
      branchOptsPanel.classList.add("open");
      branchOptsBtn.setAttribute("aria-expanded", "true");
    };

    branchOptsBtn.addEventListener("click", function (e) {
      e.preventDefault();
      e.stopPropagation();
      var isOpen = branchOptsPanel.classList.contains("open");
      if (isOpen) closeBranchOptionsPanel();
      else openBranchOptionsPanel();
    });

    document.addEventListener("click", function (e) {
      if (!branchOptsPanel.classList.contains("open")) return;
      if (!branchOptsWrap) return;
      if (branchOptsWrap.contains(e.target)) return;
      closeBranchOptionsPanel();
    });

    document.addEventListener("keydown", function (e) {
      if (e.key !== "Escape") return;
      if (!branchOptsPanel.classList.contains("open")) return;
      closeBranchOptionsPanel();
      branchOptsBtn.focus();
    });
  }

  // ── Sidebar toggle ──
  var sidebarToggle = document.getElementById("sidebar-toggle");
  var sidebar = document.getElementById("app-sidebar");
  if (sidebarToggle && sidebar) {
    sidebarToggle.addEventListener("click", function () {
      sidebar.classList.toggle("collapsed");
    });
  }

  // ── File input ──
  var fileInput = document.getElementById("file-input");
  if (fileInput) fileInput.addEventListener("change", handleFileSelect);

  // ── Drag and drop ──
  var uploadZone = document.getElementById("upload-zone");
  if (uploadZone) {
    uploadZone.addEventListener("dragover", function (e) {
      e.preventDefault();
      e.stopPropagation();
      uploadZone.classList.add("drag-over");
    });
    uploadZone.addEventListener("dragleave", function (e) {
      e.preventDefault();
      e.stopPropagation();
      uploadZone.classList.remove("drag-over");
    });
    uploadZone.addEventListener("drop", function (e) {
      e.preventDefault();
      e.stopPropagation();
      uploadZone.classList.remove("drag-over");
      if (e.dataTransfer.files && e.dataTransfer.files.length)
        loadFile(e.dataTransfer.files[0]);
    });
  }

  // ── Run button + keyboard shortcut ──
  var btnRun = document.getElementById("btn-run");
  if (btnRun)
    btnRun.addEventListener("click", function (e) {
      e.preventDefault();
      runPhytClust();
    });
  document.addEventListener("keydown", function (e) {
    if ((e.ctrlKey || e.metaKey) && e.key === "Enter") {
      e.preventDefault();
      runPhytClust();
    }
  });

  // ── Reset button ──
  var resetBtn = document.getElementById("btn-extra-reset");
  if (resetBtn)
    resetBtn.addEventListener("click", function (e) {
      e.preventDefault();
      resetExtraParams();
    });

  // ── Mode selector ──
  document.querySelectorAll("#mode-selector .mode-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document
        .querySelectorAll("#mode-selector .mode-btn")
        .forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      var mode = btn.dataset.mode;
      var paramK = document.getElementById("param-k");
      var paramTopN = document.getElementById("param-topn");
      var paramBins = document.getElementById("param-bins");
      if (paramK) paramK.style.display = mode === "k" ? "" : "none";
      if (paramTopN) paramTopN.style.display = mode === "k" ? "none" : "";
      if (paramBins)
        paramBins.style.display = mode === "resolution" ? "" : "none";
      if (extraResolutionEl) extraResolutionEl.checked = mode === "resolution";
    });
  });

  // ── Collapsible sections ──
  wireCollapse("tree-opts-toggle", "tree-opts-body");
  wireCollapse("outlier-opts-toggle", "outlier-opts-body");
  wireCollapse("support-opts-toggle", "support-opts-body");
  wireCollapse("peak-opts-toggle", "peak-opts-body");

  // ── Tabs ──
  document.querySelectorAll(".tab-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      switchTab(btn.dataset.tab);
    });
  });
  // Ensure toolbars/panels match the active default tab after any hot reload.
  switchTab("viewer");

  // ── Help sidebar ──
  var btnHelp = document.getElementById("btn-help");
  var helpSidebar = document.getElementById("help-sidebar");
  if (btnHelp && helpSidebar)
    btnHelp.addEventListener("click", function (e) {
      e.stopPropagation();
      helpSidebar.classList.toggle("open");
    });
  var btnHelpClose = document.getElementById("btn-help-close");
  if (btnHelpClose && helpSidebar)
    btnHelpClose.addEventListener("click", function () {
      helpSidebar.classList.remove("open");
    });
  // Close help when clicking outside
  if (helpSidebar) {
    document.addEventListener("click", function (e) {
      if (
        helpSidebar.classList.contains("open") &&
        !helpSidebar.contains(e.target) &&
        e.target !== btnHelp
      ) {
        helpSidebar.classList.remove("open");
      }
    });
    helpSidebar.addEventListener("click", function (e) {
      e.stopPropagation();
    });
  }

  // ── About modal ──
  var btnAbout = document.getElementById("btn-about");
  var aboutModal = document.getElementById("aboutModal");
  var aboutBg = document.getElementById("aboutBackdrop");
  function openAbout() {
    if (aboutModal) aboutModal.classList.add("show");
    if (aboutBg) aboutBg.classList.add("show");
  }
  function closeAbout() {
    if (aboutModal) aboutModal.classList.remove("show");
    if (aboutBg) aboutBg.classList.remove("show");
  }
  if (btnAbout) btnAbout.addEventListener("click", openAbout);
  var btnAboutClose1 = document.getElementById("btn-about-close");
  var btnAboutClose2 = document.getElementById("btn-about-close2");
  if (btnAboutClose1) btnAboutClose1.addEventListener("click", closeAbout);
  if (btnAboutClose2) btnAboutClose2.addEventListener("click", closeAbout);
  if (aboutBg) aboutBg.addEventListener("click", closeAbout);

  // ── Help load example ──
  var helpExampleBtn = document.getElementById("btn-help-load-example");
  if (helpExampleBtn) {
    helpExampleBtn.addEventListener("click", function (e) {
      e.preventDefault();
      newickEl.value = EXAMPLE_NEWICK;
      if (helpSidebar) helpSidebar.classList.remove("open");
      runPhytClust();
    });
  }

  // ── Save dropdown ──
  var btnSave = document.getElementById("btn-save");
  var saveDropdown = document.getElementById("save-dropdown");
  if (btnSave && saveDropdown) {
    btnSave.addEventListener("click", function (e) {
      e.stopPropagation();
      saveDropdown.classList.toggle("show");
    });
    document.addEventListener("click", function () {
      saveDropdown.classList.remove("show");
    });
    saveDropdown.addEventListener("click", function (e) {
      e.stopPropagation();
    });
    saveDropdown.querySelectorAll(".dropdown-item").forEach(function (item) {
      item.addEventListener("click", function () {
        saveDropdown.classList.remove("show");
        var type = item.dataset.saveType;
        if (type === "tsv") exportTSV();
        if (type === "tree_plot")
          exportSvgFromEl("#tree_display", "phytclust_tree.svg");
        if (type === "tree_png") {
          var dpiStr = prompt(
            "Enter DPI for publication-ready PNG (default 300):",
            "300",
          );
          if (dpiStr === null) return;
          var dpi = parseInt(dpiStr, 10);
          if (isNaN(dpi) || dpi < 72) dpi = 300;
          exportPngFromEl("#tree_display", "phytclust_tree.png", dpi);
        }
        if (type === "k_plot")
          exportSvgFromEl("#optimalk_plot", "phytclust_scores.svg");
        if (type === "save_server") saveToServer();
      });
    });
  }

  // ── Copy dropdown (tree) ──
  var btnCopyTree = document.getElementById("btn-copy-tree");
  var copyDropdown = document.getElementById("copy-dropdown");
  if (btnCopyTree && copyDropdown) {
    btnCopyTree.addEventListener("click", function (e) {
      e.stopPropagation();
      copyDropdown.classList.toggle("show");
    });
    document.addEventListener("click", function () {
      copyDropdown.classList.remove("show");
    });
    copyDropdown.addEventListener("click", function (e) {
      e.stopPropagation();
    });
    copyDropdown.querySelectorAll(".dropdown-item").forEach(function (item) {
      item.addEventListener("click", function () {
        copyDropdown.classList.remove("show");
        var fmt = item.dataset.copyFmt;
        var dpiEl = document.getElementById("copy-dpi");
        var dpi = dpiEl ? parseInt(dpiEl.value, 10) : 300;
        if (isNaN(dpi) || dpi < 72) dpi = 300;
        if (fmt === "svg") {
          copySvgToClipboard("#tree_display");
        } else if (fmt === "png") {
          copyRasterToClipboard("#tree_display", dpi, "image/png");
        } else if (fmt === "jpg") {
          copyRasterToClipboard("#tree_display", dpi, "image/jpeg");
        }
      });
    });
  }
  // ── Copy button for scores plot (keep simple) ──
  var btnCopyMaxk = document.getElementById("btn-copy-maxk");
  if (btnCopyMaxk)
    btnCopyMaxk.addEventListener("click", function () {
      copySvgToClipboard("#optimalk_plot");
    });

  // ── Cluster selector ──
  var clusterSelect = document.getElementById("cluster-select");
  if (clusterSelect)
    clusterSelect.addEventListener("change", function () {
      switchCluster(parseInt(this.value, 10));
    });
  var clusterToggle = document.getElementById("cluster-view-toggle");
  if (clusterToggle)
    clusterToggle.addEventListener("click", toggleClusterViewMode);

  // ── Cluster bulk editor ──
  var editClustersBtn = document.getElementById("btn-edit-clusters");
  var clusterEditorBg = document.getElementById("clusterEditorBackdrop");
  var clusterEditorClose = document.getElementById("btn-cluster-editor-close");
  var clusterEditorCancel = document.getElementById(
    "btn-cluster-editor-cancel",
  );
  var clusterEditorSave = document.getElementById("btn-cluster-editor-save");
  if (editClustersBtn)
    editClustersBtn.addEventListener("click", function () {
      if (!getCurrentClusterIds().length) {
        showToast("No clusters available yet.", "info", 1600);
        return;
      }
      openClusterEditor();
    });
  if (clusterEditorClose)
    clusterEditorClose.addEventListener("click", closeClusterEditor);
  if (clusterEditorCancel)
    clusterEditorCancel.addEventListener("click", closeClusterEditor);
  if (clusterEditorSave)
    clusterEditorSave.addEventListener("click", saveClusterEditor);
  if (clusterEditorBg)
    clusterEditorBg.addEventListener("click", closeClusterEditor);
  var clusterFilterInput = document.getElementById("cluster-filter-text");
  var clusterFilterOutliersCb = document.getElementById(
    "cluster-filter-show-outliers",
  );
  if (clusterFilterInput)
    clusterFilterInput.addEventListener("input", renderClusterEditorRows);
  if (clusterFilterOutliersCb)
    clusterFilterOutliersCb.addEventListener("change", renderClusterEditorRows);

  // ── Cluster prev/next buttons ──
  var clusterPrev = document.getElementById("cluster-prev");
  var clusterNext = document.getElementById("cluster-next");
  if (clusterPrev)
    clusterPrev.addEventListener("click", function () {
      cycleCluster(-1);
    });
  if (clusterNext)
    clusterNext.addEventListener("click", function () {
      cycleCluster(1);
    });

  // ── Arrow key cycling through clusters ──
  document.addEventListener("keydown", function (e) {
    if (
      e.target.tagName === "INPUT" ||
      e.target.tagName === "TEXTAREA" ||
      e.target.tagName === "SELECT"
    )
      return;

    var isArrow =
      e.key === "ArrowLeft" ||
      e.key === "ArrowRight" ||
      e.key === "ArrowUp" ||
      e.key === "ArrowDown";

    if (isArrow && COLOR_MODE === "boxes" && SELECTED_CLUSTER_IDS.size > 0) {
      e.preventDefault();
      var step = 4;
      var resize = !!e.shiftKey;
      SELECTED_CLUSTER_IDS.forEach(function (cid) {
        var adj = getBoxAdjust(cid);
        if (resize) {
          if (e.key === "ArrowLeft") adj.padX = Math.max(0, adj.padX - step);
          if (e.key === "ArrowRight") adj.padX += step;
          if (e.key === "ArrowUp") adj.padY += step;
          if (e.key === "ArrowDown") adj.padY = Math.max(0, adj.padY - step);
        } else {
          if (e.key === "ArrowLeft") adj.dx -= step;
          if (e.key === "ArrowRight") adj.dx += step;
          if (e.key === "ArrowUp") adj.dy -= step;
          if (e.key === "ArrowDown") adj.dy += step;
        }
      });
      drawTree();
      return;
    }

    if (e.key === "Escape" && SELECTED_CLUSTER_IDS.size > 0) {
      SELECTED_CLUSTER_IDS.clear();
      drawTree();
      showToast("Cluster focus cleared", "info", 1500);
      return;
    }

    if (e.key === "ArrowLeft") {
      e.preventDefault();
      cycleCluster(-1);
    }
    if (e.key === "ArrowRight") {
      e.preventDefault();
      cycleCluster(1);
    }
  });

  // ── Compare add bar ──
  var compareAddBtn = document.getElementById("compare-add-bar");
  if (compareAddBtn) {
    compareAddBtn.addEventListener("click", function () {
      COMPARE_CONFIGS.push(makeCompareConfig());
      renderCompareConfigList();
      drawComparison();
    });
  }

  // ── Layout mode ──
  var layoutSelect = document.getElementById("layout-mode");
  if (layoutSelect)
    layoutSelect.addEventListener("change", function () {
      CURRENT_LAYOUT_MODE = this.value;
      drawTree();
    });

  // ── Color mode ──
  var colorModeSelect = document.getElementById("color-mode");
  if (colorModeSelect)
    colorModeSelect.addEventListener("change", function () {
      COLOR_MODE = this.value || "bars";
      drawTree();
    });

  var branchWidthInput = document.getElementById("branch-width");
  if (branchWidthInput) {
    branchWidthInput.value = BRANCH_STROKE_WIDTH;
    branchWidthInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 0.5 && v <= 6) {
        BRANCH_STROKE_WIDTH = v;
        drawTree();
        if (document.getElementById("compare").classList.contains("active"))
          drawComparison();
      }
    });
  }

  var branchColorInput = document.getElementById("branch-color");
  var branchColorResetBtn = document.getElementById("branch-color-reset");
  if (branchColorInput) {
    const tc = getThemeColors();
    branchColorInput.value = tc.branch || "#000000";
    branchColorInput.addEventListener("input", function () {
      BRANCH_COLOR_OVERRIDE = this.value || null;
      drawTree();
      if (document.getElementById("compare").classList.contains("active"))
        drawComparison();
    });
  }
  if (branchColorResetBtn) {
    branchColorResetBtn.addEventListener("click", function () {
      BRANCH_COLOR_OVERRIDE = null;
      const style = getComputedStyle(document.documentElement);
      const themeBranch =
        style.getPropertyValue("--pc-tree-branch").trim() || "#000000";
      if (branchColorInput) branchColorInput.value = themeBranch;
      drawTree();
      if (document.getElementById("compare").classList.contains("active"))
        drawComparison();
    });
  }

  // ── Palette selector ──
  let CUSTOM_PALETTE = null;
  function normalizeHexColor(s) {
    const v = (s || "").trim();
    if (/^#[0-9a-fA-F]{6}$/.test(v)) return v;
    if (/^[0-9a-fA-F]{6}$/.test(v)) return "#" + v;
    return null;
  }
  function updatePalettePreview(colors) {
    const box = document.getElementById("palette-preview");
    if (box)
      box.innerHTML = colors
        .map((c) => `<span title="${c}" style="background:${c};"></span>`)
        .join("");
  }

  function applyPalette(type) {
    const palettes = {
      default: [
        "#b84b4b",
        "#849060",
        "#3d7c74",
        "#6e3f8a",
        "#ceb94b",
        "#3f648a",
        "#3f408a",
        "#da63aa",
      ],
      pastel: [
        "#ffb3ba",
        "#ffdfba",
        "#ffffba",
        "#baffc9",
        "#bae1ff",
        "#d7baff",
        "#ffcce6",
        "#c2f0c2",
      ],
      vivid: [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#ffff33",
        "#a65628",
        "#f781bf",
      ],
      dark: [
        "#4e79a7",
        "#59a14f",
        "#e15759",
        "#b07aa1",
        "#edc948",
        "#76b7b2",
        "#ff9da7",
        "#9c755f",
      ],
    };
    if (type === "custom" && CUSTOM_PALETTE && CUSTOM_PALETTE.length)
      BASE_COLORS.splice(0, BASE_COLORS.length, ...CUSTOM_PALETTE);
    else if (palettes[type])
      BASE_COLORS.splice(0, BASE_COLORS.length, ...palettes[type]);
    updatePalettePreview(BASE_COLORS);
    if (Object.keys(CURRENT_CLUSTERS || {}).length > 0) {
      CLUSTER_COLORS = generateClusterColors(
        Math.max(...Object.values(CURRENT_CLUSTERS)) + 1,
      );
      drawTree();
    }
  }

  var paletteSelect = document.getElementById("palette-select");
  if (paletteSelect) {
    paletteSelect.addEventListener("change", function () {
      if (this.value === "custom") {
        var initial =
          CUSTOM_PALETTE && CUSTOM_PALETTE.length
            ? CUSTOM_PALETTE.join(",")
            : BASE_COLORS.join(",");
        var raw = prompt("Enter hex colors separated by commas:", initial);
        if (raw == null) {
          this.value = "default";
          applyPalette("default");
          return;
        }
        var parsed = raw.split(",").map(normalizeHexColor).filter(Boolean);
        if (parsed.length < 2) {
          alert("Please provide at least 2 valid hex colors.");
          this.value = "default";
          applyPalette("default");
          return;
        }
        CUSTOM_PALETTE = parsed;
        applyPalette("custom");
        return;
      }
      applyPalette(this.value);
    });
  }
  updatePalettePreview(BASE_COLORS);

  // ── Shuffle colors button ──
  var btnShuffle = document.getElementById("btn-shuffle-colors");
  if (btnShuffle) {
    btnShuffle.addEventListener("click", function () {
      var shuffled = shuffle(BASE_COLORS.slice());
      BASE_COLORS.splice(0, BASE_COLORS.length);
      for (var i = 0; i < shuffled.length; i++) BASE_COLORS.push(shuffled[i]);
      updatePalettePreview(BASE_COLORS);
      if (Object.keys(CURRENT_CLUSTERS || {}).length > 0) {
        CLUSTER_COLORS = generateClusterColors(
          Math.max.apply(null, Object.values(CURRENT_CLUSTERS).map(Number)) + 1,
        );
        drawTree();
      }
    });
  }

  // ── Leaf/internal name toggles ──
  var internalCb = document.getElementById("show-internal-names");
  if (internalCb) {
    internalCb.checked = SHOW_INTERNAL_NAMES;
    internalCb.addEventListener("change", function () {
      SHOW_INTERNAL_NAMES = this.checked;
      drawTree();
    });
  }
  var leafCb = document.getElementById("show-leaf-names");
  if (leafCb) {
    leafCb.checked = SHOW_LEAF_NAMES;
    leafCb.addEventListener("change", function () {
      SHOW_LEAF_NAMES = this.checked;
      drawTree();
    });
  }

  var boxLabelCb = document.getElementById("show-box-labels");
  if (boxLabelCb) {
    boxLabelCb.checked = SHOW_BOX_LABELS;
    boxLabelCb.addEventListener("change", function () {
      SHOW_BOX_LABELS = this.checked;
      drawTree();
    });
  }
  var outlierBoxesCb = document.getElementById("show-outlier-boxes");
  if (outlierBoxesCb) {
    outlierBoxesCb.checked = SHOW_OUTLIER_BOXES;
    outlierBoxesCb.addEventListener("change", function () {
      SHOW_OUTLIER_BOXES = this.checked;
      drawTree();
    });
  }
  var boxAlphaSlider = document.getElementById("box-alpha");
  if (boxAlphaSlider) {
    boxAlphaSlider.value = BOX_ALPHA;
    boxAlphaSlider.addEventListener("input", function () {
      BOX_ALPHA = parseFloat(this.value);
      drawTree();
    });
  }
  var boxPadVSlider = document.getElementById("box-pad-v");
  if (boxPadVSlider) {
    boxPadVSlider.value = BOX_PAD_V;
    boxPadVSlider.addEventListener("input", function () {
      BOX_PAD_V = parseFloat(this.value);
      drawTree();
    });
  }
  var boxPadHSlider = document.getElementById("box-pad-h");
  if (boxPadHSlider) {
    boxPadHSlider.value = BOX_PAD_H;
    boxPadHSlider.addEventListener("input", function () {
      BOX_PAD_H = parseFloat(this.value);
      drawTree();
    });
  }
  var boxRadiusSlider = document.getElementById("box-radius");
  if (boxRadiusSlider) {
    boxRadiusSlider.value = BOX_CORNER_RADIUS;
    boxRadiusSlider.addEventListener("input", function () {
      BOX_CORNER_RADIUS = parseFloat(this.value);
      drawTree();
    });
  }

  // ── Node/label sizes ──
  var nodeSizeInput = document.getElementById("node-size");
  if (nodeSizeInput) {
    nodeSizeInput.value = LEAF_NODE_RADIUS;
    nodeSizeInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 1 && v <= 10) {
        LEAF_NODE_RADIUS = v;
        drawTree();
      }
    });
  }
  var labelSizeInput = document.getElementById("label-size");
  if (labelSizeInput) {
    labelSizeInput.value = LABEL_FONT_SIZE;
    labelSizeInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 6 && v <= 36) {
        LABEL_FONT_SIZE = v;
        drawTree();
      }
    });
  }
  var axisSizeInput = document.getElementById("axis-font-size");
  if (axisSizeInput) {
    axisSizeInput.value = AXIS_FONT_SIZE;
    axisSizeInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 6 && v <= 18) {
        AXIS_FONT_SIZE = v;
        drawTree();
      }
    });
  }

  // ── Width/Height sliders ──
  var widthSlider = document.getElementById("tree-width-scale");
  var heightSlider = document.getElementById("tree-height-scale");
  if (widthSlider)
    widthSlider.addEventListener("input", function () {
      TREE_WIDTH_SCALE = parseFloat(this.value);
      if (!isNaN(TREE_WIDTH_SCALE)) drawTree();
    });
  if (heightSlider)
    heightSlider.addEventListener("input", function () {
      TREE_HEIGHT_SCALE = parseFloat(this.value);
      if (!isNaN(TREE_HEIGHT_SCALE)) drawTree();
    });

  // ── Reset view ──
  var btnResetView = document.getElementById("btn-reset-view");
  if (btnResetView) {
    btnResetView.addEventListener("click", function (e) {
      e.preventDefault();
      if (LAST_TREE_SVG && LAST_TREE_ZOOM)
        LAST_TREE_SVG.transition()
          .duration(150)
          .call(LAST_TREE_ZOOM.transform, d3.zoomIdentity);
      if (widthSlider) widthSlider.value = "1.0";
      if (heightSlider) heightSlider.value = "1.0";
      TREE_WIDTH_SCALE = 1.0;
      TREE_HEIGHT_SCALE = 1.0;
      clearAllCollapsedFlags(NEWICK_RAW_TREE);
      drawTree();
    });
  }

  // ── Search ──
  var searchInput = document.getElementById("node-search");
  if (searchInput) {
    let searchTimer = null;
    searchInput.addEventListener("input", function () {
      clearTimeout(searchTimer);
      searchTimer = setTimeout(function () {
        SEARCH_TERM = searchInput.value.trim();
        if (NEWICK_RAW_TREE) drawTree();
      }, 200);
    });
  }

  // ── Context menu buttons ──
  (function initContextMenuButtons() {
    function btn(id, fn) {
      var el = document.getElementById(id);
      if (el)
        el.addEventListener("click", function () {
          fn();
          hideNodeContextMenu();
        });
    }
    btn("ctx-rename", function () {
      if (!CTX_TARGET_DATA) return;
      var c = getNodeCustom(CTX_TARGET_DATA);
      var val = prompt("Rename node:", nodeDisplayName(CTX_TARGET_DATA));
      if (val != null) {
        c.renamedTo = val;
        drawTree();
      }
    });
    btn("ctx-highlight", function () {
      if (!CTX_TARGET_DATA) return;
      var c = getNodeCustom(CTX_TARGET_DATA);
      c.highlighted = !c.highlighted;
      drawTree();
    });
    // Size slider in context menu
    var ctxSlider = document.getElementById("ctx-size-slider");
    var ctxSliderVal = document.getElementById("ctx-size-value");
    if (ctxSlider) {
      ctxSlider.addEventListener("input", function () {
        if (!CTX_TARGET_DATA) return;
        var c = getNodeCustom(CTX_TARGET_DATA);
        c.radiusScale = parseFloat(ctxSlider.value);
        if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);
        drawTree();
      });
      // Prevent menu from closing when interacting with slider
      ctxSlider.addEventListener("click", function (e) {
        e.stopPropagation();
      });
    }
    btn("ctx-copy-subtree", function () {
      if (!CTX_TARGET_DATA) return;
      var names = [];
      collectLeafNamesData(CTX_TARGET_DATA, names);
      navigator.clipboard.writeText(names.join("\n")).then(
        () => showToast(names.length + " leaf names copied", "success", 2000),
        () => showToast("Copy failed", "danger"),
      );
    });
    btn("ctx-original-size", function () {
      if (!CTX_TARGET_DATA) return;
      var c = getNodeCustom(CTX_TARGET_DATA);
      c.lockOriginalSize = !c.lockOriginalSize;
      if (c.lockOriginalSize) c.radiusScale = 1;
      showToast(
        c.lockOriginalSize
          ? "Node locked to original size."
          : "Node returned to auto collapsed sizing.",
        "info",
        1600,
      );
      drawTree();
    });
    btn("ctx-select-cluster", function () {
      if (!CTX_TARGET_DATA || !CURRENT_CLUSTERS) return;
      var name = CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? CURRENT_CLUSTERS[name] : null;
      if (cid == null) {
        showToast("No cluster for this node", "info");
        return;
      }
      var nCid = Number(cid);
      if (SELECTED_CLUSTER_IDS.has(nCid)) {
        SELECTED_CLUSTER_IDS.delete(nCid);
        if (SELECTED_CLUSTER_IDS.size === 0) {
          showToast("Cluster focus cleared", "info", 1800);
        } else {
          showToast(
            "Focused clusters: " +
              Array.from(SELECTED_CLUSTER_IDS)
                .sort((a, b) => a - b)
                .join(", "),
            "info",
            1800,
          );
        }
      } else {
        SELECTED_CLUSTER_IDS.add(nCid);
        showToast(
          "Focused clusters: " +
            Array.from(SELECTED_CLUSTER_IDS)
              .sort((a, b) => a - b)
              .join(", "),
          "success",
          1800,
        );
      }
      drawTree();
    });
    btn("ctx-rename-cluster", function () {
      if (!CTX_TARGET_DATA || !CURRENT_CLUSTERS) return;
      var name = CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? CURRENT_CLUSTERS[name] : null;
      if (cid == null) {
        showToast("No cluster for this node", "info");
        return;
      }
      var current = BOX_LABEL_MAP[cid] || "C" + cid;
      var val = prompt("Label for cluster " + cid + ":", current);
      if (val != null) {
        BOX_LABEL_MAP[cid] = val;
        drawTree();
      }
    });
    btn("ctx-adjust-box", function () {
      if (!CTX_TARGET_DATA || !CURRENT_CLUSTERS) return;
      var name = CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? CURRENT_CLUSTERS[name] : null;
      if (cid == null) {
        showToast("No cluster for this node", "info");
        return;
      }
      var adj = getBoxAdjust(cid);
      var padStr = prompt(
        "Box padding for cluster " + cid + " (padX, padY in px):",
        adj.padX + ", " + adj.padY,
      );
      if (padStr != null) {
        var parts = padStr.split(",").map(function (s) {
          return parseFloat(s.trim());
        });
        if (parts.length >= 1 && !isNaN(parts[0])) adj.padX = parts[0];
        if (parts.length >= 2 && !isNaN(parts[1])) adj.padY = parts[1];
        drawTree();
      }
    });
    btn("ctx-reset", function () {
      if (!CTX_TARGET_DATA) return;
      NODE_CUSTOM.delete(CTX_TARGET_DATA);
      drawTree();
    });
  })();

  // ── Mini scores panel toggle + drag ──
  var miniPanel = document.getElementById("mini-scores-panel");
  var miniScoresToggle = document.getElementById("mini-scores-toggle");
  var miniScoresBody = document.getElementById("mini-scores-body");
  if (miniScoresToggle && miniScoresBody) {
    miniScoresToggle.addEventListener("click", function (e) {
      e.stopPropagation();
      miniScoresBody.classList.toggle("collapsed");
      miniScoresToggle.innerHTML = miniScoresBody.classList.contains(
        "collapsed",
      )
        ? "&#x25B2;"
        : "&#x25BC;";
    });
  }
  // Draggable mini-scores panel
  if (miniPanel) {
    var miniHeader = miniPanel.querySelector(".mini-scores-header");
    if (miniHeader) {
      var dragState = {
        dragging: false,
        startX: 0,
        startY: 0,
        origLeft: 0,
        origTop: 0,
      };
      miniHeader.addEventListener("mousedown", function (e) {
        if (e.target === miniScoresToggle) return; // don't drag when clicking toggle
        dragState.dragging = true;
        var rect = miniPanel.getBoundingClientRect();
        dragState.startX = e.clientX;
        dragState.startY = e.clientY;
        dragState.origLeft = rect.left;
        dragState.origTop = rect.top;
        e.preventDefault();
      });
      document.addEventListener("mousemove", function (e) {
        if (!dragState.dragging) return;
        var dx = e.clientX - dragState.startX;
        var dy = e.clientY - dragState.startY;
        miniPanel.style.position = "fixed";
        miniPanel.style.left = dragState.origLeft + dx + "px";
        miniPanel.style.top = dragState.origTop + dy + "px";
        miniPanel.style.right = "auto";
        miniPanel.style.bottom = "auto";
      });
      document.addEventListener("mouseup", function () {
        dragState.dragging = false;
      });
    }
  }

  // ── Prefill example and auto-run ──
  if (newickEl) newickEl.value = EXAMPLE_NEWICK;
  updateClusterEditorAvailability();
  runPhytClust();
});
