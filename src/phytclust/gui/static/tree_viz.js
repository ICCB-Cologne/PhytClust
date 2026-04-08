/* ============================================================
   PhytClust – tree_viz.js
   ============================================================ */

// ── Globals ──
let HIER_CART = null;
let HIER_CIRC = null;

const EXAMPLE_NEWICK = "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

const newickEl  = document.getElementById("newick-input");
const resultEl  = document.getElementById("result");
const statusEl  = document.getElementById("status-message");
const statusDot = document.getElementById("status-dot");
const treeHost  = document.getElementById("tree_display");

const BASE_COLORS = [
    "#b84b4b", "#849060", "#3d7c74", "#6e3f8a",
    "#ceb94b", "#3f648a", "#3f408a", "#da63aa"
];

let SHOW_LEAF_NAMES     = true;
let SHOW_INTERNAL_NAMES = false;
let LEAF_NODE_RADIUS    = 3.0;
let INTERNAL_NODE_RADIUS = 1.8;
let LABEL_FONT_SIZE     = 9;
let TREE_WIDTH_SCALE    = 1.0;
let TREE_HEIGHT_SCALE   = 1.0;
let CLUSTER_COLORS      = [];
let CURRENT_CLUSTERS    = {};
let NEWICK_RAW_TREE     = null;
let CURRENT_LAYOUT_MODE = 'rectangular';
let COLOR_MODE          = 'bars';
const LABEL_PAD         = 6;

let isRunning          = false;
let latestOptimalKData = null;
let latestApiData      = null;
let LAST_TREE_SVG      = null;
let LAST_TREE_ZOOM     = null;
let SEARCH_TERM        = "";
let SHOW_BOX_LABELS    = true;
const BOX_LABEL_MAP    = {};  // cid -> custom name, editable via context menu

// Per-node customizations: keyed by data node reference
const NODE_CUSTOM = new WeakMap();

// Element references for params
const extraKEl          = document.getElementById("extra-k");
const extraOutgroupEl   = document.getElementById("extra-outgroup");
const extraRootTaxonEl  = document.getElementById("extra-root-taxon");
const extraTopNEl       = document.getElementById("extra-topn");
const extraResolutionEl = document.getElementById("extra-resolution");
const extraBinsEl       = document.getElementById("extra-bins");
const extraMaxKEl       = document.getElementById("extra-maxk");
const extraMaxKLimitEl  = document.getElementById("extra-maxklimit");
const extraLambdaEl     = document.getElementById("extra-lambda");
const extraMinClusterEl = document.getElementById("extra-min-cluster-size");
const extraOutlierEl    = document.getElementById("extra-outlier");

// D3 tooltip
const d3Tooltip = d3.select("body").append("div").attr("class", "d3-tooltip");


/* ─────────────────────────────────────
   Dark Mode
   ───────────────────────────────────── */
function initTheme() {
    const saved = localStorage.getItem("phytclust-theme");
    if (saved) document.documentElement.setAttribute("data-theme", saved);
}
initTheme();

function toggleTheme() {
    const html = document.documentElement;
    const next = html.getAttribute("data-theme") === "dark" ? "light" : "dark";
    html.setAttribute("data-theme", next);
    localStorage.setItem("phytclust-theme", next);
    // Redraw tree so branch colors update
    if (NEWICK_RAW_TREE) drawTree();
    if (latestOptimalKData) drawOptimalK(latestOptimalKData);
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
        setTimeout(function () { toast.remove(); }, 300);
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
        if (type === "info")    statusDot.classList.add("running");
        if (type === "success") statusDot.classList.add("ready");
        if (type === "danger")  statusDot.classList.add("error");
    }
    if (type === "danger")  showToast(message, "danger", 6000);
    else if (type === "success") showToast(message, "success", 3000);
}


/* ─────────────────────────────────────
   Utility Functions
   ───────────────────────────────────── */
function escapeRegExp(str) { return str.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }

function isOutgroupInNewick(newick, outgroup) {
    if (!outgroup) return true;
    var name = escapeRegExp(outgroup);
    return new RegExp("[\\(,]\\s*'?" + name + "'?\\s*(?:[:,\\)])").test(newick);
}

function estimateLeafCount(newick) {
    var inQuote = null, commas = 0;
    for (var i = 0; i < newick.length; i++) {
        var ch = newick[i];
        if (inQuote) { if (ch === inQuote) inQuote = null; continue; }
        if (ch === '"' || ch === "'") { inQuote = ch; continue; }
        if (ch === ',') commas++;
    }
    return Math.max(0, commas + 1);
}

function resetExtraParams() {
    ["extra-k", "extra-outgroup", "extra-root-taxon", "extra-topn",
     "extra-bins", "extra-maxk", "extra-maxklimit", "extra-lambda",
     "extra-min-cluster-size", "extra-min-prominence", "extra-outlier-threshold",
     "extra-min-support", "extra-support-weight"
    ].forEach(function(id) {
        var el = document.getElementById(id);
        if (el) el.value = "";
    });
    ["extra-outlier", "extra-optimize-polytomies"].forEach(function(id) {
        var el = document.getElementById(id);
        if (el) el.checked = true;
    });
    ["extra-resolution", "extra-compute-all", "extra-outlier-prefer-fewer",
     "extra-no-split-zero", "extra-use-support", "extra-relative-prom"
    ].forEach(function(id) {
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
    let g = ((num >> 8) & 0xFF) + factor * 255;
    let b = (num & 0xFF) + factor * 255;
    return `rgb(${Math.min(255,Math.max(0,r))|0}, ${Math.min(255,Math.max(0,g))|0}, ${Math.min(255,Math.max(0,b))|0})`;
}

function withAlpha(color, alpha) {
    if (color.startsWith("rgb")) return color.replace("rgb", "rgba").replace(")", `, ${alpha})`);
    const num = parseInt(color.slice(1), 16);
    return `rgba(${num >> 16}, ${(num >> 8) & 0xFF}, ${num & 0xFF}, ${alpha})`;
}

function generateClusterColors(nClusters) {
    let palette = shuffle(BASE_COLORS);
    if (nClusters <= palette.length) return palette.slice(0, nClusters);
    let colors = [];
    const repeats = Math.ceil(nClusters / palette.length);
    for (let r = 0; r < repeats; r++) {
        const factor = (r - 1) * 0.18;
        const alpha  = 1 - r * 0.15;
        palette.forEach(hex => {
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
    return {
        branch:   style.getPropertyValue('--pc-tree-branch').trim()   || '#94a3b8',
        label:    style.getPropertyValue('--pc-tree-label').trim()    || '#334155',
        internal: style.getPropertyValue('--pc-tree-internal').trim() || '#64748b',
        text:     style.getPropertyValue('--pc-text').trim()          || '#241b3d',
        muted:    style.getPropertyValue('--pc-text-muted').trim()    || '#9a91b2',
        bg:       style.getPropertyValue('--pc-bg').trim()            || '#f7f5fa',
        border:   style.getPropertyValue('--pc-border').trim()        || '#e4ddf0',
        accent:   style.getPropertyValue('--pc-accent').trim()        || '#6b61ac',
    };
}


/* ─────────────────────────────────────
   Tree Layout Computation
   ───────────────────────────────────── */
function accumulateBranchLength(node, length) {
    length = length || 0;
    node._bl = length;
    if (node.children) {
        for (const c of node.children) accumulateBranchLength(c, length + (c.length || 0));
    }
}

function computeLayouts() {
    if (!NEWICK_RAW_TREE) return;
    HIER_CART = d3.hierarchy(NEWICK_RAW_TREE, getVisibleChildren);
    HIER_CIRC = d3.hierarchy(NEWICK_RAW_TREE, getVisibleChildren);

    var width  = treeHost.clientWidth  || 800;
    var height = treeHost.clientHeight || 500;
    var margin = { top: 20, right: 80, bottom: 20, left: 80 };
    var innerW = width - margin.left - margin.right;
    var innerH = height - margin.top  - margin.bottom;
    var radius = Math.min(innerW, innerH) / 2;

    var maxBl = d3.max(HIER_CART.descendants(), d => d.data._bl || 0) || 1;
    var blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);
    var blToR = d3.scaleLinear().domain([0, maxBl]).range([0, radius]);

    d3.cluster().size([innerH, 1])(HIER_CART);
    HIER_CART.each(d => { d._x = d.x; d._y = blToX(d.data._bl || 0); });

    d3.cluster().size([2 * Math.PI, 1])(HIER_CIRC);
    HIER_CIRC.each(d => { d._angle = d.x; d._radius = blToR(d.data._bl || 0); });
}


/* ─────────────────────────────────────
   Node Context Menu
   ───────────────────────────────────── */
let CTX_TARGET_DATA = null;

function getNodeCustom(dataNode) {
    if (!NODE_CUSTOM.has(dataNode)) NODE_CUSTOM.set(dataNode, { highlighted: false, radiusScale: 1, renamedTo: null });
    return NODE_CUSTOM.get(dataNode);
}

function nodeDisplayName(dataNode) {
    const c = NODE_CUSTOM.has(dataNode) ? NODE_CUSTOM.get(dataNode) : null;
    return (c && c.renamedTo != null) ? c.renamedTo : (dataNode.name || "");
}

function customNodeRadius(d) {
    const base = (d.data && d.data.children && d.data.children.length) ? INTERNAL_NODE_RADIUS : LEAF_NODE_RADIUS;
    const c = NODE_CUSTOM.has(d.data) ? NODE_CUSTOM.get(d.data) : null;
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
    var c = NODE_CUSTOM.has(d.data) ? NODE_CUSTOM.get(d.data) : { radiusScale: 1 };
    if (ctxSlider) ctxSlider.value = c.radiusScale;
    if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);
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
    if (!kids) { if (node.name) out.push(node.name); return; }
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
function clearTree() { treeHost.innerHTML = ""; }

function clearAllCollapsedFlags(node) {
    if (!node) return;
    if (node._collapsed) delete node._collapsed;
    if (node.children && node.children.length) {
        for (const c of node.children) clearAllCollapsedFlags(c);
    }
}

function drawTree() {
    clearTree();
    if (!NEWICK_RAW_TREE) return;
    computeLayouts();

    const layoutMode = CURRENT_LAYOUT_MODE || "rectangular";
    const colorMode  = COLOR_MODE || "bars";
    const tc = getThemeColors();

    const container = d3.select("#tree_display");
    const width  = treeHost.clientWidth  || 800;
    const height = treeHost.clientHeight || 500;
    const margin = { top: 20, right: 80, bottom: 20, left: 80 };

    const svg = container.append("svg")
        .attr("width", width).attr("height", height)
        .style("background", "transparent");

    const zoomLayer = svg.append("g");
    const g = zoomLayer.append("g")
        .attr("transform",
            layoutMode === "circular"
                ? `translate(${width/2},${height/2})`
                : `translate(${margin.left},${margin.top})`
        );

    const zoom = d3.zoom().scaleExtent([0.3, 8])
        .on("zoom", (event) => { zoomLayer.attr("transform", event.transform); });
    svg.call(zoom);
    LAST_TREE_SVG  = svg;
    LAST_TREE_ZOOM = zoom;

    function hasOriginalChildren(d) { return !!(d && d.data && Array.isArray(d.data.children) && d.data.children.length > 0); }
    function canCollapse(d) { return hasOriginalChildren(d); }
    function isVisibleLeaf(d) { return !(d && d.children && d.children.length); }
    function nodeRadius(d) { return customNodeRadius(d); }
    function isHighlighted(d) {
        if (isSearchMatch(d.data)) return true;
        return NODE_CUSTOM.has(d.data) && NODE_CUSTOM.get(d.data).highlighted;
    }

    function highlightColor() { return "#f59e0b"; }

    function nodeFill(d) {
        const name = d.data && d.data.name;
        if (colorMode === "bars") return tc.label;
        if (name && !hasOriginalChildren(d)) {
            const cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[name] : null;
            if (cid == null || !CLUSTER_COLORS.length) return tc.label;
            return CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
        }
        const rep = representativeClusterIdFromData(d.data);
        if (rep != null && CLUSTER_COLORS.length) return CLUSTER_COLORS[rep % CLUSTER_COLORS.length];
        return tc.label;
    }

    // ── CIRCULAR LAYOUT ──
    if (layoutMode === "circular") {
        const root = HIER_CIRC;
        const maxRadiusPx = Math.min(width, height) / 2 - 30;
        const maxR = d3.max(root.descendants(), d => d._radius || 0) || 1;
        const rScale = d3.scaleLinear().domain([0, maxR]).range([0, maxRadiusPx]);

        function polarToXY(angle, radius) {
            const a = angle - Math.PI / 2;
            const rr = rScale(radius);
            return [(rr * Math.cos(a)) * TREE_WIDTH_SCALE, (rr * Math.sin(a)) * TREE_HEIGHT_SCALE];
        }
        function nodeXY(d) { return polarToXY(d._angle || 0, d._radius || 0); }

        function radialElbowPath(link) {
            const s = link.source, t = link.target;
            const sa = s._angle || 0, sr = s._radius || 0;
            const ta = t._angle || 0, tr = t._radius || 0;
            const p0 = polarToXY(sa, sr), p1 = polarToXY(ta, sr), p2 = polarToXY(ta, tr);
            const delta = ta - sa, sweep = delta >= 0 ? 1 : 0;
            const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
            const arcR = rScale(sr);
            return `M${p0[0]},${p0[1]}A${arcR},${arcR} 0 ${largeArc},${sweep} ${p1[0]},${p1[1]}L${p2[0]},${p2[1]}`;
        }

        g.append("g").selectAll(".tree-link")
            .data(root.links()).enter().append("path")
            .attr("class", "tree-link").attr("fill", "none")
            .attr("stroke", tc.branch).attr("stroke-opacity", 0.6).attr("stroke-width", 1.1)
            .attr("d", radialElbowPath);

        const node = g.append("g").selectAll(".tree-node")
            .data(root.descendants()).enter().append("g")
            .attr("class", d => "tree-node" + (isSearchMatch(d.data) ? " node-search-match" : ""))
            .attr("transform", d => { const p = nodeXY(d); return `translate(${p[0]},${p[1]})`; });

        node.append("circle")
            .attr("r", nodeRadius).attr("fill", nodeFill)
            .attr("stroke", d => isHighlighted(d) ? highlightColor() : "none")
            .attr("stroke-width", d => isHighlighted(d) ? 2.5 : 0)
            .style("cursor", "pointer")
            .on("click", function (event, d) { if (!canCollapse(d)) return; d.data._collapsed = !d.data._collapsed; drawTree(); })
            .on("contextmenu", showNodeContextMenu);

        node.filter(d => isVisibleLeaf(d)).append("text")
            .attr("dy", "0.32em").attr("font-size", LABEL_FONT_SIZE + "px")
            .attr("fill", d => isHighlighted(d) ? highlightColor() : tc.label)
            .attr("font-weight", d => isHighlighted(d) ? "bold" : "normal")
            .attr("transform", function (d) {
                const a = d._angle || 0;
                let deg = (a * 180 / Math.PI) - 90;
                const left = (deg > 90 || deg < -90);
                return `rotate(${deg}) translate(${nodeRadius(d) + LABEL_PAD},0) rotate(${left ? 180 : 0})`;
            })
            .style("text-anchor", d => { const deg = ((d._angle || 0) * 180 / Math.PI) - 90; return (deg > 90 || deg < -90) ? "end" : "start"; })
            .text(d => SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : "");

        node.filter(d => !isVisibleLeaf(d)).append("text")
            .attr("dy", "0.32em").attr("font-size", (LABEL_FONT_SIZE - 1) + "px")
            .attr("fill", d => isHighlighted(d) ? highlightColor() : tc.internal)
            .attr("transform", d => { const a = d._angle || 0; let deg = (a * 180 / Math.PI) - 90; return `rotate(${deg}) translate(0,${-(nodeRadius(d) + 5)})`; })
            .style("text-anchor", "middle")
            .text(d => SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : "");

        // Circular side bars
        if (colorMode === "bars" && CURRENT_CLUSTERS && CLUSTER_COLORS.length) {
            const visLeaves = root.descendants().filter(d => isVisibleLeaf(d)).sort((a, b) => (a._angle || 0) - (b._angle || 0));
            if (visLeaves.length >= 2) {
                const angles = visLeaves.map(d => d._angle || 0);
                const boundaries = new Array(visLeaves.length + 1);
                for (let i = 1; i < visLeaves.length; i++) boundaries[i] = (angles[i - 1] + angles[i]) / 2;
                boundaries[0] = angles[0] - (boundaries[1] - angles[0]);
                boundaries[visLeaves.length] = angles[visLeaves.length - 1] + (angles[visLeaves.length - 1] - boundaries[visLeaves.length - 1]);

                function visLeafCid(d) {
                    const nm = d.data && d.data.name;
                    const isOrigLeaf = nm && !(d.data.children && d.data.children.length);
                    if (isOrigLeaf) return CURRENT_CLUSTERS[nm];
                    return representativeClusterIdFromData(d.data);
                }

                const leafLabelChars = visLeaves.reduce((mx, d) => Math.max(mx, (d.data.name || "").length), 0);
                const estRadialLabelWidth = SHOW_LEAF_NAMES ? leafLabelChars * LABEL_FONT_SIZE * 0.55 + LABEL_PAD + LEAF_NODE_RADIUS : 10;
                const ringInner = maxRadiusPx + estRadialLabelWidth + 6;
                const ringOuter = ringInner + 14;
                const arc = d3.arc().innerRadius(ringInner).outerRadius(ringOuter);

                let runStart = 0, runCid = visLeafCid(visLeaves[0]);
                for (let i = 1; i <= visLeaves.length; i++) {
                    const cid = (i < visLeaves.length) ? visLeafCid(visLeaves[i]) : Symbol("END");
                    if (cid !== runCid) {
                        if (runCid != null) {
                            g.append("path").attr("d", arc.startAngle(boundaries[runStart]).endAngle(boundaries[i]))
                                .attr("fill", CLUSTER_COLORS[runCid % CLUSTER_COLORS.length]).attr("opacity", 0.75);
                        }
                        runStart = i;
                        runCid = (i < visLeaves.length) ? cid : null;
                    }
                }
            }
        }
        return;
    }

    // ── CARTESIAN (rectangular / cladogram) ──
    const root = HIER_CART;
    const allNodes = root.descendants();
    const maxYCart = d3.max(allNodes, d => d._y || 0) || 0;
    const leafNames = allNodes.filter(d => !d.children || !d.children.length).map(d => (d.data.name || "").length);
    const maxLabelChars = d3.max(leafNames) || 0;
    const estLabelWidth = SHOW_LEAF_NAMES ? maxLabelChars * LABEL_FONT_SIZE * 0.6 + LABEL_PAD + LEAF_NODE_RADIUS : 20;
    const labelColumnX = (maxYCart * TREE_WIDTH_SCALE) + estLabelWidth + 12;

    g.append("g").selectAll(".tree-link")
        .data(root.links()).enter().append("path")
        .attr("class", "tree-link").attr("fill", "none")
        .attr("stroke", tc.branch).attr("stroke-opacity", 0.6).attr("stroke-width", 1.2)
        .attr("d", function (d) {
            if (layoutMode === "rectangular") {
                return "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
                       "V" + (d.target._x * TREE_HEIGHT_SCALE) +
                       "H" + (d.target._y * TREE_WIDTH_SCALE);
            }
            return "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
                   "L" + (d.target._y * TREE_WIDTH_SCALE) + "," + (d.target._x * TREE_HEIGHT_SCALE);
        });

    const node = g.append("g").selectAll(".tree-node")
        .data(root.descendants()).enter().append("g")
        .attr("class", d => "tree-node" + (isSearchMatch(d.data) ? " node-search-match" : ""))
        .attr("transform", d => `translate(${(d._y || 0) * TREE_WIDTH_SCALE},${(d._x || 0) * TREE_HEIGHT_SCALE})`);

    node.append("circle")
        .attr("r", nodeRadius).attr("fill", nodeFill)
        .attr("stroke", d => isHighlighted(d) ? highlightColor() : "none")
        .attr("stroke-width", d => isHighlighted(d) ? 2.5 : 0)
        .style("cursor", "pointer")
        .on("click", function (event, d) { if (!canCollapse(d)) return; d.data._collapsed = !d.data._collapsed; drawTree(); })
        .on("contextmenu", showNodeContextMenu);

    // Leaf labels
    node.filter(d => !(d.children && d.children.length)).append("text")
        .attr("dy", 3)
        .attr("x", d => nodeRadius(d) + LABEL_PAD)
        .style("text-anchor", "start")
        .style("font-size", LABEL_FONT_SIZE + "px")
        .attr("fill", d => isHighlighted(d) ? highlightColor() : tc.label)
        .attr("font-weight", d => isHighlighted(d) ? "bold" : "normal")
        .text(d => SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : "");

    // Internal node labels — offset below-right to avoid branch overlap
    node.filter(d => !!(d.children && d.children.length)).append("text")
        .attr("dy", d => nodeRadius(d) + 12)
        .attr("x", d => nodeRadius(d) + 3)
        .style("text-anchor", "start")
        .style("font-size", (LABEL_FONT_SIZE - 1) + "px")
        .style("font-style", "italic")
        .attr("fill", d => isHighlighted(d) ? highlightColor() : tc.internal)
        .attr("font-weight", d => isHighlighted(d) ? "bold" : "normal")
        .text(d => SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : "");

    // Cluster side bars
    if (colorMode === "bars") {
        const leafNodes = allNodes.filter(d => !d.children || !d.children.length).sort((a, b) => (a._x || 0) - (b._x || 0));
        const barX = labelColumnX, barW = 20;
        const ys = leafNodes.map(d => (d._x || 0) * TREE_HEIGHT_SCALE);

        if (leafNodes.length === 1) {
            const nm = leafNodes[0].data && leafNodes[0].data.name;
            const cid = nm && CURRENT_CLUSTERS ? CURRENT_CLUSTERS[nm] : null;
            if (cid != null && CLUSTER_COLORS.length) {
                g.append("rect").attr("x", barX).attr("y", ys[0] - 6).attr("width", barW).attr("height", 12).attr("rx", 3)
                    .attr("fill", CLUSTER_COLORS[cid % CLUSTER_COLORS.length]).attr("opacity", 0.75);
            }
        } else if (leafNodes.length > 1) {
            const boundaries = new Array(leafNodes.length + 1);
            for (let i = 1; i < leafNodes.length; i++) boundaries[i] = (ys[i - 1] + ys[i]) / 2;
            boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
            boundaries[leafNodes.length] = ys[leafNodes.length - 1] + (ys[leafNodes.length - 1] - boundaries[leafNodes.length - 1]);

            const barsG = g.append("g").attr("class", "cluster-bars");
            function leafCid(i) {
                const nm = leafNodes[i].data && leafNodes[i].data.name;
                return nm ? (CURRENT_CLUSTERS ? CURRENT_CLUSTERS[nm] : null) : null;
            }

            let runStart = 0, runCid = leafCid(0);
            for (let i = 1; i <= leafNodes.length; i++) {
                const cid = (i < leafNodes.length) ? leafCid(i) : Symbol("END");
                if (cid !== runCid) {
                    if (runCid != null && CLUSTER_COLORS.length) {
                        const y0 = boundaries[runStart], y1 = boundaries[i];
                        barsG.append("rect").attr("x", barX).attr("y", y0)
                            .attr("width", barW).attr("height", Math.max(1, y1 - y0)).attr("rx", 3)
                            .attr("fill", CLUSTER_COLORS[runCid % CLUSTER_COLORS.length]).attr("opacity", 0.75);
                    }
                    runStart = i;
                    runCid = (i < leafNodes.length) ? cid : null;
                }
            }
        }
    }

    // Cluster boxes (MRCA mode) — draw coloured rectangles from MRCA to leaf tips
    if (colorMode === "boxes" && CURRENT_CLUSTERS && CLUSTER_COLORS.length) {
        // Group leaves by cluster id
        const clusterGroups = {};
        allNodes.forEach(d => {
            if (d.children && d.children.length) return;
            const nm = d.data && d.data.name;
            const cid = nm ? (CURRENT_CLUSTERS[nm] != null ? CURRENT_CLUSTERS[nm] : null) : null;
            if (cid == null || cid < 0) return;
            if (!clusterGroups[cid]) clusterGroups[cid] = [];
            clusterGroups[cid].push(d);
        });
        // Find MRCA for each cluster group
        function findMRCA(nodes) {
            if (!nodes.length) return null;
            if (nodes.length === 1) return nodes[0];
            // Collect ancestors for each node
            function ancestors(n) { const a = []; let cur = n; while (cur) { a.push(cur); cur = cur.parent; } return a; }
            let common = ancestors(nodes[0]);
            for (let i = 1; i < nodes.length; i++) {
                const anc = new Set(ancestors(nodes[i]));
                common = common.filter(n => anc.has(n));
            }
            return common[0]; // first (deepest) common ancestor
        }
        // Estimate rightmost extent including leaf labels
        var estLabelPx = SHOW_LEAF_NAMES ? (maxLabelChars * LABEL_FONT_SIZE * 0.6 + LABEL_PAD + LEAF_NODE_RADIUS + 16) : 16;
        var boxPadV = 8;   // vertical padding around topmost/bottommost leaf
        var boxPadL = 6;   // left padding before MRCA
        var boxesG = g.append("g").attr("class", "cluster-boxes");
        for (var cid of Object.keys(clusterGroups)) {
            var leaves = clusterGroups[cid];
            var mrca = findMRCA(leaves);
            if (!mrca) continue;
            var ys = leaves.map(function(d) { return (d._x || 0) * TREE_HEIGHT_SCALE; });
            var minY = Math.min.apply(null, ys) - boxPadV;
            var maxY = Math.max.apply(null, ys) + boxPadV;
            var mrcaX = (mrca._y || 0) * TREE_WIDTH_SCALE - boxPadL;
            var maxLeafX = Math.max.apply(null, leaves.map(function(d) { return (d._y || 0) * TREE_WIDTH_SCALE; }));
            var boxRight = maxLeafX + estLabelPx;
            var boxW = boxRight - mrcaX;
            var color = CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
            boxesG.append("rect")
                .attr("x", mrcaX).attr("y", minY)
                .attr("width", Math.max(boxW, 8)).attr("height", Math.max(maxY - minY, 4))
                .attr("rx", 6).attr("ry", 6)
                .attr("fill", color)
                .attr("opacity", 0.12)
                .attr("stroke", color)
                .attr("stroke-width", 1.5)
                .attr("stroke-opacity", 0.45);
            // Cluster label on right edge of box
            if (SHOW_BOX_LABELS) {
                boxesG.append("text")
                    .attr("x", boxRight + 4)
                    .attr("y", (minY + maxY) / 2)
                    .attr("dy", "0.35em")
                    .attr("font-size", (LABEL_FONT_SIZE - 1) + "px")
                    .attr("font-weight", "600")
                    .attr("fill", color)
                    .text(BOX_LABEL_MAP[cid] || ("C" + cid));
            }
        }
    }

    // Branch-length axis
    const maxX  = d3.max(allNodes, d => d._x || 0) || 0;
    const maxY  = d3.max(allNodes, d => d._y || 0) || 0;
    const maxBl = d3.max(allNodes, d => d.data._bl || 0) || 1;
    const axisY = (maxX * TREE_HEIGHT_SCALE) + 20;
    const blScale = d3.scaleLinear().domain([0, maxBl]).range([0, maxY * TREE_WIDTH_SCALE]);

    const axisG = g.append("g").attr("class", "branch-length-axis").attr("transform", `translate(0, ${axisY})`)
        .call(d3.axisBottom(blScale).ticks(5));
    axisG.selectAll("text").attr("fill", tc.internal).style("font-size", "10px");
    axisG.selectAll("line").attr("stroke", tc.branch);
    axisG.selectAll("path").attr("stroke", tc.branch);
    axisG.append("text").attr("x", (maxY * TREE_WIDTH_SCALE) / 2).attr("y", 30)
        .attr("text-anchor", "middle").attr("font-size", 10).attr("fill", tc.internal).text("Branch length");

    svg.attr("height", Math.max(height, axisY + 50 + margin.bottom));
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
    const width  = hostEl.clientWidth  || 400;
    const height = hostEl.clientHeight || 400;
    const margin = { top: 16, right: 60, bottom: 16, left: 60 };
    const innerW = width - margin.left - margin.right;
    const innerH = height - margin.top  - margin.bottom;

    const maxBl = d3.max(hier.descendants(), d => d.data._bl || 0) || 1;
    const blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);

    d3.cluster().size([innerH, 1])(hier);
    hier.each(d => { d._x = d.x; d._y = blToX(d.data._bl || 0); });

    const nClusters = clusters ? Math.max(0, ...Object.values(clusters)) + 1 : 0;
    const colors = nClusters > 0 ? generateClusterColors(nClusters) : [];

    const svg = d3.select(hostEl).append("svg").attr("width", width).attr("height", height).style("background", "transparent");
    const zoomLayer = svg.append("g");
    const g = zoomLayer.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
    svg.call(d3.zoom().scaleExtent([0.3, 8]).on("zoom", (event) => { zoomLayer.attr("transform", event.transform); }));

    g.append("g").selectAll(".tree-link").data(hier.links()).enter().append("path")
        .attr("class", "tree-link").attr("fill", "none")
        .attr("stroke", tc.branch).attr("stroke-opacity", 0.5).attr("stroke-width", 1)
        .attr("d", d => "M" + d.source._y + "," + d.source._x + "V" + d.target._x + "H" + d.target._y);

    const allNodes = hier.descendants();
    const node = g.append("g").selectAll(".tree-node").data(allNodes).enter().append("g")
        .attr("class", "tree-node").attr("transform", d => `translate(${d._y},${d._x})`);

    node.append("circle")
        .attr("r", d => (d.children && d.children.length) ? 1.5 : 2.5)
        .attr("fill", function(d) {
            const name = d.data && d.data.name;
            if (name && !(d.children && d.children.length) && clusters) {
                const cid = clusters[name];
                if (cid != null && colors.length) return colors[cid % colors.length];
            }
            return tc.label;
        });

    node.filter(d => !(d.children && d.children.length)).append("text")
        .attr("dy", 3).attr("x", 6).style("text-anchor", "start").style("font-size", "8px")
        .attr("fill", tc.label).text(d => d.data.name || "");

    // Side bars
    if (clusters && colors.length) {
        const leaves = allNodes.filter(d => !(d.children && d.children.length)).sort((a, b) => a._x - b._x);
        if (leaves.length > 1) {
            const estLabel = d3.max(leaves, d => (d.data.name || "").length) * 5 + 12;
            const barX = (d3.max(allNodes, d => d._y) || 0) + estLabel;
            const ys = leaves.map(d => d._x);
            const boundaries = new Array(leaves.length + 1);
            for (let i = 1; i < leaves.length; i++) boundaries[i] = (ys[i - 1] + ys[i]) / 2;
            boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
            boundaries[leaves.length] = ys[leaves.length - 1] + (ys[leaves.length - 1] - boundaries[leaves.length - 1]);

            let runStart = 0, runCid = clusters[leaves[0].data.name];
            for (let i = 1; i <= leaves.length; i++) {
                const cid = (i < leaves.length) ? clusters[leaves[i].data.name] : Symbol("END");
                if (cid !== runCid) {
                    if (runCid != null) {
                        g.append("rect").attr("x", barX).attr("y", boundaries[runStart])
                            .attr("width", 14).attr("height", Math.max(1, boundaries[i] - boundaries[runStart]))
                            .attr("rx", 2).attr("fill", colors[runCid % colors.length]).attr("opacity", 0.7);
                    }
                    runStart = i; runCid = (i < leaves.length) ? cid : null;
                }
            }
        }
    }
}


/* ─────────────────────────────────────
   API Call
   ───────────────────────────────────── */
async function apiPostJson(url, payload) {
    const res = await fetch(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) });
    const text = await res.text();
    let data = null;
    try { data = text ? JSON.parse(text) : null; } catch { /* ignore */ }
    if (!res.ok) { throw new Error((data && data.detail) ? data.detail : `Request failed (${res.status})`); }
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
    if (!newickText) { showToast("Please upload or paste a Newick tree.", "danger", 4000); return; }

    var numSamples = estimateLeafCount(newickText);
    var mode = getCurrentMode();

    // Read params
    var kVal = readIntParam("extra-k");
    var outgroupVal = (extraOutgroupEl ? (extraOutgroupEl.value || "").trim() : "") || null;
    var rootTaxonVal = (extraRootTaxonEl ? (extraRootTaxonEl.value || "").trim() : "") || null;
    var topNVal = readIntParam("extra-topn");
    var binsVal = readIntParam("extra-bins");
    var maxKVal = readIntParam("extra-maxk");
    var maxKLimitVal = readFloatParam("extra-maxklimit");
    var lambdaVal = readFloatParam("extra-lambda");
    var minClusterVal = readIntParam("extra-min-cluster-size");

    if (outgroupVal && !isOutgroupInNewick(newickText, outgroupVal)) {
        showStatus("Outgroup not found in Newick.", "danger"); return;
    }

    if (extraResolutionEl) extraResolutionEl.checked = (mode === "resolution");

    // Build payload
    var payload = { newick: newickText, mode: mode };
    if (mode === "k") {
        if (kVal === null) { showToast("Please enter a value for k.", "danger"); return; }
        payload.k = kVal;
    }
    if (outgroupVal)           payload.outgroup = outgroupVal;
    if (rootTaxonVal)          payload.root_taxon = rootTaxonVal;
    if (topNVal !== null)      payload.top_n = topNVal;
    if (binsVal !== null)      payload.num_bins = binsVal;
    if (maxKVal !== null)      payload.max_k = maxKVal;
    if (maxKLimitVal !== null)  payload.max_k_limit = maxKLimitVal;
    if (lambdaVal !== null)     payload.lambda_weight = lambdaVal;
    if (minClusterVal !== null) payload.min_cluster_size = minClusterVal;
    if (mode === "resolution")  payload.by_resolution = true;

    // New params
    if (readCheckParam("extra-compute-all"))  payload.compute_all_clusters = true;
    if (readCheckParam("extra-use-support"))  payload.use_branch_support = true;
    var minSup = readFloatParam("extra-min-support");
    if (minSup !== null) payload.min_support = minSup;
    var supW = readFloatParam("extra-support-weight");
    if (supW !== null) payload.support_weight = supW;

    // Outlier config
    var outlierThresh = readIntParam("extra-outlier-threshold");
    if (outlierThresh !== null) payload.outlier_size_threshold = outlierThresh;
    if (readCheckParam("extra-outlier-prefer-fewer")) payload.outlier_prefer_fewer = true;
    var ratioMode = readSelectParam("extra-outlier-ratio-mode");
    if (ratioMode && ratioMode !== "exp") payload.outlier_ratio_mode = ratioMode;

    // Polytomy config
    payload.optimize_polytomies = readCheckParam("extra-optimize-polytomies");
    if (readCheckParam("extra-no-split-zero")) payload.no_split_zero_length = true;

    // Peak config
    var rankMode = readSelectParam("extra-ranking-mode");
    if (rankMode) payload.ranking_mode = rankMode;
    var minProm = readFloatParam("extra-min-prominence");
    if (minProm !== null) payload.min_prominence = minProm;
    if (readCheckParam("extra-relative-prom")) payload.use_relative_prominence = true;

    showStatus("Running PhytClust...", "info");
    resultEl.textContent = "Running PhytClust...";
    clearTree();

    var runBtn = document.getElementById("btn-run");
    try {
        isRunning = true;
        if (runBtn) { runBtn.disabled = true; runBtn.innerHTML = '<span class="spinner"></span> Running...'; }

        const t0 = performance.now();
        const data = await apiPostJson("/api/run", payload);
        const dt = (performance.now() - t0) / 1000;

        latestApiData = data;
        showStatus(`Finished in ${dt.toFixed(2)}s`, "success");
        resultEl.textContent = JSON.stringify(data, null, 2);

        // Update leaf count label
        var lcLabel = document.getElementById("leaf-count-label");
        if (lcLabel && data.newick) lcLabel.textContent = estimateLeafCount(data.newick) + " leaves";

        if (data.newick) {
            NEWICK_RAW_TREE = parseNewick(data.newick);
            accumulateBranchLength(NEWICK_RAW_TREE);
            computeLayouts();
            populateClusterSelector(data);
            populateCompareSelectors(data);

            const clusterMap = (data.clusters && data.clusters.length > 0) ? data.clusters[0] : {};
            CURRENT_CLUSTERS = clusterMap;
            if (Object.keys(clusterMap).length > 0) {
                CLUSTER_COLORS = generateClusterColors(Math.max(...Object.values(clusterMap)) + 1);
            } else {
                CLUSTER_COLORS = [];
                showStatus("No clusters found.", "danger");
            }
            drawTree();
        } else {
            clearTree();
            showStatus("No Newick tree returned by API.", "danger");
        }

        // Draw mini scores panel on tree page
        try { drawMiniScores(data); } catch (e) { console.warn("Mini scores error:", e); }

        // Draw optimal k plot
        try {
            if (data && data.scores) {
                var plotHost = document.getElementById('optimalk_plot');
                if (plotHost) plotHost.innerHTML = '';
                latestOptimalKData = data;
                drawOptimalK(latestOptimalKData);
            } else {
                var plotHost2 = document.getElementById('optimalk_plot');
                if (plotHost2) plotHost2.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);font-size:14px;">Score plot unavailable in fixed k mode.</div>';
            }
        } catch (plotErr) { console.warn('Plot rendering error:', plotErr); }

    } catch (e) {
        console.error(e);
        showStatus("Error: " + e.message, "danger");
        resultEl.textContent = "Error: " + e;
        latestOptimalKData = null; latestApiData = null;
    } finally {
        isRunning = false;
        if (runBtn) { runBtn.disabled = false; runBtn.innerHTML = '<svg width="14" height="14" viewBox="0 0 16 16" fill="currentColor"><polygon points="3,1 13,8 3,15"/></svg> Run PhytClust'; }
    }
}


/* ─────────────────────────────────────
   Cluster Selector (multiple results)
   ───────────────────────────────────── */
let CLUSTER_VIEW_MODE = "peaks";

function populateClusterSelector(data) {
    const controls = document.getElementById("cluster-select-controls");
    const selectEl = document.getElementById("cluster-select");
    const labelEl  = document.getElementById("cluster-select-label");
    const toggleEl = document.getElementById("cluster-view-toggle");
    if (!controls || !selectEl) return;

    const hasPeaks = (data.clusters || []).length > 1;
    const hasAll   = !!(data.all_clusters && data.all_clusters.length > 1);

    if (!hasPeaks && !hasAll) { controls.classList.remove("visible"); if (toggleEl) toggleEl.style.display = "none"; return; }
    if (toggleEl) toggleEl.style.display = hasAll ? "inline-flex" : "none";

    if (hasAll && CLUSTER_VIEW_MODE === "all") {
        _fillClusterSelect(selectEl, labelEl, data.all_clusters, data.all_ks, false);
    } else {
        CLUSTER_VIEW_MODE = "peaks";
        _fillClusterSelect(selectEl, labelEl, data.clusters, data.ks, true);
    }
    controls.classList.add("visible");
}

function _fillClusterSelect(selectEl, labelEl, clusters, ks, showRank) {
    selectEl.innerHTML = "";
    for (let i = 0; i < clusters.length; i++) {
        const opt = document.createElement("option");
        opt.value = i;
        const kLabel = ks && ks[i] != null ? ks[i] : (i + 1);
        opt.textContent = showRank ? `k = ${kLabel} (rank ${i + 1})` : `k = ${kLabel}`;
        selectEl.appendChild(opt);
    }
    if (labelEl) labelEl.textContent = `1 / ${clusters.length}`;
}

function switchCluster(idx) {
    if (!latestApiData) return;
    let clusters, ks;
    if (CLUSTER_VIEW_MODE === "all" && latestApiData.all_clusters) { clusters = latestApiData.all_clusters; ks = latestApiData.all_ks; }
    else { clusters = latestApiData.clusters || []; ks = latestApiData.ks || []; }
    if (idx < 0 || idx >= clusters.length) return;
    CURRENT_CLUSTERS = clusters[idx];
    if (Object.keys(CURRENT_CLUSTERS).length > 0) {
        CLUSTER_COLORS = generateClusterColors(Math.max(...Object.values(CURRENT_CLUSTERS)) + 1);
    } else { CLUSTER_COLORS = []; }
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
    CLUSTER_VIEW_MODE = (CLUSTER_VIEW_MODE === "peaks") ? "all" : "peaks";
    populateClusterSelector(latestApiData);
    switchCluster(0);
    const toggleEl = document.getElementById("cluster-view-toggle");
    if (toggleEl) toggleEl.textContent = CLUSTER_VIEW_MODE === "all" ? "Show peaks only" : "Show all k";
}


/* ─────────────────────────────────────
   Comparison Mode
   ───────────────────────────────────── */
function populateCompareSelectors(data) {
    const selA = document.getElementById("compare-k-a");
    const selB = document.getElementById("compare-k-b");
    if (!selA || !selB) return;

    const ks = (data.all_ks && data.all_ks.length) ? data.all_ks : (data.ks || []);
    selA.innerHTML = ""; selB.innerHTML = "";

    ks.forEach(function(k) {
        const optA = document.createElement("option"); optA.value = k; optA.textContent = "k = " + k;
        const optB = document.createElement("option"); optB.value = k; optB.textContent = "k = " + k;
        selA.appendChild(optA); selB.appendChild(optB);
    });

    if (ks.length >= 2) { selA.value = ks[0]; selB.value = ks[Math.min(1, ks.length - 1)]; }
}

function drawComparison() {
    if (!latestApiData || !NEWICK_RAW_TREE) return;

    const selA = document.getElementById("compare-k-a");
    const selB = document.getElementById("compare-k-b");
    const leftHost  = document.getElementById("compare-tree-left");
    const rightHost = document.getElementById("compare-tree-right");
    const leftLabel  = document.getElementById("compare-left-label");
    const rightLabel = document.getElementById("compare-right-label");

    if (!selA || !selB || !leftHost || !rightHost) return;

    const kA = parseInt(selA.value, 10);
    const kB = parseInt(selB.value, 10);

    function getClusters(k) {
        if (latestApiData.all_clusters && latestApiData.all_ks) {
            const idx = latestApiData.all_ks.indexOf(k);
            if (idx >= 0) return latestApiData.all_clusters[idx];
        }
        if (latestApiData.ks && latestApiData.clusters) {
            const idx = latestApiData.ks.indexOf(k);
            if (idx >= 0) return latestApiData.clusters[idx];
        }
        return {};
    }

    if (leftLabel) leftLabel.textContent = "k = " + kA;
    if (rightLabel) rightLabel.textContent = "k = " + kB;

    drawTreeInto(leftHost, getClusters(kA));
    drawTreeInto(rightHost, getClusters(kB));
}


/* ─────────────────────────────────────
   Optimal k Plot
   ───────────────────────────────────── */
function drawOptimalK(data) {
    if (!data) data = latestOptimalKData;
    var plotEl = document.getElementById('optimalk_plot');
    if (!plotEl) return;

    var scores = data && data.scores ? data.scores : [];
    var peaks  = data && data.peaks  ? data.peaks  : [];
    const tc = getThemeColors();

    if (!Array.isArray(scores) || scores.length === 0) {
        plotEl.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);">No scores available to plot.</div>';
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
        return;
    }

    plotEl.innerHTML = "";
    var dataPoints = [];
    for (let i = 0; i < scores.length - 1; i++) dataPoints.push({ k: i + 2, score: scores[i + 1] });

    var peakPoints = peaks
        .map(k => { var idx = k - 2; return (idx >= 0 && idx < dataPoints.length) ? { k, score: dataPoints[idx].score } : null; })
        .filter(d => d);

    var width  = plotEl.clientWidth  || 700;
    var height = plotEl.clientHeight || 420;
    var margin = { top: 24, right: 24, bottom: 48, left: 58 };
    var innerWidth  = width  - margin.left - margin.right;
    var innerHeight = height - margin.top  - margin.bottom;

    var svg = d3.select(plotEl).append("svg").attr("width", width).attr("height", height).style("background", "transparent");
    var g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    var axisModeEl = document.getElementById('axis-mode');
    var defaultMode = (scores.length > 50) ? 'log' : 'normal';
    var mode = defaultMode;
    if (axisModeEl) {
        if (!axisModeEl.__initialized) { axisModeEl.value = defaultMode; axisModeEl.__initialized = true; }
        if (axisModeEl.__userOverride) mode = axisModeEl.value || defaultMode;
        else axisModeEl.value = defaultMode;
    }

    var xScale, xAxis;
    if (mode === 'log') {
        xScale = d3.scaleLog().domain([2, d3.max(dataPoints, d => d.k)]).range([0, innerWidth]);
        xAxis  = d3.axisBottom(xScale).tickValues(dataPoints.map(d => d.k)).tickFormat(d3.format("d"));
    } else {
        xScale = d3.scaleBand().domain(dataPoints.map(d => d.k)).range([0, innerWidth]).padding(0.2);
        var nPts = dataPoints.length;
        var dtick = Math.max(1, Math.ceil(nPts / 10));
        xAxis = d3.axisBottom(xScale).tickValues(dataPoints.map(d => d.k).filter((d, i) => i % dtick === 0)).tickFormat(d3.format("d"));
    }

    const yScale = d3.scaleLinear()
        .domain([d3.min(dataPoints, d => d.score), d3.max(dataPoints, d => d.score)])
        .nice().range([innerHeight, 0]);

    // Grid
    g.append("g").attr("class", "grid").attr("transform", `translate(0,${innerHeight})`)
        .call(d3.axisBottom(xScale).tickSize(-innerHeight).tickFormat(""))
        .selectAll("line").attr("stroke", tc.border).attr("stroke-dasharray", "2,2");
    g.append("g").attr("class", "grid")
        .call(d3.axisLeft(yScale).tickSize(-innerWidth).tickFormat(""))
        .selectAll("line").attr("stroke", tc.border).attr("stroke-dasharray", "2,2");
    g.selectAll(".grid .domain").remove();

    // Axes
    var xAxisG = g.append("g").attr("transform", `translate(0,${innerHeight})`).call(xAxis);
    var yAxisG = g.append("g").call(d3.axisLeft(yScale).ticks(6));
    [xAxisG, yAxisG].forEach(ag => {
        ag.selectAll("text").attr("fill", tc.internal).style("font-size", "11px");
        ag.selectAll("line").attr("stroke", tc.branch);
        ag.selectAll("path").attr("stroke", tc.branch);
    });

    g.append("text").attr("x", innerWidth / 2).attr("y", innerHeight + 38)
        .attr("text-anchor", "middle").attr("font-size", 12).attr("fill", tc.internal).text("k (number of clusters)");
    g.append("text").attr("transform", "rotate(-90)").attr("x", -innerHeight / 2).attr("y", -42)
        .attr("text-anchor", "middle").attr("font-size", 12).attr("fill", tc.internal).text("CalBow Score");

    // Line + points
    const lineColor = tc.accent;
    const line = d3.line()
        .x(d => mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2)
        .y(d => yScale(d.score)).curve(d3.curveMonotoneX);
    g.append("path").datum(dataPoints).attr("fill", "none").attr("stroke", lineColor).attr("stroke-width", 2).attr("d", line);

    g.selectAll(".score-point").data(dataPoints).enter().append("circle").attr("class", "score-point")
        .attr("cx", d => mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2)
        .attr("cy", d => yScale(d.score))
        .attr("r", 3.5).attr("fill", lineColor).attr("stroke", tc.bg).attr("stroke-width", 1.5)
        .on("mouseover", (event, d) => d3Tooltip.style("opacity", 1).html(`k = ${d.k}<br/>score = ${d.score.toFixed(4)}`))
        .on("mousemove", event => d3Tooltip.style("left", (event.pageX + 12) + "px").style("top", (event.pageY + 12) + "px"))
        .on("mouseout", () => d3Tooltip.style("opacity", 0));

    // Peaks
    g.selectAll(".peak-point").data(peakPoints).enter().append("path").attr("class", "peak-point")
        .attr("transform", d => `translate(${mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2},${yScale(d.score)})`)
        .attr("d", d3.symbol().type(d3.symbolDiamond).size(100))
        .attr("fill", "#ef4444").attr("stroke", tc.bg).attr("stroke-width", 1.5)
        .on("mouseover", (event, d) => d3Tooltip.style("opacity", 1).html(`<strong>Optimal k = ${d.k}</strong><br/>score = ${d.score.toFixed(4)}`))
        .on("mousemove", event => d3Tooltip.style("left", (event.pageX + 12) + "px").style("top", (event.pageY + 12) + "px"))
        .on("mouseout", () => d3Tooltip.style("opacity", 0));

    // Legend
    const legend = g.append("g").attr("transform", `translate(${innerWidth - 130}, 8)`);
    legend.append("rect").attr("x", -8).attr("y", -8).attr("width", 140).attr("height", 44)
        .attr("fill", tc.bg).attr("rx", 6).attr("stroke", tc.border).attr("opacity", 0.9);
    legend.append("line").attr("x1", 0).attr("y1", 6).attr("x2", 18).attr("y2", 6).attr("stroke", lineColor).attr("stroke-width", 2);
    legend.append("circle").attr("cx", 9).attr("cy", 6).attr("r", 3).attr("fill", lineColor);
    legend.append("text").attr("x", 24).attr("y", 10).attr("font-size", 11).attr("fill", tc.internal).text("Score");
    legend.append("path").attr("transform", "translate(9,26)").attr("d", d3.symbol().type(d3.symbolDiamond).size(80)).attr("fill", "#ef4444");
    legend.append("text").attr("x", 24).attr("y", 30).attr("font-size", 11).attr("fill", tc.internal).text("Optimal k");

    if (axisModeEl && !axisModeEl.__wired) {
        axisModeEl.__wired = true;
        axisModeEl.addEventListener('change', function () { axisModeEl.__userOverride = true; drawOptimalK(data); });
    }

    const okSvgNode = document.querySelector("#optimalk_plot svg");
    window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = okSvgNode ? new XMLSerializer().serializeToString(okSvgNode) : null;
}


/* ─────────────────────────────────────
   Mini Scores Panel (overlay on tree page)
   ───────────────────────────────────── */
function drawMiniScores(data) {
    var panel = document.getElementById("mini-scores-panel");
    var plotEl = document.getElementById("mini-scores-plot");
    if (!panel || !plotEl) return;

    var scores = data && data.scores ? data.scores : [];
    var peaks  = data && data.peaks  ? data.peaks  : [];
    if (!Array.isArray(scores) || scores.length === 0) { panel.classList.remove("visible"); return; }

    panel.classList.add("visible");
    plotEl.innerHTML = "";
    var tc = getThemeColors();

    var dataPoints = [];
    for (var i = 0; i < scores.length - 1; i++) dataPoints.push({ k: i + 2, score: scores[i + 1] });
    if (!dataPoints.length) return;

    var width = 252, height = 124;
    var margin = { top: 8, right: 8, bottom: 20, left: 36 };
    var innerW = width - margin.left - margin.right;
    var innerH = height - margin.top - margin.bottom;

    var svg = d3.select(plotEl).append("svg").attr("width", width).attr("height", height).style("background", "transparent");
    var g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var xScale = d3.scaleLinear().domain([d3.min(dataPoints, d => d.k), d3.max(dataPoints, d => d.k)]).range([0, innerW]);
    var yScale = d3.scaleLinear().domain([d3.min(dataPoints, d => d.score), d3.max(dataPoints, d => d.score)]).nice().range([innerH, 0]);

    // Axes (minimal)
    g.append("g").attr("transform", "translate(0," + innerH + ")").call(d3.axisBottom(xScale).ticks(4).tickFormat(d3.format("d")))
        .selectAll("text").attr("fill", tc.internal).style("font-size", "9px");
    g.append("g").call(d3.axisLeft(yScale).ticks(3))
        .selectAll("text").attr("fill", tc.internal).style("font-size", "9px");
    g.selectAll(".domain").attr("stroke", tc.branch);
    g.selectAll(".tick line").attr("stroke", tc.branch);

    // Line
    var line = d3.line().x(d => xScale(d.k)).y(d => yScale(d.score)).curve(d3.curveMonotoneX);
    g.append("path").datum(dataPoints).attr("fill", "none").attr("stroke", tc.accent).attr("stroke-width", 1.5).attr("d", line);

    // Clickable points
    g.selectAll(".mini-point").data(dataPoints).enter().append("circle")
        .attr("cx", d => xScale(d.k)).attr("cy", d => yScale(d.score))
        .attr("r", d => peaks.indexOf(d.k) >= 0 ? 4 : 2.5)
        .attr("fill", d => peaks.indexOf(d.k) >= 0 ? "#ef4444" : tc.accent)
        .attr("stroke", tc.bg).attr("stroke-width", 1)
        .style("cursor", "pointer")
        .on("click", function (event, d) { switchToK(d.k); })
        .on("mouseover", function (event, d) { d3Tooltip.style("opacity", 1).html("k=" + d.k + " score=" + d.score.toFixed(3)); })
        .on("mousemove", function (event) { d3Tooltip.style("left", (event.pageX + 10) + "px").style("top", (event.pageY + 10) + "px"); })
        .on("mouseout", function () { d3Tooltip.style("opacity", 0); });
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
    if (latestApiData.ks && latestApiData.clusters) {
        var idx2 = latestApiData.ks.indexOf(k);
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
    showToast("k=" + k + " not cached. Enable 'Compute all clusters' to access any k.", "info", 3000);
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
    var fileBadge  = document.getElementById("file-badge");
    if (fileNameEl && fileBadge) { fileBadge.textContent = file.name; fileNameEl.style.display = "flex"; }
    showStatus("Loading tree...", "info");
    var reader = new FileReader();
    reader.onload = function (e) { newickEl.value = (e.target.result || "").trim(); showStatus("Tree loaded", "success"); };
    reader.readAsText(file);
}


/* ─────────────────────────────────────
   Save / Export
   ───────────────────────────────────── */
function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a"); a.href = url; a.download = filename;
    document.body.appendChild(a); a.click();
    setTimeout(() => { document.body.removeChild(a); URL.revokeObjectURL(url); }, 100);
}

function exportSvgFromEl(selector, filename) {
    const svgNode = document.querySelector(selector + " svg");
    if (!svgNode) { showToast("No visualization to export.", "danger"); return; }
    downloadBlob(new Blob([new XMLSerializer().serializeToString(svgNode)], { type: "image/svg+xml" }), filename);
    showToast("Downloaded " + filename, "success", 2000);
}

async function exportTSV() {
    try {
        const res = await fetch("/api/export_tsv", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ top_n: 1, outlier: true }) });
        if (!res.ok) { const err = await res.json().catch(() => ({})); throw new Error(err.detail || "Export failed"); }
        downloadBlob(new Blob([await res.text()], { type: "text/tab-separated-values" }), "phytclust_results.tsv");
        showToast("Downloaded TSV", "success", 2000);
    } catch (e) { showToast("TSV export failed: " + e.message, "danger"); }
}

async function saveToServer() {
    var el = document.getElementById("output-dir");
    const dir = (el && el.value ? el.value : "results");
    try {
        await apiPostJson("/api/save", { results_dir: dir, top_n: 1, outlier: true });
        showToast("Saved to " + dir, "success", 3000);
    } catch (e) { showToast("Save failed: " + e.message, "danger"); }
}

function exportPngFromEl(selector, filename, dpi) {
    const svgNode = document.querySelector(selector + " svg");
    if (!svgNode) { showToast("No visualization to export.", "danger"); return; }
    dpi = dpi || 300;
    const scale = dpi / 96;  // browsers render at 96 DPI
    const svgData = new XMLSerializer().serializeToString(svgNode);
    const svgW = parseFloat(svgNode.getAttribute("width")) || svgNode.clientWidth || 800;
    const svgH = parseFloat(svgNode.getAttribute("height")) || svgNode.clientHeight || 600;
    const canvas = document.createElement("canvas");
    canvas.width = Math.ceil(svgW * scale);
    canvas.height = Math.ceil(svgH * scale);
    const ctx = canvas.getContext("2d");
    // Fill with white background for publication
    const bgColor = getComputedStyle(document.documentElement).getPropertyValue("--pc-bg").trim() || "#ffffff";
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
            if (!pngBlob) { showToast("PNG export failed.", "danger"); return; }
            downloadBlob(pngBlob, filename);
            showToast("Downloaded " + filename + " at " + dpi + " DPI", "success", 3000);
        }, "image/png");
    };
    img.onerror = function () { URL.revokeObjectURL(url); showToast("PNG rendering failed.", "danger"); };
    img.src = url;
}

function copySvgToClipboard(selector) {
    const svgNode = document.querySelector(selector + " svg");
    if (!svgNode) { showToast("Nothing to copy.", "danger"); return; }
    navigator.clipboard.writeText(new XMLSerializer().serializeToString(svgNode)).then(
        () => showToast("SVG copied to clipboard", "success", 2000),
        () => showToast("Copy failed", "danger")
    );
}


/* ─────────────────────────────────────
   Tab System
   ───────────────────────────────────── */
function switchTab(tabName) {
    document.querySelectorAll(".tab-btn").forEach(btn => btn.classList.toggle("active", btn.dataset.tab === tabName));
    document.querySelectorAll(".viz-panel").forEach(panel => panel.classList.toggle("active", panel.id === tabName));

    var treeToolbar    = document.getElementById("tree-toolbar");
    var optkToolbar    = document.getElementById("optk-toolbar");
    var compareToolbar = document.getElementById("compare-toolbar");
    if (treeToolbar)    treeToolbar.style.display    = (tabName === "viewer") ? "" : "none";
    if (optkToolbar)    optkToolbar.style.display    = (tabName === "viewer-optimal-k") ? "" : "none";
    if (compareToolbar) compareToolbar.style.display = (tabName === "compare") ? "" : "none";

    if (tabName === "viewer-optimal-k" && latestOptimalKData) drawOptimalK(latestOptimalKData);
    if (tabName === "compare") drawComparison();
}


/* ─────────────────────────────────────
   Collapsible section helper
   ───────────────────────────────────── */
function wireCollapse(toggleId, bodyId) {
    var toggle = document.getElementById(toggleId);
    var body   = document.getElementById(bodyId);
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

    // ── Theme toggle ──
    var btnTheme = document.getElementById("btn-theme");
    if (btnTheme) btnTheme.addEventListener("click", toggleTheme);

    // ── Sidebar toggle ──
    var sidebarToggle = document.getElementById("sidebar-toggle");
    var sidebar       = document.getElementById("app-sidebar");
    if (sidebarToggle && sidebar) {
        sidebarToggle.addEventListener("click", function () { sidebar.classList.toggle("collapsed"); });
    }

    // ── File input ──
    var fileInput = document.getElementById("file-input");
    if (fileInput) fileInput.addEventListener("change", handleFileSelect);

    // ── Drag and drop ──
    var uploadZone = document.getElementById("upload-zone");
    if (uploadZone) {
        uploadZone.addEventListener("dragover", function (e) { e.preventDefault(); e.stopPropagation(); uploadZone.classList.add("drag-over"); });
        uploadZone.addEventListener("dragleave", function (e) { e.preventDefault(); e.stopPropagation(); uploadZone.classList.remove("drag-over"); });
        uploadZone.addEventListener("drop", function (e) {
            e.preventDefault(); e.stopPropagation(); uploadZone.classList.remove("drag-over");
            if (e.dataTransfer.files && e.dataTransfer.files.length) loadFile(e.dataTransfer.files[0]);
        });
    }

    // ── Run button + keyboard shortcut ──
    var btnRun = document.getElementById("btn-run");
    if (btnRun) btnRun.addEventListener("click", function (e) { e.preventDefault(); runPhytClust(); });
    document.addEventListener("keydown", function (e) {
        if ((e.ctrlKey || e.metaKey) && e.key === "Enter") { e.preventDefault(); runPhytClust(); }
    });

    // ── Reset button ──
    var resetBtn = document.getElementById("btn-extra-reset");
    if (resetBtn) resetBtn.addEventListener("click", function (e) { e.preventDefault(); resetExtraParams(); });

    // ── Mode selector ──
    document.querySelectorAll("#mode-selector .mode-btn").forEach(function (btn) {
        btn.addEventListener("click", function () {
            document.querySelectorAll("#mode-selector .mode-btn").forEach(b => b.classList.remove("active"));
            btn.classList.add("active");
            var mode = btn.dataset.mode;
            var paramK    = document.getElementById("param-k");
            var paramTopN = document.getElementById("param-topn");
            var paramBins = document.getElementById("param-bins");
            if (paramK)    paramK.style.display    = (mode === "k") ? "" : "none";
            if (paramTopN) paramTopN.style.display  = (mode === "k") ? "none" : "";
            if (paramBins) paramBins.style.display  = (mode === "resolution") ? "" : "none";
            if (extraResolutionEl) extraResolutionEl.checked = (mode === "resolution");
        });
    });

    // ── Collapsible sections ──
    wireCollapse("tree-opts-toggle", "tree-opts-body");
    wireCollapse("outlier-opts-toggle", "outlier-opts-body");
    wireCollapse("support-opts-toggle", "support-opts-body");
    wireCollapse("peak-opts-toggle", "peak-opts-body");

    // ── Tabs ──
    document.querySelectorAll(".tab-btn").forEach(function (btn) {
        btn.addEventListener("click", function () { switchTab(btn.dataset.tab); });
    });

    // ── Help sidebar ──
    var btnHelp = document.getElementById("btn-help");
    var helpSidebar = document.getElementById("help-sidebar");
    if (btnHelp && helpSidebar) btnHelp.addEventListener("click", function (e) { e.stopPropagation(); helpSidebar.classList.toggle("open"); });
    var btnHelpClose = document.getElementById("btn-help-close");
    if (btnHelpClose && helpSidebar) btnHelpClose.addEventListener("click", function () { helpSidebar.classList.remove("open"); });
    // Close help when clicking outside
    if (helpSidebar) {
        document.addEventListener("click", function (e) {
            if (helpSidebar.classList.contains("open") && !helpSidebar.contains(e.target) && e.target !== btnHelp) {
                helpSidebar.classList.remove("open");
            }
        });
        helpSidebar.addEventListener("click", function (e) { e.stopPropagation(); });
    }

    // ── About modal ──
    var btnAbout = document.getElementById("btn-about");
    var aboutModal = document.getElementById("aboutModal");
    var aboutBg = document.getElementById("aboutBackdrop");
    function openAbout()  { if (aboutModal) aboutModal.classList.add("show"); if (aboutBg) aboutBg.classList.add("show"); }
    function closeAbout() { if (aboutModal) aboutModal.classList.remove("show"); if (aboutBg) aboutBg.classList.remove("show"); }
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
            e.preventDefault(); newickEl.value = EXAMPLE_NEWICK;
            if (helpSidebar) helpSidebar.classList.remove("open");
            runPhytClust();
        });
    }

    // ── Save dropdown ──
    var btnSave = document.getElementById("btn-save");
    var saveDropdown = document.getElementById("save-dropdown");
    if (btnSave && saveDropdown) {
        btnSave.addEventListener("click", function (e) { e.stopPropagation(); saveDropdown.classList.toggle("show"); });
        document.addEventListener("click", function () { saveDropdown.classList.remove("show"); });
        saveDropdown.addEventListener("click", function (e) { e.stopPropagation(); });
        saveDropdown.querySelectorAll(".dropdown-item").forEach(function (item) {
            item.addEventListener("click", function () {
                saveDropdown.classList.remove("show");
                var type = item.dataset.saveType;
                if (type === "tsv") exportTSV();
                if (type === "tree_plot") exportSvgFromEl("#tree_display", "phytclust_tree.svg");
                if (type === "tree_png") {
                    var dpiStr = prompt("Enter DPI for publication-ready PNG (default 300):", "300");
                    if (dpiStr === null) return;
                    var dpi = parseInt(dpiStr, 10);
                    if (isNaN(dpi) || dpi < 72) dpi = 300;
                    exportPngFromEl("#tree_display", "phytclust_tree.png", dpi);
                }
                if (type === "k_plot") exportSvgFromEl("#optimalk_plot", "phytclust_scores.svg");
                if (type === "save_server") saveToServer();
            });
        });
    }

    // ── Copy buttons ──
    var btnCopyTree = document.getElementById("btn-copy-tree");
    var btnCopyMaxk = document.getElementById("btn-copy-maxk");
    if (btnCopyTree) btnCopyTree.addEventListener("click", function () { copySvgToClipboard("#tree_display"); });
    if (btnCopyMaxk) btnCopyMaxk.addEventListener("click", function () { copySvgToClipboard("#optimalk_plot"); });

    // ── Cluster selector ──
    var clusterSelect = document.getElementById("cluster-select");
    if (clusterSelect) clusterSelect.addEventListener("change", function () { switchCluster(parseInt(this.value, 10)); });
    var clusterToggle = document.getElementById("cluster-view-toggle");
    if (clusterToggle) clusterToggle.addEventListener("click", toggleClusterViewMode);

    // ── Cluster prev/next buttons ──
    var clusterPrev = document.getElementById("cluster-prev");
    var clusterNext = document.getElementById("cluster-next");
    if (clusterPrev) clusterPrev.addEventListener("click", function () { cycleCluster(-1); });
    if (clusterNext) clusterNext.addEventListener("click", function () { cycleCluster(1); });

    // ── Arrow key cycling through clusters ──
    document.addEventListener("keydown", function (e) {
        if (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA" || e.target.tagName === "SELECT") return;
        if (e.key === "ArrowLeft") { e.preventDefault(); cycleCluster(-1); }
        if (e.key === "ArrowRight") { e.preventDefault(); cycleCluster(1); }
    });

    // ── Compare selectors ──
    var compareA = document.getElementById("compare-k-a");
    var compareB = document.getElementById("compare-k-b");
    if (compareA) compareA.addEventListener("change", drawComparison);
    if (compareB) compareB.addEventListener("change", drawComparison);

    // ── Layout mode ──
    var layoutSelect = document.getElementById("layout-mode");
    if (layoutSelect) layoutSelect.addEventListener("change", function () { CURRENT_LAYOUT_MODE = this.value; drawTree(); });

    // ── Color mode ──
    var colorModeSelect = document.getElementById("color-mode");
    if (colorModeSelect) colorModeSelect.addEventListener("change", function () { COLOR_MODE = this.value || "bars"; drawTree(); });

    // ── Palette selector ──
    let CUSTOM_PALETTE = null;
    function normalizeHexColor(s) { const v = (s || "").trim(); if (/^#[0-9a-fA-F]{6}$/.test(v)) return v; if (/^[0-9a-fA-F]{6}$/.test(v)) return "#" + v; return null; }
    function updatePalettePreview(colors) { const box = document.getElementById("palette-preview"); if (box) box.innerHTML = colors.map(c => `<span title="${c}" style="background:${c};"></span>`).join(""); }

    function applyPalette(type) {
        const palettes = {
            default: ["#b84b4b", "#849060", "#3d7c74", "#6e3f8a", "#ceb94b", "#3f648a", "#3f408a", "#da63aa"],
            pastel:  ["#ffb3ba", "#ffdfba", "#ffffba", "#baffc9", "#bae1ff", "#d7baff", "#ffcce6", "#c2f0c2"],
            vivid:   ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"],
            dark:    ["#4e79a7", "#59a14f", "#e15759", "#b07aa1", "#edc948", "#76b7b2", "#ff9da7", "#9c755f"]
        };
        if (type === "custom" && CUSTOM_PALETTE && CUSTOM_PALETTE.length) BASE_COLORS.splice(0, BASE_COLORS.length, ...CUSTOM_PALETTE);
        else if (palettes[type]) BASE_COLORS.splice(0, BASE_COLORS.length, ...palettes[type]);
        updatePalettePreview(BASE_COLORS);
        if (Object.keys(CURRENT_CLUSTERS || {}).length > 0) {
            CLUSTER_COLORS = generateClusterColors(Math.max(...Object.values(CURRENT_CLUSTERS)) + 1);
            drawTree();
        }
    }

    var paletteSelect = document.getElementById("palette-select");
    if (paletteSelect) {
        paletteSelect.addEventListener("change", function () {
            if (this.value === "custom") {
                var initial = (CUSTOM_PALETTE && CUSTOM_PALETTE.length) ? CUSTOM_PALETTE.join(",") : BASE_COLORS.join(",");
                var raw = prompt("Enter hex colors separated by commas:", initial);
                if (raw == null) { this.value = "default"; applyPalette("default"); return; }
                var parsed = raw.split(",").map(normalizeHexColor).filter(Boolean);
                if (parsed.length < 2) { alert("Please provide at least 2 valid hex colors."); this.value = "default"; applyPalette("default"); return; }
                CUSTOM_PALETTE = parsed;
                applyPalette("custom"); return;
            }
            applyPalette(this.value);
        });
    }
    updatePalettePreview(BASE_COLORS);

    // ── Leaf/internal name toggles ──
    var internalCb = document.getElementById("show-internal-names");
    if (internalCb) { internalCb.checked = SHOW_INTERNAL_NAMES; internalCb.addEventListener("change", function () { SHOW_INTERNAL_NAMES = this.checked; drawTree(); }); }
    var leafCb = document.getElementById("show-leaf-names");
    if (leafCb) { leafCb.checked = SHOW_LEAF_NAMES; leafCb.addEventListener("change", function () { SHOW_LEAF_NAMES = this.checked; drawTree(); }); }

    var boxLabelCb = document.getElementById("show-box-labels");
    if (boxLabelCb) { boxLabelCb.checked = SHOW_BOX_LABELS; boxLabelCb.addEventListener("change", function () { SHOW_BOX_LABELS = this.checked; drawTree(); }); }

    // ── Node/label sizes ──
    var nodeSizeInput = document.getElementById("node-size");
    if (nodeSizeInput) { nodeSizeInput.value = LEAF_NODE_RADIUS; nodeSizeInput.addEventListener("input", function () { var v = parseFloat(this.value); if (!isNaN(v) && v >= 1 && v <= 10) { LEAF_NODE_RADIUS = v; drawTree(); } }); }
    var labelSizeInput = document.getElementById("label-size");
    if (labelSizeInput) { labelSizeInput.value = LABEL_FONT_SIZE; labelSizeInput.addEventListener("input", function () { var v = parseFloat(this.value); if (!isNaN(v) && v >= 6 && v <= 18) { LABEL_FONT_SIZE = v; drawTree(); } }); }

    // ── Width/Height sliders ──
    var widthSlider  = document.getElementById("tree-width-scale");
    var heightSlider = document.getElementById("tree-height-scale");
    if (widthSlider) widthSlider.addEventListener("input", function () { TREE_WIDTH_SCALE = parseFloat(this.value); if (!isNaN(TREE_WIDTH_SCALE)) drawTree(); });
    if (heightSlider) heightSlider.addEventListener("input", function () { TREE_HEIGHT_SCALE = parseFloat(this.value); if (!isNaN(TREE_HEIGHT_SCALE)) drawTree(); });

    // ── Reset view ──
    var btnResetView = document.getElementById("btn-reset-view");
    if (btnResetView) {
        btnResetView.addEventListener("click", function (e) {
            e.preventDefault();
            if (LAST_TREE_SVG && LAST_TREE_ZOOM) LAST_TREE_SVG.transition().duration(150).call(LAST_TREE_ZOOM.transform, d3.zoomIdentity);
            if (widthSlider) widthSlider.value = "1.0";
            if (heightSlider) heightSlider.value = "1.0";
            TREE_WIDTH_SCALE = 1.0; TREE_HEIGHT_SCALE = 1.0;
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
            if (el) el.addEventListener("click", function () { fn(); hideNodeContextMenu(); });
        }
        btn("ctx-rename", function () {
            if (!CTX_TARGET_DATA) return;
            var c = getNodeCustom(CTX_TARGET_DATA);
            var val = prompt("Rename node:", nodeDisplayName(CTX_TARGET_DATA));
            if (val != null) { c.renamedTo = val; drawTree(); }
        });
        btn("ctx-highlight", function () {
            if (!CTX_TARGET_DATA) return;
            var c = getNodeCustom(CTX_TARGET_DATA);
            c.highlighted = !c.highlighted; drawTree();
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
            ctxSlider.addEventListener("click", function (e) { e.stopPropagation(); });
        }
        btn("ctx-copy-subtree", function () {
            if (!CTX_TARGET_DATA) return;
            var names = [];
            collectLeafNamesData(CTX_TARGET_DATA, names);
            navigator.clipboard.writeText(names.join("\n")).then(
                () => showToast(names.length + " leaf names copied", "success", 2000),
                () => showToast("Copy failed", "danger")
            );
        });
        btn("ctx-select-cluster", function () {
            if (!CTX_TARGET_DATA || !CURRENT_CLUSTERS) return;
            var name = CTX_TARGET_DATA.name;
            if (!name) { var leaves = []; collectLeafNamesData(CTX_TARGET_DATA, leaves); name = leaves[0]; }
            var cid = name ? CURRENT_CLUSTERS[name] : null;
            if (cid == null) { showToast("No cluster for this node", "info"); return; }
            var members = Object.entries(CURRENT_CLUSTERS).filter(([k, v]) => v === cid).map(([k]) => k);
            navigator.clipboard.writeText(members.join("\n")).then(
                () => showToast(members.length + " members of cluster " + cid + " copied", "success", 2000),
                () => showToast("Copy failed", "danger")
            );
        });
        btn("ctx-rename-cluster", function () {
            if (!CTX_TARGET_DATA || !CURRENT_CLUSTERS) return;
            var name = CTX_TARGET_DATA.name;
            if (!name) { var leaves = []; collectLeafNamesData(CTX_TARGET_DATA, leaves); name = leaves[0]; }
            var cid = name ? CURRENT_CLUSTERS[name] : null;
            if (cid == null) { showToast("No cluster for this node", "info"); return; }
            var current = BOX_LABEL_MAP[cid] || ("C" + cid);
            var val = prompt("Label for cluster " + cid + ":", current);
            if (val != null) { BOX_LABEL_MAP[cid] = val; drawTree(); }
        });
        btn("ctx-reset", function () {
            if (!CTX_TARGET_DATA) return;
            NODE_CUSTOM.delete(CTX_TARGET_DATA); drawTree();
        });
    })();

    // ── Mini scores panel toggle ──
    var miniScoresToggle = document.getElementById("mini-scores-toggle");
    var miniScoresBody = document.getElementById("mini-scores-body");
    if (miniScoresToggle && miniScoresBody) {
        miniScoresToggle.addEventListener("click", function (e) {
            e.stopPropagation();
            miniScoresBody.classList.toggle("collapsed");
            miniScoresToggle.innerHTML = miniScoresBody.classList.contains("collapsed") ? "&#x25B2;" : "&#x25BC;";
        });
    }

    // ── Prefill example and auto-run ──
    if (newickEl) newickEl.value = EXAMPLE_NEWICK;
    runPhytClust();
});
