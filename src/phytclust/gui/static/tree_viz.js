/* ============================================================
   PhytClust – tree_viz.js (2026 Modern Rewrite)
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

let isRunning         = false;
let latestOptimalKData = null;
let latestApiData      = null;
let LAST_TREE_SVG      = null;
let LAST_TREE_ZOOM     = null;

// Per-node customizations: keyed by data node reference
// { highlighted: bool, radiusScale: number, renamedTo: string }
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
    // Also show toast for important messages
    if (type === "danger") {
        showToast(message, "danger", 6000);
    } else if (type === "success") {
        showToast(message, "success", 3000);
    }
}


/* ─────────────────────────────────────
   Utility Functions
   ───────────────────────────────────── */
function escapeRegExp(str) {
    return str.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

function isOutgroupInNewick(newick, outgroup) {
    if (!outgroup) return true;
    var name = escapeRegExp(outgroup);
    var quoted   = new RegExp("[\\(,]\\s*'" + name + "'\\s*(?:[:,\\)])");
    var unquoted = new RegExp("[\\(,]\\s*"  + name + "\\s*(?:[:,\\)])");
    return quoted.test(newick) || unquoted.test(newick);
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
    if (extraKEl)          extraKEl.value = "";
    if (extraOutgroupEl)   extraOutgroupEl.value = "";
    if (extraRootTaxonEl)  extraRootTaxonEl.value = "";
    if (extraTopNEl)       extraTopNEl.value = "";
    if (extraResolutionEl) extraResolutionEl.checked = false;
    if (extraBinsEl)       extraBinsEl.value = "";
    if (extraMaxKEl)       extraMaxKEl.value = "";
    if (extraMaxKLimitEl)  extraMaxKLimitEl.value = "";
    if (extraLambdaEl)     extraLambdaEl.value = "";
    if (extraMinClusterEl) extraMinClusterEl.value = "";
    if (extraOutlierEl)    extraOutlierEl.checked = true;
    var computeAllEl = document.getElementById("extra-compute-all");
    if (computeAllEl) computeAllEl.checked = false;
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
    let r = (num >> 16)       + factor * 255;
    let g = ((num >> 8) & 0xFF) + factor * 255;
    let b = (num & 0xFF)       + factor * 255;
    r = Math.min(255, Math.max(0, r));
    g = Math.min(255, Math.max(0, g));
    b = Math.min(255, Math.max(0, b));
    return `rgb(${r | 0}, ${g | 0}, ${b | 0})`;
}

function withAlpha(color, alpha) {
    if (color.startsWith("rgb")) {
        return color.replace("rgb", "rgba").replace(")", `, ${alpha})`);
    }
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
   Tree Layout Computation
   ───────────────────────────────────── */
function accumulateBranchLength(node, length) {
    length = length || 0;
    node._bl = length;
    if (node.children) {
        for (const c of node.children) {
            accumulateBranchLength(c, length + (c.length || 0));
        }
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

    // Cartesian
    d3.cluster().size([innerH, 1])(HIER_CART);
    HIER_CART.each(d => { d._x = d.x; d._y = blToX(d.data._bl || 0); });

    // Circular
    d3.cluster().size([2 * Math.PI, 1])(HIER_CIRC);
    HIER_CIRC.each(d => { d._angle = d.x; d._radius = blToR(d.data._bl || 0); });
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

function drawTree() {
    clearTree();
    if (!NEWICK_RAW_TREE) return;

    computeLayouts();

    const layoutMode = CURRENT_LAYOUT_MODE || "rectangular";
    const colorMode  = COLOR_MODE || "bars";

    const container = d3.select("#tree_display");
    const width  = treeHost.clientWidth  || 800;
    const height = treeHost.clientHeight || 500;
    const margin = { top: 20, right: 80, bottom: 20, left: 80 };

    const svg = container.append("svg")
        .attr("width", width)
        .attr("height", height)
        .style("background", "transparent");

    const zoomLayer = svg.append("g");
    const g = zoomLayer.append("g")
        .attr("transform",
            layoutMode === "circular"
                ? `translate(${width / 2},${height / 2})`
                : `translate(${margin.left},${margin.top})`
        );

    const zoom = d3.zoom()
        .scaleExtent([0.3, 8])
        .on("zoom", (event) => { zoomLayer.attr("transform", event.transform); });

    svg.call(zoom);
    LAST_TREE_SVG  = svg;
    LAST_TREE_ZOOM = zoom;

    function hasOriginalChildren(d) {
        return !!(d && d.data && Array.isArray(d.data.children) && d.data.children.length > 0);
    }
    function canCollapse(d)   { return hasOriginalChildren(d); }
    function isVisibleLeaf(d) { return !(d && d.children && d.children.length); }
    function nodeRadius(d)    { return customNodeRadius(d); }
    function isHighlighted(d) {
        return NODE_CUSTOM.has(d.data) && NODE_CUSTOM.get(d.data).highlighted;
    }

    function nodeFill(d) {
        const name = d.data && d.data.name;
        if (colorMode === "bars") return "#334155";

        if (name && !hasOriginalChildren(d)) {
            const cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[name] : null;
            if (cid == null || !CLUSTER_COLORS.length) return "#334155";
            return CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
        }
        const rep = representativeClusterIdFromData(d.data);
        if (rep != null && CLUSTER_COLORS.length) return CLUSTER_COLORS[rep % CLUSTER_COLORS.length];
        return "#334155";
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
            const delta = ta - sa;
            const sweep = delta >= 0 ? 1 : 0;
            const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
            const arcR = rScale(sr);
            return `M${p0[0]},${p0[1]}A${arcR},${arcR} 0 ${largeArc},${sweep} ${p1[0]},${p1[1]}L${p2[0]},${p2[1]}`;
        }

        // Links
        g.append("g").selectAll(".tree-link")
            .data(root.links()).enter().append("path")
            .attr("class", "tree-link").attr("fill", "none")
            .attr("stroke", "#64748b").attr("stroke-opacity", 0.6).attr("stroke-width", 1.1)
            .attr("d", radialElbowPath);

        // Nodes
        const node = g.append("g").selectAll(".tree-node")
            .data(root.descendants()).enter().append("g")
            .attr("class", "tree-node")
            .attr("transform", d => { const p = nodeXY(d); return `translate(${p[0]},${p[1]})`; });

        node.append("circle")
            .attr("r", nodeRadius).attr("fill", nodeFill)
            .attr("stroke", d => isHighlighted(d) ? "#f59e0b" : "none")
            .attr("stroke-width", d => isHighlighted(d) ? 2.5 : 0)
            .style("cursor", "pointer")
            .on("click", function (event, d) {
                if (!canCollapse(d)) return;
                d.data._collapsed = !d.data._collapsed;
                drawTree();
            })
            .on("contextmenu", showNodeContextMenu);

        // Leaf labels — radially outward
        node.filter(d => isVisibleLeaf(d))
            .append("text")
            .attr("dy", "0.32em").attr("font-size", LABEL_FONT_SIZE + "px")
            .attr("fill", d => isHighlighted(d) ? "#f59e0b" : "#334155")
            .attr("font-weight", d => isHighlighted(d) ? "bold" : "normal")
            .attr("transform", function (d) {
                const a = d._angle || 0;
                let deg = (a * 180 / Math.PI) - 90;
                const left = (deg > 90 || deg < -90);
                const r = nodeRadius(d) + LABEL_PAD;
                const flip = left ? 180 : 0;
                return `rotate(${deg}) translate(${r},0) rotate(${flip})`;
            })
            .style("text-anchor", function (d) {
                const deg = ((d._angle || 0) * 180 / Math.PI) - 90;
                return (deg > 90 || deg < -90) ? "end" : "start";
            })
            .text(d => SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : "");

        // Internal labels — offset perpendicular (above the node), smaller
        node.filter(d => !isVisibleLeaf(d))
            .append("text")
            .attr("dy", "0.32em").attr("font-size", (LABEL_FONT_SIZE - 1) + "px")
            .attr("fill", d => isHighlighted(d) ? "#f59e0b" : "#64748b")
            .attr("transform", function (d) {
                const a = d._angle || 0;
                let deg = (a * 180 / Math.PI) - 90;
                const offset = -(nodeRadius(d) + 5);
                return `rotate(${deg}) translate(0,${offset})`;
            })
            .style("text-anchor", "middle")
            .text(d => SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : "");

        // Circular side bars
        if (colorMode === "bars" && CURRENT_CLUSTERS && CLUSTER_COLORS.length) {
            const visLeaves = root.descendants()
                .filter(d => isVisibleLeaf(d))
                .sort((a, b) => (a._angle || 0) - (b._angle || 0));

            if (visLeaves.length >= 2) {
                const angles = visLeaves.map(d => d._angle || 0);
                const boundaries = new Array(visLeaves.length + 1);
                for (let i = 1; i < visLeaves.length; i++) {
                    boundaries[i] = (angles[i - 1] + angles[i]) / 2;
                }
                boundaries[0] = angles[0] - (boundaries[1] - angles[0]);
                boundaries[visLeaves.length] = angles[visLeaves.length - 1] + (angles[visLeaves.length - 1] - boundaries[visLeaves.length - 1]);

                function visLeafCid(d) {
                    const nm = d.data && d.data.name;
                    const isOrigLeaf = nm && !(d.data.children && d.data.children.length);
                    if (isOrigLeaf) return CURRENT_CLUSTERS[nm];
                    return representativeClusterIdFromData(d.data);
                }

                // Push ring past leaf labels
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
                            g.append("path")
                                .attr("d", arc.startAngle(boundaries[runStart]).endAngle(boundaries[i]))
                                .attr("fill", CLUSTER_COLORS[runCid % CLUSTER_COLORS.length])
                                .attr("opacity", 0.75);
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
    // Estimate max leaf label width to avoid bar overlap
    const leafNames = allNodes
        .filter(d => !d.children || !d.children.length)
        .map(d => (d.data.name || "").length);
    const maxLabelChars = d3.max(leafNames) || 0;
    const estLabelWidth = SHOW_LEAF_NAMES ? maxLabelChars * LABEL_FONT_SIZE * 0.6 + LABEL_PAD + LEAF_NODE_RADIUS : 20;
    const labelColumnX = (maxYCart * TREE_WIDTH_SCALE) + estLabelWidth + 12;

    // Links
    g.append("g").selectAll(".tree-link")
        .data(root.links()).enter().append("path")
        .attr("class", "tree-link").attr("fill", "none")
        .attr("stroke", "#64748b").attr("stroke-opacity", 0.6).attr("stroke-width", 1.2)
        .attr("d", function (d) {
            if (layoutMode === "rectangular") {
                return "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
                       "V" + (d.target._x * TREE_HEIGHT_SCALE) +
                       "H" + (d.target._y * TREE_WIDTH_SCALE);
            }
            return "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
                   "L" + (d.target._y * TREE_WIDTH_SCALE) + "," + (d.target._x * TREE_HEIGHT_SCALE);
        });

    // Nodes
    const node = g.append("g").selectAll(".tree-node")
        .data(root.descendants()).enter().append("g")
        .attr("class", "tree-node")
        .attr("transform", d => `translate(${(d._y || 0) * TREE_WIDTH_SCALE},${(d._x || 0) * TREE_HEIGHT_SCALE})`);

    node.append("circle")
        .attr("r", nodeRadius).attr("fill", nodeFill)
        .attr("stroke", d => isHighlighted(d) ? "#f59e0b" : "none")
        .attr("stroke-width", d => isHighlighted(d) ? 2.5 : 0)
        .style("cursor", "pointer")
        .on("click", function (event, d) {
            if (!canCollapse(d)) return;
            d.data._collapsed = !d.data._collapsed;
            drawTree();
        })
        .on("contextmenu", showNodeContextMenu);

    node.append("text")
        .attr("dy", d => (d.children && d.children.length) ? -(nodeRadius(d) + 4) : 3)
        .attr("x", d => {
            const r = nodeRadius(d);
            return (d.children && d.children.length) ? 0 : (r + LABEL_PAD);
        })
        .style("text-anchor", d => (d.children && d.children.length) ? "middle" : "start")
        .style("font-size", d => (d.children && d.children.length) ? (LABEL_FONT_SIZE - 1) + "px" : LABEL_FONT_SIZE + "px")
        .attr("fill", function (d) {
            if (isHighlighted(d)) return "#f59e0b";
            return (d.children && d.children.length) ? "#64748b" : "#334155";
        })
        .attr("font-weight", d => isHighlighted(d) ? "bold" : "normal")
        .text(function (d) {
            const isLeaf = !d.children || d.children.length === 0;
            if (isLeaf && SHOW_LEAF_NAMES)      return nodeDisplayName(d.data);
            if (!isLeaf && SHOW_INTERNAL_NAMES)  return nodeDisplayName(d.data);
            return "";
        });

    // Cluster side bars
    if (colorMode === "bars") {
        const leafNodes = allNodes
            .filter(d => !d.children || !d.children.length)
            .sort((a, b) => (a._x || 0) - (b._x || 0));

        const barX = labelColumnX;
        const barW = 20;
        const ys = leafNodes.map(d => (d._x || 0) * TREE_HEIGHT_SCALE);

        if (leafNodes.length === 1) {
            const nm = leafNodes[0].data && leafNodes[0].data.name;
            const cid = nm && CURRENT_CLUSTERS ? CURRENT_CLUSTERS[nm] : null;
            if (cid != null && CLUSTER_COLORS.length) {
                g.append("rect").attr("x", barX).attr("y", ys[0] - 6)
                    .attr("width", barW).attr("height", 12).attr("rx", 3)
                    .attr("fill", CLUSTER_COLORS[cid % CLUSTER_COLORS.length]).attr("opacity", 0.75);
            }
        } else if (leafNodes.length > 1) {
            const boundaries = new Array(leafNodes.length + 1);
            for (let i = 1; i < leafNodes.length; i++) {
                boundaries[i] = (ys[i - 1] + ys[i]) / 2;
            }
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
                        barsG.append("rect")
                            .attr("x", barX).attr("y", y0)
                            .attr("width", barW).attr("height", Math.max(1, y1 - y0))
                            .attr("rx", 3)
                            .attr("fill", CLUSTER_COLORS[runCid % CLUSTER_COLORS.length])
                            .attr("opacity", 0.75);
                    }
                    runStart = i;
                    runCid = (i < leafNodes.length) ? cid : null;
                }
            }
        }
    }

    // Branch-length axis
    const maxX  = d3.max(allNodes, d => d._x || 0) || 0;
    const maxY  = d3.max(allNodes, d => d._y || 0) || 0;
    const maxBl = d3.max(allNodes, d => d.data._bl || 0) || 1;
    const axisY = (maxX * TREE_HEIGHT_SCALE) + 20;
    const blScale = d3.scaleLinear().domain([0, maxBl]).range([0, maxY * TREE_WIDTH_SCALE]);

    const axisG = g.append("g")
        .attr("class", "branch-length-axis")
        .attr("transform", `translate(0, ${axisY})`)
        .call(d3.axisBottom(blScale).ticks(5));

    axisG.selectAll("text").attr("fill", "#64748b").style("font-size", "10px");
    axisG.selectAll("line").attr("stroke", "#94a3b8");
    axisG.selectAll("path").attr("stroke", "#94a3b8");

    axisG.append("text")
        .attr("x", (maxY * TREE_WIDTH_SCALE) / 2)
        .attr("y", 30).attr("text-anchor", "middle")
        .attr("font-size", 10).attr("fill", "#64748b")
        .text("Branch length");

    svg.attr("height", Math.max(height, axisY + 50 + margin.bottom));
}


/* ─────────────────────────────────────
   Node Context Menu
   ───────────────────────────────────── */
let CTX_TARGET_DATA = null;

function getNodeCustom(dataNode) {
    if (!NODE_CUSTOM.has(dataNode)) {
        NODE_CUSTOM.set(dataNode, { highlighted: false, radiusScale: 1, renamedTo: null });
    }
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
}

function hideNodeContextMenu() {
    const menu = document.getElementById("node-context-menu");
    if (menu) menu.style.display = "none";
    CTX_TARGET_DATA = null;
}

document.addEventListener("click", hideNodeContextMenu);
document.addEventListener("contextmenu", function (e) {
    // Hide if clicking outside tree nodes
    if (!e.target.closest(".tree-node")) hideNodeContextMenu();
});

(function initContextMenuButtons() {
    function btn(id, fn) {
        var el = document.getElementById(id);
        if (el) el.addEventListener("click", function () { fn(); hideNodeContextMenu(); });
    }
    btn("ctx-rename", function () {
        if (!CTX_TARGET_DATA) return;
        var c = getNodeCustom(CTX_TARGET_DATA);
        var current = nodeDisplayName(CTX_TARGET_DATA);
        var val = prompt("Rename node:", current);
        if (val != null) { c.renamedTo = val; drawTree(); }
    });
    btn("ctx-highlight", function () {
        if (!CTX_TARGET_DATA) return;
        var c = getNodeCustom(CTX_TARGET_DATA);
        c.highlighted = !c.highlighted;
        drawTree();
    });
    btn("ctx-grow", function () {
        if (!CTX_TARGET_DATA) return;
        var c = getNodeCustom(CTX_TARGET_DATA);
        c.radiusScale = Math.min(c.radiusScale + 0.5, 5);
        drawTree();
    });
    btn("ctx-shrink", function () {
        if (!CTX_TARGET_DATA) return;
        var c = getNodeCustom(CTX_TARGET_DATA);
        c.radiusScale = Math.max(c.radiusScale - 0.5, 0.5);
        drawTree();
    });
    btn("ctx-reset", function () {
        if (!CTX_TARGET_DATA) return;
        NODE_CUSTOM.delete(CTX_TARGET_DATA);
        drawTree();
    });
})();


/* ─────────────────────────────────────
   API Call
   ───────────────────────────────────── */
async function apiPostJson(url, payload) {
    const res = await fetch(url, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload)
    });
    const text = await res.text();
    let data = null;
    try { data = text ? JSON.parse(text) : null; } catch { /* ignore */ }
    if (!res.ok) {
        const msg = (data && data.detail) ? data.detail : `Request failed (${res.status})`;
        throw new Error(msg);
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

    // Read param values
    var kVal = null;
    if (extraKEl) {
        var kRaw = (extraKEl.value || "").trim();
        if (kRaw) { var p = parseInt(kRaw, 10); if (!isNaN(p) && p >= 1) kVal = p; }
    }

    var outgroupVal = extraOutgroupEl ? (extraOutgroupEl.value || "").trim() || null : null;
    var rootTaxonVal = extraRootTaxonEl ? (extraRootTaxonEl.value || "").trim() || null : null;
    var topNVal = null;
    if (extraTopNEl) {
        var tnRaw = (extraTopNEl.value || "").trim();
        if (tnRaw) { var tni = parseInt(tnRaw, 10); if (!isNaN(tni) && tni >= 1) topNVal = tni; }
    }
    var binsVal = null;
    if (extraBinsEl) {
        var bRaw = (extraBinsEl.value || "").trim();
        if (bRaw) { var bi = parseInt(bRaw, 10); if (!isNaN(bi) && bi >= 1) binsVal = bi; }
    }
    var maxKVal = null;
    if (extraMaxKEl) {
        var mkRaw = (extraMaxKEl.value || "").trim();
        if (mkRaw) {
            var mki = parseInt(mkRaw, 10);
            if (!isNaN(mki)) {
                if (mki < 4) { mki = 4; extraMaxKEl.value = String(mki); }
                maxKVal = mki;
            }
        }
    }
    var maxKLimitVal = null;
    if (extraMaxKLimitEl) {
        var mklRaw = (extraMaxKLimitEl.value || "").trim();
        if (mklRaw) {
            var mkli = parseFloat(mklRaw);
            if (!isNaN(mkli)) {
                var minLimit = (numSamples > 0) ? (4 / numSamples) : 0;
                if (mkli < minLimit) { mkli = minLimit; extraMaxKLimitEl.value = String(mkli); }
                maxKLimitVal = mkli;
            }
        }
    }
    var lambdaVal = null;
    if (extraLambdaEl) {
        var lRaw = (extraLambdaEl.value || "").trim();
        if (lRaw) { var li = parseFloat(lRaw); if (!isNaN(li)) lambdaVal = li; }
    }
    var minClusterVal = null;
    if (extraMinClusterEl) {
        var mcRaw = (extraMinClusterEl.value || "").trim();
        if (mcRaw) { var mci = parseInt(mcRaw, 10); if (!isNaN(mci) && mci >= 1) minClusterVal = mci; }
    }

    if (outgroupVal && !isOutgroupInNewick(newickText, outgroupVal)) {
        showStatus("Outgroup not found in Newick.", "danger");
        return;
    }

    // Sync hidden resolution checkbox
    if (extraResolutionEl) extraResolutionEl.checked = (mode === "resolution");

    // Build payload
    var payload = { newick: newickText, mode: mode };

    if (mode === "k") {
        if (kVal === null) { showToast("Please enter a value for k.", "danger"); return; }
        payload.k = kVal;
    }
    if (outgroupVal)      payload.outgroup = outgroupVal;
    if (rootTaxonVal)     payload.root_taxon = rootTaxonVal;
    if (topNVal !== null)  payload.top_n = topNVal;
    if (binsVal !== null)  payload.num_bins = binsVal;
    if (maxKVal !== null)  payload.max_k = maxKVal;
    if (maxKLimitVal !== null) payload.max_k_limit = maxKLimitVal;
    if (lambdaVal !== null) payload.lambda_weight = lambdaVal;
    if (minClusterVal !== null) payload.min_cluster_size = minClusterVal;

    var computeAllEl = document.getElementById("extra-compute-all");
    if (computeAllEl && computeAllEl.checked) payload.compute_all_clusters = true;

    if (mode === "resolution") payload.by_resolution = true;

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
        showStatus(`Finished in ${dt.toFixed(2)}s`, "success");
        resultEl.textContent = JSON.stringify(data, null, 2);

        if (data.newick) {
            NEWICK_RAW_TREE = parseNewick(data.newick);
            accumulateBranchLength(NEWICK_RAW_TREE);
            computeLayouts();

            // Handle cluster selector for multiple results
            populateClusterSelector(data);

            const clusterIdx = 0;
            const clusterMap = (data.clusters && data.clusters.length > clusterIdx)
                ? data.clusters[clusterIdx] : {};

            CURRENT_CLUSTERS = clusterMap;

            if (Object.keys(clusterMap).length > 0) {
                const maxCluster = Math.max(...Object.values(clusterMap));
                CLUSTER_COLORS = generateClusterColors(maxCluster + 1);
            } else {
                CLUSTER_COLORS = [];
                showStatus("No clusters found.", "danger");
            }
            drawTree();
        } else {
            clearTree();
            showStatus("No Newick tree returned by API.", "danger");
        }

        // Draw optimal k plot
        try {
            if (data && data.scores) {
                var plotHost = document.getElementById('optimalk_plot');
                if (plotHost) plotHost.innerHTML = '';
                latestOptimalKData = data;
                drawOptimalK(latestOptimalKData);
            } else {
                var plotHost2 = document.getElementById('optimalk_plot');
                if (plotHost2) {
                    plotHost2.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);font-size:14px;">Score plot unavailable in fixed k mode.</div>';
                }
            }
        } catch (plotErr) {
            console.warn('Plot rendering error:', plotErr);
        }

    } catch (e) {
        console.error(e);
        showStatus("Error: " + e.message, "danger");
        resultEl.textContent = "Error: " + e;
        latestOptimalKData = null;
        latestApiData = null;
        var plotHost3 = document.getElementById("optimalk_plot");
        if (plotHost3) plotHost3.innerHTML = "";
    } finally {
        isRunning = false;
        if (runBtn) {
            runBtn.disabled = false;
            runBtn.innerHTML = '<svg width="14" height="14" viewBox="0 0 16 16" fill="currentColor"><polygon points="3,1 13,8 3,15"/></svg> Run PhytClust';
        }
    }
}


/* ─────────────────────────────────────
   Cluster Selector (multiple results)
   ───────────────────────────────────── */
let CLUSTER_VIEW_MODE = "peaks"; // "peaks" or "all"

function populateClusterSelector(data) {
    const controls = document.getElementById("cluster-select-controls");
    const selectEl = document.getElementById("cluster-select");
    const labelEl  = document.getElementById("cluster-select-label");
    const toggleEl = document.getElementById("cluster-view-toggle");

    if (!controls || !selectEl) return;

    const hasPeaks = (data.clusters || []).length > 1;
    const hasAll   = !!(data.all_clusters && data.all_clusters.length > 1);

    if (!hasPeaks && !hasAll) {
        controls.classList.remove("visible");
        if (toggleEl) toggleEl.style.display = "none";
        return;
    }

    // Show toggle only if all_clusters is available
    if (toggleEl) {
        toggleEl.style.display = hasAll ? "inline-flex" : "none";
    }

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
    if (CLUSTER_VIEW_MODE === "all" && latestApiData.all_clusters) {
        clusters = latestApiData.all_clusters;
        ks = latestApiData.all_ks;
    } else {
        clusters = latestApiData.clusters || [];
        ks = latestApiData.ks || [];
    }

    if (idx < 0 || idx >= clusters.length) return;

    const clusterMap = clusters[idx];
    CURRENT_CLUSTERS = clusterMap;

    if (Object.keys(clusterMap).length > 0) {
        const maxCluster = Math.max(...Object.values(clusterMap));
        CLUSTER_COLORS = generateClusterColors(maxCluster + 1);
    } else {
        CLUSTER_COLORS = [];
    }
    drawTree();

    const labelEl = document.getElementById("cluster-select-label");
    if (labelEl) labelEl.textContent = `${idx + 1} / ${clusters.length}`;
}

function toggleClusterViewMode() {
    if (!latestApiData) return;
    CLUSTER_VIEW_MODE = (CLUSTER_VIEW_MODE === "peaks") ? "all" : "peaks";
    populateClusterSelector(latestApiData);
    switchCluster(0);

    const toggleEl = document.getElementById("cluster-view-toggle");
    if (toggleEl) {
        toggleEl.textContent = CLUSTER_VIEW_MODE === "all" ? "Show peaks only" : "Show all k";
    }
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

    if (!Array.isArray(scores) || scores.length === 0) {
        plotEl.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);">No scores available to plot.</div>';
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
        return;
    }

    plotEl.innerHTML = "";

    var dataPoints = [];
    for (let i = 0; i < scores.length - 1; i++) {
        dataPoints.push({ k: i + 2, score: scores[i + 1] });
    }

    var peakPoints = peaks
        .map(k => { var idx = k - 2; return (idx >= 0 && idx < dataPoints.length) ? { k, score: dataPoints[idx].score } : null; })
        .filter(d => d);

    var width  = plotEl.clientWidth  || 700;
    var height = plotEl.clientHeight || 420;
    var margin = { top: 24, right: 24, bottom: 48, left: 58 };
    var innerWidth  = width  - margin.left - margin.right;
    var innerHeight = height - margin.top  - margin.bottom;

    var svg = d3.select(plotEl).append("svg")
        .attr("width", width).attr("height", height)
        .style("background", "transparent");

    var g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    // Axis mode
    var axisModeEl = document.getElementById('axis-mode');
    var defaultMode = (scores.length > 50) ? 'log' : 'normal';
    var mode = defaultMode;
    if (axisModeEl) {
        if (!axisModeEl.__initialized) { axisModeEl.value = defaultMode; axisModeEl.__initialized = true; }
        if (axisModeEl.__userOverride) { mode = axisModeEl.value || defaultMode; }
        else { axisModeEl.value = defaultMode; }
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

    // Grid lines
    g.append("g").attr("class", "grid")
        .attr("transform", `translate(0,${innerHeight})`)
        .call(d3.axisBottom(xScale).tickSize(-innerHeight).tickFormat(""))
        .selectAll("line").attr("stroke", "#e2e8f0").attr("stroke-dasharray", "2,2");

    g.append("g").attr("class", "grid")
        .call(d3.axisLeft(yScale).tickSize(-innerWidth).tickFormat(""))
        .selectAll("line").attr("stroke", "#e2e8f0").attr("stroke-dasharray", "2,2");

    g.selectAll(".grid .domain").remove();

    // Axes
    var xAxisG = g.append("g").attr("transform", `translate(0,${innerHeight})`).call(xAxis);
    var yAxisG = g.append("g").call(d3.axisLeft(yScale).ticks(6));

    [xAxisG, yAxisG].forEach(ag => {
        ag.selectAll("text").attr("fill", "#64748b").style("font-size", "11px");
        ag.selectAll("line").attr("stroke", "#94a3b8");
        ag.selectAll("path").attr("stroke", "#94a3b8");
    });

    // Labels
    g.append("text").attr("x", innerWidth / 2).attr("y", innerHeight + 38)
        .attr("text-anchor", "middle").attr("font-size", 12).attr("fill", "#64748b")
        .text("k (number of clusters)");

    g.append("text").attr("transform", "rotate(-90)")
        .attr("x", -innerHeight / 2).attr("y", -42)
        .attr("text-anchor", "middle").attr("font-size", 12).attr("fill", "#64748b")
        .text("CalBow Score");

    // Line
    const line = d3.line()
        .x(d => mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2)
        .y(d => yScale(d.score))
        .curve(d3.curveMonotoneX);

    g.append("path").datum(dataPoints)
        .attr("fill", "none").attr("stroke", "#3b82f6").attr("stroke-width", 2)
        .attr("d", line);

    // Points
    g.selectAll(".score-point").data(dataPoints).enter().append("circle")
        .attr("class", "score-point")
        .attr("cx", d => mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2)
        .attr("cy", d => yScale(d.score))
        .attr("r", 3.5).attr("fill", "#3b82f6").attr("stroke", "#fff").attr("stroke-width", 1.5)
        .on("mouseover", (event, d) => d3Tooltip.style("opacity", 1).html(`k = ${d.k}<br/>score = ${d.score.toFixed(4)}`))
        .on("mousemove", event => d3Tooltip.style("left", (event.pageX + 12) + "px").style("top", (event.pageY + 12) + "px"))
        .on("mouseout", () => d3Tooltip.style("opacity", 0));

    // Peaks
    g.selectAll(".peak-point").data(peakPoints).enter().append("path")
        .attr("class", "peak-point")
        .attr("transform", d => `translate(${mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2},${yScale(d.score)})`)
        .attr("d", d3.symbol().type(d3.symbolDiamond).size(100))
        .attr("fill", "#ef4444").attr("stroke", "#fff").attr("stroke-width", 1.5)
        .on("mouseover", (event, d) => d3Tooltip.style("opacity", 1).html(`<strong>Optimal k = ${d.k}</strong><br/>score = ${d.score.toFixed(4)}`))
        .on("mousemove", event => d3Tooltip.style("left", (event.pageX + 12) + "px").style("top", (event.pageY + 12) + "px"))
        .on("mouseout", () => d3Tooltip.style("opacity", 0));

    // Legend
    const legend = g.append("g").attr("transform", `translate(${innerWidth - 130}, 8)`);
    legend.append("rect").attr("x", -8).attr("y", -8).attr("width", 140).attr("height", 44)
        .attr("fill", "rgba(255,255,255,0.85)").attr("rx", 6).attr("stroke", "#e2e8f0");
    legend.append("line").attr("x1", 0).attr("y1", 6).attr("x2", 18).attr("y2", 6).attr("stroke", "#3b82f6").attr("stroke-width", 2);
    legend.append("circle").attr("cx", 9).attr("cy", 6).attr("r", 3).attr("fill", "#3b82f6");
    legend.append("text").attr("x", 24).attr("y", 10).attr("font-size", 11).attr("fill", "#64748b").text("Score");
    legend.append("path").attr("transform", "translate(9,26)").attr("d", d3.symbol().type(d3.symbolDiamond).size(80)).attr("fill", "#ef4444");
    legend.append("text").attr("x", 24).attr("y", 30).attr("font-size", 11).attr("fill", "#64748b").text("Optimal k");

    // Wire axis mode once
    if (axisModeEl && !axisModeEl.__wired) {
        axisModeEl.__wired = true;
        axisModeEl.addEventListener('change', function () {
            axisModeEl.__userOverride = true;
            drawOptimalK(data);
        });
    }

    // Save SVG snapshot
    const okSvgNode = document.querySelector("#optimalk_plot svg");
    if (okSvgNode) {
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = new XMLSerializer().serializeToString(okSvgNode);
    } else {
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
    }
}


/* ─────────────────────────────────────
   File Handling (upload + drag & drop)
   ───────────────────────────────────── */
function handleFileSelect(evt) {
    var file = evt.target.files ? evt.target.files[0] : null;
    if (!file) return;
    loadFile(file);
}

function loadFile(file) {
    var fileNameEl = document.getElementById("file-name");
    var fileBadge  = document.getElementById("file-badge");
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
   Save / Export Handlers
   ───────────────────────────────────── */
function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { document.body.removeChild(a); URL.revokeObjectURL(url); }, 100);
}

function exportSvgFromEl(selector, filename) {
    const svgNode = document.querySelector(selector + " svg");
    if (!svgNode) { showToast("No visualization to export.", "danger"); return; }
    const serializer = new XMLSerializer();
    const svgStr = serializer.serializeToString(svgNode);
    const blob = new Blob([svgStr], { type: "image/svg+xml" });
    downloadBlob(blob, filename);
    showToast("Downloaded " + filename, "success", 2000);
}

async function exportTSV() {
    try {
        const res = await fetch("/api/export_tsv", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ top_n: 1, outlier: true })
        });
        if (!res.ok) {
            const err = await res.json().catch(() => ({}));
            throw new Error(err.detail || "Export failed");
        }
        const text = await res.text();
        const blob = new Blob([text], { type: "text/tab-separated-values" });
        downloadBlob(blob, "phytclust_results.tsv");
        showToast("Downloaded TSV", "success", 2000);
    } catch (e) {
        showToast("TSV export failed: " + e.message, "danger");
    }
}

async function saveToServer() {
    const outputDir = document.getElementById("output-dir");
    const dir = outputDir ? (outputDir.value || "results") : "results";
    try {
        const data = await apiPostJson("/api/save", { results_dir: dir, top_n: 1, outlier: true });
        showToast("Saved to " + dir, "success", 3000);
    } catch (e) {
        showToast("Save failed: " + e.message, "danger");
    }
}

function copySvgToClipboard(selector) {
    const svgNode = document.querySelector(selector + " svg");
    if (!svgNode) { showToast("Nothing to copy.", "danger"); return; }
    const svgStr = new XMLSerializer().serializeToString(svgNode);
    navigator.clipboard.writeText(svgStr).then(
        () => showToast("SVG copied to clipboard", "success", 2000),
        () => showToast("Copy failed", "danger")
    );
}


/* ─────────────────────────────────────
   Tab System (custom, no Bootstrap)
   ───────────────────────────────────── */
function switchTab(tabName) {
    // Update tab buttons
    document.querySelectorAll(".tab-btn").forEach(btn => {
        btn.classList.toggle("active", btn.dataset.tab === tabName);
    });

    // Update panels
    document.querySelectorAll(".viz-panel").forEach(panel => {
        panel.classList.toggle("active", panel.id === tabName);
    });

    // Toggle toolbars
    const treeToolbar = document.getElementById("tree-toolbar");
    const optkToolbar = document.getElementById("optk-toolbar");
    if (treeToolbar) treeToolbar.style.display = (tabName === "viewer") ? "" : "none";
    if (optkToolbar) optkToolbar.style.display = (tabName === "viewer-optimal-k") ? "" : "none";

    // Redraw optimal k chart when tab is shown
    if (tabName === "viewer-optimal-k" && latestOptimalKData) {
        drawOptimalK(latestOptimalKData);
    }
}


/* ─────────────────────────────────────
   DOM Ready – Wire Everything
   ───────────────────────────────────── */
document.addEventListener("DOMContentLoaded", function () {

    // Hide loading
    var appLoading = document.getElementById("app-loading");
    if (appLoading) appLoading.style.display = "none";

    showStatus("Ready", "success");

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
            if (e.dataTransfer.files && e.dataTransfer.files.length) {
                loadFile(e.dataTransfer.files[0]);
            }
        });
    }

    // ── Run button ──
    var btnRun = document.getElementById("btn-run");
    if (btnRun) btnRun.addEventListener("click", function (e) { e.preventDefault(); runPhytClust(); });

    // ── Reset button ──
    var resetBtn = document.getElementById("btn-extra-reset");
    if (resetBtn) resetBtn.addEventListener("click", function (e) { e.preventDefault(); resetExtraParams(); });

    // ── Mode selector ──
    document.querySelectorAll("#mode-selector .mode-btn").forEach(function (btn) {
        btn.addEventListener("click", function () {
            document.querySelectorAll("#mode-selector .mode-btn").forEach(b => b.classList.remove("active"));
            btn.classList.add("active");

            var mode = btn.dataset.mode;
            // Toggle param visibility
            var paramK    = document.getElementById("param-k");
            var paramTopN = document.getElementById("param-topn");
            var paramBins = document.getElementById("param-bins");

            if (paramK)    paramK.style.display    = (mode === "k") ? "" : "none";
            if (paramTopN) paramTopN.style.display  = (mode === "k") ? "none" : "";
            if (paramBins) paramBins.style.display  = (mode === "resolution") ? "" : "none";

            // Sync hidden checkbox
            if (extraResolutionEl) extraResolutionEl.checked = (mode === "resolution");
        });
    });

    // ── Advanced options collapse ──
    var collapseToggle = document.getElementById("extra-params-toggle");
    var collapseBody   = document.getElementById("extra-params");
    var collapseIcon   = document.getElementById("extra-params-icon");
    if (collapseToggle && collapseBody) {
        collapseToggle.addEventListener("click", function () {
            var isOpen = collapseBody.style.display !== "none";
            collapseBody.style.display = isOpen ? "none" : "";
            collapseToggle.setAttribute("aria-expanded", !isOpen);
            if (collapseIcon) collapseIcon.innerHTML = isOpen ? "&#x25B6;" : "&#x25BC;";
        });
    }

    // ── Tabs ──
    document.querySelectorAll(".tab-btn").forEach(function (btn) {
        btn.addEventListener("click", function () {
            switchTab(btn.dataset.tab);
        });
    });

    // ── Help sidebar ──
    var btnHelp = document.getElementById("btn-help");
    var helpSidebar = document.getElementById("help-sidebar");
    if (btnHelp && helpSidebar) {
        btnHelp.addEventListener("click", function () {
            helpSidebar.classList.toggle("open");
        });
    }
    var btnHelpClose = document.getElementById("btn-help-close");
    if (btnHelpClose && helpSidebar) {
        btnHelpClose.addEventListener("click", function () {
            helpSidebar.classList.remove("open");
        });
    }

    // ── About modal ──
    var btnAbout     = document.getElementById("btn-about");
    var aboutModal   = document.getElementById("aboutModal");
    var aboutBg      = document.getElementById("aboutBackdrop");
    var btnAboutClose  = document.getElementById("btn-about-close");
    var btnAboutClose2 = document.getElementById("btn-about-close2");

    function openAbout()  { if (aboutModal) aboutModal.classList.add("show"); if (aboutBg) aboutBg.classList.add("show"); }
    function closeAbout() { if (aboutModal) aboutModal.classList.remove("show"); if (aboutBg) aboutBg.classList.remove("show"); }

    if (btnAbout)       btnAbout.addEventListener("click", openAbout);
    if (btnAboutClose)  btnAboutClose.addEventListener("click", closeAbout);
    if (btnAboutClose2) btnAboutClose2.addEventListener("click", closeAbout);
    if (aboutBg)        aboutBg.addEventListener("click", closeAbout);

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
    var btnSave      = document.getElementById("btn-save");
    var saveDropdown  = document.getElementById("save-dropdown");
    if (btnSave && saveDropdown) {
        btnSave.addEventListener("click", function (e) {
            e.stopPropagation();
            saveDropdown.classList.toggle("show");
        });
        document.addEventListener("click", function () {
            saveDropdown.classList.remove("show");
        });
        saveDropdown.addEventListener("click", function (e) { e.stopPropagation(); });

        saveDropdown.querySelectorAll(".dropdown-item").forEach(function (item) {
            item.addEventListener("click", function () {
                saveDropdown.classList.remove("show");
                var type = item.dataset.saveType;
                if (type === "tsv")         exportTSV();
                if (type === "tree_plot")   exportSvgFromEl("#tree_display", "phytclust_tree.svg");
                if (type === "k_plot")      exportSvgFromEl("#optimalk_plot", "phytclust_scores.svg");
                if (type === "save_server") saveToServer();
            });
        });
    }

    // ── Copy buttons ──
    var btnCopyTree = document.getElementById("btn-copy-tree");
    if (btnCopyTree) btnCopyTree.addEventListener("click", function () { copySvgToClipboard("#tree_display"); });

    var btnCopyMaxK = document.getElementById("btn-copy-maxk");
    if (btnCopyMaxK) btnCopyMaxK.addEventListener("click", function () { copySvgToClipboard("#optimalk_plot"); });

    // ── Cluster selector ──
    var clusterSelect = document.getElementById("cluster-select");
    if (clusterSelect) {
        clusterSelect.addEventListener("change", function () {
            switchCluster(parseInt(this.value, 10));
        });
    }
    var clusterToggle = document.getElementById("cluster-view-toggle");
    if (clusterToggle) {
        clusterToggle.addEventListener("click", toggleClusterViewMode);
    }

    // ── Layout mode selector ──
    var layoutSelect = document.getElementById("layout-mode");
    if (layoutSelect) {
        layoutSelect.addEventListener("change", function () {
            CURRENT_LAYOUT_MODE = this.value;
            drawTree();
        });
    }

    // ── Color mode ──
    var colorModeSelect = document.getElementById("color-mode");
    if (colorModeSelect) {
        colorModeSelect.addEventListener("change", function () {
            COLOR_MODE = this.value || "bars";
            drawTree();
        });
    }

    // ── Palette selector ──
    let CUSTOM_PALETTE = null;

    function normalizeHexColor(s) {
        const v = (s || "").trim();
        if (!v) return null;
        if (/^#[0-9a-fA-F]{6}$/.test(v)) return v;
        if (/^[0-9a-fA-F]{6}$/.test(v)) return "#" + v;
        return null;
    }

    function updatePalettePreview(colors) {
        const box = document.getElementById("palette-preview");
        if (!box) return;
        box.innerHTML = colors
            .map(c => `<span title="${c}" style="background:${c};"></span>`)
            .join("");
    }

    function applyPalette(type) {
        const palettes = {
            default: ["#b84b4b", "#849060", "#3d7c74", "#6e3f8a", "#ceb94b", "#3f648a", "#3f408a", "#da63aa"],
            pastel:  ["#ffb3ba", "#ffdfba", "#ffffba", "#baffc9", "#bae1ff", "#d7baff", "#ffcce6", "#c2f0c2"],
            vivid:   ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"],
            dark:    ["#4e79a7", "#59a14f", "#e15759", "#b07aa1", "#edc948", "#76b7b2", "#ff9da7", "#9c755f"]
        };

        if (type === "custom" && CUSTOM_PALETTE && CUSTOM_PALETTE.length) {
            BASE_COLORS.splice(0, BASE_COLORS.length, ...CUSTOM_PALETTE);
        } else if (palettes[type]) {
            BASE_COLORS.splice(0, BASE_COLORS.length, ...palettes[type]);
        }

        updatePalettePreview(BASE_COLORS);

        if (Object.keys(CURRENT_CLUSTERS || {}).length > 0) {
            const maxCluster = Math.max(...Object.values(CURRENT_CLUSTERS));
            CLUSTER_COLORS = generateClusterColors(maxCluster + 1);
            drawTree();
        }
    }

    var paletteSelect = document.getElementById("palette-select");
    if (paletteSelect) {
        paletteSelect.addEventListener("change", function () {
            if (this.value === "custom") {
                var initial = (CUSTOM_PALETTE && CUSTOM_PALETTE.length) ? CUSTOM_PALETTE.join(",") : BASE_COLORS.join(",");
                var raw = prompt("Enter hex colors separated by commas (e.g. #ff0000,#00ff00,#0000ff):", initial);
                if (raw == null) { this.value = "default"; applyPalette("default"); return; }
                var parsed = raw.split(",").map(normalizeHexColor).filter(Boolean);
                if (parsed.length < 2) { alert("Please provide at least 2 valid hex colors."); this.value = "default"; applyPalette("default"); return; }
                CUSTOM_PALETTE = parsed;
                applyPalette("custom");
                return;
            }
            applyPalette(this.value);
        });
    }
    updatePalettePreview(BASE_COLORS);

    // ── Checkboxes: leaf/internal names ──
    var internalCb = document.getElementById("show-internal-names");
    if (internalCb) {
        internalCb.checked = SHOW_INTERNAL_NAMES;
        internalCb.addEventListener("change", function () { SHOW_INTERNAL_NAMES = this.checked; drawTree(); });
    }

    var leafCb = document.getElementById("show-leaf-names");
    if (leafCb) {
        leafCb.checked = SHOW_LEAF_NAMES;
        leafCb.addEventListener("change", function () { SHOW_LEAF_NAMES = this.checked; drawTree(); });
    }

    // ── Node/label sizes ──
    var nodeSizeInput = document.getElementById("node-size");
    if (nodeSizeInput) {
        nodeSizeInput.value = LEAF_NODE_RADIUS;
        nodeSizeInput.addEventListener("input", function () {
            var v = parseFloat(this.value);
            if (!isNaN(v) && v >= 1 && v <= 10) { LEAF_NODE_RADIUS = v; drawTree(); }
        });
    }

    var labelSizeInput = document.getElementById("label-size");
    if (labelSizeInput) {
        labelSizeInput.value = LABEL_FONT_SIZE;
        labelSizeInput.addEventListener("input", function () {
            var v = parseFloat(this.value);
            if (!isNaN(v) && v >= 6 && v <= 18) { LABEL_FONT_SIZE = v; drawTree(); }
        });
    }

    // ── Width/Height sliders ──
    var widthSlider  = document.getElementById("tree-width-scale");
    var heightSlider = document.getElementById("tree-height-scale");
    if (widthSlider) {
        widthSlider.addEventListener("input", function () {
            TREE_WIDTH_SCALE = parseFloat(this.value);
            if (!isNaN(TREE_WIDTH_SCALE)) drawTree();
        });
    }
    if (heightSlider) {
        heightSlider.addEventListener("input", function () {
            TREE_HEIGHT_SCALE = parseFloat(this.value);
            if (!isNaN(TREE_HEIGHT_SCALE)) drawTree();
        });
    }

    // ── Reset view ──
    var btnResetView = document.getElementById("btn-reset-view");
    if (btnResetView) {
        btnResetView.addEventListener("click", function (e) {
            e.preventDefault();
            if (LAST_TREE_SVG && LAST_TREE_ZOOM) {
                LAST_TREE_SVG.transition().duration(150).call(LAST_TREE_ZOOM.transform, d3.zoomIdentity);
            }
            if (widthSlider)  widthSlider.value  = "1.0";
            if (heightSlider) heightSlider.value  = "1.0";
            TREE_WIDTH_SCALE  = 1.0;
            TREE_HEIGHT_SCALE = 1.0;
            clearAllCollapsedFlags(NEWICK_RAW_TREE);
            drawTree();
        });
    }

    // ── Prefill example and auto-run ──
    if (newickEl) newickEl.value = EXAMPLE_NEWICK;
    runPhytClust();
});
