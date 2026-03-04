let HIER_CART = null;
let HIER_CIRC = null;
var EXAMPLE_NEWICK = "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

var newickEl = document.getElementById("newick-input");
var resultEl = document.getElementById("result");
var statusEl = document.getElementById("status-message");
var treeHost = document.getElementById("tree_display");

const state = {
  isRunning: false,
  latestOptimalKData: null,
  newickRawTree: null,
  currentClusters: {},
  clusterColors: [],
  layoutMode: "rectangular",
  colorMode: "bars",
  showLeafNames: true,
  showInternalNames: false,
  leafNodeRadius: 3.0,
  internalNodeRadius: 1.8,
  labelFontSize: 9,
  treeWidthScale: 1.0,
  treeHeightScale: 1.0
};

const BASE_COLORS = [
        "#b84b4b",
        "#849060",
        "#3d7c74",
        "#6e3f8a",
        "#ceb94b",
        "#3f648a",
        "#3f408a",
        "#da63aa"
    ];

    let SHOW_LEAF_NAMES = true;
    let SHOW_INTERNAL_NAMES = false;

    let LEAF_NODE_RADIUS = 3.0;
    let INTERNAL_NODE_RADIUS = 1.8;

    let LABEL_FONT_SIZE = 9;

    // Add missing globals used by drawTree() + slider handlers
    let TREE_WIDTH_SCALE = 1.0;
    let TREE_HEIGHT_SCALE = 1.0;

    let CLUSTER_COLORS = [];
    let CURRENT_CLUSTERS = {};
    let NEWICK_RAW_TREE = null;
    let CURRENT_LAYOUT_MODE = 'rectangular';
    let COLOR_MODE = 'bars';

    const LABEL_PAD = 6; // constant visual gap between node edge and label

var statusHideTimeout = null;
var isRunning = false;
var extraKEl = document.getElementById("extra-k");
var extraOutgroupEl = document.getElementById("extra-outgroup");
var extraTopNEl = document.getElementById("extra-topn");
var extraResolutionEl = document.getElementById("extra-resolution");
var extraBinsEl = document.getElementById("extra-bins");
var extraMaxKEl = document.getElementById("extra-maxk");
var extraMaxKLimitEl = document.getElementById("extra-maxklimit");
var extraLambdaEl = document.getElementById("extra-lambda");

var latestOptimalKData = null;

var d3Tooltip = d3.select("body").append("div")
    .attr("class", "d3-tooltip");

function resetExtraParams() {
    if (extraKEl) extraKEl.value = "";
    if (extraOutgroupEl) extraOutgroupEl.value = "";
    if (extraTopNEl) extraTopNEl.value = "";
    if (extraResolutionEl) extraResolutionEl.checked = false;
    if (extraBinsEl) extraBinsEl.value = "";
    if (extraMaxKEl) extraMaxKEl.value = "";
    if (extraMaxKLimitEl) extraMaxKLimitEl.value = "";
    if (extraLambdaEl) extraLambdaEl.value = "";
    showStatus("Extra parameters reset to defaults.", "info");
}

function escapeRegExp(str) {
    return str.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

function isOutgroupInNewick(newick, outgroup) {
    if (!outgroup) return true;
    var name = escapeRegExp(outgroup);
    var quoted = new RegExp("[\\(,]\\s*'" + name + "'\\s*(?:[:,\\)])");
    var unquoted = new RegExp("[\\(,]\\s*" + name + "\\s*(?:[:,\\)])");
    return quoted.test(newick) || unquoted.test(newick);
}

// Rough leaf estimate
function estimateLeafCount(newick) {
    var inQuote = null;
    var commas = 0;
    for (var i = 0; i < newick.length; i++) {
        var ch = newick[i];
        if (inQuote) {
            if (ch === inQuote) { inQuote = null; }
            continue;
        }
        if (ch === '"' || ch === "'") { inQuote = ch; continue; }
        if (ch === ',') { commas++; }
    }
    return Math.max(0, commas + 1);
}

// Bootstrap 5 / no-jQuery status helper (replace the whole function with this)
function showStatus(message, type) {
    type = type || "info";

    var el = document.getElementById("status-message");
    if (!el) return;

    if (statusHideTimeout) { clearTimeout(statusHideTimeout); statusHideTimeout = null; }

    el.className = "";
    el.classList.add("alert", "alert-" + type);
    el.textContent = message;
    el.style.display = "block";

    if (type === "success") {
        statusHideTimeout = setTimeout(function () {
            el.style.display = "none";
        }, 5000);
    }
}

function clearTree() {
    treeHost.innerHTML = "";
}

// ---- Newick parsing (no external libs) ----
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
                if (token && token !== ";") {
                    current.name = token;
                }
        }
    }
    return current;
}

function getVisibleChildren(node) {
    // If this node is collapsed, treat it as a leaf in the rendered hierarchy
    if (node && node._collapsed) return null;
    return node && node.children && node.children.length ? node.children : null;
}

async function runPhytClust() {
    if (isRunning) { return; }
    var newickText = (newickEl.value || "").trim();
    if (!newickText) {
        alert("Please upload a Newick tree file or paste a Newick string.");
        return;
    }
    var numSamples = estimateLeafCount(newickText);

    var kVal = null;
    if (extraKEl) {
        var kRaw = (extraKEl.value || "").trim();
        if (kRaw) {
            var parsed = parseInt(kRaw, 10);
            if (!isNaN(parsed) && parsed >= 1) {
                kVal = parsed;
            }
        }
    }

    var payload;
    var outgroupVal = null;
    if (extraOutgroupEl) {
        var ogRaw = (extraOutgroupEl.value || "").trim();
        if (ogRaw) { outgroupVal = ogRaw; }
    }
    var topNVal = null;
    if (extraTopNEl) {
        var tnRaw = (extraTopNEl.value || "").trim();
        if (tnRaw) { var tni = parseInt(tnRaw, 10); if (!isNaN(tni) && tni >= 1) topNVal = tni; }
    }
    var resolutionFlag = !!(extraResolutionEl && extraResolutionEl.checked);
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
                if (mki < 4) { mki = 4; if (extraMaxKEl) extraMaxKEl.value = String(mki); }
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
                if (mkli < minLimit) { mkli = minLimit; if (extraMaxKLimitEl) extraMaxKLimitEl.value = String(mkli); }
                maxKLimitVal = mkli;
            }
        }
    }
    var lambdaVal = null;
    if (extraLambdaEl) {
        var lRaw = (extraLambdaEl.value || "").trim();
        if (lRaw) { var li = parseFloat(lRaw); if (!isNaN(li)) lambdaVal = li; }
    }

    if (outgroupVal && !isOutgroupInNewick(newickText, outgroupVal)) {
        showStatus("Outgroup not found in Newick.", "danger");
        return;
    }
    if (kVal !== null) {
        payload = { newick: newickText, mode: "k", k: kVal };
        if (outgroupVal !== null) { payload.outgroup = outgroupVal; }
        if (lambdaVal !== null) { payload.lambda_weight = lambdaVal; }
        showStatus("Running PhytClust (fixed k-mode)…", "info");
        resultEl.textContent = "Running PhytClust (k-mode)…";
    } else if (resolutionFlag) {
        payload = { newick: newickText, mode: "resolution", top_n: (topNVal || 1) };
        if (outgroupVal !== null) { payload.outgroup = outgroupVal; }
        if (binsVal !== null) { payload.num_bins = binsVal; }
        if (maxKVal !== null) { payload.max_k = maxKVal; }
        if (maxKLimitVal !== null) { payload.max_k_limit = maxKLimitVal; }
        if (lambdaVal !== null) { payload.lambda_weight = lambdaVal; }
        showStatus("Running PhytClust (resolution mode)…", "info");
        resultEl.textContent = "Running PhytClust (resolution mode)…";
    } else {
        payload = { newick: newickText, mode: "global", top_n: (topNVal || 1) };
        if (outgroupVal !== null) { payload.outgroup = outgroupVal; }
        if (maxKVal !== null) { payload.max_k = maxKVal; }
        if (maxKLimitVal !== null) { payload.max_k_limit = maxKLimitVal; }
        if (lambdaVal !== null) { payload.lambda_weight = lambdaVal; }
        showStatus("Running PhytClust...", "info");
        resultEl.textContent = "Running PhytClust...";
    }
    clearTree();

    try {
        isRunning = true;
        var saveBtn = document.getElementById('btn-save');
        if (saveBtn) { saveBtn.disabled = true; }
        const t0 = performance.now();
        const data = await apiPostJson("/api/run", payload);
        const dt = (performance.now() - t0) / 1000;
        showStatus(`PhytClust finished in ${dt.toFixed(2)}s`, "success");
        resultEl.textContent = JSON.stringify(data, null, 2);

        var logEl = document.getElementById("log-output");
        if (logEl) {
            if (data && data.log) {
                logEl.textContent = data.log;
            } else {
                logEl.textContent = "(no logs)";
            }
        }

        if (data.newick) {
            NEWICK_RAW_TREE = parseNewick(data.newick);
            accumulateBranchLength(NEWICK_RAW_TREE);
            computeLayouts();
            const clusterMap = (data.clusters && data.clusters.length > 0)
                ? data.clusters[0]
                : {};

            CURRENT_CLUSTERS = clusterMap;

            // Generate cluster colors BEFORE drawing the tree
            if (Object.keys(clusterMap).length > 0) {
                const maxCluster = Math.max(...Object.values(clusterMap));
                CLUSTER_COLORS = generateClusterColors(maxCluster + 1);
                showStatus("PhytClust finished.", "success");
            } else {
                CLUSTER_COLORS = [];
                showStatus("No clusters found (likely not enough samples).", "danger");
            }

            // Now draw using updated cluster colors
            drawTree();

        } else {
            clearTree();
            showStatus("No Newick tree returned by API.", "danger");
        }

        try {
            if (data && data.scores) {
                var plotHostReady = document.getElementById('optimalk_plot');
                if (plotHostReady) {
                    plotHostReady.innerHTML = '';
                }
                latestOptimalKData = data;
                drawOptimalK(latestOptimalKData);
            } else {
                var plotHost = document.getElementById('optimalk_plot');
                if (plotHost) {
                    plotHost.innerHTML = '<div class="alert alert-info" style="width:100%; margin:0;">Optimal number of clusters plot is unavailable in fixed k mode or when no scores are returned.</div>';
                }
            }
        } catch (plotErr) {
            console.warn('Plot rendering error:', plotErr);
        }
    } catch (e) {
        console.error(e);
        showStatus("PhytClust failed with Error: " + e, "danger");
        resultEl.textContent = "Error: " + e;
        latestOptimalKData = null;
        let plotHost = document.getElementById("optimalk_plot");
        if (plotHost) plotHost.innerHTML = "";
    } finally {
        isRunning = false;
        var saveBtn2 = document.getElementById('btn-save');
        if (saveBtn2) { saveBtn2.disabled = false; }
    }
}

// Shuffle Fisher–Yates
function shuffle(arr) {
    let a = arr.slice();
    for (let i = a.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [a[i], a[j]] = [a[j], a[i]];
    }
    return a;
}

// Lighten or darken by factor (-0.5 to 0.5)
function adjustLight(hex, factor) {
    const num = parseInt(hex.slice(1), 16);
    let r = (num >> 16) + factor * 255;
    let g = ((num >> 8) & 0xFF) + factor * 255;
    let b = (num & 0xFF) + factor * 255;
    r = Math.min(255, Math.max(0, r));
    g = Math.min(255, Math.max(0, g));
    b = Math.min(255, Math.max(0, b));
    return `rgb(${r | 0}, ${g | 0}, ${b | 0})`;
}

// Add alpha to color
function withAlpha(color, alpha) {
        if (color.startsWith("rgb")) {
            return color.replace("rgb", "rgba").replace(")", `, ${alpha})`);
        }
        const num = parseInt(color.slice(1), 16);
        const r = num >> 16;
        const g = (num >> 8) & 0xFF;
        const b = num & 0xFF;
        return `rgba(${r}, ${g}, ${b}, ${alpha})`;
    }

function generateClusterColors(nClusters) {
    let palette = shuffle(BASE_COLORS);   // randomize ordering
    let colors = [];

    if (nClusters <= palette.length) {
        return palette.slice(0, nClusters); // simple case
    }

    const repeats = Math.ceil(nClusters / palette.length);

    for (let r = 0; r < repeats; r++) {
        const factor = (r - 1) * 0.18;  // lighten/darken steps
        const alpha = 1 - r * 0.15;

        palette.forEach(hex => {
            const adjusted = adjustLight(hex, factor);
            const finalCol = alpha < 1 ? withAlpha(adjusted, alpha) : adjusted;
            colors.push(finalCol);
        });
    }

    return colors.slice(0, nClusters);
}

function accumulateBranchLength(node, length = 0) {
    node._bl = length;
    if (node.children) {
        for (const c of node.children) {
            accumulateBranchLength(c, length + (c.length || 0));
        }
    }
}
function computeLayouts() {
    if (!NEWICK_RAW_TREE) return;

    // Build hierarchy while respecting collapsed nodes
    HIER_CART = d3.hierarchy(NEWICK_RAW_TREE, getVisibleChildren);
    HIER_CIRC = d3.hierarchy(NEWICK_RAW_TREE, getVisibleChildren);

    var width  = treeHost.clientWidth  || 800;
    var height = treeHost.clientHeight || 500;

    var margin = { top:20, right:80, bottom:20, left:80 };
    var innerW = width - margin.left - margin.right;
    var innerH = height - margin.top - margin.bottom;
    var radius = Math.min(innerW, innerH) / 2;

    var maxBl = d3.max(HIER_CART.descendants(), d => d.data._bl || 0) || 1;

    var blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);
    var blToR = d3.scaleLinear().domain([0, maxBl]).range([0, radius]);

    // CARTESIAN
    var cartLayout = d3.cluster().size([innerH, 1]);
    cartLayout(HIER_CART);

    HIER_CART.each(d => {
        d._x = d.x;
        d._y = blToX(d.data._bl || 0);
    });

    // CIRCULAR
    var circLayout = d3.cluster().size([2 * Math.PI, 1]);
    circLayout(HIER_CIRC);

    HIER_CIRC.each(d => {
        d._angle = d.x;
        d._radius = blToR(d.data._bl || 0);
    });
}

let LAST_TREE_SVG = null;
let LAST_TREE_ZOOM = null;

// Clear collapse flags in the backing Newick tree (data)
function clearAllCollapsedFlags(node) {
    if (!node) return;
    if (node._collapsed) delete node._collapsed;
    if (node.children && node.children.length) {
        for (const c of node.children) clearAllCollapsedFlags(c);
    }
}

// Collect leaf names under a *data* node (respects original structure, not visible d3 children)
function collectLeafNamesData(node, out) {
    if (!node) return;
    const kids = node.children && node.children.length ? node.children : null;
    if (!kids) {
        if (node.name) out.push(node.name);
        return;
    }
    for (const c of kids) collectLeafNamesData(c, out);
}

// If a subtree is uniform-cluster, return that cluster id; else null ("mixed"/unknown)
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
        .attr("height", height);

    const zoomLayer = svg.append("g");
    const g = zoomLayer.append("g")
        .attr("transform",
            layoutMode === "circular"
                ? `translate(${width / 2},${height / 2})`
                : `translate(${margin.left},${margin.top})`
        );

    const zoom = d3.zoom()
        .scaleExtent([0.5, 5])
        .on("zoom", (event) => {
            zoomLayer.attr("transform", event.transform);
        });

    svg.call(zoom);
    LAST_TREE_SVG = svg;
    LAST_TREE_ZOOM = zoom;

    function hasOriginalChildren(d) {
        return !!(d && d.data && Array.isArray(d.data.children) && d.data.children.length > 0);
    }

    function canCollapse(d) {
        // collapse is only meaningful if the original data node has children
        return hasOriginalChildren(d);
    }

    // A "visible leaf" is a node with no visible children in the hierarchy (which includes collapsed internal nodes)
    function isVisibleLeaf(d) {
        return !(d && d.children && d.children.length);
    }

    function nodeRadius(d) {
        return hasOriginalChildren(d) ? INTERNAL_NODE_RADIUS : LEAF_NODE_RADIUS;
    }

    function nodeFill(d) {
        const name = d.data && d.data.name;

        // In bar mode, circles are black (bars carry color)
        if (colorMode === "bars") return "black";

        // Leaves: color by cluster
        if (name && !hasOriginalChildren(d)) {
            const cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[name] : null;
            if (cid == null || !CLUSTER_COLORS.length) return "black";
            return CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
        }

        // Internal/collapsed nodes: inherit if subtree is uniform
        const rep = representativeClusterIdFromData(d.data);
        if (rep != null && CLUSTER_COLORS.length) return CLUSTER_COLORS[rep % CLUSTER_COLORS.length];

        return "black";
    }

    // -----------------------------
    // CIRCULAR
    // -----------------------------
    if (layoutMode === "circular") {
        const root = HIER_CIRC;

        const maxRadiusPx = Math.min(width, height) / 2 - 30;
        const maxR = d3.max(root.descendants(), d => d._radius || 0) || 1;

        const rScale = d3.scaleLinear()
            .domain([0, maxR])
            .range([0, maxRadiusPx]);

        function polarToXY(angle, radius) {
            const a = angle - Math.PI / 2;
            const rr = rScale(radius);
            return [
                (rr * Math.cos(a)) * TREE_WIDTH_SCALE,
                (rr * Math.sin(a)) * TREE_HEIGHT_SCALE
            ];
        }

        function nodeXY(d) {
            return polarToXY(d._angle || 0, d._radius || 0);
        }

        function radialElbowPath(link) {
            const s = link.source;
            const t = link.target;

            const sa = s._angle || 0;
            const sr = s._radius || 0;
            const ta = t._angle || 0;
            const tr = t._radius || 0;

            const p0 = polarToXY(sa, sr);
            const p1 = polarToXY(sa, tr);
            const p2 = polarToXY(ta, tr);

            const delta = (ta - sa);
            const sweep = delta >= 0 ? 1 : 0;
            const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
            const arcR = rScale(tr);

            return `M${p0[0]},${p0[1]}L${p1[0]},${p1[1]}A${arcR},${arcR} 0 ${largeArc},${sweep} ${p2[0]},${p2[1]}`;
        }

        // Links
        g.append("g")
            .selectAll(".tree-link")
            .data(root.links())
            .enter()
            .append("path")
            .attr("class", "tree-link")
            .attr("fill", "none")
            .attr("stroke", "#111")
            .attr("stroke-opacity", 0.9)
            .attr("stroke-width", 1.1)
            .attr("d", radialElbowPath);

        // Nodes
        const node = g.append("g")
            .selectAll(".tree-node")
            .data(root.descendants())
            .enter()
            .append("g")
            .attr("class", "tree-node")
            .attr("transform", d => {
                const p = nodeXY(d);
                return `translate(${p[0]},${p[1]})`;
            });

        node.append("circle")
            .attr("r", nodeRadius)
            .attr("fill", nodeFill)
             .style("cursor", d => canCollapse(d) ? "pointer" : "default")
            .on("click", function (event, d) {
                if (!canCollapse(d)) return;
                d.data._collapsed = !d.data._collapsed;
                drawTree(); // recomputeLayouts() will rebuild hierarchy with getVisibleChildren()
            });

        // Labels: always pushed outward, flipped on left half so text reads outward
        node.append("text")
            .attr("dy", "0.32em")
            .attr("font-size", LABEL_FONT_SIZE + "px")
            .attr("transform", function (d) {
                const a = (d._angle || 0);

                // radial direction (outward)
                let deg = (a * 180 / Math.PI) - 90;

                // If on left half, flip text to keep it readable but still outward
                const left = (deg > 90 || deg < -90);

                // For radial text, we rotate by (deg + 90) so baseline points outward (along radius)
                // Then translate outward by node radius + pad
                const radialDeg = deg + 90;

                const r = nodeRadius(d) + LABEL_PAD;

                // On left side, rotate an extra 180 so text isn't upside down
                const flip = left ? 180 : 0;

                return `rotate(${radialDeg}) translate(${r},0) rotate(${flip})`;
            })
            .style("text-anchor", function (d) {
                const deg = ((d._angle || 0) * 180 / Math.PI) - 90;
                const left = (deg > 90 || deg < -90);
                return left ? "end" : "start";
            })
            .text(d => {
                // show label for visible leaves, plus internal names if enabled
                if (isVisibleLeaf(d) && SHOW_LEAF_NAMES) return d.data.name || "";
                if (!isVisibleLeaf(d) && SHOW_INTERNAL_NAMES) return d.data.name || "";
                return "";
            });

        // Circular "side bars" (ring) when colorMode=bars
        if (colorMode === "bars" && CURRENT_CLUSTERS && CLUSTER_COLORS.length) {
            // Use visible leaves (including collapsed internal nodes as leaves)
            const visLeaves = root.descendants()
                .filter(d => isVisibleLeaf(d))
                .sort((a, b) => (a._angle || 0) - (b._angle || 0));

            if (visLeaves.length >= 2) {
                const angles = visLeaves.map(d => d._angle || 0);

                // Boundaries tile the circle (midpoints). Make them wrap correctly.
                const boundaries = new Array(visLeaves.length + 1);
                for (let i = 1; i < visLeaves.length; i++) {
                    boundaries[i] = (angles[i - 1] + angles[i]) / 2;
                }
                // extrapolate ends
                boundaries[0] = angles[0] - (boundaries[1] - angles[0]);
                boundaries[visLeaves.length] = angles[visLeaves.length - 1] + (angles[visLeaves.length - 1] - boundaries[visLeaves.length - 1]);

                function visibleLeafClusterId(d) {
                    // If it's an original leaf, use its cluster
                    const nm = d.data && d.data.name;
                    const isOrigLeaf = nm && !(d.data.children && d.data.children.length);
                    if (isOrigLeaf) return CURRENT_CLUSTERS[nm];

                    // collapsed internal node: inherit if uniform
                    return representativeClusterIdFromData(d.data);
                }

                const ringInner = maxRadiusPx + 10;
                const ringOuter = ringInner + 14;

                const arc = d3.arc().innerRadius(ringInner).outerRadius(ringOuter);

                // Contiguous runs by cluster id (null -> skip)
                let runStart = 0;
                let runCid = visibleLeafClusterId(visLeaves[0]);

                for (let i = 1; i <= visLeaves.length; i++) {
                    const cid = (i < visLeaves.length) ? visibleLeafClusterId(visLeaves[i]) : Symbol("END");
                    if (cid !== runCid) {
                        if (runCid != null) {
                            const col = CLUSTER_COLORS[runCid % CLUSTER_COLORS.length];
                            g.append("path")
                                .attr("d", arc.startAngle(boundaries[runStart]).endAngle(boundaries[i]))
                                .attr("fill", col)
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

    // -----------------------------
    // CARTESIAN (rectangular/cladogram)
    // -----------------------------
    const root = HIER_CART;
    const allNodes = root.descendants();
    const maxYCart = d3.max(allNodes, d => d._y || 0) || 0;

    const labelWidth  = 24;
    const labelColumnX = (maxYCart * TREE_WIDTH_SCALE) + 40;

    // Links
    g.append("g")
        .selectAll(".tree-link")
        .data(root.links())
        .enter()
        .append("path")
        .attr("class", "tree-link")
        .attr("fill", "none")
        .attr("stroke", "black")
        .attr("stroke-width", 1.2)
        .attr("d", function (d) {
            if (layoutMode === "rectangular") {
                return (
                    "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
                    "V" + (d.target._x * TREE_HEIGHT_SCALE) +
                    "H" + (d.target._y * TREE_WIDTH_SCALE)
                );
            }
            return (
                "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
                "L" + (d.target._y * TREE_WIDTH_SCALE) + "," + (d.target._x * TREE_HEIGHT_SCALE)
            );
        });

    // Nodes
    const node = g.append("g")
        .selectAll(".tree-node")
        .data(root.descendants())
        .enter()
        .append("g")
        .attr("class", "tree-node")
        .attr("transform", function (d) {
            const x = (d._y || 0) * TREE_WIDTH_SCALE;
            const y = (d._x || 0) * TREE_HEIGHT_SCALE;
            return `translate(${x},${y})`;
        });

    node.append("circle")
        .attr("r", nodeRadius)
        .attr("fill", nodeFill)
        .style("cursor", d => canCollapse(d) ? "pointer" : "default")
        .on("click", function (event, d) {
            if (!canCollapse(d)) return;
            d.data._collapsed = !d.data._collapsed;
            drawTree(); // recomputeLayouts() will rebuild hierarchy with getVisibleChildren()
        });

    node.append("text")
        .attr("dy", 3)
        .attr("x", function (d) {
            const r = nodeRadius(d);
            return (d.children && d.children.length) ? -(r + LABEL_PAD) : (r + LABEL_PAD);
        })
        .style("text-anchor", function (d) {
            return (d.children && d.children.length) ? "end" : "start";
        })
        .style("font-size", LABEL_FONT_SIZE + "px")
        .text(function (d) {
            const isLeaf = !d.children || d.children.length === 0;
            if (isLeaf && SHOW_LEAF_NAMES)       return d.data.name || "";
            if (!isLeaf && SHOW_INTERNAL_NAMES)  return d.data.name || "";
            return "";
        });

    // Cluster side bars (cartesian only) — contiguous runs that TOUCH (no gaps)
    if (colorMode === "bars") {
        const leafNodes = allNodes
            .filter(d => !d.children || !d.children.length)
            .sort((a, b) => (a._x || 0) - (b._x || 0));

        const barX = labelColumnX;
        const barW = 24;

        // Estimate band boundaries based on midpoints between adjacent leaves
        // so bars perfectly tile the vertical axis without gaps.
        const ys = leafNodes.map(d => (d._x || 0) * TREE_HEIGHT_SCALE);

        // Edge case: 0/1 leaves
        if (leafNodes.length === 1) {
            const nm = leafNodes[0].data && leafNodes[0].data.name;
            const cid = nm && CURRENT_CLUSTERS ? CURRENT_CLUSTERS[nm] : null;
            if (cid != null && CLUSTER_COLORS.length) {
                const y = ys[0];
                barsG = g.append("g").attr("class", "cluster-bars");
                barsG.append("rect")
                    .attr("x", barX)
                    .attr("y", y - 6)
                    .attr("width", barW)
                    .attr("height", 12)
                    .attr("fill", CLUSTER_COLORS[cid % CLUSTER_COLORS.length])
                    .attr("opacity", 0.7);
            }
        } else if (leafNodes.length > 1) {
            // boundaries[i] is the top edge of leaf i's "cell"
            const boundaries = new Array(leafNodes.length + 1);

            // interior boundaries are midpoints
            for (let i = 1; i < leafNodes.length; i++) {
                boundaries[i] = (ys[i - 1] + ys[i]) / 2;
            }

            // extrapolate first/last boundary using nearest spacing
            boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
            boundaries[leafNodes.length] = ys[leafNodes.length - 1] + (ys[leafNodes.length - 1] - boundaries[leafNodes.length - 1]);

            const barsG = g.append("g").attr("class", "cluster-bars");

            // Build contiguous runs by cluster id
            let runStart = 0;
            let runCid = null;

            function leafCid(i) {
                const nm = leafNodes[i].data && leafNodes[i].data.name;
                if (!nm) return null;
                return (CURRENT_CLUSTERS ? CURRENT_CLUSTERS[nm] : null);
            }

            runCid = leafCid(0);

            for (let i = 1; i <= leafNodes.length; i++) {
                const cid = (i < leafNodes.length) ? leafCid(i) : Symbol("END");

                const changed = (cid !== runCid);
                if (changed) {
                    if (runCid != null && CLUSTER_COLORS.length) {
                        const y0 = boundaries[runStart];
                        const y1 = boundaries[i];
                        const col = CLUSTER_COLORS[runCid % CLUSTER_COLORS.length];

                        barsG.append("rect")
                            .attr("x", barX)
                            .attr("y", y0)
                            .attr("width", barW)
                            .attr("height", Math.max(1, y1 - y0))
                            .attr("fill", col)
                            .attr("opacity", 0.7);
                    }
                    runStart = i;
                    runCid = (i < leafNodes.length) ? cid : null;
                }
            }
        }
    }

    // branch-length axis
    const maxX = d3.max(allNodes, d => d._x || 0) || 0;
    const maxY = d3.max(allNodes, d => d._y || 0) || 0;
    const maxBl = d3.max(allNodes, d => d.data._bl || 0) || 1;

    const axisY = (maxX * TREE_HEIGHT_SCALE) + 20;
    const blScale = d3.scaleLinear()
        .domain([0, maxBl])
        .range([0, maxY * TREE_WIDTH_SCALE]);

    g.selectAll(".branch-length-axis").remove();
    const axisG = g.append("g")
        .attr("class", "branch-length-axis")
        .attr("transform", `translate(0, ${axisY})`)
        .call(d3.axisBottom(blScale).ticks(5));

    axisG.append("text")
        .attr("x", (maxY * TREE_WIDTH_SCALE) / 2)
        .attr("y", 30)
        .attr("text-anchor", "middle")
        .attr("font-size", 10)
        .text("Branch length");

    svg.attr("height", Math.max(height, axisY + 50 + margin.bottom));
}

// Clean apiPostJson (your current copy got spliced into other code)
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

// ---- D3 Optimal-k plot ----
function drawOptimalK(data) {
    if (!data) data = latestOptimalKData;
    var plotEl = document.getElementById('optimalk_plot');
    if (!plotEl) return;

    var scores = data && data.scores ? data.scores : [];
    var peaks = data && data.peaks ? data.peaks : [];

    if (!Array.isArray(scores) || scores.length === 0) {
        plotEl.innerHTML = '<div class="alert alert-warning">No scores available to plot.</div>';
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
        return;
    }

    plotEl.innerHTML = "";

    // Build data points
    var dataPoints = [];
    for (let i = 0; i < scores.length - 1; i++) {
        dataPoints.push({ k: i + 2, score: scores[i + 1] });
    }

    var peakPoints = peaks
        .map(k => {
            var idx = k - 2;
            if (idx >= 0 && idx < dataPoints.length) return { k, score: dataPoints[idx].score };
            return null;
        })
        .filter(d => d);

    var width = plotEl.clientWidth || 700;
    var height = plotEl.clientHeight || 420;
    var margin = { top: 20, right: 20, bottom: 45, left: 55 };
    var innerWidth = width - margin.left - margin.right;
    var innerHeight = height - margin.top - margin.bottom;

    var svg = d3.select(plotEl).append("svg")
        .attr("width", width)
        .attr("height", height);

    var g = svg.append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    // Axis mode handling with user override
    var axisModeEl = document.getElementById('axis-mode');
    var defaultMode = (scores.length > 50) ? 'log' : 'normal';
    var mode = defaultMode;
    if (axisModeEl) {
        if (!axisModeEl.__initialized) {
            axisModeEl.value = defaultMode;
            axisModeEl.__initialized = true;
        }
        if (axisModeEl.__userOverride) {
            mode = axisModeEl.value || defaultMode;
        } else {
            axisModeEl.value = defaultMode;
        }
    }

    // X scale
    var xScale, xAxis;
    if (mode === 'log') {
        xScale = d3.scaleLog()
            .domain([2, d3.max(dataPoints, d => d.k)])
            .range([0, innerWidth]);

        // Custom tick array for log mode
        var xVals = dataPoints.map(d => d.k);
        xAxis = d3.axisBottom(xScale)
            .tickValues(xVals) // always integer k
            .tickFormat(d3.format("d"));
    } else {
        xScale = d3.scaleBand()
            .domain(dataPoints.map(d => d.k))
            .range([0, innerWidth])
            .padding(0.2);

        // Limit x ticks to at most 10 for normal mode
        var nPoints = dataPoints.length;
        var dtick = Math.max(1, Math.ceil(nPoints / 10));
        xAxis = d3.axisBottom(xScale)
            .tickValues(dataPoints.map(d => d.k).filter((d, i) => i % dtick === 0))
            .tickFormat(d3.format("d"));
    }

    const yScale = d3.scaleLinear()
        .domain([d3.min(dataPoints, d => d.score), d3.max(dataPoints, d => d.score)])
        .nice()
        .range([innerHeight, 0]);

    g.append("g")
        .attr("transform", `translate(0,${innerHeight})`)
        .call(xAxis);

    g.append("g").call(d3.axisLeft(yScale).ticks(6));

    // Labels
    g.append("text")
        .attr("x", innerWidth / 2)
        .attr("y", innerHeight + 35)
        .attr("text-anchor", "middle")
        .attr("font-size", 12)
        .text("k (number of clusters)");

    g.append("text")
        .attr("transform", "rotate(-90)")
        .attr("x", -innerHeight / 2)
        .attr("y", -40)
        .attr("text-anchor", "middle")
        .attr("font-size", 12)
        .text("Score");

    // Line generator
    const line = d3.line()
        .x(d => mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth()/2)
        .y(d => yScale(d.score));

    g.append("path")
        .datum(dataPoints)
        .attr("fill", "none")
        .attr("stroke", "#1f77b4")
        .attr("stroke-width", 2)
        .attr("d", line);

    // Points
    g.selectAll(".score-point")
        .data(dataPoints)
        .enter().append("circle")
        .attr("class", "score-point")
        .attr("cx", d => mode === 'log' ? xScale(d.k) : xScale(d.k) + xScale.bandwidth()/2)
        .attr("cy", d => yScale(d.score))
        .attr("r", 4)
        .attr("fill", "#1f77b4")
        .on("mouseover", (event,d) => d3Tooltip.style("opacity",1).html(`k=${d.k}<br/>score=${d.score.toFixed(4)}`))
        .on("mousemove", event => d3Tooltip.style("left",(event.pageX+10)+"px").style("top",(event.pageY+10)+"px"))
        .on("mouseout", () => d3Tooltip.style("opacity",0));

    // Peaks
    g.selectAll(".peak-point")
        .data(peakPoints)
        .enter().append("path")
        .attr("class", "peak-point")
        .attr("transform", d => `translate(${mode==='log'?xScale(d.k):xScale(d.k)+xScale.bandwidth()/2},${yScale(d.score)})`)
        .attr("d", d3.symbol().type(d3.symbolStar).size(120))
        .attr("fill", "#d62728")
        .on("mouseover", (event,d) => d3Tooltip.style("opacity",1).html(`<strong>Optimal k</strong><br/>k=${d.k}<br/>score=${d.score.toFixed(4)}`))
        .on("mousemove", event => d3Tooltip.style("left",(event.pageX+10)+"px").style("top",(event.pageY+10)+"px"))
        .on("mouseout", () => d3Tooltip.style("opacity",0));

    // Legend
    const legend = g.append("g").attr("transform", `translate(${innerWidth - 120},10)`);
    legend.append("line").attr("x1",0).attr("y1",0).attr("x2",20).attr("y2",0).attr("stroke","#1f77b4").attr("stroke-width",2);
    legend.append("text").attr("x",26).attr("y",4).attr("font-size",11).text("Score");
    legend.append("path").attr("transform","translate(10,20)").attr("d",d3.symbol().type(d3.symbolStar).size(120)).attr("fill","#d62728");
    legend.append("text").attr("x",26).attr("y",24).attr("font-size",11).text("Optimal k");

    // Wire up axis mode change ONCE
    if (axisModeEl && !axisModeEl.__wired) {
        axisModeEl.__wired = true;
        axisModeEl.addEventListener('change', function () {
            axisModeEl.__userOverride = true;
            drawOptimalK(data);
        });
    }

    // Save latest SVG snapshot for export
    const okSvgNode = document.querySelector("#optimalk_plot svg");
    if (okSvgNode) {
        const serializer = new XMLSerializer();
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = serializer.serializeToString(okSvgNode);
    } else {
        window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
    }
}

function handleFileSelect(evt) {
    var file = evt.target.files[0];
    if (!file) return;

    // Update UI label (mirrors what index.html used to do)
    var fileNameSpan = document.getElementById("file-name");
    if (fileNameSpan) {
        fileNameSpan.innerHTML =
            '<span style="color: #555;">Chosen file:</span> ' +
            '<span style="color: #000;">' + file.name + '</span>';
    }

    showStatus("Tree is loading", "info");
    var reader = new FileReader();
    reader.onload = function (e) {
        var text = e.target.result || "";
        newickEl.value = text.trim();
        showStatus("Tree successfully loaded", "success");
    };
    reader.readAsText(file);
}

// On DOM ready (Bootstrap 5 / no jQuery)
document.addEventListener("DOMContentLoaded", function () {
    // Help sidebar toggle
    const btnHelp = document.getElementById("btn-help");
    if (btnHelp) {
        btnHelp.addEventListener("click", function () {
            const sb = document.getElementById("help-sidebar");
            if (!sb) return;
            if (sb.style.display === "none" || sb.style.display === "") sb.style.display = "block";
            else sb.style.display = "none";
        });
    }

    const btnHelpClose = document.getElementById("btn-help-close");
    if (btnHelpClose) {
        btnHelpClose.addEventListener("click", function () {
            const sb = document.getElementById("help-sidebar");
            if (sb) sb.style.display = "none";
        });
    }

    // Hide loader (no jQuery)
    const appLoading = document.getElementById("app-loading");
    if (appLoading) appLoading.style.display = "none";

    showStatus("Application loaded. Example tree is prefilled.", "info");

    // Wire controls
    const fileInput = document.getElementById("file-input");
    if (fileInput) fileInput.addEventListener("change", handleFileSelect);

    const btnRun = document.getElementById("btn-run");
    if (btnRun) btnRun.addEventListener("click", function (e) { e.preventDefault(); runPhytClust(); });

    const resetBtn = document.getElementById("btn-extra-reset");
    if (resetBtn) resetBtn.addEventListener("click", function (e) { e.preventDefault(); resetExtraParams(); });

    const helpExampleBtn = document.getElementById("btn-help-load-example");
    if (helpExampleBtn) {
        helpExampleBtn.addEventListener("click", function (e) {
            e.preventDefault();
            showStatus("tree is loading", "info");
            newickEl.value = EXAMPLE_NEWICK;
            showStatus("tree was loaded", "success");
            runPhytClust();
        });
    }

    // Typing feedback when Newick changes
    let loadTimer = null;
    if (newickEl) {
        newickEl.addEventListener("input", function () {
            if (loadTimer) clearTimeout(loadTimer);
            showStatus("tree is loading", "info");
            loadTimer = setTimeout(function () {
                showStatus("tree was loaded", "success");
            }, 250);
        });
    }

    // Extra params icon +/-
    // (Bootstrap 5 fires show.bs.collapse / hide.bs.collapse on the collapse element)
    const extra = document.getElementById("extra-params");
    const icon = document.getElementById("extra-params-icon");
    if (extra && icon) {
        extra.addEventListener("show.bs.collapse", () => { icon.textContent = "−"; });
        extra.addEventListener("hide.bs.collapse", () => { icon.textContent = "+"; });
    }

    // Prefill example and auto-run once for convenience
    if (newickEl) newickEl.value = EXAMPLE_NEWICK;
    runPhytClust();

    // Re-draw optimal k chart when tab is shown (to adapt to width)
    // Bootstrap 5: use the tab button id, and shown.bs.tab event
    const optimalKTabBtn = document.getElementById("viewer-optimal-k-tab");
    if (optimalKTabBtn) {
        optimalKTabBtn.addEventListener("shown.bs.tab", function () {
            if (latestOptimalKData) drawOptimalK(latestOptimalKData);
        });
    }

    // Layout mode selector
    const layoutSelect = document.getElementById("layout-mode");
    if (layoutSelect) {
        layoutSelect.addEventListener("change", function () {
            CURRENT_LAYOUT_MODE = this.value;
            drawTree();
        });
    }

    // Color mode selector
    const colorModeSelect = document.getElementById("color-mode");
    if (colorModeSelect) {
        colorModeSelect.addEventListener("change", function () {
            COLOR_MODE = this.value || "nodes";
            drawTree();
        });
    }

    // ----------------------
    // Palette Selector Logic
    // ----------------------
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
            .map(c => `<span title="${c}" style="display:inline-block;width:12px;height:12px;background:${c};border:1px solid #555;margin-right:2px;"></span>`)
            .join("");
    }

    function applyPalette(type) {
        if (type === "default") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#b84b4b", "#849060", "#3d7c74", "#6e3f8a",
                "#ceb94b", "#3f648a", "#3f408a", "#da63aa"
            );
        } else if (type === "pastel") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#ffb3ba", "#ffdfba", "#ffffba", "#baffc9",
                "#bae1ff", "#d7baff", "#ffcce6", "#c2f0c2"
            );
        } else if (type === "vivid") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                "#ff7f00", "#ffff33", "#a65628", "#f781bf"
            );
        } else if (type === "dark") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#4e79a7", "#59a14f", "#e15759", "#b07aa1",
                "#edc948", "#76b7b2", "#ff9da7", "#9c755f"
            );
        } else if (type === "custom") {
            if (CUSTOM_PALETTE && CUSTOM_PALETTE.length) {
                BASE_COLORS.splice(0, BASE_COLORS.length, ...CUSTOM_PALETTE);
            }
        }

        updatePalettePreview(BASE_COLORS);

        // Recalculate cluster colors + redraw
        if (Object.keys(CURRENT_CLUSTERS || {}).length > 0) {
            const maxCluster = Math.max(...Object.values(CURRENT_CLUSTERS));
            CLUSTER_COLORS = generateClusterColors(maxCluster + 1);
            drawTree();
        }
    }

    const paletteSelect = document.getElementById("palette-select");
    if (paletteSelect) {
        paletteSelect.addEventListener("change", function () {
            if (this.value === "custom") {
                const initial = (CUSTOM_PALETTE && CUSTOM_PALETTE.length)
                    ? CUSTOM_PALETTE.join(",")
                    : BASE_COLORS.join(",");
                const raw = prompt(
                    "Enter custom palette as comma-separated hex colors (e.g. #ff0000,#00ff00,#0000ff):",
                    initial
                );
                if (raw == null) {
                    // user cancelled: revert to default option visually
                    this.value = "default";
                    applyPalette("default");
                    return;
                }

                const parsed = raw.split(",")
                    .map(x => normalizeHexColor(x))
                    .filter(Boolean);

                if (parsed.length < 2) {
                    alert("Please provide at least 2 valid hex colors (like #RRGGBB).");
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

    // leaf/internal name toggles
    const internalNamesCheckbox = document.getElementById("show-internal-names");
    if (internalNamesCheckbox) {
        internalNamesCheckbox.checked = SHOW_INTERNAL_NAMES;
        internalNamesCheckbox.addEventListener("change", function () {
            SHOW_INTERNAL_NAMES = this.checked;
            drawTree();
        });
    }

    const leafNamesCheckbox = document.getElementById("show-leaf-names");
    if (leafNamesCheckbox) {
        leafNamesCheckbox.checked = SHOW_LEAF_NAMES;
        leafNamesCheckbox.addEventListener("change", function () {
            SHOW_LEAF_NAMES = this.checked;
            drawTree();
        });
    }

    // node/label sizes
    const nodeSizeInput = document.getElementById("node-size");
    if (nodeSizeInput) {
        nodeSizeInput.value = LEAF_NODE_RADIUS;
        nodeSizeInput.addEventListener("input", function () {
            const v = parseFloat(this.value);
            if (!isNaN(v) && v >= 1 && v <= 10) {
                LEAF_NODE_RADIUS = v;
                drawTree();
            }
        });
    }

    const labelSizeInput = document.getElementById("label-size");
    if (labelSizeInput) {
        labelSizeInput.value = LABEL_FONT_SIZE;
        labelSizeInput.addEventListener("input", function () {
            const v = parseFloat(this.value);
            if (!isNaN(v) && v >= 6 && v <= 14) {
                LABEL_FONT_SIZE = v;
                drawTree();
            }
        });
    }

    // width/height scale sliders
    const widthScaleInput = document.getElementById("tree-width-scale");
    const heightScaleInput = document.getElementById("tree-height-scale");

    if (widthScaleInput) {
        widthScaleInput.addEventListener("input", function () {
            TREE_WIDTH_SCALE = parseFloat(this.value);
            if (!isNaN(TREE_WIDTH_SCALE)) drawTree();
        });
    }

    if (heightScaleInput) {
        heightScaleInput.addEventListener("input", function () {
            TREE_HEIGHT_SCALE = parseFloat(this.value);
            if (!isNaN(TREE_HEIGHT_SCALE)) drawTree();
        });
    }

    // Fix reset view handler (replace your current broken block with this)
    const btnResetView = document.getElementById("btn-reset-view");
    if (btnResetView) {
        btnResetView.addEventListener("click", function (e) {
            e.preventDefault();

            if (LAST_TREE_SVG && LAST_TREE_ZOOM) {
                LAST_TREE_SVG
                    .transition()
                    .duration(150)
                    .call(LAST_TREE_ZOOM.transform, d3.zoomIdentity);
            }

            const w = document.getElementById("tree-width-scale");
            const h = document.getElementById("tree-height-scale");
            if (w) w.value = "1.0";
            if (h) h.value = "1.0";
            TREE_WIDTH_SCALE = 1.0;
            TREE_HEIGHT_SCALE = 1.0;

            // Reset collapsed nodes
            clearAllCollapsedFlags(NEWICK_RAW_TREE);

            drawTree();
        });
    }
});
