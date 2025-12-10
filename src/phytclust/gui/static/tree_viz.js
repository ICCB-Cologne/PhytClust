let HIER_CART = null;
let HIER_CIRC = null;
var EXAMPLE_NEWICK = "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

var newickEl = document.getElementById("newick-input");
var resultEl = document.getElementById("result");
var statusEl = document.getElementById("status-message");
var treeHost = document.getElementById("tree_display");

let TREE_WIDTH_SCALE = 1.0;
let TREE_HEIGHT_SCALE = 1.0;

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


    let CLUSTER_COLORS = [];
    let CURRENT_CLUSTERS = {};
    let NEWICK_RAW_TREE = null;
    let CURRENT_LAYOUT_MODE = 'rectangular';
    let COLOR_MODE = 'bars';



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

// Simple global tooltip shared by tree + plot
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

function showStatus(message, type) {
    var el = $("#status-message");
    type = type || "info";
    if (statusHideTimeout) { clearTimeout(statusHideTimeout); statusHideTimeout = null; }
    el.removeClass();
    el.addClass("alert alert-" + type);
    el.text(message);
    el.show();
    if (type === "success") {
        statusHideTimeout = setTimeout(function () { el.hide(); }, 5000);
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
        const res = await fetch("/api/run", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload)
        });
        const data = await res.json();
        resultEl.textContent = JSON.stringify(data, null, 2);

        if (!res.ok) {
            showStatus("Error running PhytClust", "danger");
        }

        var logEl = document.getElementById("log-output");
        if (logEl) {
            if (data && data.log) {
                logEl.textContent = data.log;
            } else {
                logEl.textContent = "(no logs)";
            }
        }

        if (data.newick) {
            // Parse tree
            NEWICK_RAW_TREE = parseNewick(data.newick);
            accumulateBranchLength(NEWICK_RAW_TREE);
            computeLayouts();
            // Extract cluster map (may be empty)
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

    // Build SINGLE hierarchy once – all coordinates attach to these objects
    HIER_CART = d3.hierarchy(NEWICK_RAW_TREE, d => d.children);
    HIER_CIRC = d3.hierarchy(NEWICK_RAW_TREE, d => d.children);

    var width  = treeHost.clientWidth  || 800;
    var height = treeHost.clientHeight || 500;

    var margin = { top:20, right:80, bottom:20, left:80 };
    var innerW = width - margin.left - margin.right;
    var innerH = height - margin.top - margin.bottom;
    var radius = Math.min(innerW, innerH) / 2;

    // max branch length
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


function drawTree() {
    clearTree();
    if (!NEWICK_RAW_TREE) return;

    var layoutMode = CURRENT_LAYOUT_MODE || 'rectangular';
    var colorMode = COLOR_MODE || 'bars';
    var container = d3.select("#tree_display");
    var width = treeHost.clientWidth || 800;
    var height = treeHost.clientHeight || 500;
    var margin = { top: 20, right: 80, bottom: 20, left: 80 };

    var svg = container.append("svg")
        .attr("width", width)
        .attr("height", height);

    var zoomLayer = svg.append("g");
    var g = zoomLayer.append("g");

    g.attr("transform",
        (layoutMode === "circular")
            ? `translate(${width/2},${height/2})`
            : `translate(${margin.left},${margin.top})`
    );

    // if (layoutMode === "circular") {
    //     g.attr("transform", "translate(" + (width / 2) + "," + (height / 2) + ")");
    // } else {
    //     g.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    // }

    var zoom = d3.zoom()
        .scaleExtent([0.5, 5])
        .on("zoom", function (event) {
            zoomLayer.attr("transform", event.transform);
        });
    svg.call(zoom);

    var root = (layoutMode === "circular" ? HIER_CIRC : HIER_CART);


    function nodeRadius(d) {
        const isLeaf = !d.children || !d.children.length;
        return isLeaf ? LEAF_NODE_RADIUS : INTERNAL_NODE_RADIUS;
    }
    function nodeFill(d) {
        var name = d.data.name;

        // Internal / collapsed nodes: always black
        if (
            !name ||
            (d.children && d.children.length) ||
            (d.data.children && d.data.children.length && d.data._collapsed)
        ) {
            return "black";
        }

        if (colorMode === "bars") {
            return "black";
        }

        var cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[name] : null;
        if (cid == null || !CLUSTER_COLORS.length) return "black";

        return CLUSTER_COLORS[cid % CLUSTER_COLORS.length];
    }


    if (layoutMode === "circular") {

        function radialPoint(d) {
            var angle = (d._angle || 0) - Math.PI / 2;
            var r = d._radius || 0;
            return [
                (r * Math.cos(angle)) * TREE_WIDTH_SCALE,
                (r * Math.sin(angle)) * TREE_HEIGHT_SCALE
            ];
        }

        var link = g.append("g")
            .selectAll(".tree-link")
            .data(root.links())
            .enter().append("path")
            .attr("class", "tree-link")
            .attr("d", function (d) {
                var s = radialPoint(d.source);
                var t = radialPoint(d.target);
                return "M" + s[0] + "," + s[1] +
                "L" + t[0] + "," + t[1];
            })
            .attr("fill", "none")
            .attr("stroke", d => d.target.color || "black")
            .attr("stroke-width", 1.2);

        var node = g.append("g")
            .selectAll(".tree-node")
            .data(root.descendants())
            .enter().append("g")
            .attr("class", "tree-node")
            .attr("transform", function (d) {
                var p = radialPoint(d);
                return "translate(" + p[0] + "," + p[1] + ")";
            });

        node.append("circle")
            .attr("r", nodeRadius)
            .attr("fill", nodeFill)
            .style("cursor", function (d) {
                return (d.data.children && d.data.children.length) ? "pointer" : "default";
            })
            .on("click", function (event, d) {
                if (d.data.children && d.data.children.length) {
                    d.data._collapsed = !d.data._collapsed;
                    drawTree();
                }
            })
            .on("mouseover", function (event, d) {
                var name = d.data.name || "(internal)";
                var cid = (CURRENT_CLUSTERS && d.data.name) ? CURRENT_CLUSTERS[d.data.name] : null;
                var clusterStr = (cid == null ? "none" : cid.toString());
                var bl = (typeof d.data.length === "number")
                    ? d.data.length.toFixed(4)
                    : "n/a";

                d3Tooltip
                    .style("opacity", 1)
                    .html(
                        "<strong>" + name + "</strong><br/>" +
                        "Depth: " + d.depth + "<br/>" +
                        "Branch length: " + bl + "<br/>" +
                        "Cluster: " + clusterStr
                    );
            })
            .on("mousemove", function (event) {
                d3Tooltip
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY + 10) + "px");
            })
            .on("mouseout", function () {
                d3Tooltip.style("opacity", 0);
            });

        node.append("text")
        .attr("dy", 3)
        .attr("x", function (d) {
            return (d.children && d.children.length) ? -8 : 8;
        })
        .style("text-anchor", function (d) {
            return (d.children && d.children.length) ? "end" : "start";
        })
        .style("font-size", LABEL_FONT_SIZE + "px")
        .text(function (d) {
            const isLeaf = !d.children || d.children.length === 0;

            if (isLeaf && SHOW_LEAF_NAMES) {
                return d.data.name || "";
            }
            if (!isLeaf && SHOW_INTERNAL_NAMES) {
                return d.data.name || "";
            }
            return "";
        });
        return;
    }

    // ========== CLADOGRAM / RECTANGULAR (cartesian) ==============
    // Precompute extents for cartesian layouts
    var allNodes = root.descendants();
    var maxYCart = d3.max(allNodes, d => d._y || 0) || 0;
    var labelColumnX = maxYCart + 30;   // x-position of coloured bars
    var labelWidth = 24;
    var labelHeight = 10;

    var link = g.append("g")
    .selectAll(".tree-link")
    .data(root.links())
    .enter().append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", "black")
    .attr("stroke-width", 1.2)
    .attr("d", function (d) {

        if (layoutMode === "rectangular") {
            return "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
            "V" + (d.target._x * TREE_HEIGHT_SCALE) +
            "H" + (d.target._y * TREE_WIDTH_SCALE);
        } else {
            return "M" + (d.source._y * TREE_WIDTH_SCALE) + "," + (d.source._x * TREE_HEIGHT_SCALE) +
            "L" + (d.target._y * TREE_WIDTH_SCALE) + "," + (d.target._x * TREE_HEIGHT_SCALE);
        }
    });


    var node = g.append("g")
        .selectAll(".tree-node")
        .data(root.descendants())
        .enter().append("g")
        .attr("class", "tree-node")
        .attr("transform", function (d) {
            const x = (d._y || 0) * TREE_WIDTH_SCALE;
            const y = (d._x || 0) * TREE_HEIGHT_SCALE;
            return `translate(${x},${y})`;
        });

    node.append("circle")
        .attr("r", nodeRadius)
        .attr("fill", nodeFill)
        .style("cursor", function (d) {
            return (d.data.children && d.data.children.length) ? "pointer" : "default";
        })
        .on("click", function (event, d) {
            if (d.data.children && d.data.children.length) {
                d.data._collapsed = !d.data._collapsed;
                drawTree();
            }
        })
        .on("mouseover", function (event, d) {
            var name = d.data.name || "(internal)";
            var cid = (CURRENT_CLUSTERS && d.data.name) ? CURRENT_CLUSTERS[d.data.name] : null;
            var clusterStr = (cid == null ? "none" : cid.toString());
            var bl = (typeof d.data.length === "number")
            ? d.data.length.toFixed(4)
            : "n/a";

            d3Tooltip
                .style("opacity", 1)
                .html(
                    "<strong>" + name + "</strong><br/>" +
                    "Depth: " + d.depth + "<br/>" +
                    "Branch length: " + bl + "<br/>" +
                    "Cluster: " + clusterStr
                );
        })
        .on("mousemove", function (event) {
            d3Tooltip
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY + 10) + "px");
        })
        .on("mouseout", function () {
            d3Tooltip.style("opacity", 0);
        });

    node.append("text")
    .attr("dy", 3)
    .attr("x", function (d) {
        return (d.children && d.children.length) ? -8 : 8;
    })
    .style("text-anchor", function (d) {
        return (d.children && d.children.length) ? "end" : "start";
    })
    .style("font-size", LABEL_FONT_SIZE + "px")
    .text(function (d) {
        const isLeaf = !d.children || d.children.length === 0;

        if (isLeaf && SHOW_LEAF_NAMES) {
            return d.data.name || "";
        }
        if (!isLeaf && SHOW_INTERNAL_NAMES) {
            return d.data.name || "";
        }
        return "";
    });


    if (colorMode === "bars") {

        // All leaves, sorted top-to-bottom
        const leafNodes = allNodes
            .filter(d => !d.children || !d.children.length)
            .sort((a, b) => a._x - b._x)

        const labelsG = g.append("g").attr("class", "cluster-bars");

        // Compute average spacing between consecutive leaves (for extrapolation at ends)
        const spacings = [];
        for (let i = 0; i < leafNodes.length - 1; i++) {
            const y1 = leafNodes[i]._x || 0;
            const y2 = leafNodes[i + 1]._x || 0;
            spacings.push(y2 - y1);
        }
        const meanSpacing = spacings.length
            ? spacings.reduce((a, b) => a + b, 0) / spacings.length
            : 12;

        // Map: cluster id -> {firstIndex, lastIndex}
        const clusterInfo = new Map();
        leafNodes.forEach((d, idx) => {
            const cid = CURRENT_CLUSTERS ? CURRENT_CLUSTERS[d.data.name] : null;
            if (cid == null) return;

            if (!clusterInfo.has(cid)) {
                clusterInfo.set(cid, { firstIndex: idx, lastIndex: idx });
            } else {
                clusterInfo.get(cid).lastIndex = idx;
            }
        });
        const usedLabelYs = [];
        function findNonOverlappingY(y0) {
            let y = y0;
            const minGap = 10; // px between labels
            let safety = 0;
            while (usedLabelYs.some(u => Math.abs(u - y) < minGap) && safety < 50) {
                y += minGap;
                safety++;
            }
            usedLabelYs.push(y);
            return y;
        }
        // Draw one contiguous bar per cluster
        clusterInfo.forEach((info, cid) => {
            const first = leafNodes[info.firstIndex];
            const last  = leafNodes[info.lastIndex];

            const yFirst = first._x || 0;
            const yLast  = last._x || 0;

            // Top boundary: midpoint with previous cluster block (or extrapolated)
            let top;
            if (info.firstIndex === 0) {
                top = yFirst - meanSpacing / 2;
            } else {
                const prevLeaf = leafNodes[info.firstIndex - 1];
                const yPrev = prevLeaf._x || 0;
                top = (yPrev + yFirst) / 2;
            }

            // Bottom boundary: midpoint with next cluster block (or extrapolated)
            let bottom;
            if (info.lastIndex === leafNodes.length - 1) {
                bottom = yLast + meanSpacing / 2;
            } else {
                const nextLeaf = leafNodes[info.lastIndex + 1];
                const yNext = nextLeaf._x || 0;
                bottom = (yLast + yNext) / 2;
            }
            const labelY = findNonOverlappingY((top + bottom) / 2);

            const color = (cid == null || !CLUSTER_COLORS.length)
                ? "lightgrey"
                : CLUSTER_COLORS[cid % CLUSTER_COLORS.length];

            labelsG.append("rect")
            .attr("x", labelColumnX * TREE_WIDTH_SCALE)
            .attr("y", top * TREE_HEIGHT_SCALE)
            .attr("width", labelWidth)
            .attr("height", (bottom - top) * TREE_HEIGHT_SCALE)
            .attr("fill", color);

                const totalLeaves = leafNodes.length;
                const clusterSize = info.lastIndex - info.firstIndex + 1;

                // Only label clusters ≥ 10% of total leaves
                if (clusterSize / totalLeaves >= 0.10) {
                    labelsG.append("text")
                        .attr("x", (labelColumnX + labelWidth + 4) * TREE_WIDTH_SCALE)
                        .attr("y", labelY * TREE_HEIGHT_SCALE)
                        .attr("dominant-baseline", "middle")
                        .attr("font-size", 10)
                        .text("C" + cid);
                }
        });
        // --------------------------------------------------
        //  BRANCH-LENGTH AXIS (for cladogram/rectangular)
        // --------------------------------------------------
        if (layoutMode !== "circular") {

            // compute horizontal scale (branch length → y position)
            const allNodes = root.descendants();
            const maxY = d3.max(allNodes, d => d._y || 0) || 0;

            const blVals = allNodes.map(d => d.data._bl || 0);
            const maxBl = d3.max(blVals) || 1;

            const blScale = d3.scaleLinear()
                .domain([0, maxBl])
                .range([0, maxY]);

            const axisGroup = g.append("g")
                .attr("class", "branch-length-axis")
                .attr("transform", `translate(0, ${height - margin.bottom - margin.top - 5})`)
                .call(d3.axisBottom(blScale).ticks(5));

            g.append("text")
                .attr("x", maxY / 2)
                .attr("y", height - margin.bottom - margin.top + 15)
                .attr("text-anchor", "middle")
                .attr("font-size", 10)
                .text("Branch length");
        }
    }
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
        return;
    }

    plotEl.innerHTML = ""; // clear previous

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

    // Wire up axis mode change
    if (axisModeEl && !axisModeEl.__wired) {
        axisModeEl.__wired = true;
        axisModeEl.addEventListener('change', function() {
            axisModeEl.__userOverride = true;
            drawOptimalK(data);
        });
    }
}

function handleFileSelect(evt) {
    var file = evt.target.files[0];
    if (!file) return;
    showStatus("Tree is loading", "info");
    var reader = new FileReader();
    reader.onload = function (e) {
        var text = e.target.result || "";
        newickEl.value = text.trim();
        showStatus("Tree successfully loaded", "success");
    };
    reader.readAsText(file);
}

// On DOM ready
$(function () {
    // Help sidebar toggle
    document.getElementById('btn-help').addEventListener('click', function () {
        var sb = document.getElementById('help-sidebar');
        if (sb.style.display === 'none' || sb.style.display === '') { sb.style.display = 'block'; }
        else { sb.style.display = 'none'; }
    });
    document.getElementById('btn-help-close').addEventListener('click', function () {
        document.getElementById('help-sidebar').style.display = 'none';
    });
    $("#app-loading").hide();
    showStatus("Application loaded. Example tree is prefilled.", "info");

    // Wire controls
    document.getElementById("file-input").addEventListener("change", handleFileSelect);
    document.getElementById("btn-run").addEventListener("click", function (e) { e.preventDefault(); runPhytClust(); });

    document.getElementById("btn-save").addEventListener("click", async function (e) {
        e.preventDefault();
        if (isRunning) { showStatus('Cannot save while PhytClust is running.', 'warning'); return; }
        var fname = prompt("Enter filename (tsv):", "phytclust_results.tsv");
        if (!fname) { fname = "phytclust_results.tsv"; }

        if (window.showDirectoryPicker) {
            try {
                const dirHandle = await window.showDirectoryPicker();
                const res = await fetch('/api/export_tsv', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        filename: fname,
                        top_n: (extraTopNEl && extraTopNEl.value ? parseInt(extraTopNEl.value, 10) || 1 : 1),
                        outlier: true,
                        output_all: false
                    })
                });
                if (!res.ok) {
                    const errData = await res.json().catch(() => ({}));
                    const msg = (errData && errData.detail) ? errData.detail : 'Export failed';
                    showStatus(msg, 'danger');
                    return;
                }
                const text = await res.text();
                const fileHandle = await dirHandle.getFileHandle(fname, { create: true });
                const writable = await fileHandle.createWritable();
                await writable.write(text);
                await writable.close();
                showStatus('Saved ' + fname + ' to selected directory', 'success');
            } catch (err) {
                if (err && err.name === 'AbortError') { return; }
                try {
                    var dirInput = document.getElementById('output-dir');
                    var dir = dirInput && dirInput.value ? dirInput.value.trim() : 'results';
                    const res2 = await fetch('/api/save', {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({
                            results_dir: dir,
                            filename: fname,
                            top_n: (extraTopNEl && extraTopNEl.value ? parseInt(extraTopNEl.value, 10) || 1 : 1),
                            outlier: true,
                            output_all: false
                        })
                    });
                    const data2 = await res2.json();
                    if (!res2.ok) {
                        var msg2 = (data2 && data2.detail) ? data2.detail : 'Save failed';
                        showStatus(msg2, 'danger');
                    } else {
                        showStatus('Saved to ' + data2.results_dir + '/' + data2.filename, 'success');
                    }
                } catch (err2) {
                    showStatus('Save error: ' + err2, 'danger');
                }
            }
        } else {
            try {
                var dirInputF = document.getElementById('output-dir');
                var dirF = dirInputF && dirInputF.value ? dirInputF.value.trim() : 'results';
                const res3 = await fetch('/api/save', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        results_dir: dirF,
                        filename: fname,
                        top_n: (extraTopNEl && extraTopNEl.value ? parseInt(extraTopNEl.value, 10) || 1 : 1),
                        outlier: true,
                        output_all: false
                    })
                });
                const data3 = await res3.json();
                if (!res3.ok) {
                    var msg3 = (data3 && data3.detail) ? data3.detail : 'Save failed';
                    showStatus(msg3, 'danger');
                } else {
                    showStatus('Saved to ' + data3.results_dir + '/' + data3.filename, 'success');
                }
            } catch (err3) {
                showStatus('Save error: ' + err3, 'danger');
            }
        }
    });

    var resetBtn = document.getElementById("btn-extra-reset");
    if (resetBtn) {
        resetBtn.addEventListener("click", function (e) { e.preventDefault(); resetExtraParams(); });
    }

    var helpExampleBtn = document.getElementById('btn-help-load-example');
    if (helpExampleBtn) {
        helpExampleBtn.addEventListener('click', function (e) {
            e.preventDefault();
            showStatus("tree is loading", "info");
            newickEl.value = EXAMPLE_NEWICK;
            showStatus("tree was loaded", "success");
            runPhytClust();
        });
    }

    var loadTimer = null;
    newickEl.addEventListener('input', function () {
        if (loadTimer) { clearTimeout(loadTimer); }
        showStatus("tree is loading", "info");
        loadTimer = setTimeout(function () {
            showStatus("tree was loaded", "success");
        }, 250);
    });

    $('#extra-params').on('shown.bs.collapse', function () {
        $('#extra-params-icon').removeClass('glyphicon-plus').addClass('glyphicon-minus');
    });
    $('#extra-params').on('hidden.bs.collapse', function () {
        $('#extra-params-icon').removeClass('glyphicon-minus').addClass('glyphicon-plus');
    });

    // Prefill example and auto-run once for convenience
    newickEl.value = EXAMPLE_NEWICK;
    runPhytClust();

    // Re-draw optimal k chart when tab is shown (to adapt to width)
    $('a[href="#viewer-optimal-k"]').on('shown.bs.tab', function () {
        if (latestOptimalKData) {
            drawOptimalK(latestOptimalKData);
        }
    });
    var layoutSelect = document.getElementById('layout-mode');
    if (layoutSelect) {
        layoutSelect.addEventListener('change', function () {
            CURRENT_LAYOUT_MODE = this.value;
            drawTree(); // re-render current tree with new layout
        });
    const colorModeSelect = document.getElementById("color-mode");
    if (colorModeSelect) {
        colorModeSelect.addEventListener("change", function () {
            COLOR_MODE = this.value || 'nodes';
            drawTree();   // re-render tree with new colouring mode
        });
    }

    // ----------------------
    // Palette Selector Logic
    // ----------------------

    let CUSTOM_PALETTE = null;

    // Show the colors in the small preview bar
    function updatePalettePreview(colors) {
        const box = document.getElementById("palette-preview");
        if (!box) return;
        box.innerHTML = colors
            .map(c => `<span style="display:inline-block;width:12px;height:12px;background:${c};
        border:1px solid #555;margin-right:2px;"></span>`)
            .join("");
    }

    // Apply a chosen palette to BASE_COLORS
    function applyPalette(type) {
        if (type === "default") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#b84b4b", "#849060", "#3d7c74", "#6e3f8a",
                "#ceb94b", "#3f648a", "#3f408a", "#da63aa"
            );
        }

        if (type === "pastel") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#ffb3ba", "#ffdfba", "#ffffba", "#baffc9",
                "#bae1ff", "#d7baff", "#ffcce6", "#c2f0c2"
            );
        }

        if (type === "vivid") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                "#ff7f00", "#ffff33", "#a65628", "#f781bf"
            );
        }

        if (type === "dark") {
            BASE_COLORS.splice(0, BASE_COLORS.length,
                "#4e79a7",
                "#59a14f",
                "#e15759",
                "#b07aa1",
                "#edc948",
                "#76b7b2",
                "#ff9da7",
                "#9c755f"

            );
        }

        if (type === "custom" && CUSTOM_PALETTE) {
            BASE_COLORS.splice(0, BASE_COLORS.length, ...CUSTOM_PALETTE);
        }

        updatePalettePreview(BASE_COLORS);

        // Recalculate cluster colors + redraw
        if (Object.keys(CURRENT_CLUSTERS).length > 0) {
            const maxCluster = Math.max(...Object.values(CURRENT_CLUSTERS));
            CLUSTER_COLORS = generateClusterColors(maxCluster + 1);
            drawTree();
        }
    }

    // Dropdown change listener
    const paletteSelect = document.getElementById("palette-select");
    if (paletteSelect) {
        paletteSelect.addEventListener("change", function () {
            if (this.value === "custom") {
                alert("Custom palette not implemented yet.");
                return;
            }
            applyPalette(this.value);
        });
    }

    updatePalettePreview(BASE_COLORS);

    const internalNamesCheckbox = document.getElementById("show-internal-names");
    if (internalNamesCheckbox) {
        internalNamesCheckbox.checked = SHOW_INTERNAL_NAMES;
        internalNamesCheckbox.addEventListener("change", function () {
            SHOW_INTERNAL_NAMES = this.checked;
            drawTree();   // redraw tree
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
    const nodeSizeInput = document.getElementById("node-size");
    if (nodeSizeInput) {
        nodeSizeInput.value = LEAF_NODE_RADIUS;
        nodeSizeInput.addEventListener("input", function () {
            const v = parseFloat(this.value);
            if (!isNaN(v) && v >= 1 && v <= 10) {
                LEAF_NODE_RADIUS = v;  // only terminal nodes
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
    const widthScaleInput  = document.getElementById("tree-width-scale");
    const heightScaleInput = document.getElementById("tree-height-scale");

    widthScaleInput.addEventListener("input", function () {
        TREE_WIDTH_SCALE = parseFloat(this.value);
        drawTree();
    });

    heightScaleInput.addEventListener("input", function () {
        TREE_HEIGHT_SCALE = parseFloat(this.value);
        drawTree();
    });
    }
});
