/* ============================================================
   PhytClust – tree/draw.js
   Tree layout computation, drawing (rectangular + circular),
   drawTreeInto (comparison), drawComparisonBarsInto,
   context menu helpers, and node collapse logic.
   ============================================================ */

import { state, LABEL_PAD, SELECTED_CLUSTER_IDS, BOX_LABEL_MAP, BOX_ADJUST_MAP, NODE_CUSTOM, hasClusterFocus, getBoxAdjust } from "../state.js";
import { treeHost } from "../dom.js";
import { getVisibleChildren, getAllChildren, collectLeafNamesData } from "../utils.js";
import { generateClusterColors, getThemeColors } from "../colors.js";

/* ─────────────────────────────────────
   Tree Layout Computation
   ───────────────────────────────────── */
export function accumulateBranchLength(node, length) {
  length = length || 0;
  node._bl = length;
  if (node.children) {
    for (const c of node.children)
      accumulateBranchLength(c, length + (c.length || 0));
  }
}

export function computeLayouts() {
  if (!state.NEWICK_RAW_TREE) return;
  state.HIER_CART = d3.hierarchy(state.NEWICK_RAW_TREE, getAllChildren);
  state.HIER_CIRC = d3.hierarchy(state.NEWICK_RAW_TREE, getAllChildren);

  var width = treeHost.clientWidth || 800;
  var height = treeHost.clientHeight || 500;
  var margin = { top: 20, right: 80, bottom: 20, left: 80 };
  var innerW = width - margin.left - margin.right;
  var innerH = height - margin.top - margin.bottom;
  var radius = Math.min(innerW, innerH) / 2;

  var maxBl = d3.max(state.HIER_CART.descendants(), (d) => d.data._bl || 0) || 1;
  var blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);
  var blToR = d3.scaleLinear().domain([0, maxBl]).range([0, radius]);

  d3.cluster().size([innerH, 1])(state.HIER_CART);
  state.HIER_CART.each((d) => {
    d._x = d.x;
    d._y = blToX(d.data._bl || 0);
  });

  d3.cluster().size([2 * Math.PI, 1])(state.HIER_CIRC);
  state.HIER_CIRC.each((d) => {
    d._angle = d.x;
    d._radius = blToR(d.data._bl || 0);
  });
}

/* ─────────────────────────────────────
   Node Context Menu
   ───────────────────────────────────── */
export function getNodeCustom(dataNode) {
  if (!NODE_CUSTOM.has(dataNode))
    NODE_CUSTOM.set(dataNode, {
      highlighted: false,
      radiusScale: 1,
      renamedTo: null,
      lockOriginalSize: false,
    });
  return NODE_CUSTOM.get(dataNode);
}

export function nodeDisplayName(dataNode) {
  const c = NODE_CUSTOM.has(dataNode) ? NODE_CUSTOM.get(dataNode) : null;
  return c && c.renamedTo != null ? c.renamedTo : dataNode.name || "";
}

export function countLeafDescendantsData(node) {
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
      ? state.INTERNAL_NODE_RADIUS
      : state.LEAF_NODE_RADIUS;
  const c = NODE_CUSTOM.has(d.data) ? NODE_CUSTOM.get(d.data) : null;

  if (c && c.lockOriginalSize) {
    return base;
  }

  if (
    d.data &&
    d.data._collapsed &&
    d.data.children &&
    d.data.children.length
  ) {
    var leafCount = countLeafDescendantsData(d.data);
    if (leafCount > 1) {
      var autoCollapsed =
        state.INTERNAL_NODE_RADIUS + state.LEAF_NODE_RADIUS * 0.55 * Math.log2(leafCount);
      base = Math.max(base, Math.min(14, autoCollapsed));
    }
  }

  return base * (c ? c.radiusScale : 1);
}

export function representativeClusterIdFromData(dataNode) {
  if (!dataNode || !state.CURRENT_CLUSTERS) return null;
  const leaves = [];
  collectLeafNamesData(dataNode, leaves);
  if (!leaves.length) return null;
  let cid = null;
  for (const nm of leaves) {
    const v = state.CURRENT_CLUSTERS[nm];
    if (v == null) return null;
    if (cid == null) cid = v;
    else if (cid !== v) return null;
  }
  return cid;
}

export function showNodeContextMenu(event, d) {
  event.preventDefault();
  event.stopPropagation();
  state.CTX_TARGET_DATA = d.data;
  const menu = document.getElementById("node-context-menu");
  if (!menu) return;
  menu.style.display = "block";
  menu.style.left = event.clientX + "px";
  menu.style.top = event.clientY + "px";
  var ctxSlider = document.getElementById("ctx-size-slider");
  var ctxSliderVal = document.getElementById("ctx-size-value");
  var c = NODE_CUSTOM.has(d.data)
    ? NODE_CUSTOM.get(d.data)
    : { radiusScale: 1 };
  if (ctxSlider) ctxSlider.value = c.radiusScale;
  if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);

  var isLeaf = !(d.data.children && d.data.children.length);
  var hasClusters =
    state.CURRENT_CLUSTERS && Object.keys(state.CURRENT_CLUSTERS).length > 0;
  var isBoxMode = state.COLOR_MODE === "boxes";

  var copyBtn = document.getElementById("ctx-copy-subtree");
  var originalSizeBtn = document.getElementById("ctx-original-size");
  var clusterBtn = document.getElementById("ctx-select-cluster");
  var renameClusterBtn = document.getElementById("ctx-rename-cluster");
  var adjustBoxBtn = document.getElementById("ctx-adjust-box");
  var c = getNodeCustom(d.data);

  if (copyBtn) copyBtn.style.display = isLeaf ? "none" : "";
  if (originalSizeBtn) {
    originalSizeBtn.style.display = isLeaf ? "none" : "";
    originalSizeBtn.textContent = c.lockOriginalSize
      ? "Use auto collapsed size"
      : "Keep original size";
  }
  if (clusterBtn) clusterBtn.style.display = hasClusters ? "" : "none";
  if (renameClusterBtn)
    renameClusterBtn.style.display = isBoxMode && hasClusters ? "" : "none";
  if (adjustBoxBtn)
    adjustBoxBtn.style.display = isBoxMode && hasClusters ? "" : "none";
}

export function hideNodeContextMenu() {
  const menu = document.getElementById("node-context-menu");
  if (menu) menu.style.display = "none";
  state.CTX_TARGET_DATA = null;
}

/* ─────────────────────────────────────
   Search highlighting helper
   ───────────────────────────────────── */
function isSearchMatch(dataNode) {
  if (!state.SEARCH_TERM) return false;
  const name = nodeDisplayName(dataNode).toLowerCase();
  return name.includes(state.SEARCH_TERM.toLowerCase());
}

/* ─────────────────────────────────────
   Collapse helpers
   ───────────────────────────────────── */
export function clearTree() {
  treeHost.innerHTML = "";
}

export function clearAllCollapsedFlags(node) {
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
  clearDescendantCollapsedFlags(dataNode);
}

/* ─────────────────────────────────────
   Tree Drawing
   ───────────────────────────────────── */
export function drawTree() {
  clearTree();
  if (!state.NEWICK_RAW_TREE) return;
  computeLayouts();

  const layoutMode = state.CURRENT_LAYOUT_MODE || "rectangular";
  const colorMode = state.COLOR_MODE || "bars";
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
      state.LAST_TREE_TRANSFORM = event.transform;
    });
  const previousTransform = state.LAST_TREE_TRANSFORM || d3.zoomIdentity;
  state.LAST_TREE_SVG = svg;
  state.LAST_TREE_ZOOM = zoom;

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
    if (cid == null || !state.CLUSTER_COLORS.length) return null;
    if (Number(cid) < 0) return tc.internal;
    return state.CLUSTER_COLORS[cid % state.CLUSTER_COLORS.length];
  }

  function nodeFill(d) {
    const name = d.data && d.data.name;
    if (colorMode === "bars") return tc.label;
    if (name && !hasOriginalChildren(d)) {
      const cid = state.CURRENT_CLUSTERS ? state.CURRENT_CLUSTERS[name] : null;
      return clusterColor(cid) || tc.label;
    }
    const rep = representativeClusterIdFromData(d.data);
    return clusterColor(rep) || tc.label;
  }

  function nodeClusterId(d) {
    const name = d.data && d.data.name;
    if (name && !hasOriginalChildren(d)) {
      const cid = state.CURRENT_CLUSTERS ? state.CURRENT_CLUSTERS[name] : null;
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
    const root = state.HIER_CIRC;
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
        rr * Math.cos(a) * state.TREE_WIDTH_SCALE,
        rr * Math.sin(a) * state.TREE_HEIGHT_SCALE,
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
      .attr("stroke-width", state.BRANCH_STROKE_WIDTH)
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
      .attr("font-size", state.LABEL_FONT_SIZE + "px")
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
      .text((d) => (state.SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : ""));

    node
      .filter((d) => !isVisibleLeaf(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", state.LABEL_FONT_SIZE - 1 + "px")
      .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.internal))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
      .attr("transform", (d) => {
        const a = d._angle || 0;
        let deg = (a * 180) / Math.PI - 90;
        return `rotate(${deg}) translate(0,${-(nodeRadius(d) + 5)})`;
      })
      .style("text-anchor", "middle")
      .text((d) => (state.SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : ""));

    node
      .filter((d) => d.data && d.data._collapsed && hasOriginalChildren(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", Math.max(8, state.LABEL_FONT_SIZE - 2) + "px")
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
    if (state.CURRENT_CLUSTERS && state.CLUSTER_COLORS.length) {
      const visLeaves = root
        .descendants()
        .filter((d) => isNodeVisible(d))
        .filter((d) => isVisibleLeaf(d))
        .sort((a, b) => (a._angle || 0) - (b._angle || 0));

      function visLeafCid(d) {
        const nm = d.data && d.data.name;
        const isOrigLeaf = nm && !(d.data.children && d.data.children.length);
        if (isOrigLeaf) return state.CURRENT_CLUSTERS[nm];
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
        const estRadialLabelWidth = state.SHOW_LEAF_NAMES
          ? leafLabelChars * state.LABEL_FONT_SIZE * 0.55 +
            LABEL_PAD +
            state.LEAF_NODE_RADIUS
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
          if (!state.SHOW_OUTLIER_BOXES && Number(cid) < 0) return;
          if (!circClusterGroups[cid]) circClusterGroups[cid] = [];
          circClusterGroups[cid].push(d);
        });
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

          var angPad = state.BOX_PAD_V * 0.01;
          var minAngle = leafAngles[0] - angPad;
          var maxAngle = leafAngles[leafAngles.length - 1] + angPad;
          if (cLeaves.length === 1) {
            var spread = 0.03;
            minAngle -= spread;
            maxAngle += spread;
          }

          var innerR = rScale(mrcaR) - state.BOX_PAD_H;
          var outerR =
            rScale(maxLeafR) + state.BOX_PAD_H + (state.SHOW_LEAF_NAMES ? 40 : 5);
          if (innerR < 0) innerR = 0;

          var isOutlierWedge = Number(cid) < 0;
          var wedgeColor = isOutlierWedge
            ? tc.internal
            : clusterColor(cid) || tc.internal;
          var wedgeArc = d3
            .arc()
            .innerRadius(innerR * state.TREE_WIDTH_SCALE)
            .outerRadius(outerR * state.TREE_WIDTH_SCALE)
            .startAngle(minAngle)
            .endAngle(maxAngle)
            .cornerRadius(state.BOX_CORNER_RADIUS);

          var wedge = boxesG
            .append("path")
            .attr("d", wedgeArc)
            .attr("fill", isOutlierWedge ? "none" : wedgeColor)
            .attr("opacity", hasClusterFocus(cid) ? state.BOX_ALPHA : 0.02)
            .attr("stroke", wedgeColor)
            .attr("stroke-width", 1.5)
            .attr(
              "stroke-opacity",
              hasClusterFocus(cid) ? Math.min(state.BOX_ALPHA * 3.5, 0.8) : 0.08,
            );
          if (isOutlierWedge) wedge.attr("stroke-dasharray", "5,3");

          if (state.SHOW_BOX_LABELS) {
            var midAngle = (minAngle + maxAngle) / 2;
            var labelR = outerR * state.TREE_WIDTH_SCALE + 4;
            var lx = labelR * Math.cos(midAngle - Math.PI / 2);
            var ly = labelR * Math.sin(midAngle - Math.PI / 2);
            var deg = (midAngle * 180) / Math.PI - 90;
            var flip = deg > 90 || deg < -90;
            boxesG
              .append("text")
              .attr("x", lx)
              .attr("y", ly)
              .attr("dy", "0.35em")
              .attr("font-size", state.CLUSTER_LABEL_FONT_SIZE + "px")
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
  const root = state.HIER_CART;
  const allNodes = root.descendants().filter((d) => isNodeVisible(d));
  const allLinks = root
    .links()
    .filter((l) => isNodeVisible(l.source) && isNodeVisible(l.target));
  const maxYCart = d3.max(allNodes, (d) => d._y || 0) || 0;
  const leafNames = allNodes
    .filter((d) => !d.children || !d.children.length)
    .map((d) => (d.data.name || "").length);
  const maxLabelChars = d3.max(leafNames) || 0;
  const estLabelWidth = state.SHOW_LEAF_NAMES
    ? maxLabelChars * state.LABEL_FONT_SIZE * 0.6 + LABEL_PAD + state.LEAF_NODE_RADIUS
    : 20;
  const labelColumnX = maxYCart * state.TREE_WIDTH_SCALE + estLabelWidth + 12;

  g.append("g")
    .selectAll(".tree-link")
    .data(allLinks)
    .enter()
    .append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", tc.branch)
    .attr("stroke-opacity", 1)
    .attr("stroke-width", state.BRANCH_STROKE_WIDTH)
    .attr("d", function (d) {
      if (layoutMode === "rectangular") {
        return (
          "M" +
          d.source._y * state.TREE_WIDTH_SCALE +
          "," +
          d.source._x * state.TREE_HEIGHT_SCALE +
          "V" +
          d.target._x * state.TREE_HEIGHT_SCALE +
          "H" +
          d.target._y * state.TREE_WIDTH_SCALE
        );
      }
      return (
        "M" +
        d.source._y * state.TREE_WIDTH_SCALE +
        "," +
        d.source._x * state.TREE_HEIGHT_SCALE +
        "L" +
        d.target._y * state.TREE_WIDTH_SCALE +
        "," +
        d.target._x * state.TREE_HEIGHT_SCALE
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
        `translate(${(d._y || 0) * state.TREE_WIDTH_SCALE},${(d._x || 0) * state.TREE_HEIGHT_SCALE})`,
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
    .style("font-size", state.LABEL_FONT_SIZE + "px")
    .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.label))
    .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
    .text((d) => (state.SHOW_LEAF_NAMES ? nodeDisplayName(d.data) : ""));

  // Internal node labels
  node
    .filter((d) => !!(d.children && d.children.length))
    .append("text")
    .attr("dy", (d) => nodeRadius(d) + 12)
    .attr("x", (d) => nodeRadius(d) + 3)
    .style("text-anchor", "start")
    .style("font-size", state.LABEL_FONT_SIZE - 1 + "px")
    .style("font-style", "italic")
    .attr("fill", (d) => (isHighlighted(d) ? highlightColor() : tc.internal))
    .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
    .text((d) => (state.SHOW_INTERNAL_NAMES ? nodeDisplayName(d.data) : ""));

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
    .style("font-size", Math.max(8, state.LABEL_FONT_SIZE - 2) + "px")
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
    const ys = leafNodes.map((d) => (d._x || 0) * state.TREE_HEIGHT_SCALE);

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

  // Cluster boxes (MRCA mode)
  if (colorMode === "boxes" && state.CURRENT_CLUSTERS && state.CLUSTER_COLORS.length) {
    const clusterGroups = {};
    allNodes.forEach((d) => {
      if (!isVisibleLeaf(d)) return;
      const cid = nodeClusterId(d);
      if (cid == null) return;
      if (!state.SHOW_OUTLIER_BOXES && Number(cid) < 0) return;
      if (!clusterGroups[cid]) clusterGroups[cid] = [];
      clusterGroups[cid].push(d);
    });
    function findMRCA(nodes) {
      if (!nodes.length) return null;
      if (nodes.length === 1) return nodes[0];
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
      return common[0];
    }
    var estLabelPx = state.SHOW_LEAF_NAMES
      ? maxLabelChars * state.LABEL_FONT_SIZE * 0.6 +
        LABEL_PAD +
        state.LEAF_NODE_RADIUS +
        16
      : 16;
    var boxPadV = state.BOX_PAD_V;
    var boxPadL = state.BOX_PAD_H;
    var boxesG = g.append("g").attr("class", "cluster-boxes");
    for (var cid of Object.keys(clusterGroups)) {
      var leaves = clusterGroups[cid];
      var mrca = findMRCA(leaves);
      if (!mrca) continue;
      var ys = leaves.map(function (d) {
        return (d._x || 0) * state.TREE_HEIGHT_SCALE;
      });
      var minY = Math.min.apply(null, ys) - boxPadV;
      var maxY = Math.max.apply(null, ys) + boxPadV;
      var mrcaX = (mrca._y || 0) * state.TREE_WIDTH_SCALE - boxPadL;
      var maxLeafX = Math.max.apply(
        null,
        leaves.map(function (d) {
          return (d._y || 0) * state.TREE_WIDTH_SCALE;
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
        : state.CLUSTER_COLORS[cid % state.CLUSTER_COLORS.length];
      var rect = boxesG
        .append("rect")
        .attr("x", mrcaX)
        .attr("y", minY)
        .attr("width", Math.max(boxW, 8))
        .attr("height", Math.max(maxY - minY, 4))
        .attr("rx", state.BOX_CORNER_RADIUS)
        .attr("ry", state.BOX_CORNER_RADIUS)
        .attr("fill", isOutlier ? "none" : color)
        .attr("opacity", hasClusterFocus(cid) ? state.BOX_ALPHA : 0.02)
        .attr("stroke", color)
        .attr("stroke-width", isOutlier ? 1.5 : 1.5)
        .attr(
          "stroke-opacity",
          hasClusterFocus(cid) ? Math.min(state.BOX_ALPHA * 3.5, 0.8) : 0.08,
        );
      if (isOutlier) {
        rect.attr("stroke-dasharray", "5,3");
      }
      if (state.SHOW_BOX_LABELS) {
        boxesG
          .append("text")
          .attr("x", boxRight + 4)
          .attr("y", (minY + maxY) / 2)
          .attr("dy", "0.35em")
          .attr("font-size", state.CLUSTER_LABEL_FONT_SIZE + "px")
          .attr("font-weight", "600")
          .attr("fill", color)
          .attr("font-style", isOutlier ? "italic" : "normal")
          .attr("opacity", hasClusterFocus(cid) ? 1 : 0.25)
          .text(BOX_LABEL_MAP[cid] || (isOutlier ? "outlier" : "C" + cid));
      }
    }
    boxesG.attr("pointer-events", "none");
    boxesG.lower();
  }

  // Branch-length axis
  var axisMaxX = d3.max(allNodes, (d) => d._x || 0) || 0;
  var axisMaxY = d3.max(allNodes, (d) => d._y || 0) || 0;
  var axisMaxBl = d3.max(allNodes, (d) => d.data._bl || 0) || 1;
  var axisYPos = axisMaxX * state.TREE_HEIGHT_SCALE + 20;
  var blScale = d3
    .scaleLinear()
    .domain([0, axisMaxBl])
    .range([0, axisMaxY * state.TREE_WIDTH_SCALE]);

  var axisG = g
    .append("g")
    .attr("class", "branch-length-axis")
    .attr("transform", "translate(0, " + axisYPos + ")")
    .call(d3.axisBottom(blScale).ticks(5));
  axisG
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", state.AXIS_FONT_SIZE + "px");
  axisG.selectAll("line").attr("stroke", tc.branch);
  axisG.selectAll("path").attr("stroke", tc.branch);
  axisG
    .append("text")
    .attr("x", (axisMaxY * state.TREE_WIDTH_SCALE) / 2)
    .attr("y", state.AXIS_FONT_SIZE + 20)
    .attr("text-anchor", "middle")
    .attr("font-size", state.AXIS_FONT_SIZE)
    .attr("fill", tc.internal)
    .text("Branch length");

  svg.attr("height", Math.max(height, axisYPos + 50 + margin.bottom));
  applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform);
}

/* ─────────────────────────────────────
   Draw tree into a specific container
   (used for comparison mode)
   ───────────────────────────────────── */
export function drawTreeInto(hostEl, clusters) {
  hostEl.innerHTML = "";
  if (!state.NEWICK_RAW_TREE) return;

  const tc = getThemeColors();
  const hier = d3.hierarchy(state.NEWICK_RAW_TREE, getVisibleChildren);
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
    .attr("stroke-width", state.BRANCH_STROKE_WIDTH)
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
export function drawComparisonBarsInto(hostEl, comparisons) {
  hostEl.innerHTML = "";
  if (!state.NEWICK_RAW_TREE) return;

  const tc = getThemeColors();
  const hier = d3.hierarchy(state.NEWICK_RAW_TREE, getVisibleChildren);
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
  const estLabelWidth = state.SHOW_LEAF_NAMES
    ? maxLabelChars * state.LABEL_FONT_SIZE * 0.6 + LABEL_PAD + state.LEAF_NODE_RADIUS
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
    .attr("stroke-width", state.BRANCH_STROKE_WIDTH)
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
      d.children && d.children.length ? state.INTERNAL_NODE_RADIUS : state.LEAF_NODE_RADIUS,
    )
    .attr("fill", tc.label)
    .attr("stroke", "none");

  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", state.LEAF_NODE_RADIUS + LABEL_PAD)
    .style("text-anchor", "start")
    .style("font-size", state.LABEL_FONT_SIZE + "px")
    .attr("fill", tc.label)
    .text((d) => (state.SHOW_LEAF_NAMES ? d.data.name || "" : ""));

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

  // Branch-length axis
  const axisYPos = innerH + 20;
  const axisG = g
    .append("g")
    .attr("class", "branch-length-axis")
    .attr("transform", "translate(0, " + axisYPos + ")")
    .call(d3.axisBottom(blToX).ticks(5));
  axisG
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", state.AXIS_FONT_SIZE + "px");
  axisG.selectAll("line").attr("stroke", tc.branch);
  axisG.selectAll("path").attr("stroke", tc.branch);
  axisG
    .append("text")
    .attr("x", innerW / 2)
    .attr("y", state.AXIS_FONT_SIZE + 20)
    .attr("text-anchor", "middle")
    .attr("font-size", state.AXIS_FONT_SIZE)
    .attr("fill", tc.internal)
    .text("Branch length");
}
