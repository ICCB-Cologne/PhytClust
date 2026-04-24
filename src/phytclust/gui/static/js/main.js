/* ============================================================
   PhytClust – main.js
   Entry point: DOM ready — wire everything.
   ============================================================ */

import { state, EXAMPLE_NEWICK, SELECTED_CLUSTER_IDS, BOX_LABEL_MAP, NODE_CUSTOM, getBoxAdjust } from "./state.js";
import { newickEl, extraResolutionEl } from "./dom.js";
import { resetExtraParams, collectLeafNamesData } from "./utils.js";
import { BASE_COLORS, shuffle, generateClusterColors, getThemeColors, initTheme } from "./colors.js";
import { showToast } from "./ui/toast.js";
import { showStatus } from "./ui/status.js";
import { switchTab } from "./ui/tabs.js";
import { wireCollapse } from "./ui/collapsibles.js";
import { handleFileSelect, loadFile } from "./ui/file.js";
import { exportSvgFromEl, exportPngFromEl, exportTSV, saveToServer, copySvgToClipboard, copyRasterToClipboard } from "./ui/export.js";
import { runPhytClust } from "./api.js";
import { drawTree, clearAllCollapsedFlags, hideNodeContextMenu, showNodeContextMenu, getNodeCustom, nodeDisplayName, clearTree } from "./tree/draw.js";
import { drawComparison, makeCompareConfig, renderCompareConfigList } from "./views/compare.js";
import { drawOptimalK } from "./views/optimal_k.js";
import {
  getCurrentClusterIds,
  updateClusterEditorAvailability,
  openClusterEditor,
  closeClusterEditor,
  saveClusterEditor,
  renderClusterEditorRows,
  populateClusterSelector,
  switchCluster,
  cycleCluster,
  toggleClusterViewMode,
} from "./views/cluster_editor.js";

// Initialize theme
initTheme();

// Wire context menu dismiss
document.addEventListener("click", hideNodeContextMenu);
document.addEventListener("contextmenu", function (e) {
  if (!e.target.closest(".tree-node")) hideNodeContextMenu();
});

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
      showToast("Parameters reset to defaults.", "info", 2000);
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
  wireCollapse("peak-opts-toggle", "peak-opts-body");
  wireCollapse("outlier-opts-toggle", "outlier-opts-body");
  wireCollapse("polytomy-opts-toggle", "polytomy-opts-body");
  wireCollapse("support-opts-toggle", "support-opts-body");
  wireCollapse("advanced-opts-toggle", "advanced-opts-body");

  // ── Tabs ──
  document.querySelectorAll(".tab-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      switchTab(btn.dataset.tab);
    });
  });
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

    if (isArrow && state.COLOR_MODE === "boxes" && SELECTED_CLUSTER_IDS.size > 0) {
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
      state.COMPARE_CONFIGS.push(makeCompareConfig());
      renderCompareConfigList();
      drawComparison();
    });
  }

  // ── Layout mode ──
  var layoutSelect = document.getElementById("layout-mode");
  if (layoutSelect)
    layoutSelect.addEventListener("change", function () {
      state.CURRENT_LAYOUT_MODE = this.value;
      drawTree();
    });

  // ── Color mode ──
  var colorModeSelect = document.getElementById("color-mode");
  if (colorModeSelect)
    colorModeSelect.addEventListener("change", function () {
      state.COLOR_MODE = this.value || "bars";
      drawTree();
    });

  var branchWidthInput = document.getElementById("branch-width");
  if (branchWidthInput) {
    branchWidthInput.value = state.BRANCH_STROKE_WIDTH;
    branchWidthInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 0.5 && v <= 6) {
        state.BRANCH_STROKE_WIDTH = v;
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
      state.BRANCH_COLOR_OVERRIDE = this.value || null;
      drawTree();
      if (document.getElementById("compare").classList.contains("active"))
        drawComparison();
    });
  }
  if (branchColorResetBtn) {
    branchColorResetBtn.addEventListener("click", function () {
      state.BRANCH_COLOR_OVERRIDE = null;
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
    if (Object.keys(state.CURRENT_CLUSTERS || {}).length > 0) {
      state.CLUSTER_COLORS = generateClusterColors(
        Math.max(...Object.values(state.CURRENT_CLUSTERS)) + 1,
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
      if (Object.keys(state.CURRENT_CLUSTERS || {}).length > 0) {
        state.CLUSTER_COLORS = generateClusterColors(
          Math.max.apply(null, Object.values(state.CURRENT_CLUSTERS).map(Number)) + 1,
        );
        drawTree();
      }
    });
  }

  // ── Leaf/internal name toggles ──
  var internalCb = document.getElementById("show-internal-names");
  if (internalCb) {
    internalCb.checked = state.SHOW_INTERNAL_NAMES;
    internalCb.addEventListener("change", function () {
      state.SHOW_INTERNAL_NAMES = this.checked;
      drawTree();
    });
  }
  var leafCb = document.getElementById("show-leaf-names");
  if (leafCb) {
    leafCb.checked = state.SHOW_LEAF_NAMES;
    leafCb.addEventListener("change", function () {
      state.SHOW_LEAF_NAMES = this.checked;
      drawTree();
    });
  }

  var boxLabelCb = document.getElementById("show-box-labels");
  if (boxLabelCb) {
    boxLabelCb.checked = state.SHOW_BOX_LABELS;
    boxLabelCb.addEventListener("change", function () {
      state.SHOW_BOX_LABELS = this.checked;
      drawTree();
    });
  }
  var outlierBoxesCb = document.getElementById("show-outlier-boxes");
  if (outlierBoxesCb) {
    outlierBoxesCb.checked = state.SHOW_OUTLIER_BOXES;
    outlierBoxesCb.addEventListener("change", function () {
      state.SHOW_OUTLIER_BOXES = this.checked;
      drawTree();
    });
  }
  var boxAlphaSlider = document.getElementById("box-alpha");
  if (boxAlphaSlider) {
    boxAlphaSlider.value = state.BOX_ALPHA;
    boxAlphaSlider.addEventListener("input", function () {
      state.BOX_ALPHA = parseFloat(this.value);
      drawTree();
    });
  }
  var boxPadVSlider = document.getElementById("box-pad-v");
  if (boxPadVSlider) {
    boxPadVSlider.value = state.BOX_PAD_V;
    boxPadVSlider.addEventListener("input", function () {
      state.BOX_PAD_V = parseFloat(this.value);
      drawTree();
    });
  }
  var boxPadHSlider = document.getElementById("box-pad-h");
  if (boxPadHSlider) {
    boxPadHSlider.value = state.BOX_PAD_H;
    boxPadHSlider.addEventListener("input", function () {
      state.BOX_PAD_H = parseFloat(this.value);
      drawTree();
    });
  }
  var boxRadiusSlider = document.getElementById("box-radius");
  if (boxRadiusSlider) {
    boxRadiusSlider.value = state.BOX_CORNER_RADIUS;
    boxRadiusSlider.addEventListener("input", function () {
      state.BOX_CORNER_RADIUS = parseFloat(this.value);
      drawTree();
    });
  }

  // ── Node/label sizes ──
  var nodeSizeInput = document.getElementById("node-size");
  if (nodeSizeInput) {
    nodeSizeInput.value = state.LEAF_NODE_RADIUS;
    nodeSizeInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 1 && v <= 10) {
        state.LEAF_NODE_RADIUS = v;
        drawTree();
      }
    });
  }
  var labelSizeInput = document.getElementById("label-size");
  if (labelSizeInput) {
    labelSizeInput.value = state.LABEL_FONT_SIZE;
    labelSizeInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 6 && v <= 36) {
        state.LABEL_FONT_SIZE = v;
        drawTree();
      }
    });
  }
  var axisSizeInput = document.getElementById("axis-font-size");
  if (axisSizeInput) {
    axisSizeInput.value = state.AXIS_FONT_SIZE;
    axisSizeInput.addEventListener("input", function () {
      var v = parseFloat(this.value);
      if (!isNaN(v) && v >= 6 && v <= 18) {
        state.AXIS_FONT_SIZE = v;
        drawTree();
      }
    });
  }

  // ── Width/Height sliders ──
  var widthSlider = document.getElementById("tree-width-scale");
  var heightSlider = document.getElementById("tree-height-scale");
  if (widthSlider)
    widthSlider.addEventListener("input", function () {
      state.TREE_WIDTH_SCALE = parseFloat(this.value);
      if (!isNaN(state.TREE_WIDTH_SCALE)) drawTree();
    });
  if (heightSlider)
    heightSlider.addEventListener("input", function () {
      state.TREE_HEIGHT_SCALE = parseFloat(this.value);
      if (!isNaN(state.TREE_HEIGHT_SCALE)) drawTree();
    });

  // ── Reset view ──
  var btnResetView = document.getElementById("btn-reset-view");
  if (btnResetView) {
    btnResetView.addEventListener("click", function (e) {
      e.preventDefault();
      if (state.LAST_TREE_SVG && state.LAST_TREE_ZOOM)
        state.LAST_TREE_SVG.transition()
          .duration(150)
          .call(state.LAST_TREE_ZOOM.transform, d3.zoomIdentity);
      if (widthSlider) widthSlider.value = "1.0";
      if (heightSlider) heightSlider.value = "1.0";
      state.TREE_WIDTH_SCALE = 1.0;
      state.TREE_HEIGHT_SCALE = 1.0;
      clearAllCollapsedFlags(state.NEWICK_RAW_TREE);
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
        state.SEARCH_TERM = searchInput.value.trim();
        if (state.NEWICK_RAW_TREE) drawTree();
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
      if (!state.CTX_TARGET_DATA) return;
      var c = getNodeCustom(state.CTX_TARGET_DATA);
      var val = prompt("Rename node:", nodeDisplayName(state.CTX_TARGET_DATA));
      if (val != null) {
        c.renamedTo = val;
        drawTree();
      }
    });
    btn("ctx-highlight", function () {
      if (!state.CTX_TARGET_DATA) return;
      var c = getNodeCustom(state.CTX_TARGET_DATA);
      c.highlighted = !c.highlighted;
      drawTree();
    });
    var ctxSlider = document.getElementById("ctx-size-slider");
    var ctxSliderVal = document.getElementById("ctx-size-value");
    if (ctxSlider) {
      ctxSlider.addEventListener("input", function () {
        if (!state.CTX_TARGET_DATA) return;
        var c = getNodeCustom(state.CTX_TARGET_DATA);
        c.radiusScale = parseFloat(ctxSlider.value);
        if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);
        drawTree();
      });
      ctxSlider.addEventListener("click", function (e) {
        e.stopPropagation();
      });
    }
    btn("ctx-copy-subtree", function () {
      if (!state.CTX_TARGET_DATA) return;
      var names = [];
      collectLeafNamesData(state.CTX_TARGET_DATA, names);
      navigator.clipboard.writeText(names.join("\n")).then(
        () => showToast(names.length + " leaf names copied", "success", 2000),
        () => showToast("Copy failed", "danger"),
      );
    });
    btn("ctx-original-size", function () {
      if (!state.CTX_TARGET_DATA) return;
      var c = getNodeCustom(state.CTX_TARGET_DATA);
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
      if (!state.CTX_TARGET_DATA || !state.CURRENT_CLUSTERS) return;
      var name = state.CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? state.CURRENT_CLUSTERS[name] : null;
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
      if (!state.CTX_TARGET_DATA || !state.CURRENT_CLUSTERS) return;
      var name = state.CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? state.CURRENT_CLUSTERS[name] : null;
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
      if (!state.CTX_TARGET_DATA || !state.CURRENT_CLUSTERS) return;
      var name = state.CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? state.CURRENT_CLUSTERS[name] : null;
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
      if (!state.CTX_TARGET_DATA) return;
      NODE_CUSTOM.delete(state.CTX_TARGET_DATA);
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
        if (e.target === miniScoresToggle) return;
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
