/* ============================================================
   PhytClust – state.js
   All mutable application state in one place.
   ============================================================ */

export const state = {
  HIER_CART: null,
  HIER_CIRC: null,

  SHOW_LEAF_NAMES: true,
  SHOW_INTERNAL_NAMES: false,
  LEAF_NODE_RADIUS: 3.0,
  INTERNAL_NODE_RADIUS: 1.8,
  LABEL_FONT_SIZE: 9,
  AXIS_FONT_SIZE: 10,
  TREE_WIDTH_SCALE: 1.0,
  TREE_HEIGHT_SCALE: 1.0,
  BRANCH_STROKE_WIDTH: 1.2,
  BRANCH_COLOR_OVERRIDE: null,
  CLUSTER_COLORS: [],
  CURRENT_CLUSTERS: {},
  NEWICK_RAW_TREE: null,
  CURRENT_LAYOUT_MODE: "rectangular",
  COLOR_MODE: "bars",

  isRunning: false,
  latestOptimalKData: null,
  latestApiData: null,
  latestRunId: null,
  LAST_TREE_SVG: null,
  LAST_TREE_ZOOM: null,
  LAST_TREE_TRANSFORM: d3.zoomIdentity,
  SEARCH_TERM: "",
  SHOW_BOX_LABELS: true,
  SHOW_OUTLIER_BOXES: true,
  BOX_ALPHA: 0.12,
  BOX_PAD_V: 8,
  BOX_PAD_H: 6,
  BOX_CORNER_RADIUS: 6,
  CLUSTER_LABEL_FONT_SIZE: 8,

  CLUSTER_VIEW_MODE: "peaks",
  COMPARE_CONFIGS: [],
  COMPARE_CONFIG_NEXT_ID: 1,

  CTX_TARGET_DATA: null,
};

export const LABEL_PAD = 6;

export const EXAMPLE_NEWICK =
  "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

export const SELECTED_CLUSTER_IDS = new Set();
export const BOX_LABEL_MAP = {}; // cid -> custom name, editable via context menu
export const BOX_ADJUST_MAP = {}; // cid -> {dx, dy, padX, padY} in px

// Per-node customizations: keyed by data node reference
export const NODE_CUSTOM = new WeakMap();

export function hasClusterFocus(cid) {
  if (cid == null) return SELECTED_CLUSTER_IDS.size === 0;
  return (
    SELECTED_CLUSTER_IDS.size === 0 || SELECTED_CLUSTER_IDS.has(Number(cid))
  );
}

export function getBoxAdjust(cid) {
  const key = String(cid);
  if (!BOX_ADJUST_MAP[key]) {
    BOX_ADJUST_MAP[key] = { dx: 0, dy: 0, padX: 0, padY: 0 };
  }
  return BOX_ADJUST_MAP[key];
}
