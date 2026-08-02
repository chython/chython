// smiles-drawer 2.4.1 restricts its package `exports` to the bundled dist, so the
// layout internals are imported directly from the package `src/` by path.
import DrawerBase from '../node_modules/smiles-drawer/src/DrawerBase.js';


// Lay out a molecule from a pre-built smiles-drawer parse tree (see clean2d/README
// for the tree schema). The tree is constructed directly from the chython molecule on
// the Python side, so no SMILES string round-trip or positional atom re-mapping is
// needed. Heavy-atom `idx` is assigned by DrawerBase in tree-DFS creation order, which
// matches the order the Python builder reports, so positions come back index-aligned.
function clean2d(tree) {
    const drawer = new DrawerBase({});
    drawer.initDraw(tree, 'light', false);
    drawer.processGraph();

    const g = drawer.graph;
    const xy = [];
    for (let i = 0; i < g.atomIdxToVertexId.length; i++) {
        const position = g.vertices[g.atomIdxToVertexId[i]].position;
        xy.push([position.x, position.y]);
    }
    return xy;
}

export { clean2d };
