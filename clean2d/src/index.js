const { Parser, Drawer } = require('smiles-drawer');

const options = {
      debug: false,
      atomVisualization: 'default'
    };
const drawer = new Drawer(options);


function clean2d(smiles) {
    let parsed = Parser.parse(smiles);
    drawer.initDraw(parsed, 'light', false);
    drawer.processGraph();

    let vertices = drawer.graph.vertices;
    let xy = Array();
    for (let i = 0; i < vertices.length; i++) {
      let position = vertices[i].position;
      xy.push([position.x, position.y]);
    }
    return xy;
}

export { clean2d };
