const path = require('path');

module.exports = {
  entry: path.resolve(__dirname, "src/index.js"),
  mode: "production",
  output: {
    filename: 'clean2d.js',
    path: path.resolve(__dirname, 'dist'),
    library: "$",
    libraryTarget: "umd"
  }
};
