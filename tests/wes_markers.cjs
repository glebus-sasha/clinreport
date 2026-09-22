const assert = require('node:assert/strict');
const fs = require('node:fs');

// Exercise the actual canvas callback with visible and hidden WES-positive nodes.
const source = fs.readFileSync('R/server/drug_network_plot.R', 'utf8');
const callback = source.match(/afterDrawing = JS\(\s*"([\s\S]*?)"\s*\)/);
assert.ok(callback, 'WES drawing callback exists');
const draw = new Function(`return (${callback[1]});`)();
const arcs = [];
const ctx = {
  save() {}, restore() {}, beginPath() {}, fill() {}, stroke() {},
  arc(x, y) { arcs.push([x, y]); }
};
const visible = { x: 1, y: 2, options: { mutated: true, hidden: false } };
const hidden = { x: 3, y: 4, options: { mutated: true, hidden: true } };
const unmutated = { x: 5, y: 6, options: { mutated: false } };
const network = { body: { nodes: { visible, hidden, unmutated } } };
draw.call(network, ctx);
assert.deepEqual(arcs, [[1, 2]]);
hidden.options.hidden = false;
arcs.length = 0;
draw.call(network, ctx);
assert.deepEqual(arcs, [[1, 2], [3, 4]]);
visible.options.hidden = true;
hidden.options.hidden = true;
arcs.length = 0;
draw.call(network, ctx);
assert.deepEqual(arcs, []);
console.log('PASS: WES markers follow node visibility on each redraw');
