var expect = require('expect.js');

var LinalgModule = require('../index');

describe('dcopy', function() {
  it('calculate y := x', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;

    linalg.dcopy(4, x.byteOffset, 1, y.byteOffset, 1);

    expect(y[0]).to.be(1);
    expect(y[1]).to.be(2);
    expect(y[2]).to.be(3);
    expect(y[3]).to.be(4);
  });
});
