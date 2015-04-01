var expect = require('expect.js'),
    linalgModule = require('../index');

describe('daxpy', function() {
  it('calculate y := a x + y', function() {
    var heap = new ArrayBuffer(80),
        x = new Float64Array(heap, 16, 4),
        y = new Float64Array(heap, 48, 4),
        linalg = linalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    y[0] = 5;
    y[1] = 6;
    y[2] = 7;
    y[3] = 8;

    linalg.daxpy(4, 1.5, x.byteOffset, 1, y.byteOffset, 1);

    expect(y[0]).to.be(6.5);
    expect(y[1]).to.be(9);
    expect(y[2]).to.be(11.5);
    expect(y[3]).to.be(14);
  });
});
