var expect = require('expect.js');

var LinalgModule = require('./index');

describe('daxpy', function() {
  it('calculate y := a x + y', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    y[0] = 5;
    y[1] = 6;
    y[2] = 7;
    y[3] = 8;

    linalg.daxpy(4, 1.5, 0, 1, 4, 1);

    expect(y[0]).to.be(6.5);
    expect(y[1]).to.be(9);
    expect(y[2]).to.be(11.5);
    expect(y[3]).to.be(14);
  });
});

describe('dasum', function() {
  it('return |x_1| + ... + |x_n|', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.dasum(4, 0, 1);

    expect(result).to.be(8);
  });
});

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

    linalg.dcopy(4, 0, 1, 4, 1);

    expect(y[0]).to.be(1);
    expect(y[1]).to.be(2);
    expect(y[2]).to.be(3);
    expect(y[3]).to.be(4);
  });
});
