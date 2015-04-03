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

    linalg.daxpy(x.length, 1.5, x.byteOffset, 1, y.byteOffset, 1);

    expect(y[0]).to.be(6.5);
    expect(y[1]).to.be(9);
    expect(y[2]).to.be(11.5);
    expect(y[3]).to.be(14);
  });

  it('calculate y := a x + y', function() {
    var heap = new ArrayBuffer(128),
        x = new Float64Array(heap, 16, 7),
        y = new Float64Array(heap, 72, 7),
        linalg = linalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    x[4] = 5;
    x[5] = 6;
    x[6] = 7;
    y[0] = 8;
    y[1] = 9;
    y[2] = 10;
    y[3] = 11;
    y[4] = 12;
    y[5] = 13;
    y[6] = 14;

    linalg.daxpy(x.length, 1.5, x.byteOffset, 1, y.byteOffset, 1);

    expect(y[0]).to.be(9.5);
    expect(y[1]).to.be(12);
    expect(y[2]).to.be(14.5);
    expect(y[3]).to.be(17);
    expect(y[4]).to.be(19.5);
    expect(y[5]).to.be(22);
    expect(y[6]).to.be(24.5);
  });
});
