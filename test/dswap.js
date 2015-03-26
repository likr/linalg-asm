var expect = require('expect.js');

var LinalgModule = require('../index');

describe('dswap', function() {
  it('calculates x, y := y, x', function() {
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

    var value = linalg.dswap(4, 0, 1, 4, 1);

    expect(x[0]).to.be(5);
    expect(x[1]).to.be(6);
    expect(x[2]).to.be(7);
    expect(x[3]).to.be(8);
    expect(y[0]).to.be(1);
    expect(y[1]).to.be(2);
    expect(y[2]).to.be(3);
    expect(y[3]).to.be(4);
  });

  it('calculates x, y := y, x with offset', function() {
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

    var value = linalg.dswap(3, 1, 1, 4, 1);

    expect(x[0]).to.be(1);
    expect(x[1]).to.be(5);
    expect(x[2]).to.be(6);
    expect(x[3]).to.be(7);
    expect(y[0]).to.be(2);
    expect(y[1]).to.be(3);
    expect(y[2]).to.be(4);
    expect(y[3]).to.be(8);
  });
});
