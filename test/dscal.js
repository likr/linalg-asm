var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dscal', function() {
  it('computes x := a x', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = linalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;

    linalg.dscal(4, 2.5, x.byteOffset, 1);

    expect(x[0]).to.be(2.5);
    expect(x[1]).to.be(5);
    expect(x[2]).to.be(7.5);
    expect(x[3]).to.be(10);
  });
});

