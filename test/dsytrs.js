var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dsytrs', function() {
  it('solves U D U^t x = B', function() {
    var heap = new ArrayBuffer(136),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        ipiv = new Uint32Array(heap, 120, 3),
        linalg = linalgModule(global, null, heap);

    a.set([0.5, 0.5, 1 / 3, 0, 2 / 3, 2 / 3, 0, 0, 3]);
    b.set([6, 4, 11, 10, 14, 15]);
    ipiv.set([0, 1, 2]);
    linalg.dsytrs(0, 3, 2, a.byteOffset, 3, ipiv.byteOffset, b.byteOffset, 2);

    expect(+b[0].toFixed(6)).to.be(1);
    expect(+b[1].toFixed(6)).to.be(-2);
    expect(+b[2].toFixed(6)).to.be(2);
    expect(+b[3].toFixed(6)).to.be(1);
    expect(+b[4].toFixed(6)).to.be(3);
    expect(+b[5].toFixed(6)).to.be(5);
  });
});
