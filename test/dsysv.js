var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dsysv', function() {
  it('solves A x = B where A = U D U^t', function() {
    var heap = new ArrayBuffer(136),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        ipiv = new Uint32Array(heap, 120, 3),
        linalg = linalgModule(global, null, heap);

    a.set([1, 1, 1, 1, 2, 2, 1, 2, 3]);
    b.set([6, 4, 11, 10, 14, 15]);
    linalg.dsysv(0, 3, 2, a.byteOffset, 3, ipiv.byteOffset, b.byteOffset, 2);

    expect(+b[0].toFixed(6)).to.be(1);
    expect(+b[1].toFixed(6)).to.be(-2);
    expect(+b[2].toFixed(6)).to.be(2);
    expect(+b[3].toFixed(6)).to.be(1);
    expect(+b[4].toFixed(6)).to.be(3);
    expect(+b[5].toFixed(6)).to.be(5);
  });
});
