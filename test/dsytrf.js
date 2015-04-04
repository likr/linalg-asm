var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dsytrf', function() {
  it('computes the Cholesky factorization A = U^t D U', function() {
    var heap = new ArrayBuffer(88),
        a = new Float64Array(heap, 0, 9),
        ipiv = new Uint32Array(heap, 72, 3),
        linalg = linalgModule(global, null, heap);

    a.set([157, -177, -36, 0, 226, 45, 0, 0, 9]);
    linalg.dsytrf(0, 3, a.byteOffset, 3, ipiv.byteOffset);

    expect(a[0]).to.be(4);
    expect(a[1]).to.be(3);
    expect(a[2]).to.be(-4);
    expect(a[4]).to.be(1);
    expect(a[5]).to.be(5);
    expect(a[8]).to.be(9);
  });
});
