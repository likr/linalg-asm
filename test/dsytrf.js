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

  it('computes the Cholesky factorization A = L^t D L', function() {
    var heap = new ArrayBuffer(88),
        a = new Float64Array(heap, 0, 9),
        ipiv = new Uint32Array(heap, 72, 3),
        linalg = linalgModule(global, null, heap);

    a.set([1, 1, 2, 1, 13, 8, 2, 8, 9]);
    linalg.dsytrf(1, 3, a.byteOffset, 3, ipiv.byteOffset);

    expect(a[0]).to.be(1);
    expect(a[3]).to.be(1);
    expect(a[4]).to.be(12);
    expect(a[6]).to.be(2);
    expect(a[7]).to.be(1 / 2);
    expect(a[8]).to.be(2);
  });

  it('computes the Cholesky factorization A = L^t D L', function() {
    var heap = new ArrayBuffer(88),
        a = new Float64Array(heap, 0, 9),
        ipiv = new Uint32Array(heap, 72, 3),
        linalg = linalgModule(global, null, heap);

    a.set([10, -10, -4, -10, 1, -14, -4, -14, -2]);
    linalg.dsytrf(1, 3, a.byteOffset, 3, ipiv.byteOffset);

    console.log(a, ipiv);
    expect(a[0]).to.be(10);
    expect(a[3]).to.be(-1);
    expect(a[4]).to.be(-9);
    expect(a[6]).to.be(-2 / 5);
    expect(a[7]).to.be(2);
    expect(a[8]).to.be(162 / 5);
  });
});
