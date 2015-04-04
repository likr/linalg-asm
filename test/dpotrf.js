var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dgetrf', function() {
  it('computes the Cholesky factorization A = L L^t', function() {
    var heap = new ArrayBuffer(72),
        a = new Float64Array(heap, 0, 9),
        linalg = linalgModule(global, null, heap);

    a.set([1, 1, 1, 1, 2, 2, 1, 2, 3]);
    linalg.dpotrf(1, 3, a.byteOffset, 3);

    expect(a[0]).to.be(1);
    expect(a[3]).to.be(1);
    expect(a[4]).to.be(1);
    expect(a[6]).to.be(1);
    expect(a[7]).to.be(1);
    expect(a[8]).to.be(1);
  });

  it('computes the Cholesky factorization A = U^t U', function() {
    var heap = new ArrayBuffer(72),
        a = new Float64Array(heap, 0, 9),
        linalg = linalgModule(global, null, heap);

    a.set([1, 1, 1, 1, 2, 2, 1, 2, 3]);
    linalg.dpotrf(0, 3, a.byteOffset, 3);

    expect(a[0]).to.be(1);
    expect(a[1]).to.be(1);
    expect(a[2]).to.be(1);
    expect(a[4]).to.be(1);
    expect(a[5]).to.be(1);
    expect(a[8]).to.be(1);
  });
});
