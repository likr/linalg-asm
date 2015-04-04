var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dpotrs', function() {
  it('solves A x = B where A = L L^t', function() {
    var heap = new ArrayBuffer(120),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = linalgModule(global, null, heap);

    a.set([1, 0, 0, 1, 1, 0, 1, 1, 1]);
    b.set([6, 4, 11, 10, 14, 15]);
    linalg.dpotrs(1, 3, 2, a.byteOffset, 3, b.byteOffset, 2);

    expect(b[0]).to.be(1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be(2);
    expect(b[3]).to.be(1);
    expect(b[4]).to.be(3);
    expect(b[5]).to.be(5);
  });

  it('solves A x = B where A = U^t U', function() {
    var heap = new ArrayBuffer(120),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = linalgModule(global, null, heap);

    a.set([1, 1, 1, 0, 1, 1, 0, 0, 1]);
    b.set([6, 4, 11, 10, 14, 15]);
    linalg.dpotrs(0, 3, 2, a.byteOffset, 3, b.byteOffset, 2);

    expect(b[0]).to.be(1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be(2);
    expect(b[3]).to.be(1);
    expect(b[4]).to.be(3);
    expect(b[5]).to.be(5);
  });
});
