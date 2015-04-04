var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dsyr', function() {
  it('A = a x x^t + A', function() {
    var heap = new ArrayBuffer(96),
        a = new Float64Array(heap, 0, 9),
        x = new Float64Array(heap, 72, 3),
        linalg = linalgModule(global, null, heap);

    a.set([1, 1, 1, 1, 1, 1, 1, 1, 1]);
    x.set([1, 2, 3]);
    linalg.dsyr(0, 3, 0.5, x.byteOffset, 1, a.byteOffset, 3);

    expect(a[0]).to.be(1.5);
    expect(a[1]).to.be(2);
    expect(a[2]).to.be(2.5);
    expect(a[4]).to.be(3);
    expect(a[5]).to.be(4);
    expect(a[8]).to.be(5.5);
  });
});
