var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dgetri', function() {
  it('computes A ^ -1', function() {
    var heap = new ArrayBuffer(136),
        a = new Float64Array(heap, 0, 9),
        ipiv = new Int32Array(heap, 72, 9),
        work = new Float64Array(heap, 112, 3),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 2;
    a[2] = 5;
    a[3] = 1;
    a[4] = -1;
    a[5] = 1;
    a[6] = 0;
    a[7] = 1;
    a[8] = 2;

    linalg.dgetrf(3, 3, a.byteOffset, 3, ipiv.byteOffset);
    linalg.dgetri(3, a.byteOffset, 3, ipiv.byteOffset, work.byteOffset, 3);

    expect(+a[0].toPrecision(5)).to.be(1.5);
    expect(+a[1].toPrecision(5)).to.be(-0.5);
    expect(+a[2].toPrecision(5)).to.be(-3.5);
    expect(+a[3].toPrecision(5)).to.be( 1);
    expect(+a[4].toPrecision(5)).to.be(-1);
    expect(+a[5].toPrecision(5)).to.be(-2);
    expect(+a[6].toPrecision(5)).to.be(-0.5);
    expect(+a[7].toPrecision(5)).to.be(0.5);
    expect(+a[8].toPrecision(5)).to.be(1.5);
  });
});
