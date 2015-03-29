var expect = require('expect.js'),
    LinalgModule = require('../index');

describe('dtrsm', function() {
  it('solves A X = B', function() {
    var heap = new ArrayBuffer(160),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = LinalgModule(global, null, heap);

    a[0] =   1;
    a[1] =   4;
    a[2] =   8;
    a[3] =   0;
    a[4] =   5;
    a[5] =   2;
    a[6] =   0;
    a[7] =   0;
    a[8] = - 3;
    b[0] =  33;
    b[1] =  18;
    b[2] =  16;
    b[3] =   9;
    b[4] = - 9;
    b[5] = - 6;

    linalg.dtrsm(0, 0, 0, 0, 4, 2, 1, a.byteOffset, 3, b.byteOffset, 3);

    expect(b[0]).to.be( 1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be( 2);
    expect(b[3]).to.be( 1);
    expect(b[4]).to.be( 3);
    expect(b[5]).to.be( 2);
  });
});
