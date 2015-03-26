var expect = require('expect.js'),
    LinalgModule = require('../index');

describe('dtrsv', function() {
  it('solves Lx = b where L is lower triangular matrix', function() {
    var heap = new ArrayBuffer(160),
        a = new Float64Array(heap, 0, 16),
        b = new Float64Array(heap, 128, 4),
        linalg = LinalgModule(global, null, heap);

    a[ 0] =   1                   ;
    a[ 1] =   0                   ;
    a[ 2] =   0                   ;
    a[ 3] =   0                   ;
    a[ 4] =   0.34285714285714286 ;
    a[ 5] =   1                   ;
    a[ 6] =   0                   ;
    a[ 7] =   0                   ;
    a[ 8] =   0.30095238095238097 ;
    a[ 9] = - 0.46311796377875664 ;
    a[10] =   1                   ;
    a[11] =   0                   ;
    a[12] = - 0.21142857142857144 ;
    a[13] = - 0.3298825256975037  ;
    a[14] =   0.004723367663983699;
    a[15] =   1                   ;
    b[ 0] =  24.35;
    b[ 1] =   9.52;
    b[ 2] =   0.77;
    b[ 3] = - 6.22;

    linalg.dtrsv(1, 0, 0, 4, 0, 4, 128, 1);

    expect(b[0]).to.be( 24.35              );
    expect(b[1]).to.be(  1.1714285714285708);
    expect(b[2]).to.be(- 6.01568086147822  );
    expect(b[3]).to.be(- 0.6568661974392578);
  });

  it('solves Lx = b where L is unit lower triangular matrix', function() {
    var heap = new ArrayBuffer(160),
        a = new Float64Array(heap, 0, 16),
        b = new Float64Array(heap, 128, 4),
        linalg = LinalgModule(global, null, heap);

    a[ 0] =   0                   ;
    a[ 1] =   0                   ;
    a[ 2] =   0                   ;
    a[ 3] =   0                   ;
    a[ 4] =   0.34285714285714286 ;
    a[ 5] =   0                   ;
    a[ 6] =   0                   ;
    a[ 7] =   0                   ;
    a[ 8] =   0.30095238095238097 ;
    a[ 9] = - 0.46311796377875664 ;
    a[10] =   0                   ;
    a[11] =   0                   ;
    a[12] = - 0.21142857142857144 ;
    a[13] = - 0.3298825256975037  ;
    a[14] =   0.004723367663983699;
    a[15] =   0                   ;
    b[ 0] =  24.35;
    b[ 1] =   9.52;
    b[ 2] =   0.77;
    b[ 3] = - 6.22;

    linalg.dtrsv(1, 0, 1, 4, 0, 4, 128, 1);

    expect(b[0]).to.be( 24.35              );
    expect(b[1]).to.be(  1.1714285714285708);
    expect(b[2]).to.be(- 6.01568086147822  );
    expect(b[3]).to.be(- 0.6568661974392578);
  });

  it('solves Ux = b where U is upper triangular matrix', function() {
    var heap = new ArrayBuffer(160),
        a = new Float64Array(heap, 0, 16),
        b = new Float64Array(heap, 128, 4),
        linalg = LinalgModule(global, null, heap);

    a[ 0] =   5.25               ;
    a[ 1] = - 2.95               ;
    a[ 2] = - 0.95               ;
    a[ 3] = - 3.8                ;
    a[ 4] =   0                  ;
    a[ 5] =   3.8914285714285715 ;
    a[ 6] =   2.3757142857142854 ;
    a[ 7] =   0.4128571428571427 ;
    a[ 8] =   0                  ;
    a[ 9] =   0                  ;
    a[10] = - 1.5138592755751348 ;
    a[11] =   0.29482060695056267;
    a[12] =   0                  ;
    a[13] =   0                  ;
    a[14] =   0                  ;
    a[15] =   0.13137323948785168;
    b[ 0] =  24.35               ;
    b[ 1] =   1.1714285714285708 ;
    b[ 2] = - 6.01568086147822   ;
    b[ 3] = - 0.6568661974392578 ;

    linalg.dtrsv(0, 0, 0, 4, 0, 4, 128, 1);

    expect(+b[0].toPrecision(5)).to.be( 1);
    expect(+b[1].toPrecision(5)).to.be(-1);
    expect(+b[2].toPrecision(5)).to.be( 3);
    expect(+b[3].toPrecision(5)).to.be(-5);
  });
});
