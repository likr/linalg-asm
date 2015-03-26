var expect = require('expect.js');

var LinalgModule = require('../index');

describe('dgetrf', function() {
  it('computes the LU factorization of a real m by n matrix', function() {
    var heap = new ArrayBuffer(144),
        a = new Float64Array(heap, 0, 16),
        ipiv = new Int32Array(heap, 128, 4),
        linalg = LinalgModule(global, null, heap);

    a[ 0] =  1.80;
    a[ 1] =  2.88;
    a[ 2] =  2.05;
    a[ 3] = -0.89;
    a[ 4] =  5.25;
    a[ 5] = -2.95;
    a[ 6] = -0.95;
    a[ 7] = -3.80;
    a[ 8] =  1.58;
    a[ 9] = -2.69;
    a[10] = -2.90;
    a[11] = -1.04;
    a[12] = -1.11;
    a[13] = -0.66;
    a[14] = -0.59;
    a[15] =  0.80;

    linalg.dgetrf(4, 4, 0, 4, 128);

    expect(a[ 0]).to.be( 5.25                );
    expect(a[ 1]).to.be(-2.95                );
    expect(a[ 2]).to.be(-0.95                );
    expect(a[ 3]).to.be(-3.8                 );
    expect(a[ 4]).to.be( 0.34285714285714286 );
    expect(a[ 5]).to.be( 3.8914285714285715  );
    expect(a[ 6]).to.be( 2.3757142857142854  );
    expect(a[ 7]).to.be( 0.4128571428571427  );
    expect(a[ 8]).to.be( 0.30095238095238097 );
    expect(a[ 9]).to.be(-0.46311796377875664 );
    expect(a[10]).to.be(-1.5138592755751348  );
    expect(a[11]).to.be( 0.29482060695056267 );
    expect(a[12]).to.be(-0.21142857142857144 );
    expect(a[13]).to.be(-0.3298825256975037  );
    expect(a[14]).to.be( 0.004723367663983699);
    expect(a[15]).to.be( 0.13137323948785168 );
    expect(ipiv[0]).to.be(1);
    expect(ipiv[1]).to.be(1);
    expect(ipiv[2]).to.be(2);
    expect(ipiv[3]).to.be(3);
  });
});
