var expect = require('expect.js');

var LinalgModule = require('./index');

describe('daxpy', function() {
  it('calculate y := a x + y', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    y[0] = 5;
    y[1] = 6;
    y[2] = 7;
    y[3] = 8;

    linalg.daxpy(4, 1.5, 0, 1, 4, 1);

    expect(y[0]).to.be(6.5);
    expect(y[1]).to.be(9);
    expect(y[2]).to.be(11.5);
    expect(y[3]).to.be(14);
  });
});

describe('dasum', function() {
  it('return |x_1| + ... + |x_n|', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.dasum(4, 0, 1);

    expect(result).to.be(8);
  });
});

describe('dcopy', function() {
  it('calculate y := x', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;

    linalg.dcopy(4, 0, 1, 4, 1);

    expect(y[0]).to.be(1);
    expect(y[1]).to.be(2);
    expect(y[2]).to.be(3);
    expect(y[3]).to.be(4);
  });
});

describe('ddot', function() {
  it('return x^t y', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    y[0] = 5;
    y[1] = 6;
    y[2] = 7;
    y[3] = 8;

    var value = linalg.ddot(4, 0, 1, 4, 1);

    expect(value).to.be(70);
  });
});

describe('dswap', function() {
  it('calculates x, y := y, x', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    y[0] = 5;
    y[1] = 6;
    y[2] = 7;
    y[3] = 8;

    var value = linalg.dswap(4, 0, 1, 4, 1);

    expect(x[0]).to.be(5);
    expect(x[1]).to.be(6);
    expect(x[2]).to.be(7);
    expect(x[3]).to.be(8);
    expect(y[0]).to.be(1);
    expect(y[1]).to.be(2);
    expect(y[2]).to.be(3);
    expect(y[3]).to.be(4);
  });

  it('calculates x, y := y, x with offset', function() {
    var heap = new ArrayBuffer(64),
        x = new Float64Array(heap, 0, 4),
        y = new Float64Array(heap, 32, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    y[0] = 5;
    y[1] = 6;
    y[2] = 7;
    y[3] = 8;

    var value = linalg.dswap(3, 1, 1, 4, 1);

    expect(x[0]).to.be(1);
    expect(x[1]).to.be(5);
    expect(x[2]).to.be(6);
    expect(x[3]).to.be(7);
    expect(y[0]).to.be(2);
    expect(y[1]).to.be(3);
    expect(y[2]).to.be(4);
    expect(y[3]).to.be(8);
  });
});

describe('idamax', function() {
  it('returns index of max abs value', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.idamax(4, 0, 1);

    expect(result).to.be(3);
  });

  it('returns index of max abs value with incX', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = LinalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.idamax(2, 1, 2);

    expect(result).to.be(1);
  });
});

describe('dgemm', function() {
  it('calculate C := a A B + b C', function() {
    var heap = new ArrayBuffer(216),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 9),
        c = new Float64Array(heap, 144, 9),
        linalg = LinalgModule(global, null, heap);

    a[0] =   1  ;
    a[1] =   8  ;
    a[2] =   3  ;
    a[3] =   2  ;
    a[4] =  10  ;
    a[5] =   8  ;
    a[6] =   9  ;
    a[7] = - 5  ;
    a[8] = - 1  ;
    b[0] =   9  ;
    b[1] =   8  ;
    b[2] =   3  ;
    b[3] =   3  ;
    b[4] =  11  ;
    b[5] =   2.3;
    b[6] = - 8  ;
    b[7] =   6  ;
    b[8] =   1  ;
    c[0] =   3  ;
    c[1] =   3  ;
    c[2] =   1.2;
    c[3] =   8  ;
    c[4] =   4  ;
    c[5] =   8  ;
    c[6] =   6  ;
    c[7] =   1  ;
    c[8] = - 2  ;

    linalg.dgemm(0, 0, 3, 3, 3, 3.0, 0, 3, 9, 3, -2.0, 18, 3);

    expect(c[0]).to.be(  21  );
    expect(c[1]).to.be( 336  );
    expect(c[2]).to.be(  70.8);
    expect(c[3]).to.be(- 64  );
    expect(c[4]).to.be( 514  );
    expect(c[5]).to.be(  95  );
    expect(c[6]).to.be( 210  );
    expect(c[7]).to.be(  31  );
    expect(c[8]).to.be(  47.5);
  });
});

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
