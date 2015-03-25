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
