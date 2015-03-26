var expect = require('expect.js'),
    LinalgModule = require('../index');

describe('dgemv', function() {
  it('computes y := a A x + b y', function() {
    var heap = new ArrayBuffer(192),
        A = new Float64Array(heap, 0, 16),
        x = new Float64Array(heap, 128, 4),
        y = new Float64Array(heap, 160, 4),
        linalg = LinalgModule(global, null, heap);

    A[ 0] =   1.80;
    A[ 1] =   2.88;
    A[ 2] =   2.05;
    A[ 3] = - 0.89;
    A[ 4] =   5.25;
    A[ 5] = - 2.95;
    A[ 6] = - 0.95;
    A[ 7] = - 3.80;
    A[ 8] =   1.58;
    A[ 9] = - 2.69;
    A[10] = - 2.90;
    A[11] = - 1.04;
    A[12] = - 1.11;
    A[13] = - 0.66;
    A[14] = - 0.59;
    A[15] =   0.80;
    x[ 0] =   1;
    x[ 1] = - 1;
    x[ 2] =   3;
    x[ 3] = - 5;

    linalg.dgemv(0, 4, 4, 1, 0, 4, 128, 1, 0, 160, 1);

    expect(+y[0].toPrecision(5)).to.be(  9.52);
    expect(+y[1].toPrecision(5)).to.be( 24.35);
    expect(+y[2].toPrecision(5)).to.be(  0.77);
    expect(+y[3].toPrecision(5)).to.be(- 6.22);
  });
});
