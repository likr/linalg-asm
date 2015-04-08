var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dgemm', function() {
  it('calculate C := a A B + b C', function() {
    var heap = new ArrayBuffer(216),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 9),
        c = new Float64Array(heap, 144, 9),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 8;
    a[2] = 3;
    a[3] = 2;
    a[4] = 10;
    a[5] = 8;
    a[6] = 9;
    a[7] = -5;
    a[8] = -1;
    b[0] = 9;
    b[1] = 8;
    b[2] = 3;
    b[3] = 3;
    b[4] = 11;
    b[5] = 2.3;
    b[6] = -8;
    b[7] = 6;
    b[8] = 1;
    c[0] = 3;
    c[1] = 3;
    c[2] = 1.2;
    c[3] = 8;
    c[4] = 4;
    c[5] = 8;
    c[6] = 6;
    c[7] = 1;
    c[8] = -2;

    linalg.dgemm(0, 0, 3, 3, 3, 3.0, a.byteOffset, 3, b.byteOffset, 3, -2.0, c.byteOffset, 3);

    expect(+c[0].toFixed(6)).to.be(21);
    expect(+c[1].toFixed(6)).to.be(336);
    expect(+c[2].toFixed(6)).to.be(70.8);
    expect(+c[3].toFixed(6)).to.be(-64);
    expect(+c[4].toFixed(6)).to.be(514);
    expect(+c[5].toFixed(6)).to.be(95);
    expect(+c[6].toFixed(6)).to.be(210);
    expect(+c[7].toFixed(6)).to.be(31);
    expect(+c[8].toFixed(6)).to.be(47.5);
  });

  it('calculate C := a A^t B + b C', function() {
    var heap = new ArrayBuffer(216),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 9),
        c = new Float64Array(heap, 144, 9),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 2;
    a[2] = 9;
    a[3] = 8;
    a[4] = 10;
    a[5] = -5;
    a[6] = 3;
    a[7] = 8;
    a[8] = -1;
    b[0] = 9;
    b[1] = 8;
    b[2] = 3;
    b[3] = 3;
    b[4] = 11;
    b[5] = 2.3;
    b[6] = -8;
    b[7] = 6;
    b[8] = 1;
    c[0] = 3;
    c[1] = 3;
    c[2] = 1.2;
    c[3] = 8;
    c[4] = 4;
    c[5] = 8;
    c[6] = 6;
    c[7] = 1;
    c[8] = -2;

    linalg.dgemm(1, 0, 3, 3, 3, 3.0, a.byteOffset, 3, b.byteOffset, 3, -2.0, c.byteOffset, 3);

    expect(+c[0].toFixed(6)).to.be(21);
    expect(+c[1].toFixed(6)).to.be(336);
    expect(+c[2].toFixed(6)).to.be(70.8);
    expect(+c[3].toFixed(6)).to.be(-64);
    expect(+c[4].toFixed(6)).to.be(514);
    expect(+c[5].toFixed(6)).to.be(95);
    expect(+c[6].toFixed(6)).to.be(210);
    expect(+c[7].toFixed(6)).to.be(31);
    expect(+c[8].toFixed(6)).to.be(47.5);
  });

  it('calculate C := a A B^t + b C', function() {
    var heap = new ArrayBuffer(216),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 9),
        c = new Float64Array(heap, 144, 9),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 8;
    a[2] = 3;
    a[3] = 2;
    a[4] = 10;
    a[5] = 8;
    a[6] = 9;
    a[7] = -5;
    a[8] = -1;
    b[0] = 9;
    b[1] = 3;
    b[2] = -8;
    b[3] = 8;
    b[4] = 11;
    b[5] = 6;
    b[6] = 3;
    b[7] = 2.3;
    b[8] = 1;
    c[0] = 3;
    c[1] = 3;
    c[2] = 1.2;
    c[3] = 8;
    c[4] = 4;
    c[5] = 8;
    c[6] = 6;
    c[7] = 1;
    c[8] = -2;

    linalg.dgemm(0, 1, 3, 3, 3, 3.0, a.byteOffset, 3, b.byteOffset, 3, -2.0, c.byteOffset, 3);

    expect(+c[0].toFixed(6)).to.be(21);
    expect(+c[1].toFixed(6)).to.be(336);
    expect(+c[2].toFixed(6)).to.be(70.8);
    expect(+c[3].toFixed(6)).to.be(-64);
    expect(+c[4].toFixed(6)).to.be(514);
    expect(+c[5].toFixed(6)).to.be(95);
    expect(+c[6].toFixed(6)).to.be(210);
    expect(+c[7].toFixed(6)).to.be(31);
    expect(+c[8].toFixed(6)).to.be(47.5);
  });

  it('calculate C := a A^t B^t + b C', function() {
    var heap = new ArrayBuffer(216),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 9),
        c = new Float64Array(heap, 144, 9),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 2;
    a[2] = 9;
    a[3] = 8;
    a[4] = 10;
    a[5] = -5;
    a[6] = 3;
    a[7] = 8;
    a[8] = -1;
    b[0] = 9;
    b[1] = 3;
    b[2] = -8;
    b[3] = 8;
    b[4] = 11;
    b[5] = 6;
    b[6] = 3;
    b[7] = 2.3;
    b[8] = 1;
    c[0] = 3;
    c[1] = 3;
    c[2] = 1.2;
    c[3] = 8;
    c[4] = 4;
    c[5] = 8;
    c[6] = 6;
    c[7] = 1;
    c[8] = -2;

    linalg.dgemm(1, 1, 3, 3, 3, 3.0, a.byteOffset, 3, b.byteOffset, 3, -2.0, c.byteOffset, 3);

    expect(+c[0].toFixed(6)).to.be(21);
    expect(+c[1].toFixed(6)).to.be(336);
    expect(+c[2].toFixed(6)).to.be(70.8);
    expect(+c[3].toFixed(6)).to.be(-64);
    expect(+c[4].toFixed(6)).to.be(514);
    expect(+c[5].toFixed(6)).to.be(95);
    expect(+c[6].toFixed(6)).to.be(210);
    expect(+c[7].toFixed(6)).to.be(31);
    expect(+c[8].toFixed(6)).to.be(47.5);
  });
});
