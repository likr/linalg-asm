var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dtrsm', function() {
  it('solves U X = B', function() {
    var heap = new ArrayBuffer(120),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 4;
    a[2] = 8;
    a[3] = 0;
    a[4] = 5;
    a[5] = 2;
    a[6] = 0;
    a[7] = 0;
    a[8] = -3;
    b[0] = 33;
    b[1] = 18;
    b[2] = 16;
    b[3] = 9;
    b[4] = -9;
    b[5] = -6;

    linalg.dtrsm(0, 0, 0, 0, 3, 2, 1, a.byteOffset, 3, b.byteOffset, 3);

    expect(b[0]).to.be(1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be(2);
    expect(b[3]).to.be(1);
    expect(b[4]).to.be(3);
    expect(b[5]).to.be(2);
  });

  it('solves L X = B', function() {
    var heap = new ArrayBuffer(120),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 0;
    a[2] = 0;
    a[3] = 4;
    a[4] = 5;
    a[5] = 0;
    a[6] = 8;
    a[7] = 2;
    a[8] = -3;
    b[0] = 1;
    b[1] = -2;
    b[2] = 14;
    b[3] = -3;
    b[4] = 3;
    b[5] = -20;

    linalg.dtrsm(0, 1, 0, 0, 3, 2, 1, a.byteOffset, 3, b.byteOffset, 3);

    expect(b[0]).to.be(1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be(2);
    expect(b[3]).to.be(1);
    expect(b[4]).to.be(3);
    expect(b[5]).to.be(2);
  });

  it('solves U X = B where U_ii = 1', function() {
    var heap = new ArrayBuffer(120),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = linalgModule(global, null, heap);

    a[0] = 0;
    a[1] = -2;
    a[2] = 3;
    a[3] = 0;
    a[4] = 0;
    a[5] = 4;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    b[0] = 6;
    b[1] = 2;
    b[2] = 14;
    b[3] = 9;
    b[4] = 3;
    b[5] = 2;

    linalg.dtrsm(0, 0, 0, 1, 3, 2, 1, a.byteOffset, 3, b.byteOffset, 3);

    expect(b[0]).to.be(1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be(2);
    expect(b[3]).to.be(1);
    expect(b[4]).to.be(3);
    expect(b[5]).to.be(2);
  });

  it('solves L X = B where L_ii = 1', function() {
    var heap = new ArrayBuffer(120),
        a = new Float64Array(heap, 0, 9),
        b = new Float64Array(heap, 72, 6),
        linalg = linalgModule(global, null, heap);

    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
    a[3] = -2;
    a[4] = 0;
    a[5] = 0;
    a[6] = 3;
    a[7] = 4;
    a[8] = 0;
    b[0] = 1;
    b[1] = -2;
    b[2] = 0;
    b[3] = 5;
    b[4] = 14;
    b[5] = 0;

    linalg.dtrsm(0, 1, 0, 1, 3, 2, 1, a.byteOffset, 3, b.byteOffset, 3);

    expect(b[0]).to.be(1);
    expect(b[1]).to.be(-2);
    expect(b[2]).to.be(2);
    expect(b[3]).to.be(1);
    expect(b[4]).to.be(3);
    expect(b[5]).to.be(2);
  });
});
