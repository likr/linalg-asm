var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dtrmv', function() {
  it('computes x := U x', function () {
    var heap = new ArrayBuffer(96),
        a = new Float64Array(heap, 0, 9),
        x = new Float64Array(heap, 72, 3),
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
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    linalg.dtrmv(0, 0, 0, 3, a.byteOffset, 3, x.byteOffset, 1);

    expect(x[0]).to.be(33);
    expect(x[1]).to.be(16);
    expect(x[2]).to.be(-9);
  });

  it('computes x := L x', function () {
    var heap = new ArrayBuffer(96),
        a = new Float64Array(heap, 0, 9),
        x = new Float64Array(heap, 72, 3),
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
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;

    linalg.dtrmv(1, 0, 0, 3, a.byteOffset, 3, x.byteOffset, 1);

    expect(x[0]).to.be(1);
    expect(x[1]).to.be(14);
    expect(x[2]).to.be(3);
  });
});
