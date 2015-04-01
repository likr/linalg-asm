var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dtrtri', function() {
  it('computes U ^ -1', function () {
    var heap = new ArrayBuffer(72),
        a = new Float64Array(heap, 0, 9),
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

    linalg.dtrtri(0, 0, 3, a.byteOffset, 3);

    expect(a[0]).to.be(1);
    expect(a[1]).to.be(-0.8);
    expect(a[2]).to.be(32 / 15);
    expect(a[3]).to.be(0);
    expect(a[4]).to.be(0.2);
    expect(a[5]).to.be(2 / 15);
    expect(a[6]).to.be(0);
    expect(a[7]).to.be(0);
    expect(a[8]).to.be(-1 / 3);
  });

  it('computes L ^ -1', function () {
    var heap = new ArrayBuffer(72),
        a = new Float64Array(heap, 0, 9),
        linalg = linalgModule(global, null, heap);

    a[0] = 1;
    a[1] = 0;
    a[2] = 0;
    a[3] = 3;
    a[4] = 4;
    a[5] = 0;
    a[6] = 5;
    a[7] = 8;
    a[8] = -1;

    linalg.dtrtri(1, 0, 3, a.byteOffset, 3);

    expect(a[0]).to.be(1);
    expect(a[1]).to.be(0);
    expect(a[2]).to.be(0);
    expect(a[3]).to.be(-0.75);
    expect(a[4]).to.be(0.25);
    expect(a[5]).to.be(0);
    expect(a[6]).to.be(-1);
    expect(a[7]).to.be(2);
    expect(a[8]).to.be(-1);
  });
});
