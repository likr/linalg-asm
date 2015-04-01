var expect = require('expect.js'),
    linalgModule = require('../index');

describe('idamax', function() {
  it('returns index of max abs value', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = linalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.idamax(4, x.byteOffset, 1);

    expect(result).to.be(3);
  });

  it('returns index of max abs value with incX', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = linalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.idamax(2, x.byteOffset, 2);

    expect(result).to.be(1);
  });

  it('returns index of max abs value with non zero offset', function() {
    var heap = new ArrayBuffer(48),
        x = new Float64Array(heap, 16, 4),
        linalg = linalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.idamax(2, x.byteOffset, 1);

    expect(result).to.be(1);
  });
});
