var expect = require('expect.js');

var LinalgModule = require('../index');

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
