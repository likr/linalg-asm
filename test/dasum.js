var expect = require('expect.js'),
    linalgModule = require('../index');

describe('dasum', function() {
  it('return |x_1| + ... + |x_n|', function() {
    var heap = new ArrayBuffer(32),
        x = new Float64Array(heap, 0, 4),
        linalg = linalgModule(global, null, heap);

    x[0] = 0.5;
    x[1] = -1.5;
    x[2] = 2.5;
    x[3] = -3.5;

    var result = linalg.dasum(4, 0, 1);

    expect(result).to.be(8);
  });
});
