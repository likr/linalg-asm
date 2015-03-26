var expect = require('expect.js');

var LinalgModule = require('../index');

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
