function LinalgModule(stdlib, foreign, heap) {
  "use asm";

  var darray = new stdlib.Float64Array(heap);

  function daxpy(n, alpha, x, incX, y, incY) {
    n = n | 0;
    alpha = +alpha;
    x = x | 0;
    incX = incX | 0;
    y = y | 0;
    incY = incY | 0;

    n = n << 3;
    x = x << 3;
    incX = incX << 3;
    y = y << 3;
    incY = incY << 3;

    var i = 0,
        j = 0;

    for (i = x, j = y; (i | 0) < (n | 0); i = i + incX, j = j + incY) {
      darray[j >> 3] = alpha * darray[i >> 3] + darray[j >> 3];
    }
  }

  return {
    daxpy: daxpy
  };
}

module.exports = LinalgModule;
