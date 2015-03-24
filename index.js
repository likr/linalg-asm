function LinalgModule(stdlib, foreign, heap) {
  "use asm";

  var abs = stdlib.Math.abs;
  var darray = new stdlib.Float64Array(heap);

  function daxpy(n, alpha, x, incX, y, incY) {
    n = n | 0;
    alpha = +alpha;
    x = x | 0;
    incX = incX | 0;
    y = y | 0;
    incY = incY | 0;

    var i = 0,
        j = 0;

    n = n << 3;
    x = x << 3;
    incX = incX << 3;
    y = y << 3;
    incY = incY << 3;

    for (i = x, j = y; (i | 0) < (n | 0); i = i + incX | 0, j = j + incY | 0) {
      darray[j >> 3] = alpha * darray[i >> 3] + +darray[j >> 3];
    }
  }

  function dasum(n, x, incX) {
    n = n | 0;
    x = x | 0;
    incX = incX | 0;

    var i = 0,
        value = 0.0;

    n = n << 3;
    x = x << 3;
    incX = incX << 3;

    for (i = x; (i | 0) < (n | 0); i = i + incX | 0) {
      value = value + +abs(darray[i >> 3]);
    }

    return value;
  }

  function dcopy(n, x, incX, y, incY) {
    n = n | 0;
    x = x | 0;
    incX = incX | 0;
    y = y | 0;
    incY = incY | 0;

    var i = 0,
        j = 0;

    n = n << 3;
    x = x << 3;
    incX = incX << 3;
    y = y << 3;
    incY = incY << 3;

    for (i = x, j = y; (i | 0) < (n | 0); i = i + incX | 0, j = j + incY | 0) {
      darray[j >> 3] = darray[i >> 3];
    }
  }

  function ddot(n, x, incX, y, incY) {
    n = n | 0;
    x = x | 0;
    incX = incX | 0;
    y = y | 0;
    incY = incY | 0;

    var i = 0,
        j = 0,
        value = 0.0;

    n = n << 3;
    x = x << 3;
    incX = incX << 3;
    y = y << 3;
    incY = incY << 3;

    for (i = x, j = y; (i | 0) < (n | 0); i = i + incX | 0, j = j + incY | 0) {
      value = value + darray[i >> 3] * darray[j >> 3];
    }

    return value;
  }

  return {
    daxpy: daxpy,
    dasum: dasum,
    dcopy: dcopy,
    ddot: ddot
  };
}

module.exports = LinalgModule;
