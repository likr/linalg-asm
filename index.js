function LinalgModule(stdlib, foreign, heap) {
  "use asm";

  var abs = stdlib.Math.abs,
      imul = stdlib.Math.imul,
      darray = new stdlib.Float64Array(heap);

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
      darray[j >> 3] = alpha * darray[i >> 3] + (+darray[j >> 3]);
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
      value = value + (+abs(darray[i >> 3]));
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

  function dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) {
    // TODO support transa, transb, lda, ldb, ldc
    transa = transa | 0;
    transb = transb | 0;
    m = m | 0;
    n = n | 0;
    k = k | 0;
    alpha = +alpha;
    a = a | 0;
    lda = lda | 0;
    b = b | 0;
    ldb = ldb | 0;
    beta = +beta;
    c = c | 0;
    ldc = ldc | 0;

    var i = 0,
        j = 0,
        l = 0,
        aIndex = 0,
        bIndex = 0,
        cIndex = 0;

    for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
      for (j = 0; (j | 0) < (n | 0); j = j + 1 | 0) {
        cIndex = c + (imul(i, n) | 0) + j << 3;
        darray[cIndex >> 3] = beta * darray[cIndex >> 3];
        for (l = 0; (l | 0) < (k | 0); l = l + 1 | 0) {
          aIndex = a + (imul(i, k) | 0) + l << 3;
          bIndex = b + (imul(l, k) | 0) + j << 3;
          darray[cIndex >> 3] = +darray[cIndex >> 3] + alpha * darray[aIndex >> 3] * darray[bIndex >> 3];
        }
      }
    }
  }

  return {
    daxpy: daxpy,
    dasum: dasum,
    dcopy: dcopy,
    ddot: ddot,
    dgemm: dgemm
  };
}

module.exports = LinalgModule;
