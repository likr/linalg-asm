function LinalgModule(stdlib, foreign, heap) {
  "use asm";

  var abs = stdlib.Math.abs,
      imul = stdlib.Math.imul,
      uiarray = new stdlib.Uint32Array(heap),
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

  function dswap(n, x, incX, y, incY) {
    n = n | 0;
    x = x | 0;
    incX = incX | 0;
    y = y | 0;
    incY = incY | 0;

    var i = 0,
        j = 0,
        tmp = 0.0;

    n = x + n << 3;
    x = x << 3;
    incX = incX << 3;
    y = y << 3;
    incY = incY << 3;

    for (i = x, j = y; (i | 0) < (n | 0); i = i + incX | 0, j = j + incY | 0) {
      tmp = +darray[i >> 3];
      darray[i >> 3] = darray[j >> 3];
      darray[j >> 3] = tmp;
    }
  }

  function idamax(n, x, incX) {
    n = n | 0;
    x = x | 0;
    incX = incX | 0;

    var i = 0,
        index = 0,
        tmp = 0.0,
        value = 0.0;

    n = (imul(n, incX) | 0) << 3;
    x = x << 3;
    incX = incX << 3;

    for (i = x; (i | 0) < (n | 0); i = i + incX | 0) {
      tmp = +abs(darray[i >> 3]);
      if (tmp > value) {
        index = (i - x | 0) / (incX | 0) |  0;
        value = tmp;
      }
    }

    return index | 0;
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

  function dgetrf(m, n, a, lda, ipiv) {
    // TODO support lda
    m = m | 0;
    n = n | 0;
    a = a | 0;
    lda = lda | 0;
    ipiv = ipiv | 0;

    var i = 0,
        j = 0,
        k = 0,
        l = 0,
        jk = 0,
        jl = 0,
        kk = 0,
        kl = 0;

    for (k = 0; (k | 0) < (n - 1 | 0); k = k + 1 | 0) {
      i = k + (idamax(n - k | 0, a + (imul(i, n) | 0) + k | 0, m) | 0) | 0;
      uiarray[ipiv + (k << 2) >> 2] = i;
      dswap(n, a + (imul(i, n) | 0) | 0, 1, a + (imul(k, n) | 0) | 0, 1);
      for (j = k + 1 | 0; (j | 0) < (n | 0); j = j + 1 | 0) {
        jk = a + (imul(j, n) | 0) + k << 3;
        kk = a + (imul(k, n) | 0) + k << 3;
        darray[jk >> 3] = darray[jk >> 3] / darray[kk >> 3];
        for (l = k + 1 | 0; (l | 0) < (n | 0); l = l + 1 | 0) {
          jl = a + (imul(j, n) | 0) + l << 3;
          kl = a + (imul(k, n) | 0) + l << 3;
          darray[jl >> 3] = darray[jl >> 3] - darray[jk >> 3] * darray[kl >> 3];
        }
      }
    }
    uiarray[ipiv + (n - 1 << 2) >> 2] = n - 1;
  }

  return {
    daxpy: daxpy,
    dasum: dasum,
    dcopy: dcopy,
    ddot: ddot,
    dswap: dswap,
    idamax: idamax,
    dgemm: dgemm,
    dgetrf: dgetrf
  };
}

module.exports = LinalgModule;
