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

  function dcopy(n, x, incx, y, incy) {
    n = n | 0;
    x = x | 0;
    incx = incx | 0;
    y = y | 0;
    incy = incy | 0;

    var i = 0,
        pxi = 0,
        pyi = 0;

    incx = incx << 3;
    incy = incy << 3;

    for (i = 0, pxi = x, pyi = y; (i | 0) < (n | 0); i = i + 1 | 0, pxi = pxi + incx | 0, pyi = pyi + incy | 0) {
      darray[pyi >> 3] = darray[pxi >> 3];
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

  function dswap(n, x, incx, y, incy) {
    n = n | 0;
    x = x | 0;
    incx = incx | 0;
    y = y | 0;
    incy = incy | 0;

    var i = 0,
        pxi = 0,
        pyi = 0,
        tmp = 0.0;

    incx = incx << 3;
    incy = incy << 3;

    for (i = 0, pxi = x, pyi = y; (i | 0) < (n | 0); i = i + 1 | 0, pxi = pxi + incx | 0, pyi = pyi + incy | 0) {
      tmp = +darray[pxi >> 3];
      darray[pxi >> 3] = darray[pyi >> 3];
      darray[pyi >> 3] = tmp;
    }
  }

  function idamax(n, x, incx) {
    n = n | 0;
    x = x | 0;
    incx = incx | 0;

    var i = 0,
        pxi = 0,
        index = 0,
        tmp = 0.0,
        value = 0.0;

    incx = incx << 3;

    for (i = 0, pxi = x; (i | 0) < (n | 0); i = i + 1 | 0, pxi = pxi + incx | 0) {
      tmp = +abs(darray[pxi >> 3]);
      if (tmp > value) {
        index = i;
        value = tmp;
      }
    }

    return index | 0;
  }

  function dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy) {
    // TODO support lda
    trans = trans | 0;
    m = m | 0;
    n = n | 0;
    alpha = +alpha;
    a = a | 0;
    lda = lda | 0;
    x = x | 0;
    incx = incx | 0;
    beta = +beta;
    y = y | 0;
    incy = incy | 0;

    var i = 0,
        j = 0,
        paij = 0,
        pxj = 0,
        pyi = 0,
        val = 0.0;

    for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
      pyi = y + ((imul(i, incx) | 0) << 3) | 0;
      val = beta * darray[pyi >> 3];
      for (j = 0; (j | 0) < (n | 0); j = j + 1 | 0) {
        paij = a + ((imul(i, n) | 0) + j << 3) | 0;
        pxj = x + ((imul(j, incy) | 0) << 3) | 0;
        val = val + alpha * darray[paij >> 3] * darray[pxj >> 3];
      }
      darray[pyi >> 3] = val;
    }
  }

  function dtrsv(uplo, trans, diag, n, a, lda, x, incx) {
    // TODO support trans, lda
    uplo = uplo | 0;
    trans = trans | 0;
    diag = diag | 0;
    n = n | 0;
    a = a | 0;
    lda = lda | 0;
    x = x | 0;
    incx = incx | 0;

    var i = 0,
        j = 0,
        paii = 0,
        paji = 0,
        pxi = 0,
        pxj = 0;

    if (uplo) {
      // lower triangular matrix
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        pxi = x + ((imul(i, incx) | 0) << 3) | 0;
        if (!diag) {
          paii = a + ((imul(i, n) | 0) + i << 3) | 0;
          darray[pxi >> 3] = darray[pxi >> 3] / darray[paii >> 3];
        }
        for (j = i + 1 | 0; (j | 0) < (n | 0); j = j + 1 | 0) {
          paji = a + ((imul(j, n) | 0) + i << 3) | 0;
          pxj = x + ((imul(j, incx) | 0) << 3) | 0;
          darray[pxj >> 3] = darray[pxj >> 3] - darray[pxi >> 3] * darray[paji >> 3];
        }
      }
    } else {
      // upper triangular matrix
      for (i = n - 1 | 0; (i | 0) >= 0; i = i - 1 | 0) {
        pxi = x + ((imul(i, incx) | 0) << 3) | 0;
        if (!diag) {
          paii = a + ((imul(i, n) | 0) + i << 3) | 0;
          darray[pxi >> 3] = darray[pxi >> 3] / darray[paii >> 3];
        }
        for (j = 0; (j | 0) < (i | 0); j = j + 1 | 0) {
          paji = a + ((imul(j, n) | 0) + i << 3) | 0;
          pxj = x + ((imul(j, incx) | 0) << 3) | 0;
          darray[pxj >> 3] = darray[pxj >> 3] - darray[pxi >> 3] * darray[paji >> 3];
        }
      }
    }
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

  function dgesv(n, nrhs, a, lda, ipiv, b, ldb) {
    n = n | 0;
    nrhs = nrhs | 0;
    a = a | 0;
    lda = lda | 0;
    ipiv = ipiv | 0;
    b = b | 0;
    ldb = ldb | 0;

    var i = 0,
        j = 0,
        tmp = 0.0,
        pipivi = 0,
        pbi = 0,
        pbj = 0;

    dgetrf(n, n, a, lda, ipiv);

    for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
      pipivi = ipiv + (i << 2) | 0;
      j = uiarray[pipivi >> 2] | 0;
      pbi = b + (i << 3) | 0;
      pbj = b + (j << 3) | 0;
      tmp = +darray[pbi >> 3];
      darray[pbi >> 3] = darray[pbj >> 3];
      darray[pbj >> 3] = tmp;
    }

    dtrsv(1, 0, 1, n, a, lda, b, 1);
    dtrsv(0, 0, 0, n, a, lda, b, 1);
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
        pipivk = 0,
        paik = 0,
        pajk = 0,
        pakk = 0,
        pajl = 0,
        pakl = 0;

    for (k = 0; (k | 0) < (n - 1 | 0); k = k + 1 | 0) {
      i = k + (idamax(n - k | 0, a + ((imul(k, n) | 0) + k << 3) | 0, m) | 0) | 0;
      pipivk = ipiv + (k << 2) | 0;
      uiarray[pipivk >> 2] = i;
      dswap(n, a + ((imul(i, n) | 0) << 3) | 0, 1, a + ((imul(k, n) | 0) << 3) | 0, 1);
      for (j = k + 1 | 0; (j | 0) < (n | 0); j = j + 1 | 0) {
        pajk = a + ((imul(j, n) | 0) + k << 3) | 0;
        pakk = a + ((imul(k, n) | 0) + k << 3) | 0;
        darray[pajk >> 3] = darray[pajk >> 3] / darray[pakk >> 3];
        for (l = k + 1 | 0; (l | 0) < (n | 0); l = l + 1 | 0) {
          pajl = a + ((imul(j, n) | 0) + l << 3) | 0;
          pakl = a + ((imul(k, n) | 0) + l << 3) | 0;
          darray[pajl >> 3] = darray[pajl >> 3] - darray[pajk >> 3] * darray[pakl >> 3];
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
    dgemv: dgemv,
    dtrsv: dtrsv,
    dgemm: dgemm,
    dgesv: dgesv,
    dgetrf: dgetrf
  };
}

module.exports = LinalgModule;
