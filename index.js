function LinalgModule(stdlib, foreign, heap) {
  'use asm';

  var abs = stdlib.Math.abs,
      imul = stdlib.Math.imul,
      sqrt = stdlib.Math.sqrt,
      uiarray = new stdlib.Uint32Array(heap),
      darray = new stdlib.Float64Array(heap);

  function daxpy(n, alpha, x, incx, y, incy) {
    n = n | 0;
    alpha = +alpha;
    x = x | 0;
    incx = incx | 0;
    y = y | 0;
    incy = incy | 0;

    var i = 0,
        nloop = 0,
        pxi = 0,
        pyi = 0,
        dx1 = 0,
        dx2 = 0,
        dx3 = 0,
        dy1 = 0,
        dy2 = 0,
        dy3 = 0;

    if (alpha == 0.0) {
      return;
    }

    nloop = n >> 2;
    incx = incx << 3;
    incy = incy << 3;
    dx1 = incx;
    dx2 = dx1 + incx | 0;
    dx3 = dx2 + incx | 0;
    dy1 = incy;
    dy2 = dy1 + incy | 0;
    dy3 = dy2 + incy | 0;

    incx = incx << 2;
    incy = incy << 2;
    for (i = 0, pxi = x, pyi = y; (i | 0) < (nloop | 0); i = i + 1 | 0, pxi = pxi + incx | 0, pyi = pyi + incy | 0) {
      darray[pyi >> 3] = +darray[pyi >> 3] + alpha * darray[pxi >> 3];
      darray[pyi + dy1 >> 3] = +darray[pyi + dy1 >> 3] + alpha * darray[pxi + dx1 >> 3];
      darray[pyi + dy2 >> 3] = +darray[pyi + dy2 >> 3] + alpha * darray[pxi + dx2 >> 3];
      darray[pyi + dy3 >> 3] = +darray[pyi + dy3 >> 3] + alpha * darray[pxi + dx3 >> 3];
    }

    i = i << 2;
    for (incx = incx >> 2, incy = incy >> 2; (i | 0) < (n | 0); i = i + 1 | 0, pxi = pxi + incx | 0, pyi = pyi + incy | 0) {
      darray[pyi >> 3] = +darray[pyi >> 3] + alpha * darray[pxi >> 3];
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

  function ddot(n, x, incx, y, incy) {
    n = n | 0;
    x = x | 0;
    incx = incx | 0;
    y = y | 0;
    incy = incy | 0;

    var i = 0,
        nloop = 0,
        pxi = 0,
        pyi = 0,
        dx1 = 0,
        dx2 = 0,
        dx3 = 0,
        dy1 = 0,
        dy2 = 0,
        dy3 = 0,
        value = 0.0;

    nloop = n >> 2;
    incx = incx << 3;
    incy = incy << 3;
    dx1 = incx;
    dx2 = dx1 + incx | 0;
    dx3 = dx2 + incx | 0;
    dy1 = incy;
    dy2 = dy1 + incy | 0;
    dy3 = dy2 + incy | 0;

    incx = incx << 2;
    incy = incy << 2;
    for (i = 0, pxi = x, pyi = y; (i | 0) < (nloop | 0); i = i + 1 | 0, pxi = pxi + incx | 0, pyi = pyi + incy | 0) {
      value = value + darray[pyi >> 3] * darray[pxi >> 3];
      value = value + darray[pyi + dy1 >> 3] * darray[pxi + dx1 >> 3];
      value = value + darray[pyi + dy2 >> 3] * darray[pxi + dx2 >> 3];
      value = value + darray[pyi + dy3 >> 3] * darray[pxi + dx3 >> 3];
    }

    i = i << 2;
    for (incx = incx >> 2, incy = incy >> 2; (i | 0) < (n | 0); i = i + 1 | 0, pxi = pxi + incx | 0, pyi = pyi + incy | 0) {
      value = value + darray[pyi >> 3] * darray[pxi >> 3];
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
        pai0 = 0,
        pyi = 0;

    if (trans) {
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        pyi = y + ((imul(i, incy) | 0) << 3) | 0;
        pai0 = a + (i << 3) | 0;
        darray[pyi >> 3] = alpha * +ddot(m, pai0, lda, x, incx) + beta * darray[pyi >> 3];
      }
    } else {
      for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
        pyi = y + ((imul(i, incy) | 0) << 3) | 0;
        pai0 = a + ((imul(i, lda) | 0) << 3) | 0;
        darray[pyi >> 3] = alpha * +ddot(n, pai0, 1, x, incx) + beta * darray[pyi >> 3];
      }
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
          darray[pxi >> 3] = +darray[pxi >> 3] / +darray[paii >> 3];
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
          darray[pxi >> 3] = +darray[pxi >> 3] / +darray[paii >> 3];
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
        pai0 = 0,
        pb0j = 0,
        pcij = 0;

    if (transa) {
      if (transb) {
        for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
          for (j = 0; (j | 0) < (n | 0); j = j + 1 | 0) {
            pai0 = a + (i << 3) | 0;
            pb0j = b + ((imul(j, ldb) | 0) << 3) | 0;
            pcij = c + ((imul(i, ldc) | 0) + j << 3) | 0;
            darray[pcij >> 3] = beta * darray[pcij >> 3] + alpha * +ddot(k, pai0, lda, pb0j, 1);
          }
        }
      } else {
        for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
          for (j = 0; (j | 0) < (n | 0); j = j + 1 | 0) {
            pai0 = a + (i << 3) | 0;
            pb0j = b + (j << 3) | 0;
            pcij = c + ((imul(i, ldc) | 0) + j << 3) | 0;
            darray[pcij >> 3] = beta * darray[pcij >> 3] + alpha * +ddot(k, pai0, lda, pb0j, ldb);
          }
        }
      }
    } else {
      if (transb) {
        for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
          for (j = 0; (j | 0) < (n | 0); j = j + 1 | 0) {
            pai0 = a + ((imul(i, lda) | 0) << 3) | 0;
            pb0j = b + ((imul(j, ldb) | 0) << 3) | 0;
            pcij = c + ((imul(i, ldc) | 0) + j << 3) | 0;
            darray[pcij >> 3] = beta * darray[pcij >> 3] + alpha * +ddot(k, pai0, 1, pb0j, 1);
          }
        }
      } else {
        for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
          for (j = 0; (j | 0) < (n | 0); j = j + 1 | 0) {
            pai0 = a + ((imul(i, lda) | 0) << 3) | 0;
            pb0j = b + (j << 3) | 0;
            pcij = c + ((imul(i, ldc) | 0) + j << 3) | 0;
            darray[pcij >> 3] = beta * darray[pcij >> 3] + alpha * +ddot(k, pai0, 1, pb0j, ldb);
          }
        }
      }
    }
  }

  function dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) {
    // TODO support side
    side = side | 0;
    uplo = uplo | 0;
    transa = transa | 0;
    diag = diag | 0;
    m = m | 0;
    n = n | 0;
    alpha = +alpha;
    a = a | 0;
    lda = lda | 0;
    b = b | 0;
    ldb = ldb | 0;

    var i = 0,
        j = 0,
        k = 0,
        paii = 0,
        paji = 0,
        pbik = 0,
        pbjk = 0;

    if (transa) {
      if (uplo) {
        // solve a U^t x = B
        for (i = m - 1 | 0; (i | 0) >= 0; i = i - 1 | 0) {
          if (!diag) {
            paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              darray[pbik >> 3] = +darray[pbik >> 3] / +darray[paii >> 3];
            }
          }
          for (j = 0; (j | 0) < (i | 0); j = j + 1 | 0) {
            paji = a + ((imul(i, lda) | 0) + j << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              pbjk = b + ((imul(j, ldb) | 0) + k << 3) | 0;
              darray[pbjk >> 3] = darray[pbjk >> 3] - darray[pbik >> 3] * darray[paji >> 3];
            }
          }
        }
      } else {
        // solve a L^t x = B
        for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
          if (!diag) {
            paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              darray[pbik >> 3] = +darray[pbik >> 3] / +darray[paii >> 3];
            }
          }
          for (j = i + 1 | 0; (j | 0) < (m | 0); j = j + 1 | 0) {
            paji = a + ((imul(i, lda) | 0) + j << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              pbjk = b + ((imul(j, ldb) | 0) + k << 3) | 0;
              darray[pbjk >> 3] = darray[pbjk >> 3] - darray[pbik >> 3] * darray[paji >> 3];
            }
          }
        }
      }
    } else {
      if (uplo) {
        // solve a L x = B
        for (i = 0; (i | 0) < (m | 0); i = i + 1 | 0) {
          if (!diag) {
            paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              darray[pbik >> 3] = +darray[pbik >> 3] / +darray[paii >> 3];
            }
          }
          for (j = i + 1 | 0; (j | 0) < (m | 0); j = j + 1 | 0) {
            paji = a + ((imul(j, lda) | 0) + i << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              pbjk = b + ((imul(j, ldb) | 0) + k << 3) | 0;
              darray[pbjk >> 3] = darray[pbjk >> 3] - darray[pbik >> 3] * darray[paji >> 3];
            }
          }
        }
      } else {
        // solve a U x = B
        for (i = m - 1 | 0; (i | 0) >= 0; i = i - 1 | 0) {
          if (!diag) {
            paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              darray[pbik >> 3] = +darray[pbik >> 3] / +darray[paii >> 3];
            }
          }
          for (j = 0; (j | 0) < (i | 0); j = j + 1 | 0) {
            paji = a + ((imul(j, lda) | 0) + i << 3) | 0;
            for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
              pbik = b + ((imul(i, ldb) | 0) + k << 3) | 0;
              pbjk = b + ((imul(j, ldb) | 0) + k << 3) | 0;
              darray[pbjk >> 3] = darray[pbjk >> 3] - darray[pbik >> 3] * darray[paji >> 3];
            }
          }
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

  function dgetri(n, a, lda, ipiv, work, lwork) {
    n = n | 0;
    a = a | 0;
    lda = lda | 0;
    ipiv = ipiv | 0;
    work = work | 0;
    lwork = lwork | 0;

    var i = 0,
        j = 0,
        pworkj = 0,
        paij = 0;

    dtrtri(0, 0, n, a, lda);

    for (i = n - 2 | 0; (i | 0) >= 0; i = i - 1 | 0) {
      for (j = i + 1 | 0; (j | 0) < (n | 0); j = j + 1 | 0) {
        pworkj = work + (j << 3) | 0;
        paij = a + ((imul(j, lda) | 0) + i << 3) | 0;
        darray[pworkj >> 3] = darray[paij >> 3];
        darray[paij >> 3] = 0.0;
      }
      dgemv(0, n, n - i - 1 | 0, -1.0, a + (i + 1 << 3) | 0, lda, work + (i + 1 << 3) | 0, 1, 1.0, a + (i << 3) | 0, lda);
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
        jn = 0,
        kn = 0,
        pipivk = 0,
        pajk = 0,
        pakk = 0,
        ajk = 0.0,
        dakk = 0.0,
        inca = 0;
    inca = n << 3;

    for (k = 0, kn = 0, pakk = a; (k | 0) < (n - 1 | 0); k = k + 1 | 0, kn = kn + n | 0, pakk = pakk + inca + 8 | 0) {
      i = k + (idamax(n - k | 0, pakk, m) | 0) | 0;
      pipivk = ipiv + (k << 2) | 0;
      uiarray[pipivk >> 2] = i;
      if ((i | 0) != (k | 0)) {
        dswap(n, a + ((imul(i, n) | 0) << 3) | 0, 1, a + (kn << 3) | 0, 1);
      }
      dakk = 1.0 / +darray[pakk >> 3];
      for (j = k + 1 | 0, jn = imul(j, n) | 0, pajk = a + (jn + k << 3) | 0; (j | 0) < (n | 0); j = j + 1 | 0, jn = jn + n | 0, pajk = pajk + inca | 0) {
        ajk = darray[pajk >> 3] = darray[pajk >> 3] * dakk;
        daxpy(n - k - 1 | 0, -ajk,
              a + (kn + k + 1 << 3) | 0, 1,
              a + (jn + k + 1 << 3) | 0, 1);
      }
    }
    uiarray[ipiv + (n - 1 << 2) >> 2] = n - 1;
  }

  function dposv(uplo, n, nrhs, a, lda, b, ldb) {
    uplo = uplo | 0;
    n = n | 0;
    nrhs = nrhs | 0;
    a = a | 0;
    lda = lda | 0;
    b = b | 0;
    ldb = ldb | 0;

    dpotrf(uplo, n, a, lda);
    dpotrs(uplo, n, nrhs, a, lda, b, ldb);
  }

  function dpotrf(uplo, n, a, lda) {
    uplo = uplo | 0;
    n = n | 0;
    a = a | 0;
    lda = lda | 0;

    var i = 0,
        paii = 0,
        pai = 0,
        aii = 0.0;

    if (uplo) {
      // lower triangular matrix
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
        pai = a + ((imul(i, lda) | 0) << 3) | 0;
        aii = +darray[paii >> 3];
        aii = darray[paii >> 3] = +sqrt(aii - +ddot(i | 0, pai, 1, pai, 1));
        if ((n - i - 1 | 0) > 0) {
          dgemv(0, n - i - 1 | 0, i, -1.0, pai + (lda << 3) | 0, 1, pai, 1, 1.0, paii + (lda << 3) | 0, lda);
          dscal(n - i - 1 | 0, 1.0 / aii, paii + (lda << 3) | 0, lda);
        }
      }
    } else {
      // upper triangular matrix
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
        pai = a + (i << 3) | 0;
        aii = +darray[paii >> 3];
        aii = darray[paii >> 3] = +sqrt(aii - +ddot(i | 0, pai, lda, pai, lda));
        if ((n - i - 1 | 0) > 0) {
          dgemv(1, i, n - i - 1 | 0, -1.0, pai + 8 | 0, lda, pai, lda, 1.0, paii + 8 | 0, 1);
          dscal(n - i - 1 | 0, 1.0 / aii, paii + 8 | 0, lda);
        }
      }
    }
  }

  function dpotrs(uplo, n, nrhs, a, lda, b, ldb) {
    uplo = uplo | 0;
    n = n | 0;
    nrhs = nrhs | 0;
    a = a | 0;
    lda = lda | 0;
    b = b | 0;
    ldb = ldb | 0;

    if (uplo) {
      // lower triangular matrix
      dtrsm(0, 1, 0, 0, n, nrhs, 1.0, a, lda, b, ldb);
      dtrsm(0, 1, 1, 0, n, nrhs, 1.0, a, lda, b, ldb);
    } else {
      dtrsm(0, 0, 1, 0, n, nrhs, 1.0, a, lda, b, ldb);
      dtrsm(0, 0, 0, 0, n, nrhs, 1.0, a, lda, b, ldb);
    }
  }

  function dscal(n, a, x, incx) {
    n = n | 0;
    a = +a;
    x = x | 0;
    incx = incx | 0;

    var i = 0,
        pxi = 0;

    incx = incx << 3;

    for (i = 0, pxi = x; (i | 0) < (n | 0); i = i + 1 | 0, pxi = pxi + incx | 0) {
      darray[pxi >> 3] = a * darray[pxi >> 3];
    }
  }

  function dsyr(uplo, n, alpha, x, incx, a, lda) {
    uplo = uplo | 0;
    n = n | 0;
    alpha = +alpha;
    x = x | 0;
    incx = incx | 0;
    a = a | 0;
    lda = lda | 0;

    var i = 0,
        j = 0,
        paij = 0,
        pxi = 0,
        pxj = 0;

    if (uplo) {
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        pxi = x + ((imul(i, incx) | 0) << 3) | 0;
        for (j = i; (j | 0) < (n | 0); j = j + 1 | 0) {
          paij = a + ((imul(j, lda) | 0) + i << 3) | 0;
          pxj = x + ((imul(j, incx) | 0) << 3) | 0;
          darray[paij >> 3] = +darray[paij >> 3] + alpha * +darray[pxi >> 3] * +darray[pxj >> 3];
        }
      }
    } else {
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        pxi = x + ((imul(i, incx) | 0) << 3) | 0;
        for (j = i; (j | 0) < (n | 0); j = j + 1 | 0) {
          paij = a + ((imul(i, lda) | 0) + j << 3) | 0;
          pxj = x + ((imul(j, incx) | 0) << 3) | 0;
          darray[paij >> 3] = +darray[paij >> 3] + alpha * +darray[pxi >> 3] * +darray[pxj >> 3];
        }
      }
    }
  }

  function dsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb) {
    uplo = uplo | 0;
    n = n | 0;
    nrhs = nrhs | 0;
    a = a | 0;
    lda = lda | 0;
    ipiv = ipiv | 0;
    b = b | 0;
    ldb = ldb | 0;

    dsytrf(uplo, n, a, lda, ipiv);
    dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
  }

  function dsytrf(uplo, n, a, lda, ipiv) {
    uplo = uplo | 0;
    n = n | 0;
    a = a | 0;
    lda = lda | 0;
    ipiv = ipiv | 0;

    var k = 0,
        kp = 0,
        imax = 0,
        jmax = 0,
        pakk = 0,
        pakpkp = 0,
        absakk = 0.0,
        colmax = 0.0,
        rowmax = 0.0,
        t = 0.0,
        r1 = 0.0;

    if (uplo) {
      for (k = 0; (k | 0) < (n | 0); k = k + 1 | 0) {
        pakk = a + ((imul(k, lda) | 0) + k << 3) | 0;
        absakk = +abs(darray[pakk >> 3]);
        if ((k | 0) < (n - 1 | 0)) {
          imax = k + 1 + (idamax(n - k - 1 | 0, a + ((imul(k + 1 | 0, lda) | 0) + k << 3) | 0, lda) | 0) | 0;
          colmax = +abs(darray[a + ((imul(imax, lda) | 0) + k << 3) >> 3]);
        } else {
          colmax = 0.0;
        }

        if (absakk >= colmax) {
          kp = k;
        } else {
          jmax = k + (idamax(imax - k | 0, a + ((imul(imax, lda) | 0) + k << 3) | 0, 1) | 0) | 0;
          rowmax = +abs(darray[a + ((imul(imax, lda) | 0) + jmax << 3) >> 3]);
          if (absakk >= colmax * colmax / rowmax) {
            kp = k;
          } else if(+abs(darray[a + ((imul(imax, lda) | 0) + imax << 3) >> 3]) < rowmax) {
            kp = k;
          } else {
            kp = imax;
          }
        }

        if ((kp | 0) != (k | 0)) {
          pakpkp = a + ((imul(kp, lda) | 0) + kp << 3) | 0;
          if ((kp | 0) < (n - 1 | 0)) {
            dswap(n - kp - 1 | 0,
                a + ((imul(kp + 1 | 0, lda) | 0) + k << 3) | 0, lda,
                a + ((imul(kp + 1 | 0, lda) | 0) + kp << 3) | 0, lda);
          }
          dswap(kp - k - 1 | 0,
              a + ((imul(k + 1 | 0, lda) | 0) + k << 3) | 0, lda,
              a + ((imul(kp, lda) | 0) + k + 1 << 3) | 0, 1);

          t = +darray[pakk >> 3];
          darray[pakk >> 3] = darray[pakpkp >> 3];
          darray[pakpkp >> 3] = t;
        }

        if ((k | 0) < (n - 1 | 0)) {
          r1 = 1.0 / +darray[pakk >> 3];
          dsyr(uplo, n - k - 1 | 0, -r1,
              a + ((imul(k + 1 | 0, lda) | 0) + k << 3) | 0, lda,
              a + ((imul(k + 1 | 0, lda) | 0) + k + 1 << 3) | 0, lda);
          dscal(n - k - 1 | 0, r1, a + ((imul(k + 1 | 0, lda) | 0) + k << 3) | 0, lda);
        }
        uiarray[ipiv + (k << 2) >> 2] = kp;
      }
    } else {
      for (k = n - 1 | 0; (k | 0) >= 0; k = k - 1 | 0) {
        pakk = a + ((imul(k, lda) | 0) + k << 3) | 0;
        absakk = +abs(darray[pakk >> 3]);
        if ((k | 0) > 1) {
          imax = idamax(k, a + (k << 3) | 0, lda) | 0;
          colmax = +abs(darray[a + ((imul(imax, lda) | 0) + k << 3) >> 3]);
        } else {
          colmax = 0.0;
        }

        if (absakk > colmax) {
          kp = k;
        } else {
          jmax = imax + (idamax(k - imax + 1 | 0, a + ((imul(imax, lda) | 0) + imax + 1 << 3) | 0, lda) | 0) | 0;
          rowmax = +abs(darray[a + ((imul(imax, lda) | 0) + jmax << 3) >> 3]);
          if (absakk >= colmax * colmax / rowmax) {
            kp = k;
          } else {
            kp = imax;
          }
        }

        if ((kp | 0) != (k | 0)) {
          pakpkp = a + ((imul(kp, lda) | 0) + kp << 3) | 0;
          dswap(kp, a + (k << 3) | 0, lda, a + (kp << 3) | 0, lda);
          dswap(k - kp | 0,
              a + ((imul(kp + 1 | 0, lda) | 0) + k << 3) | 0, lda,
              a + ((imul(kp, lda) | 0) + kp + 1 << 3) | 0, 1);

          t = +darray[pakk >> 3];
          darray[pakk >> 3] = darray[pakpkp >> 3];
          darray[pakpkp >> 3] = t;
        }

        r1 = 1.0 / +darray[pakk >> 3];
        dsyr(uplo, k, -r1, a + (k << 3) | 0, lda, a, lda);
        dscal(k, r1, a + (k << 3) | 0, lda);
        uiarray[ipiv + (k << 2) >> 2] = kp;
      }
    }
  }

  function dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb) {
    uplo = uplo | 0;
    n = n | 0;
    nrhs = nrhs | 0;
    a = a | 0;
    lda = lda | 0;
    ipiv = ipiv | 0;
    b = b | 0;
    ldb = ldb | 0;

    var i = 0,
        j = 0,
        paii = 0,
        pai = 0,
        paj = 0,
        pbi = 0,
        pbj = 0,
        pipivi = 0;

    // z := L^-1 P^t b
    for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
      pipivi = ipiv + (i << 2) | 0;
      j = uiarray[pipivi >> 2] | 0;
      if ((i | 0) != (j | 0)) {
        pai = a + ((imul(i, lda) | 0) << 3) | 0;
        paj = a + ((imul(j, lda) | 0) << 3) | 0;
        dswap(i, pai, 1, paj, 1);
        pbi = b + ((imul(i, ldb) | 0) << 3) | 0;
        pbj = b + ((imul(j, ldb) | 0) << 3) | 0;
        dswap(nrhs, pbi, 1, pbj, 1);
      }
    }
    dtrsm(0, uplo, 0, 1, n, nrhs, 1.0, a, lda, b, ldb);

    // y := L^t^-1 D^-1 y
    for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
      paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
      dscal(nrhs, 1.0 / +darray[paii >> 3], b + ((imul(i, ldb) | 0) << 3) | 0, 1);
    }
    dtrsm(0, uplo, 1, 1, n, nrhs, 1.0, a, lda, b, ldb);

    // x := P y
    for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
      pipivi = ipiv + (i << 2) | 0;
      j = uiarray[pipivi >> 2] | 0;
      if ((i | 0) != (j | 0)) {
        pbi = b + ((imul(i, ldb) | 0) << 3) | 0;
        pbj = b + ((imul(j, ldb) | 0) << 3) | 0;
        dswap(nrhs, pbi, 1, pbj, 1);
      }
    }
  }

  function dtrmv(uplo, trans, diag, n, a, lda, x, incx) {
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
        paij = 0,
        pxi = 0,
        pxj = 0,
        val = 0.0;

    if (uplo) {
      // lower triangular matrix
      for (i = n - 1 | 0; (i | 0) >= 0; i = i - 1 | 0) {
        pxi = x + ((imul(i, incx) | 0) << 3) | 0;
        paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
        val = darray[pxi >> 3] * darray[paii >> 3];
        for (j = 0; (j | 0) < (i | 0); j = j + 1 | 0) {
          pxj = x + ((imul(j, incx) | 0) << 3) | 0;
          paij = a + ((imul(i, lda) | 0) + j << 3) | 0;
          val = val + darray[pxj >> 3] * darray[paij >> 3];
        }
        darray[pxi >> 3] = val;
      }
    } else {
      // upper triangular matrix
      for (i = 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        pxi = x + ((imul(i, incx) | 0) << 3) | 0;
        paii = a + ((imul(i, lda) | 0) + i << 3) | 0;
        val = darray[pxi >> 3] * darray[paii >> 3];
        for (j = i + 1 | 0; (j | 0) < (n | 0); j = j + 1 | 0) {
          pxj = x + ((imul(j, incx) | 0) << 3) | 0;
          paij = a + ((imul(i, lda) | 0) + j << 3) | 0;
          val = val + darray[pxj >> 3] * darray[paij >> 3];
        }
        darray[pxi >> 3] = val;
      }
    }
  }

  function dtrtri(uplo, diag, n, a, lda) {
    // TODO support diag, lda
    uplo = uplo | 0;
    diag = diag | 0;
    n = n | 0;
    a = a | 0;
    lda = lda | 0;

    var i = 0,
        paii = 0,
        paij = 0,
        pajj = 0,
        aii = 0.0;

    if (uplo) {
      // lower triangular matrix
      paii = a + ((imul(n, n) | 0) - 1 << 3) | 0;
      darray[paii >> 3] = 1.0 / +darray[paii >> 3];
      for (i = n - 2 | 0; (i | 0) >= 0; i = i - 1 | 0) {
        paii = a + ((imul(i, n) | 0) + i << 3) | 0;
        paij = a + ((imul(i + 1 | 0, n) | 0) + i << 3) | 0;
        pajj = a + ((imul(i + 1 | 0, n) | 0) + i + 1 << 3) | 0;
        aii = darray[paii >> 3] = 1.0 / +darray[paii >> 3];
        dtrmv(1, 0, 0, n - i | 0, pajj, lda, paij, n);
        dscal(n - i | 0, -aii, paij, n);
      }
    } else {
      // upper triangular matrix
      darray[a >> 3] = 1.0 / +darray[a >> 3];
      for (i = 1 | 0; (i | 0) < (n | 0); i = i + 1 | 0) {
        paii = a + ((imul(i, n) | 0) + i << 3) | 0;
        paij = a + (i << 3) | 0;
        aii = darray[paii >> 3] = 1.0 / +darray[paii >> 3];
        dtrmv(0, 0, 0, i, a, lda, paij, n);
        dscal(i, -aii, paij, n);
      }
    }
  }

  return {
    dasum: dasum,
    daxpy: daxpy,
    dcopy: dcopy,
    ddot: ddot,
    dgemm: dgemm,
    dgemv: dgemv,
    dgesv: dgesv,
    dgetri: dgetri,
    dgetrf: dgetrf,
    dposv: dposv,
    dpotrf: dpotrf,
    dpotrs: dpotrs,
    dscal: dscal,
    dswap: dswap,
    dsyr: dsyr,
    dsysv: dsysv,
    dsytrf: dsytrf,
    dsytrs: dsytrs,
    dtrmv: dtrmv,
    dtrsv: dtrsv,
    dtrsm: dtrsm,
    dtrtri: dtrtri,
    idamax: idamax
  };
}

module.exports = LinalgModule;
