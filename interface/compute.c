/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 65536.0
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMM "
#else
#define ERROR_NAME "SGEMM "
#endif
#else
#define SMP_THRESHOLD_MIN 8192.0
#ifndef GEMM3M
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM "
#else
#define ERROR_NAME "CGEMM "
#endif
#else
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM3M "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM3M "
#else
#define ERROR_NAME "CGEMM3M "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

//interface for A & B gemm_compute
//N or T or Packed
static int (*gemmcompute[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_AN_BN, GEMM_AT_BN, GEMM_AP_BN, NULL,
  GEMM_AN_BT, GEMM_AT_BT, GEMM_AP_BT, NULL,
  GEMM_AN_BP, GEMM_AT_BP, GEMM_AP_BP, NULL,
  NULL, NULL, NULL, NULL,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_AN_BN, GEMM_THREAD_AT_BN, GEMM_THREAD_AP_BN, NULL,
  GEMM_THREAD_AN_BT, GEMM_THREAD_AT_BT, GEMM_THREAD_AP_BT, NULL,
  GEMM_THREAD_AN_BP, GEMM_THREAD_AT_BP, GEMM_THREAD_AP_BP, NULL,
  NULL, NULL, NULL, NULL,
#endif
};

#ifdef CBLAS

void CNAME(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB,
	   blasint m, blasint n, blasint k,
	   FLOAT *a, blasint lda, FLOAT *b, blasint ldb, FLOAT beta,
	   FLOAT *c, blasint ldc) {

  blas_arg_t args;
  int transa, transb;
  blasint nrowa, nrowb, info;

  XFLOAT *buffer;
  XFLOAT *sa, *sb;

#ifdef SMP
  double MNK;
#ifndef COMPLEX
#ifdef XDOUBLE
  int mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  int mode  =  BLAS_DOUBLE  | BLAS_REAL;
#else
  int mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#else
#ifdef XDOUBLE
  int mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  int mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#else
  int mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif

#if defined(SMP) && !defined(NO_AFFINITY) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  int nodes;
#endif

  PRINT_DEBUG_CNAME;

#if !defined(COMPLEX) && !defined(DOUBLE) && defined(USE_SGEMM_KERNEL_DIRECT)
#ifdef DYNAMIC_ARCH
 if (support_avx512() )
#endif  
  if (beta == 0 && alpha == 1.0 && order == CblasRowMajor && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
	SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
	return;
  }

#endif

#ifndef COMPLEX
  args.beta  = (void *)&beta;
#else
  args.beta  = (void *)beta;
#endif
  transa = -1;
  transb = -1;
  info   =  0;

  if (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;

    args.a = (void *)a;
    args.b = (void *)b;
    args.c = (void *)c;

    args.lda = lda;
    args.ldb = ldb;
    args.ldc = ldc;

    if (TransA == CblasNoTrans)     transa = 0;
    if (TransA == CblasTrans)       transa = 1;
    if (TransA == CblasPacked)      transa = 2;
    
    if (TransB == CblasNoTrans)     transb = 0;
    if (TransB == CblasTrans)       transb = 1;
    if (TransB == CblasPacked)      transb = 2;

    nrowa = args.m;
    if (transa & 1) nrowa = args.k;
    nrowb = args.k;
    if (transb & 1) nrowb = args.n;

    info = -1;

    if (args.ldc < args.m) info = 13;
    if (args.ldb < nrowb)  info = 10;
    if (args.lda < nrowa)  info =  8;
    if (args.k < 0)        info =  5;
    if (args.n < 0)        info =  4;
    if (args.m < 0)        info =  3;
    if (transb < 0)        info =  2;
    if (transa < 0)        info =  1;
  }

  if (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    args.a = (void *)b;
    args.b = (void *)a;
    args.c = (void *)c;

    args.lda = ldb;
    args.ldb = lda;
    args.ldc = ldc;

    if (TransB == CblasNoTrans)     transa = 0;
    if (TransB == CblasTrans)       transa = 1;
    if (TransB == CblasPacked)      transb = 2;
    
    if (TransA == CblasNoTrans)     transb = 0;
    if (TransA == CblasTrans)       transb = 1;
    if (TransA == CblasPacked)      transa = 2;

    nrowa = args.m;
    if (transa & 1) nrowa = args.k;
    nrowb = args.k;
    if (transb & 1) nrowb = args.n;

    info = -1;

    if (args.ldc < args.m) info = 13;
    if (args.ldb < nrowb)  info = 10;
    if (args.lda < nrowa)  info =  8;
    if (args.k < 0)        info =  5;
    if (args.n < 0)        info =  4;
    if (args.m < 0)        info =  3;
    if (transb < 0)        info =  2;
    if (transa < 0)        info =  1;

  }

  if (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

#endif

  if ((args.m == 0) || (args.n == 0)) return;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  
  buffer = (XFLOAT *)blas_memory_alloc(0);

  sa = (XFLOAT *)((BLASLONG)buffer +GEMM_OFFSET_A);
  sb = (XFLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
  mode |= (transa << BLAS_TRANSA_SHIFT);
  mode |= (transb << BLAS_TRANSB_SHIFT);

  MNK = (double) args.m * (double) args.n * (double) args.k;
  if ( MNK <= (SMP_THRESHOLD_MIN  * (double) GEMM_MULTITHREAD_THRESHOLD)  )
	args.nthreads = 1;
  else
	args.nthreads = num_cpu_avail(3);
  args.common = NULL;

 if (args.nthreads == 1) {
#endif

    (gemmcompute[transA | transB << 2])(&args, NULL, NULL, sa, sb, 0);

#ifdef SMP

  } else {

#ifndef USE_SIMPLE_THREADED_LEVEL3

#ifndef NO_AFFINITY
      nodes = get_num_nodes();

      if ((nodes > 1) && get_node_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, NULL, NULL, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } else {
#endif

	(gemmcompute[16 | transA | transB << 2])(&args, NULL, NULL, sa, sb, 0);

#else

	GEMM_THREAD(mode, &args, NULL, NULL, gemmcompute[16 | transA | transB << 2], sa, sb, args.nthreads);

#endif

#ifndef USE_SIMPLE_THREADED_LEVEL3
#ifndef NO_AFFINITY
      }
#endif
#endif

#endif

#ifdef SMP
  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n, 2 * args.m * args.n * args.k);

  IDEBUG_END;

  return;
}
