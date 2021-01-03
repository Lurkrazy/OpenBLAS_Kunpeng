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

//pack , so direct is not supported
#undef USE_SGEMM_KERNEL_DIRECT 

#undef USE_SIMPLE_THREADED_LEVEL3

//todo: GEMM_NN -> GEMM_NN_PACK ...
// add defines in common_macro.h for GEMM_NN_PACK ...
// only for N & T
static int (*gemmpack[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_AN_PACK, GEMM_AT_PACK,
  GEMM_BN_PACK, GEMM_BT_PACK,
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  GEMM_THREAD_AN_PACK, GEMM_THREAD_AT_PACK,
  GEMM_THREAD_BN_PACK, GEMM_THREAD_BT_PACK,
#endif
};

#ifdef CBLAS
void CNAME(enum CBLAS_ORDER order, const CBLAS_IDENTIFIER identifier, enum CBLAS_TRANSPOSE Trans,
	   blasint m, blasint n, blasint k, FLOAT alpha, FLOAT *src, blasint ld, FLOAT *dest)
{
  blas_arg_t args;
  int trans, meta;
  blasint nrow, info;

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

//#if !defined(COMPLEX) && !defined(DOUBLE) && defined(USE_SGEMM_KERNEL_DIRECT)
////#ifdef DYNAMIC_ARCH
//// if (support_avx512() )
////#endif  
//  //if (beta == 0 && alpha == 1.0 && order == CblasRowMajor && TransA == CblasNoTrans && TransB == CblasNoTrans && SGEMM_DIRECT_PERFORMANT(m,n,k)) {
//  //  SGEMM_DIRECT(m, n, k, a, lda, b, ldb, c, ldc);
//  //  return;
//  //}
//
//#endif

//#ifndef COMPLEX
  args.alpha = (void *)&alpha;
//#else
//  args.alpha = (void *)alpha;
//#endif

  trans = -1;
  info   =  0;

  if (order == CblasColMajor) {
    args.m = m;
    args.n = n;
    args.k = k;


    args.lda = ld;

    if (Trans == CblasNoTrans)     trans = 0;
    if (Trans == CblasTrans)       trans = 1;
    if (identifier == CblasAMatrix) meta = 0;
    else meta = 1;

    // always use args.a 
    if(meta == 0) {
        nrow = args.m;
        if (trans & 1) nrow = args.k;
        args.a = (void *)src;
    } else {
        nrow = args.k;
        if (trans & 1) nrow = args.n;
        args.a = (void *)src;
    }
    info = -1;

    if (args.lda < nrow)  info =  8;
    if (args.k < 0)        info =  5;
    if (args.n < 0)        info =  4;
    if (args.m < 0)        info =  3;
    if (trans < 0)        info =  1;
  }

  if (order == CblasRowMajor) {
    args.m = n;
    args.n = m;
    args.k = k;

    //???potential error
    args.lda = ld;

    if (Trans == CblasNoTrans)     trans = 0;
    if (Trans == CblasTrans)       trans = 1;
    if (identifier == CblasAMatrix) meta = 1;
    else meta = 0;

    if(meta == 0) {
        nrow = args.m;
        if (trans & 1) nrow = args.k;
        args.a = (void *)src;
    } else {
        nrow = args.k;
        if (trans & 1) nrow = args.n;
	args.a = (void *)src;
    }

    info = -1;

    if (args.lda < nrow)  info =  8;
    if (args.k < 0)        info =  5;
    if (args.n < 0)        info =  4;
    if (args.m < 0)        info =  3;
    if (trans < 0)        info =  1;

  }

  if (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    return;
  }

#endif

  if ((args.m == 0) || (args.n == 0)) return;


  IDEBUG_START;

  FUNCTION_PROFILE_START();

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

    (gemmpack[(meta<<1) | trans])(&args, NULL, NULL, dest, NULL, 0);

#ifdef SMP

  } else {

#ifndef USE_SIMPLE_THREADED_LEVEL3
/*
#ifndef NO_AFFINITY
      nodes = get_num_nodes();

      if ((nodes > 1) && get_node_equal()) {

	args.nthreads /= nodes;

	gemm_thread_mn(mode, &args, NULL, NULL, gemm[16 | (transb << 2) | transa], sa, sb, nodes);

      } else {
#endif
*/
	(gemmpack[4 | (meta<<1) | trans])(&args, NULL, NULL, sa, sb, 0);

#else
//not consider now
	GEMM_THREAD(mode, &args, NULL, NULL, gemm[(transb << 2) | transa], sa, sb, args.nthreads);

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


  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.m * args.k + args.k * args.n + args.m * args.n, 2 * args.m * args.n * args.k);

  IDEBUG_END;

  return;
}
