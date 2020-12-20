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

/* This file is a template for level 3 operation */

#ifndef BETA_OPERATION
#define BETA_OPERATION(M_FROM, M_TO, N_FROM, N_TO, BETA, C, LDC) \
	GEMM_BETA((M_TO) - (M_FROM), (N_TO - N_FROM), 0, \
		  BETA[0], NULL, 0, NULL, 0, \
		  (FLOAT *)(C) + ((M_FROM) + (N_FROM) * (LDC)) * COMPSIZE, LDC)
#endif

#ifndef ICOPY_OPERATION
#if defined(NN) || defined(NT) || defined(NC) || defined(NR) || \
    defined(RN) || defined(RT) || defined(RC) || defined(RR)
#define ICOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_ITCOPY_PACK(M, N, (IFLOAT *)(A) + ((Y) + (X) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#else
#define ICOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_INCOPY_PACK(M, N, (IFLOAT *)(A) + ((X) + (Y) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#endif
#endif

#ifndef OCOPY_OPERATION
#if defined(NN) || defined(TN) || defined(CN) || defined(RN) || \
    defined(NR) || defined(TR) || defined(CR) || defined(RR)
#define OCOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_ONCOPY_PACK(M, N, (IFLOAT *)(A) + ((X) + (Y) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#else
#define OCOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_OTCOPY_PACK(M, N, (IFLOAT *)(A) + ((Y) + (X) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#endif
#endif

#ifndef KERNEL_FUNC
#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
#define KERNEL_FUNC	GEMM_KERNEL_N
#endif
#if defined(CN) || defined(CT) || defined(RN) || defined(RT)
#define KERNEL_FUNC	GEMM_KERNEL_L
#endif
#if defined(NC) || defined(TC) || defined(NR) || defined(TR)
#define KERNEL_FUNC	GEMM_KERNEL_R
#endif
#if defined(CC) || defined(CR) || defined(RC) || defined(RR)
#define KERNEL_FUNC	GEMM_KERNEL_B
#endif
#endif

#ifndef KERNEL_OPERATION
#if !defined(XDOUBLE) || !defined(QUAD_PRECISION)
#define KERNEL_OPERATION(M, N, K, ALPHA, SA, SB, C, LDC, X, Y) \
	KERNEL_FUNC(M, N, K, ALPHA[0], SA, SB, (FLOAT *)(C) + ((X) + (Y) * LDC) * COMPSIZE, LDC)
#else
#define KERNEL_OPERATION(M, N, K, ALPHA, SA, SB, C, LDC, X, Y) \
	KERNEL_FUNC(M, N, K, ALPHA, SA, SB, (FLOAT *)(C) + ((X) + (Y) * LDC) * COMPSIZE, LDC)
#endif
#endif

#ifndef A
#define A	args -> a
#endif
#ifndef LDA
#define LDA	args -> lda
#endif
#ifndef B
#define B	args -> b
#endif
#ifndef LDB
#define LDB	args -> ldb
#endif
#ifndef C
#define C	args -> c
#endif
#ifndef LDC
#define LDC	args -> ldc
#endif
#ifndef M
#define M	args -> m
#endif
#ifndef N
#define N	args -> n
#endif
#ifndef K
#define K	args -> k
#endif

int CNAME(blas_arg_t *args, BLASLONG *range_m, BLASLONG *range_n,
		  //XFLOAT *sa, XFLOAT *sb, blasint transa, blasint transb, BLASLONG dummy){
		  void *sa, void *sb, blasint transa, blasint transb, BLASLONG dummy){
  BLASLONG k, lda, ldb, ldc;
  FLOAT *beta;
  IFLOAT *a, *b;
  FLOAT *c;
  BLASLONG m_from, m_to, n_from, n_to;

  BLASLONG ls, is, js;
  BLASLONG min_l, min_i, min_j;

  BLASLONG l1stride, gemm_p, l2size;

  k = K;

  a = (IFLOAT *)A;
  b = (IFLOAT *)B;
  c = (FLOAT *)C;

  lda = LDA;
  ldb = LDB;
  ldc = LDC;
#define IFPACKED(TRAN) (TRAN == 2)

#if IFPACKED(transa)
#define AP
#endif

#if IFPACKED(transb)
#define BP
#endif

  alpha = (FLOAT *)args -> alpha;
  beta  = (FLOAT *)args -> beta;

  m_from = 0;
  m_to   = M;

  if (range_m) {
    m_from = *(((BLASLONG *)range_m) + 0);
    m_to   = *(((BLASLONG *)range_m) + 1);
  }

  n_from = 0;
  n_to   = N;

  if (range_n) {
    n_from = *(((BLASLONG *)range_n) + 0);
    n_to   = *(((BLASLONG *)range_n) + 1);
  }

  if (beta) {
    if (beta[0] != ONE) {
	  BETA_OPERATION(m_from, m_to, n_from, n_to, beta, c, ldc);
	}
  }


  l2size = GEMM_P * GEMM_Q;

#if defined(AP) || defined(BP)
  unsigned long desta, destb;
  //init desta and destb
  desta = *sa;
  destb = *sb;
#endif


  for(js = n_from; js < n_to; js += GEMM_R){
    min_j = n_to - js;
    if (min_j > GEMM_R) min_j = GEMM_R;
    //reset desta
    desta = *sa;

    for(ls = 0; ls < k; ls += min_l){

      min_l = k - ls;

      if (min_l >= GEMM_Q * 2) {
	// gemm_p = GEMM_P;
	min_l  = GEMM_Q;
      } else {
	if (min_l > GEMM_Q) {
	  min_l = ((min_l / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
	}
	gemm_p = ((l2size / min_l + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
	while (gemm_p * min_l > l2size) gemm_p -= GEMM_UNROLL_M;
      }

      /* First, we have to move data A to L2 cache */
      min_i = m_to - m_from;
      l1stride = 1;

      if (min_i >= GEMM_P * 2) {
	min_i = GEMM_P;
      } else {
	if (min_i > GEMM_P) {
	  min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
	} else {
	  l1stride = 0;
	}
      }
#if defined(AP)
      //a packed
      desta = *(sa++);
#else
      //ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_from, sa);
      ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_from, desta);
#endif
      
#if defined(BP)
      destb = *(sb++);
#else
      OCOPY_OPERATION(min_l, min_j, b, ldb, ls, js, destb);
#endif
      
#if defined(AP)
      KERNEL_OPERATION(min_i, min_j, min_l, 1, desta, destb, c, ldc, m_from, js);
#else 
      KERNEL_OPERATION(min_i, min_j, min_l, alpha, desta, destb, c, ldc, m_from, js);
#endif


//      for(jjs = js; jjs < js + min_j; jjs += min_jj){
//	min_jj = min_j + js - jjs;
//#if defined(SKYLAKEX) || defined(COOPERLAKE)
//	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve best performance */
//	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
//#else
//        if (min_jj >= 3*GEMM_UNROLL_N) min_jj = 3*GEMM_UNROLL_N;
//        else
///*
//		if (min_jj >= 2*GEMM_UNROLL_N) min_jj = 2*GEMM_UNROLL_N;
//        	else
//*/
//          		if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
//#endif
//
//	OCOPY_OPERATION(min_l, min_jj, b, ldb, ls, jjs,
//			sb + min_l * (jjs - js) * COMPSIZE * l1stride);
//
//	KERNEL_OPERATION(min_i, min_jj, min_l, alpha,
//			 sa, sb + min_l * (jjs - js)  * COMPSIZE * l1stride, c, ldc, m_from, jjs);
//
//      }

      for(is = m_from + min_i; is < m_to; is += min_i){
	min_i = m_to - is;

	if (min_i >= GEMM_P * 2) {
	  min_i = GEMM_P;
	} else
	  if (min_i > GEMM_P) {
	    min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
	  }
#if defined(AP)
    desta = *(a++);
#else
	ICOPY_OPERATION(min_l, min_i, a, lda, ls, is, desta);
#endif

#if defined(AP)
	KERNEL_OPERATION(min_i, min_j, min_l, 1, desta, destb, c, ldc, is, js);
#else
	KERNEL_OPERATION(min_i, min_j, min_l, alpha, desta, destb, c, ldc, is, js);
#endif

      } /* end of is */
    } /* end of js */
  } /* end of ls */

  return 0;
}
