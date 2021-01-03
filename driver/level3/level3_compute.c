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
#if defined(AN)
#define ICOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER) GEMM_ITCOPY(M, N, (IFLOAT *)(A) + ((Y) + (X) * (LDA)) * COMPSIZE, LDA, BUFFER);
#elif defined(AT)
#define ICOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER) GEMM_INCOPY(M, N, (IFLOAT *)(A) + ((X) + (Y) * (LDA)) * COMPSIZE, LDA, BUFFER);
#endif
#endif

#ifndef OCOPY_OPERATION
#if defined(BN)
#define OCOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER) GEMM_ONCOPY(M, N, (IFLOAT *)(A) + ((X) + (Y) * (LDA)) * COMPSIZE, LDA, BUFFER);
#elif defined(BT)
#define OCOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER) GEMM_OTCOPY(M, N, (IFLOAT *)(A) + ((Y) + (X) * (LDA)) * COMPSIZE, LDA, BUFFER);
#endif
#endif

#ifndef KERNEL_FUNC
#define KERNEL_FUNC	GEMM_KERNEL_N
#endif

//#ifndef KERNEL_FUNC 
//#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
//#define KERNEL_FUNC	GEMM_KERNEL_N
//#endif
//#if defined(CN) || defined(CT) || defined(RN) || defined(RT)
//#define KERNEL_FUNC	GEMM_KERNEL_L
//#endif
//#if defined(NC) || defined(TC) || defined(NR) || defined(TR)
//#define KERNEL_FUNC	GEMM_KERNEL_R
//#endif
//#if defined(CC) || defined(CR) || defined(RC) || defined(RR)
//#define KERNEL_FUNC	GEMM_KERNEL_B
//#endif
//#endif

#ifndef KERNEL_OPERATION
#define KERNEL_OPERATION(M, N, K, ALPHA, SA, SB, C, LDC, X, Y) \
	KERNEL_FUNC(M, N, K, ALPHA[0], SA, SB, (FLOAT *)(C) + ((X) + (Y) * LDC) * COMPSIZE, LDC)
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

int CNAME(blas_arg_t *args, BLASLONG *range_m, BLASLONG *range_n, FLOAT *sa, FLOAT *sb, BLASLONG dummy){
  BLASLONG k, lda, ldb, ldc;
  FLOAT *alpha, *beta;
  IFLOAT *a, *b, *DESTA, *DESTB;
  FLOAT *c;
  BLASLONG m_from, m_to, n_from, n_to;

  BLASLONG ls, is, js;
  BLASLONG min_l, min_i, min_j;

  BLASLONG l1stride, gemm_p, l2size;

  k = K;

  DESTA = (IFLOAT *)A;
  DESTB = (IFLOAT *)B;
  c = (FLOAT *)C;

  lda = LDA;
  ldb = LDB;
  ldc = LDC;

  //alphaT just for initialize alpha = 1
  FLOAT alphaT[1] = {1};
  alpha = alphaT;
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

//  unsigned long long desta;
//  unsigned long long destb;
//
//#if defined(AP)
//  unsigned long long* tmpa;
//#else
////  FLOAT* tmpa;
//#endif

  for(js = n_from; js < n_to; js += GEMM_R){
    min_j = n_to - js;
    if (min_j > GEMM_R) min_j = GEMM_R;

#if defined(AP)
  //reset the pointer to packed A
  unsigned long long* pointerToEveryBlockofPackedA = (unsigned long long *)DESTA;
#endif

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
      FLOAT* desta = (FLOAT *)(*pointerToEveryBlockofPackedA);
      pointerToEveryBlockofPackedA++;
#else
      ICOPY_OPERATION(min_l, min_i, DESTA, lda, ls, m_from, sa);
      FLOAT* desta = sa;
#endif
      
#if defined(BP)
      unsigned long long* pointerToEveryBlockofPackedB = DESTB;
      FLOAT* destb = (FLOAT *)(*pointerToEveryBlockofPackedB);
      pointerToEveryBlockofPackedB++;
#else
      OCOPY_OPERATION(min_l, min_j, DESTB, ldb, ls, js, sb);
      FLOAT* destb = sb;
#endif
      
      KERNEL_OPERATION(min_i, min_j, min_l, alpha, desta, destb, c, ldc, m_from, js);

      for(is = m_from + min_i; is < m_to; is += min_i){
	min_i = m_to - is;

	if (min_i >= GEMM_P * 2) {
	  min_i = GEMM_P;
	} else
	  if (min_i > GEMM_P) {
	    min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
	  }
#if defined(AP)
      //a packed
      FLOAT* desta = (FLOAT *)(*pointerToEveryBlockofPackedA);
      pointerToEveryBlockofPackedA++;
#else
	    ICOPY_OPERATION(min_l, min_i, DESTA, lda, ls, is, sa);
      FLOAT* desta = sa;
#endif

	KERNEL_OPERATION(min_i, min_j, min_l, alpha, desta, destb, c, ldc, is, js);

      } /* end of is */
    } /* end of js */
  } /* end of ls */

  return 0;
}
