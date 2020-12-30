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

#ifndef ICOPY_OPERATION
#define ICOPY_OPERATION
#ifdef AN
#define ICOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_ITCOPY_PACK(M, N, (IFLOAT *)(A) + ((Y) + (X) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#endif
#ifdef AT
#define ICOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_INCOPY_PACK(M, N, (IFLOAT *)(A) + ((X) + (Y) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#endif
#endif

#ifndef OCOPY_OPERATION
#define OCOPY_OPERATION
#ifdef BN
#define OCOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_ONCOPY_PACK(M, N, (IFLOAT *)(A) + ((X) + (Y) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#endif
#ifdef BT
#define OCOPY_OPERATION(M, N, A, LDA, X, Y, BUFFER, ALPHA) GEMM_OTCOPY_PACK(M, N, (IFLOAT *)(A) + ((Y) + (X) * (LDA)) * COMPSIZE, LDA, BUFFER, ALPHA);
#endif
#endif

#ifndef A
#define A	args -> a
#endif
#ifndef LDA
#define LDA	args -> lda
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

#define DEST_BASE 1048576
#define MAX_THREADS 128
/*
 *  sa is the base addr of dest.
 *  the length of sa must > m*k*sizeof)FLOAT) + 1048576.
 *
 *  addr map:
 *  addr of thread0 base, thread1 base .....
 *  (thread0 base)addr of thread0 block0, thread0 block1 ....
 *  (thread1 base)addr of thread1 block0, ....
 *  ...
 *  thread0 block0
 *  thread0 block1
 *  ...
 * */

int CNAME(blas_arg_t *args, BLASLONG *range_m, BLASLONG *range_n,
		  XFLOAT *sa, XFLOAT *sb, BLASLONG dummy){
  BLASLONG k, lda;
  FLOAT *alpha;
  IFLOAT *a;
  IFLOAT *dest = (IFLOAT*)((unsigned long) sa + DEST_BASE);
  unsigned long *block = (unsigned long)((unsigned long)sa + MAX_THREADS*8);
  BLASLONG m_from, m_to, n_from, n_to;

  BLASLONG ls, is, js;
  BLASLONG min_l, min_i, min_j;

  BLASLONG l1stride, gemm_p, l2size;

  //threads = 1
  *(unsigned long*)((unsigned long)sa) = ((unsigned long)block);

  k = K;

  a = (IFLOAT *)A;

  lda = LDA;

  alpha = (FLOAT *)args -> alpha;

  m_from = 0;
  m_to   = M;

  n_from = 0;
  n_to   = N;


  l2size = GEMM_P * GEMM_Q;


#if defined(BP)
  for(js = n_from; js < n_to; js += GEMM_R){
    min_j = n_to - js;
    if (min_j > GEMM_R) min_j = GEMM_R;
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

#if defined (AP)
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
      
      ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_from, dest);
      *block = dest; dest += min_l * min_i; block++;
#endif
#if defined(BP)
//here is different from level3.c . pack for only once
	OCOPY_OPERATION(min_l, min_j, b, ldb, ls, js, dest);
      *block = dest; dest += min_l * min_i; block++;

#endif
#if defined(AP)
      for(is = m_from + min_i; is < m_to; is += min_i){
	min_i = m_to - is;

	if (min_i >= GEMM_P * 2) {
	  min_i = GEMM_P;
	} else
	  if (min_i > GEMM_P) {
	    min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
	  }
	ICOPY_OPERATION(min_l, min_i, a, lda, ls, is, dest);
	*block = dest; dest += min_l * min_i; block++;
      } /* end of is */
#endif
    } /* end of js */
#if defined(BP)
  } /* end of ls */
#endif

  return 0;
}
