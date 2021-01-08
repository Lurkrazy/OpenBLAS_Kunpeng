/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "bench.h"
#include "cblas.h"
#include <float.h>
#undef GEMM

#ifdef DOUBLE
#define GEMM   BLASFUNC(dgemm)
#define MAX_NUM 100.0
//#define MAX_NUM DBL_MAX
#else
#define GEMM   BLASFUNC(sgemm)
#define MAX_NUM 100.0
//#define MAX_NUM FLT_MAX
#endif

int main(int argc, char *argv[]){

  IFLOAT *a, *b;
  FLOAT *c, *sc, *sd;
  FLOAT alpha[] = {2.0, 1.0, 0.0};
  //FLOAT alpha[] = {1.0, 1.0, 0.0};
  FLOAT beta [] = {2.0, 0.0};
  char transa = 'N';
  char transb = 'N';
  blasint m, n, k, i, j, lda, ldb, ldc;
  int loops = 2000;
  int has_param_m = 0;
  int has_param_n = 0;
  int has_param_k = 0;
  char *p;
  int error_flag = 0;

  int from =   1;
  int to   = 1500;
  int step =   1;

  double time1, timeg;

  argc--;argv++;

  if (argc > 0) { from = atol(*argv);            argc--; argv++; }
  if (argc > 0) { to   = MAX(atol(*argv), from); argc--; argv++; }
  if (argc > 0) { step = atol(*argv);            argc--; argv++; }

  if ((p = getenv("OPENBLAS_TRANS"))) {
    transa=*p;
    transb=*p;
  }
  if ((p = getenv("OPENBLAS_TRANSA"))) {
    transa=*p;
  }
  if ((p = getenv("OPENBLAS_TRANSB"))) {
    transb=*p;
  }
  TOUPPER(transa);
  TOUPPER(transb);

  fprintf(stderr, "From : %3d  To : %3d Step=%d : Transa=%c : Transb=%c\n", from, to, step, transa, transb);

  p = getenv("OPENBLAS_LOOPS");
  if ( p != NULL ) {
    loops = atoi(p);
  }

  if ((p = getenv("OPENBLAS_PARAM_M"))) {
    m = atoi(p);
    has_param_m=1;
  } else {
    m = to;
  }
  if ((p = getenv("OPENBLAS_PARAM_N"))) {
    n = atoi(p);
    has_param_n=1;
  } else {
    n = to;
  }
  if ((p = getenv("OPENBLAS_PARAM_K"))) {
    k = atoi(p);
    has_param_k=1;
  } else {
    k = to;
  }

  if (( a = (IFLOAT *)malloc(sizeof(IFLOAT) * m * k * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( b = (IFLOAT *)malloc(sizeof(IFLOAT) * k * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( c = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( sc = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( sd = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
#ifdef __linux
  srandom(getpid());
#endif

//#define REAL_RAND 1

#ifdef REAL_RAND
  for (i = 0; i < m * k * COMPSIZE; i++) {
    a[i] = ((IFLOAT) rand() / (IFLOAT) RAND_MAX) - 0.5;
  }
  for (i = 0; i < k * n * COMPSIZE; i++) {
    b[i] = ((IFLOAT) rand() / (IFLOAT) RAND_MAX) - 0.5;
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    c[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    sc[i] = c[i];
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    sd[i] = c[i];
  }
#else
  srand(0x216);
  //generate floating-points in [0.5, RAND_MAX)
  for (i = 0; i < m * k * COMPSIZE; i++) {
    a[i] = ((IFLOAT) rand() / (IFLOAT) RAND_MAX) * (MAX_NUM);
  }
  for (i = 0; i < k * n * COMPSIZE; i++) {
    b[i] = ((IFLOAT) rand() / (IFLOAT) RAND_MAX) * (MAX_NUM);
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    c[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) * (MAX_NUM);
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    sc[i] = c[i];
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    sd[i] = c[i];
  }
#endif


  fprintf(stderr, "          SIZE                   Flops             Time\n");

  FLOAT *dest;
#define PACK 1
#define TESTA 1
#define TESTB 1

  for (i = from; i <= to; i += step) {
   
    error_flag = 0;
  //test real random
  //printf("\na[0] = %f, b[0] = %f\n", a[0], b[0], c[0]);
    timeg=0;

    if (!has_param_m) { m = i; }
    if (!has_param_n) { n = i; }
    if (!has_param_k) { k = i; }

    if (transa == 'N') { lda = m; }
    else { lda = k; }
    if (transb == 'N') { ldb = k; }
    else { ldb = n; }
    ldc = m;

#ifdef PACK 
    if (( dest = (FLOAT *)malloc(sizeof(FLOAT) * m * n * k + 1048576 + 1024)) == NULL) {
        fprintf(stderr,"Out of Memory!!\n");exit(1);
    }
#endif

    GEMM (&transa, &transb, &m, &n, &k, alpha, a, &lda, b, &ldb, beta, c, &ldc);
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha[0], a, lda, b, ldb, beta[0], c, ldc);


#ifdef PACK
#ifdef DOUBLE 
    //pack and compute
    //test A packed 
    #ifdef TESTA
    cblas_dgemm_pack(CblasColMajor, CblasAMatrix, CblasNoTrans, m, n, k, alpha[0], a, lda, dest);
    cblas_dgemm_compute(CblasColMajor, CblasPacked, CblasNoTrans, m, n, k, dest, lda, b, ldb, beta[0], sc, ldc);
    #endif
    //test B packed 
    #ifdef TESTB
    cblas_dgemm_pack(CblasColMajor, CblasBMatrix, CblasNoTrans, m, n, k, alpha[0], b, ldb, dest);
    cblas_dgemm_compute(CblasColMajor, CblasNoTrans, CblasPacked, m, n, k, a, lda, dest, ldb, beta[0], sd, ldc);
    #endif
#else
    //pack and compute
    //test A packed 
    #ifdef TESTA
    cblas_sgemm_pack(CblasColMajor, CblasAMatrix, CblasNoTrans, m, n, k, alpha[0], a, lda, dest);
    cblas_sgemm_compute(CblasColMajor, CblasPacked, CblasNoTrans, m, n, k, dest, lda, b, ldb, beta[0], sc, ldc);
    #endif
    //test B packed 
    #ifdef TESTB
    cblas_sgemm_pack(CblasColMajor, CblasBMatrix, CblasNoTrans, m, n, k, alpha[0], b, ldb, dest);
    cblas_sgemm_compute(CblasColMajor, CblasNoTrans, CblasPacked, m, n, k, a, lda, dest, ldb, beta[0], sd, ldc);
    #endif
#endif
#else
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha[0], a, lda, b, ldb, beta[0], sc, ldc);
#endif
    
#ifdef PACK
    //correctness test
    //test A packed    
    #ifdef TESTA
    for(j = 0; j < m * n; j++)
    {
        #ifdef DOUBLE
        //fprintf(stderr, "c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
            //fprintf(stderr, "For case %d, m = %d, n = %d, k = %d\n", i, m, n, k);
        //if(fabs(c[j] - sc[j])/fabs(c[j]) <= 1e-15)
        //{
        //    printf("there is a test passed the test, c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
        ////    exit(-1);
        //}
        //if((c[j] - sc[j])/c[j] >= 1e-5)
        if(fabs(c[j] - sc[j])/fabs(c[j]) >= 1e-15)
        {
            //fprintf(stderr, "Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n\n", j, c[j], j, sc[j]);
            //if(j > m*256)
            //if(j < 512)
            printf("Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
            error_flag = 1;
            //exit(-1);
        }
        #else
        if(fabs(c[j] - sc[j])/fabs(c[j]) >= 1e-6)
        {
            //fprintf(stderr, "Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n\n", j, c[j], j, sc[j]);
            //if(j > m*256)
            //if(j < 512)
            printf("Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
            error_flag = 1;
            //exit(-1);
        }
        #endif
    }
    #endif
    //test B packed    
    #ifdef TESTB
    for(j = 0; j < m * n; j++)
    {
        #ifdef DOUBLE
        //fprintf(stderr, "c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
            //fprintf(stderr, "For case %d, m = %d, n = %d, k = %d\n", i, m, n, k);
        //if(fabs(c[j] - sc[j])/fabs(c[j]) <= 1e-15)
        //{
        //    printf("there is a test passed the test, c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
        ////    exit(-1);
        //}
        //if((c[j] - sc[j])/c[j] >= 1e-5)
        if(fabs(c[j] - sc[j])/fabs(c[j]) >= 1e-15)
        {
            //fprintf(stderr, "Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n\n", j, c[j], j, sc[j]);
            //if(j > m*256)
            //if(j < 512)
            printf("Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
            error_flag = 1;
            //exit(-1);
        }
        #else
        if(fabs(c[j] - sc[j])/fabs(c[j]) >= 1e-6)
        {
            //fprintf(stderr, "Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n\n", j, c[j], j, sc[j]);
            //if(j > m*256)
            //if(j < 512)
            printf("Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
            error_flag = 1;
            //exit(-1);
        }
        #endif
    }
    #endif
#else
    for(j = 0; j < m * n; j++)
    {
        //fprintf(stderr, "c[%d] = %.17f, sc[%d] = %.17f\n", j, c[j], j, sc[j]);
        if((c[j] - sc[j])/c[j] >= 1e-5)
        {
            fprintf(stderr, "For case %d, m = %d, n = %d, k = %d\n", i, m, n, k);
            fprintf(stderr, "Fail to pass the test for A packed, c[%d] = %.17f, sc[%d] = %.17f\n\n", j, c[j], j, sc[j]);
            error_flag = 1;
            //exit(-1);
        }
    }
#endif    

    //fprintf(stderr, "FOR M=%4d, N=%4d, K=%4d: ", (int)m, (int)n, (int)k);
    printf("FOR M=%4d, N=%4d, K=%4d: ", (int)m, (int)n, (int)k);
    if(error_flag) 
        printf("\nTest failed!\n***********************************\n\n");
    else
        printf("\nTest passed, done well!\n***********************************\n\n");

//    begin();
//
//    for (j=0; j<loops; j++) {
//#ifdef PACK
//#ifdef DOUBLE 
//    //pack and compute                                                                                          
//    cblas_dgemm_pack(CblasColMajor, CblasAMatrix, CblasNoTrans, m, n, k, alpha[0], a, lda, dest);                  
//    cblas_dgemm_compute(CblasColMajor, CblasPacked, CblasNoTrans, m, n, k, dest, lda, b, ldb, beta[0], sc, ldc);
//#else
//    //pack and compute                                                                                          
//    cblas_sgemm_pack(CblasColMajor, CblasAMatrix, CblasNoTrans, m, n, k, alpha[0], a, lda, dest);                  
//    cblas_sgemm_compute(CblasColMajor, CblasPacked, CblasNoTrans, m, n, k, dest, lda, b, ldb, beta[0], sc, ldc);
//#endif 
//#else 
//    GEMM (&transa, &transb, &m, &n, &k, alpha, a, &lda, b, &ldb, beta, c, &ldc);
//#endif 
//    }
//
//    end();
//    time1 = getsec();
//
//    timeg = time1/loops;
//    fprintf(stderr,
//	    " %10.2f MFlops %10.6f sec\n",
//	    COMPSIZE * COMPSIZE * 2. * (double)k * (double)m * (double)n / timeg * 1.e-6, time1);
//    
  }
  
  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
