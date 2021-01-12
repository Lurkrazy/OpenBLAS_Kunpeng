/*******************************************************************************
Copyright (c) 2015, The OpenBLAS Project
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
*******************************************************************************/

#define ASSEMBLER
#include "common.h"

/*                   X0          X1          X2          s0        X3        x4       x5           x6  */
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc) */

#define origM		x0
#define origN		x1
#define origK		x2
#define origPA		x3
#define origPB		x4
#define pC		x5
#define LDC		x6
#define temp		x7
#define counterL	x8
#define counterI	x9
#define counterJ	x10
#define pB		x11
#define pCRow0		x12
#define pCRow1		x13
#define pCRow2		x14
#define pCRow3		x15
#define pA		x16
#define alpha		w17
#define pCRow4		x18
#define pCRow5		x19
#define pCRow6		x20
#define pCRow7		x21

#define alpha0		s3
#define alphaV0		v3.s[0]

#define A_PRE_SIZE	5120	
#define B_PRE_SIZE	448
#define C_PRE_SIZE	320


// 00 origM
// 01 origN
// 02 origK
// 03 origPA
// 04 origPB
// 05 pC
// 06 origLDC -> LDC
// 07 offset
// 08 counterL
// 09 counterI
// 10 counterJ
// 11 pB
// 12 pCRow0
// 13 pCRow1
// 14 pCRow2
// 15 pA
// 16 temp
// 17
// 18 must save
// 19 must save
// 20 must save
// 21 must save
// 22 must save
// 23 must save
// 24 must save
// 25 must save
// 26 must save
// 27 must save
// 28 must save
// 29 frame
// 30 link
// 31 sp

//v00 ALPHA -> pA0_00, pA0_01, pA0_02, pA0_03
//v01 pA0_04, pA0_05, pA0_06, pA0_07
//v02 pA0_08, pA0_09, pA0_10, pA0_11
//v03 pA0_12, pA0_13, pA0_14, pA0_15
//v04 pA1_00, pA1_01, pA1_02, pA1_03
//v05 pA1_04, pA1_05, pA1_06, pA1_07
//v06 pA1_08, pA1_09, pA1_10, pA1_11
//v07 pA1_12, pA1_13, pA1_14, pA1_15
//v08 must save pB00
//v09 must save pB01
//v10 must save pB02
//v11 must save pB03
//v12 must save pB10
//v13 must save pB11
//v14 must save pB12
//v15 must save pB13
//v16 must save C00, C01, C02, C03
//v17 must save C04, C05, C06, C07
//v18 C08, C09, C10, C11
//v19 C12, C13, C14, C15
//v20 C16, C17, C18, C19
//v21 C20, C21, C22, C23
//v22 C24, C25, C26, C27
//v23 C28, C29, C30, C31
//v24 C32, C33, C34, C35
//v25 C36, C37, C38, C39
//v26 C40, C41, C42, C43
//v27 C44, C45, C46, C47
//v28 C48, C49, C50, C51
//v29 C52, C53, C54, C55
//v30 C56, C57, C58, C59
//v31 C60, C61, C62, C63

/*******************************************************************************
* Macro definitions
*******************************************************************************/

.macro INIT12x8
	fmov		s8,  wzr
	fmov		s9,  wzr
	fmov		s10, s8
	fmov		s11, s9
	fmov		s12, s8
	fmov		s13, s9
	fmov		s14, s8
	fmov		s15, s9
	fmov		s16, wzr
	fmov		s17, wzr
	fmov		s18, s16
	fmov		s19, s17
	fmov		s20, wzr
	fmov		s21, s16
	fmov		s22, s17
	fmov		s23, s18
	fmov		s24, wzr
	fmov		s25, s16
	fmov		s26, s17
	fmov		s27, s18
	fmov		s28, wzr
	fmov		s29, s16
	fmov		s30, s17
	fmov		s31, s18
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL12x8_I
	ldp	q0, q1, [pA], #32

	ldp	q4, q5, [pB], #32

	fmul	v8.4s, v0.4s, v4.s[0]
	fmul	v11.4s, v0.4s, v4.s[1]
	ldr	q6, [pA] 
	add	pA, pA, #16
	fmul	v14.4s, v0.4s, v4.s[2]
	fmul	v17.4s, v0.4s, v4.s[3]
	fmul	v20.4s, v0.4s, v5.s[0]
	fmul	v23.4s, v0.4s, v5.s[1]
	fmul	v26.4s, v0.4s, v5.s[2]
	fmul	v29.4s, v0.4s, v5.s[3]


	fmul	v9.4s, v1.4s, v4.s[0]
	fmul	v12.4s, v1.4s, v4.s[1]
	fmul	v15.4s, v1.4s, v4.s[2]
	fmul	v18.4s, v1.4s, v4.s[3]
	fmul	v21.4s, v1.4s, v5.s[0]
	fmul	v24.4s, v1.4s, v5.s[1]
	fmul	v27.4s, v1.4s, v5.s[2]
	fmul	v30.4s, v1.4s, v5.s[3]
	ldp	q0, q1, [pA], #32

	ldp	q2, q3, [pB], #32

	fmul	v10.4s, v6.4s, v4.s[0]
	fmul	v13.4s, v6.4s, v4.s[1]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	fmul	v16.4s, v6.4s, v4.s[2]
	fmul	v19.4s, v6.4s, v4.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
	fmul	v22.4s, v6.4s, v5.s[0]
	fmul	v25.4s, v6.4s, v5.s[1]
	fmul	v28.4s, v6.4s, v5.s[2]
	fmul	v31.4s, v6.4s, v5.s[3]
	ldr	q7, [pA]
	add	pA, pA, #16

.endm

.macro KERNEL12x8_M1
	fmla	v8.4s, v0.4s, v4.s[0]
	fmla	v11.4s, v0.4s, v4.s[1]
	ldr	q7, [pA, #32] 
	fmla	v14.4s, v0.4s, v4.s[2]
	fmla	v17.4s, v0.4s, v4.s[3]
	fmla	v20.4s, v0.4s, v5.s[0]
	fmla	v23.4s, v0.4s, v5.s[1]
	ldp	q2, q3, [pB], #32
	fmla	v26.4s, v0.4s, v5.s[2]
	fmla	v29.4s, v0.4s, v5.s[3]


	fmla	v9.4s, v1.4s, v4.s[0]
	fmla	v12.4s, v1.4s, v4.s[1]
	fmla	v15.4s, v1.4s, v4.s[2]
	fmla	v18.4s, v1.4s, v4.s[3]
	fmla	v21.4s, v1.4s, v5.s[0]
	fmla	v24.4s, v1.4s, v5.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v5.s[2]
	fmla	v30.4s, v1.4s, v5.s[3]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v6.4s, v4.s[0]
	fmla	v13.4s, v6.4s, v4.s[1]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	fmla	v16.4s, v6.4s, v4.s[2]
	fmla	v19.4s, v6.4s, v4.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
	fmla	v22.4s, v6.4s, v5.s[0]
	fmla	v25.4s, v6.4s, v5.s[1]
	fmla	v28.4s, v6.4s, v5.s[2]
	fmla	v31.4s, v6.4s, v5.s[3]
.endm

.macro KERNEL12x8_M2
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	ldr	q6, [pA, #32] 
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	ldp	q4, q5, [pB], #32
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	fmla	v21.4s, v1.4s, v3.s[0]
	fmla	v24.4s, v1.4s, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v3.s[2]
	fmla	v30.4s, v1.4s, v3.s[3]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v7.4s, v2.s[0]
	fmla	v13.4s, v7.4s, v2.s[1]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	fmla	v16.4s, v7.4s, v2.s[2]
	fmla	v19.4s, v7.4s, v2.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
	fmla	v22.4s, v7.4s, v3.s[0]
	fmla	v25.4s, v7.4s, v3.s[1]
	fmla	v28.4s, v7.4s, v3.s[2]
	fmla	v31.4s, v7.4s, v3.s[3]
.endm

.macro KERNEL12x8_E
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	fmla	v21.4s, v1.4s, v3.s[0]
	fmla	v24.4s, v1.4s, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v3.s[2]
	fmla	v30.4s, v1.4s, v3.s[3]
	fmla	v10.4s, v7.4s, v2.s[0]
	fmla	v13.4s, v7.4s, v2.s[1]
	fmla	v16.4s, v7.4s, v2.s[2]
	fmla	v19.4s, v7.4s, v2.s[3]
	fmla	v22.4s, v7.4s, v3.s[0]
	fmla	v25.4s, v7.4s, v3.s[1]
	fmla	v28.4s, v7.4s, v3.s[2]
	fmla	v31.4s, v7.4s, v3.s[3]
.endm

.macro KERNEL12x8_SUB
	ldp	q0, q1, [pA], #32
	ldp	q2, q3, [pB], #32

	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]

	ldr	q6, [pA]
	add	pA, pA, #16


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	fmla	v21.4s, v1.4s, v3.s[0]
	fmla	v24.4s, v1.4s, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v3.s[2]
	fmla	v30.4s, v1.4s, v3.s[3]
	fmla	v10.4s, v6.4s, v2.s[0]
	fmla	v13.4s, v6.4s, v2.s[1]
	fmla	v16.4s, v6.4s, v2.s[2]
	fmla	v19.4s, v6.4s, v2.s[3]
	fmla	v22.4s, v6.4s, v3.s[0]
	fmla	v25.4s, v6.4s, v3.s[1]
	fmla	v28.4s, v6.4s, v3.s[2]
	fmla	v31.4s, v6.4s, v3.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE12x8

        fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow0]
	ldr	q2, [pCRow0, #32]
	fmla	v0.4s, v8.4s, alphaV0
	fmla	v1.4s, v9.4s, alphaV0
	fmla	v2.4s, v10.4s, alphaV0
	stp 	q0, q1, [pCRow0]
	str	q2, [pCRow0, #32]
	add	pCRow0, pCRow0, #48

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow1]
	ldr	q2, [pCRow1, #32] 
	fmla	v0.4s, v11.4s, alphaV0
	fmla	v1.4s, v12.4s, alphaV0
	fmla	v2.4s, v13.4s, alphaV0
	stp 	q0, q1, [pCRow1]
	str	q2, [pCRow1, #32]
	add	pCRow1, pCRow1, #48

	prfm	PLDL2KEEP, [pCRow2, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow2]
	ldr	q2, [pCRow2, #32]
	fmla	v0.4s, v14.4s, alphaV0
	fmla	v1.4s, v15.4s, alphaV0
	fmla	v2.4s, v16.4s, alphaV0
	stp 	q0, q1, [pCRow2]
	str	q2, [pCRow2, #32]
	add	pCRow2, pCRow2, #48

	prfm	PLDL2KEEP, [pCRow3, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow3]
	ldr	q2, [pCRow3, #32]
	fmla	v0.4s, v17.4s, alphaV0
	fmla	v1.4s, v18.4s, alphaV0
	fmla	v2.4s, v19.4s, alphaV0
	stp 	q0, q1, [pCRow3]
	str	q2, [pCRow3, #32]
	add	pCRow3, pCRow3, #48

	prfm	PLDL2KEEP, [pCRow4, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow4]
	ldr	q2, [pCRow4, #32]
	fmla	v0.4s, v20.4s, alphaV0
	fmla	v1.4s, v21.4s, alphaV0
	fmla	v2.4s, v22.4s, alphaV0
	stp 	q0, q1, [pCRow4]
	str	q2, [pCRow4, #32]
	add	pCRow4, pCRow4, #48

	prfm	PLDL2KEEP, [pCRow5, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow5]
	ldr	q2, [pCRow5, #32] 
	fmla	v0.4s, v23.4s, alphaV0
	fmla	v1.4s, v24.4s, alphaV0
	fmla	v2.4s, v25.4s, alphaV0
	stp 	q0, q1, [pCRow5]
	str	q2, [pCRow5, #32]
	add	pCRow5, pCRow5, #48

	prfm	PLDL2KEEP, [pCRow6, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow6]
	ldr	q2, [pCRow6, #32]
	fmla	v0.4s, v26.4s, alphaV0
	fmla	v1.4s, v27.4s, alphaV0
	fmla	v2.4s, v28.4s, alphaV0
	stp 	q0, q1, [pCRow6]
	str	q2, [pCRow6, #32]
	add	pCRow6, pCRow6, #48

	prfm	PLDL2KEEP, [pCRow7, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow7]
	ldr	q2, [pCRow7, #32]
	fmla	v0.4s, v29.4s, alphaV0
	fmla	v1.4s, v30.4s, alphaV0
	fmla	v2.4s, v31.4s, alphaV0
	stp 	q0, q1, [pCRow7]
	str	q2, [pCRow7, #32]
	add	pCRow7, pCRow7, #48

.endm

/******************************************************************************/
.macro INIT8x8
	fmov		s8,  wzr
	fmov		s9,  wzr
	fmov		s11, s9
	fmov		s12, s8
	fmov		s14, s8
	fmov		s15, s9
	fmov		s17, wzr
	fmov		s18, s11
	fmov		s20, wzr
	fmov		s21, s11
	fmov		s23, s12
	fmov		s24, wzr
	fmov		s26, s12
	fmov		s27, s14
	fmov		s29, s15
	fmov		s30, s17
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL8x8_I
	ldp	q0, q1, [pA], #32

	ldp	q4, q5, [pB], #32

	fmul	v8.4s, v0.4s, v4.s[0]
	fmul	v11.4s, v0.4s, v4.s[1]
	fmul	v14.4s, v0.4s, v4.s[2]
	fmul	v17.4s, v0.4s, v4.s[3]
	fmul	v20.4s, v0.4s, v5.s[0]
	fmul	v23.4s, v0.4s, v5.s[1]
	fmul	v26.4s, v0.4s, v5.s[2]
	fmul	v29.4s, v0.4s, v5.s[3]


	fmul	v9.4s, v1.4s, v4.s[0]
	fmul	v12.4s, v1.4s, v4.s[1]
	fmul	v15.4s, v1.4s, v4.s[2]
	fmul	v18.4s, v1.4s, v4.s[3]
	fmul	v21.4s, v1.4s, v5.s[0]
	fmul	v24.4s, v1.4s, v5.s[1]
	fmul	v27.4s, v1.4s, v5.s[2]
	fmul	v30.4s, v1.4s, v5.s[3]
	ldp	q0, q1, [pA], #32

	ldp	q2, q3, [pB], #32

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]

.endm

.macro KERNEL8x8_M1
	fmla	v8.4s, v0.4s, v4.s[0]
	fmla	v11.4s, v0.4s, v4.s[1]
	fmla	v14.4s, v0.4s, v4.s[2]
	fmla	v17.4s, v0.4s, v4.s[3]
	fmla	v20.4s, v0.4s, v5.s[0]
	fmla	v23.4s, v0.4s, v5.s[1]
	ldp	q2, q3, [pB], #32
	fmla	v26.4s, v0.4s, v5.s[2]
	fmla	v29.4s, v0.4s, v5.s[3]


	fmla	v9.4s, v1.4s, v4.s[0]
	fmla	v12.4s, v1.4s, v4.s[1]
	fmla	v15.4s, v1.4s, v4.s[2]
	fmla	v18.4s, v1.4s, v4.s[3]
	fmla	v21.4s, v1.4s, v5.s[0]
	fmla	v24.4s, v1.4s, v5.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v5.s[2]
	fmla	v30.4s, v1.4s, v5.s[3]

	ldp	q0, q1, [pA], #32

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL8x8_M2
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	ldp	q4, q5, [pB], #32
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	fmla	v21.4s, v1.4s, v3.s[0]
	fmla	v24.4s, v1.4s, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v3.s[2]
	fmla	v30.4s, v1.4s, v3.s[3]

	ldp	q0, q1, [pA], #32

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL8x8_E
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	fmla	v21.4s, v1.4s, v3.s[0]
	fmla	v24.4s, v1.4s, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v3.s[2]
	fmla	v30.4s, v1.4s, v3.s[3]
.endm

.macro KERNEL8x8_SUB
	ldp	q0, q1, [pA], #32
	ldp	q2, q3, [pB], #32

	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]

	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	fmla	v21.4s, v1.4s, v3.s[0]
	fmla	v24.4s, v1.4s, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v27.4s, v1.4s, v3.s[2]
	fmla	v30.4s, v1.4s, v3.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE8x8
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow0]
	fmla	v0.4s, v8.4s, alphaV0
	fmla	v1.4s, v9.4s, alphaV0
	stp 	q0, q1, [pCRow0]
	add	pCRow0, pCRow0, #32

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow1]
	fmla	v0.4s, v11.4s, alphaV0
	fmla	v1.4s, v12.4s, alphaV0
	stp 	q0, q1, [pCRow1]
	
	add	pCRow1, pCRow1, #32

	prfm	PLDL2KEEP, [pCRow2, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow2]
	fmla	v0.4s, v14.4s, alphaV0
	fmla	v1.4s, v15.4s, alphaV0
	stp 	q0, q1, [pCRow2]

	add	pCRow2, pCRow2, #32

	prfm	PLDL2KEEP, [pCRow3, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow3]
	fmla	v0.4s, v17.4s, alphaV0
	fmla	v1.4s, v18.4s, alphaV0
	stp 	q0, q1, [pCRow3]

	add	pCRow3, pCRow3, #32

	prfm	PLDL2KEEP, [pCRow4, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow4]
	fmla	v0.4s, v20.4s, alphaV0
	fmla	v1.4s, v21.4s, alphaV0
	stp 	q0, q1, [pCRow4]
	add	pCRow4, pCRow4, #32

	prfm	PLDL2KEEP, [pCRow5, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow5]
	fmla	v0.4s, v23.4s, alphaV0
	fmla	v1.4s, v24.4s, alphaV0
	stp 	q0, q1, [pCRow5]
	
	add	pCRow5, pCRow5, #32

	prfm	PLDL2KEEP, [pCRow6, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow6]
	fmla	v0.4s, v26.4s, alphaV0
	fmla	v1.4s, v27.4s, alphaV0
	stp 	q0, q1, [pCRow6]

	add	pCRow6, pCRow6, #32

	prfm	PLDL2KEEP, [pCRow7, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow7]
	fmla	v0.4s, v29.4s, alphaV0
	fmla	v1.4s, v30.4s, alphaV0
	stp 	q0, q1, [pCRow7]

	add	pCRow7, pCRow7, #32

.endm

/******************************************************************************/

.macro INIT4x8
	fmov		s8,  wzr
	fmov		s11, s8
	fmov		s14, s8
	fmov		s17, wzr
	fmov		s20, wzr
	fmov		s23, s11
	fmov		s26, s14
	fmov		s29, s17
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL4x8_I
	ldr	q0, [pA]
	add	pA, pA, #16
	ldp	q4, q5, [pB], #32

	fmul	v8.4s, v0.4s, v4.s[0]
	fmul	v11.4s, v0.4s, v4.s[1]
	fmul	v14.4s, v0.4s, v4.s[2]
	fmul	v17.4s, v0.4s, v4.s[3]
	fmul	v20.4s, v0.4s, v5.s[0]
	fmul	v23.4s, v0.4s, v5.s[1]
	fmul	v26.4s, v0.4s, v5.s[2]
	fmul	v29.4s, v0.4s, v5.s[3]


	ldr	q0, [pA]
	add	pA, pA, #16

	ldp	q2, q3, [pB], #32

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]

.endm

.macro KERNEL4x8_M1
	fmla	v8.4s, v0.4s, v4.s[0]
	fmla	v11.4s, v0.4s, v4.s[1]
	fmla	v14.4s, v0.4s, v4.s[2]
	fmla	v17.4s, v0.4s, v4.s[3]
	fmla	v20.4s, v0.4s, v5.s[0]
	fmla	v23.4s, v0.4s, v5.s[1]
	ldp	q2, q3, [pB], #32
	fmla	v26.4s, v0.4s, v5.s[2]
	fmla	v29.4s, v0.4s, v5.s[3]


	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldr	q0, [pA]
	add	pA, pA, #16

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL4x8_M2
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	ldp	q4, q5, [pB], #32
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]


	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldr	q0, [pA]
	add	pA, pA, #16

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL4x8_E
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]


	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
.endm

.macro KERNEL4x8_SUB
	ldr	q0, [pA]
	add	pA, pA, #16
	ldp	q2, q3, [pB], #32

	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	fmla	v20.4s, v0.4s, v3.s[0]
	fmla	v23.4s, v0.4s, v3.s[1]
	fmla	v26.4s, v0.4s, v3.s[2]
	fmla	v29.4s, v0.4s, v3.s[3]

	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE4x8
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]
	ldr	q0,  [pCRow0]
	fmla	v0.4s, v8.4s, alphaV0
	str 	q0, [pCRow0]
	add	pCRow0, pCRow0, #16

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldr	q0, [pCRow1]
	fmla	v0.4s, v11.4s, alphaV0
	str 	q0, [pCRow1]
	add	pCRow1, pCRow1, #16

	prfm	PLDL2KEEP, [pCRow2, #C_PRE_SIZE]
	ldr	q0, [pCRow2]
	fmla	v0.4s, v14.4s, alphaV0
	str 	q0, [pCRow2]
	add	pCRow2, pCRow2, #16

	prfm	PLDL2KEEP, [pCRow3, #C_PRE_SIZE]
	ldr	q0, [pCRow3]
	fmla	v0.4s, v17.4s, alphaV0
	str 	q0, [pCRow3]
	add	pCRow3, pCRow3, #16

	prfm	PLDL2KEEP, [pCRow4, #C_PRE_SIZE]
	ldr	q0, [pCRow4]
	fmla	v0.4s, v20.4s, alphaV0
	str 	q0, [pCRow4]
	add	pCRow4, pCRow4, #16

	prfm	PLDL2KEEP, [pCRow5, #C_PRE_SIZE]
	ldr	q0, [pCRow5]
	fmla	v0.4s, v23.4s, alphaV0
	str 	q0, [pCRow5]
	add	pCRow5, pCRow5, #16

	prfm	PLDL2KEEP, [pCRow6, #C_PRE_SIZE]
	ldr	q0, [pCRow6]
	fmla	v0.4s, v26.4s, alphaV0
	str 	q0, [pCRow6]
	add	pCRow6, pCRow6, #16

	prfm	PLDL2KEEP, [pCRow7, #C_PRE_SIZE]
	ldr	q0, [pCRow7]
	fmla	v0.4s, v29.4s, alphaV0
	str 	q0, [pCRow7]
	add	pCRow7, pCRow7, #16

.endm

/******************************************************************************/
.macro INIT2x8
	fmov		s8,  wzr
	fmov		s9,  wzr
	fmov		s11, s9
	fmov		s12, s8
	fmov		s14, s8
	fmov		s15, s9
	fmov		s17, wzr
	fmov		s18, s11
	fmov		s20, wzr
	fmov		s21, s11
	fmov		s23, s12
	fmov		s24, wzr
	fmov		s26, s12
	fmov		s27, s14
	fmov		s29, s15
	fmov		s30, s17
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL2x8_I
	ldp	s0, s1, [pA], #8

	ldp	q4, q5, [pB], #32

	fmul	s8, s0, v4.s[0]
	fmul	s11, s0, v4.s[1]
	fmul	s14, s0, v4.s[2]
	fmul	s17, s0, v4.s[3]
	fmul	s20, s0, v5.s[0]
	fmul	s23, s0, v5.s[1]
	fmul	s26, s0, v5.s[2]
	fmul	s29, s0, v5.s[3]


	fmul	s9, s1, v4.s[0]
	fmul	s12, s1, v4.s[1]
	fmul	s15, s1, v4.s[2]
	fmul	s18, s1, v4.s[3]
	fmul	s21, s1, v5.s[0]
	fmul	s24, s1, v5.s[1]
	fmul	s27, s1, v5.s[2]
	fmul	s30, s1, v5.s[3]
	ldp	s0, s1, [pA], #8

	ldp	q2, q3, [pB], #32

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]

.endm

.macro KERNEL2x8_M1
	fmla	s8, s0, v4.s[0]
	fmla	s11, s0, v4.s[1]
	fmla	s14, s0, v4.s[2]
	fmla	s17, s0, v4.s[3]
	fmla	s20, s0, v5.s[0]
	fmla	s23, s0, v5.s[1]
	ldp	q2, q3, [pB], #32
	fmla	s26, s0, v5.s[2]
	fmla	s29, s0, v5.s[3]


	fmla	s9, s1, v4.s[0]
	fmla	s12, s1, v4.s[1]
	fmla	s15, s1, v4.s[2]
	fmla	s18, s1, v4.s[3]
	fmla	s21, s1, v5.s[0]
	fmla	s24, s1, v5.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	s27, s1, v5.s[2]
	fmla	s30, s1, v5.s[3]

	ldp	s0, s1, [pA], #8

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL2x8_M2
	fmla	s8, s0, v2.s[0]
	fmla	s11, s0, v2.s[1]
	fmla	s14, s0, v2.s[2]
	fmla	s17, s0, v2.s[3]
	fmla	s20, s0, v3.s[0]
	fmla	s23, s0, v3.s[1]
	ldp	q4, q5, [pB], #32
	fmla	s26, s0, v3.s[2]
	fmla	s29, s0, v3.s[3]


	fmla	s9, s1, v2.s[0]
	fmla	s12, s1, v2.s[1]
	fmla	s15, s1, v2.s[2]
	fmla	s18, s1, v2.s[3]
	fmla	s21, s1, v3.s[0]
	fmla	s24, s1, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	s27, s1, v3.s[2]
	fmla	s30, s1, v3.s[3]

	ldp	s0, s1, [pA], #8

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL2x8_E
	fmla	s8, s0, v2.s[0]
	fmla	s11, s0, v2.s[1]
	fmla	s14, s0, v2.s[2]
	fmla	s17, s0, v2.s[3]
	fmla	s20, s0, v3.s[0]
	fmla	s23, s0, v3.s[1]
	fmla	s26, s0, v3.s[2]
	fmla	s29, s0, v3.s[3]


	fmla	s9, s1, v2.s[0]
	fmla	s12, s1, v2.s[1]
	fmla	s15, s1, v2.s[2]
	fmla	s18, s1, v2.s[3]
	fmla	s21, s1, v3.s[0]
	fmla	s24, s1, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	s27, s1, v3.s[2]
	fmla	s30, s1, v3.s[3]
.endm

.macro KERNEL2x8_SUB
	ldp	s0, s1, [pA], #8
	ldp	q2, q3, [pB], #32

	fmla	s8, s0, v2.s[0]
	fmla	s11, s0, v2.s[1]
	fmla	s14, s0, v2.s[2]
	fmla	s17, s0, v2.s[3]
	fmla	s20, s0, v3.s[0]
	fmla	s23, s0, v3.s[1]
	fmla	s26, s0, v3.s[2]
	fmla	s29, s0, v3.s[3]

	fmla	s9, s1, v2.s[0]
	fmla	s12, s1, v2.s[1]
	fmla	s15, s1, v2.s[2]
	fmla	s18, s1, v2.s[3]
	fmla	s21, s1, v3.s[0]
	fmla	s24, s1, v3.s[1]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	s27, s1, v3.s[2]
	fmla	s30, s1, v3.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE2x8
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldp	s0, s1, [pCRow0]
	fmla	s0, s8, alphaV0
	fmla	s1, s9, alphaV0
	stp 	s0, s1, [pCRow0]
	add	pCRow0, pCRow0, #8

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldp	s0, s1, [pCRow1]
	fmla	s0, s11, alphaV0
	fmla	s1, s12, alphaV0
	stp 	s0, s1, [pCRow1]
	
	add	pCRow1, pCRow1, #8

	prfm	PLDL2KEEP, [pCRow2, #C_PRE_SIZE]

	ldp	s0, s1, [pCRow2]
	fmla	s0, s14, alphaV0
	fmla	s1, s15, alphaV0
	stp 	s0, s1, [pCRow2]

	add	pCRow2, pCRow2, #8

	prfm	PLDL2KEEP, [pCRow3, #C_PRE_SIZE]

	ldp	s0, s1, [pCRow3]
	fmla	s0, s17, alphaV0
	fmla	s1, s18, alphaV0
	stp 	s0, s1, [pCRow3]

	add	pCRow3, pCRow3, #8

	prfm	PLDL2KEEP, [pCRow4, #C_PRE_SIZE]

	ldp	s0, s1, [pCRow4]
	fmla	s0, s20, alphaV0
	fmla	s1, s21, alphaV0
	stp 	s0, s1, [pCRow4]
	add	pCRow4, pCRow4, #8

	prfm	PLDL2KEEP, [pCRow5, #C_PRE_SIZE]
	ldp	s0, s1, [pCRow5]
	fmla	s0, s23, alphaV0
	fmla	s1, s24, alphaV0
	stp 	s0, s1, [pCRow5]
	
	add	pCRow5, pCRow5, #8

	prfm	PLDL2KEEP, [pCRow6, #C_PRE_SIZE]

	ldp	s0, s1, [pCRow6]
	fmla	s0, s26, alphaV0
	fmla	s1, s27, alphaV0
	stp 	s0, s1, [pCRow6]

	add	pCRow6, pCRow6, #8

	prfm	PLDL2KEEP, [pCRow7, #C_PRE_SIZE]

	ldp	s0, s1, [pCRow7]
	fmla	s0, s29, alphaV0
	fmla	s1, s30, alphaV0
	stp 	s0, s1, [pCRow7]

	add	pCRow7, pCRow7, #8

.endm

/******************************************************************************/
.macro INIT1x8
	fmov		s8,  wzr
	fmov		s11, s9
	fmov		s14, s8
	fmov		s17, wzr
	fmov		s20, wzr
	fmov		s23, s12
	fmov		s26, s12
	fmov		s29, s15
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL1x8_I
	ldr	s0, [pA]
	add	pA, pA, #4

	ldp	q4, q5, [pB], #32

	fmul	s8, s0, v4.s[0]
	fmul	s11, s0, v4.s[1]
	fmul	s14, s0, v4.s[2]
	fmul	s17, s0, v4.s[3]
	fmul	s20, s0, v5.s[0]
	fmul	s23, s0, v5.s[1]
	fmul	s26, s0, v5.s[2]
	fmul	s29, s0, v5.s[3]

	ldr	s0, [pA]
	add	pA, pA, #4

	ldp	q2, q3, [pB], #32

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]

.endm

.macro KERNEL1x8_M1
	fmla	s8, s0, v4.s[0]
	fmla	s11, s0, v4.s[1]
	fmla	s14, s0, v4.s[2]
	fmla	s17, s0, v4.s[3]
	fmla	s20, s0, v5.s[0]
	fmla	s23, s0, v5.s[1]
	ldp	q2, q3, [pB], #32
	fmla	s26, s0, v5.s[2]
	fmla	s29, s0, v5.s[3]


	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldr	s0, [pA]
	add	pA, pA, #4

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL1x8_M2
	fmla	s8, s0, v2.s[0]
	fmla	s11, s0, v2.s[1]
	fmla	s14, s0, v2.s[2]
	fmla	s17, s0, v2.s[3]
	fmla	s20, s0, v3.s[0]
	fmla	s23, s0, v3.s[1]
	ldp	q4, q5, [pB], #32
	fmla	s26, s0, v3.s[2]
	fmla	s29, s0, v3.s[3]


	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	ldr	s0, [pA]
	add	pA, pA, #4

	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL1x8_E
	fmla	s8, s0, v2.s[0]
	fmla	s11, s0, v2.s[1]
	fmla	s14, s0, v2.s[2]
	fmla	s17, s0, v2.s[3]
	fmla	s20, s0, v3.s[0]
	fmla	s23, s0, v3.s[1]
	fmla	s26, s0, v3.s[2]
	fmla	s29, s0, v3.s[3]


	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
.endm

.macro KERNEL1x8_SUB
	ldr	s0, [pA]
	add	pA, pA, #4
	ldp	q2, q3, [pB], #32

	fmla	s8, s0, v2.s[0]
	fmla	s11, s0, v2.s[1]
	fmla	s14, s0, v2.s[2]
	fmla	s17, s0, v2.s[3]
	fmla	s20, s0, v3.s[0]
	fmla	s23, s0, v3.s[1]
	fmla	s26, s0, v3.s[2]
	fmla	s29, s0, v3.s[3]

	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE1x8
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldr	s0,  [pCRow0]
	fmla	s0, s8, alphaV0
	str 	s0, [pCRow0]
	add	pCRow0, pCRow0, #4

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldr	s0, [pCRow1]
	fmla	s0, s11, alphaV0
	str 	s0, [pCRow1]
	
	add	pCRow1, pCRow1, #4

	prfm	PLDL2KEEP, [pCRow2, #C_PRE_SIZE]

	ldr	s0, [pCRow2]
	fmla	s0, s14, alphaV0
	str 	s0, [pCRow2]

	add	pCRow2, pCRow2, #4

	prfm	PLDL2KEEP, [pCRow3, #C_PRE_SIZE]

	ldr	s0, [pCRow3]
	fmla	s0, s17, alphaV0
	str 	s0, [pCRow3]

	add	pCRow3, pCRow3, #4

	prfm	PLDL2KEEP, [pCRow4, #C_PRE_SIZE]

	ldr	s0, [pCRow4]
	fmla	s0, s20, alphaV0
	str 	s0, [pCRow4]
	add	pCRow4, pCRow4, #4

	prfm	PLDL2KEEP, [pCRow5, #C_PRE_SIZE]
	ldr	s0, [pCRow5]
	fmla	s0, s23, alphaV0
	str 	s0, [pCRow5]
	
	add	pCRow5, pCRow5, #4

	prfm	PLDL2KEEP, [pCRow6, #C_PRE_SIZE]

	ldr	s0, [pCRow6]
	fmla	s0, s26, alphaV0
	str 	s0, [pCRow6]

	add	pCRow6, pCRow6, #4

	prfm	PLDL2KEEP, [pCRow7, #C_PRE_SIZE]

	ldr	s0, [pCRow7]
	fmla	s0, s29, alphaV0
	str 	s0, [pCRow7]

	add	pCRow7, pCRow7, #4

.endm

/******************************************************************************/



.macro INIT12x4
	fmov		s8,  wzr
	fmov		s9,  wzr
	fmov		s10, s8
	fmov		s11, s9
	fmov		s12, s8
	fmov		s13, s9
	fmov		s14, s8
	fmov		s15, s9
	fmov		s16, wzr
	fmov		s17, wzr
	fmov		s18, s16
	fmov		s19, s17
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL12x4_I
	ldp	q0, q1, [pA], #32

	ldr	q4, [pB]
	add	pB, pB, #16

	fmul	v8.4s, v0.4s, v4.s[0]
	fmul	v11.4s, v0.4s, v4.s[1]
	ldr	q6, [pA] 
	add	pA, pA, #16
	fmul	v14.4s, v0.4s, v4.s[2]
	fmul	v17.4s, v0.4s, v4.s[3]


	fmul	v9.4s, v1.4s, v4.s[0]
	fmul	v12.4s, v1.4s, v4.s[1]
	fmul	v15.4s, v1.4s, v4.s[2]
	fmul	v18.4s, v1.4s, v4.s[3]
	ldp	q0, q1, [pA], #32

	ldr	q2, [pB]
	add	pB, pB, #16

	fmul	v10.4s, v6.4s, v4.s[0]
	fmul	v13.4s, v6.4s, v4.s[1]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	fmul	v16.4s, v6.4s, v4.s[2]
	fmul	v19.4s, v6.4s, v4.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
	ldr	q7, [pA]
	add	pA, pA, #16

.endm

.macro KERNEL12x4_M1
	fmla	v8.4s, v0.4s, v4.s[0]
	fmla	v11.4s, v0.4s, v4.s[1]
	ldr	q7, [pA, #32] 
	fmla	v14.4s, v0.4s, v4.s[2]
	fmla	v17.4s, v0.4s, v4.s[3]
	ldr	q2, [pB]
	add	pB, pB, #16


	fmla	v9.4s, v1.4s, v4.s[0]
	fmla	v12.4s, v1.4s, v4.s[1]
	fmla	v15.4s, v1.4s, v4.s[2]
	fmla	v18.4s, v1.4s, v4.s[3]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v6.4s, v4.s[0]
	fmla	v13.4s, v6.4s, v4.s[1]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	fmla	v16.4s, v6.4s, v4.s[2]
	fmla	v19.4s, v6.4s, v4.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL12x4_M2
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	ldr	q6, [pA, #32] 
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]
	ldr	q4, [pB]
	add	pB, pB, #16


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v7.4s, v2.s[0]
	fmla	v13.4s, v7.4s, v2.s[1]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	fmla	v16.4s, v7.4s, v2.s[2]
	fmla	v19.4s, v7.4s, v2.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL12x4_E
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v10.4s, v7.4s, v2.s[0]
	fmla	v13.4s, v7.4s, v2.s[1]
	fmla	v16.4s, v7.4s, v2.s[2]
	fmla	v19.4s, v7.4s, v2.s[3]
.endm

.macro KERNEL12x4_SUB
	ldp	q0, q1, [pA], #32
	ldr	q2, [pB]
	add	pB, pB, #16

	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v2.s[1]
	fmla	v14.4s, v0.4s, v2.s[2]
	fmla	v17.4s, v0.4s, v2.s[3]

	ldr	q6, [pA]
	add	pA, pA, #16


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v2.s[1]
	fmla	v15.4s, v1.4s, v2.s[2]
	fmla	v18.4s, v1.4s, v2.s[3]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v10.4s, v6.4s, v2.s[0]
	fmla	v13.4s, v6.4s, v2.s[1]
	fmla	v16.4s, v6.4s, v2.s[2]
	fmla	v19.4s, v6.4s, v2.s[3]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE12x4
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow0]
	ldr	q2, [pCRow0, #32]
	fmla	v0.4s, v8.4s, alphaV0
	fmla	v1.4s, v9.4s, alphaV0
	fmla	v2.4s, v10.4s, alphaV0
	stp 	q0, q1, [pCRow0]
	str	q2, [pCRow0, #32]
	add	pCRow0, pCRow0, #48

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow1]
	ldr	q2, [pCRow1, #32] 
	fmla	v0.4s, v11.4s, alphaV0
	fmla	v1.4s, v12.4s, alphaV0
	fmla	v2.4s, v13.4s, alphaV0
	stp 	q0, q1, [pCRow1]
	str	q2, [pCRow1, #32]
	
	add	pCRow1, pCRow1, #48

	prfm	PLDL2KEEP, [pCRow2, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow2]
	ldr	q2, [pCRow2, #32]
	fmla	v0.4s, v14.4s, alphaV0
	fmla	v1.4s, v15.4s, alphaV0
	fmla	v2.4s, v16.4s, alphaV0
	stp 	q0, q1, [pCRow2]
	str	q2, [pCRow2, #32]

	add	pCRow2, pCRow2, #48

	prfm	PLDL2KEEP, [pCRow3, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow3]
	ldr	q2, [pCRow3, #32]
	fmla	v0.4s, v17.4s, alphaV0
	fmla	v1.4s, v18.4s, alphaV0
	fmla	v2.4s, v19.4s, alphaV0
	stp 	q0, q1, [pCRow3]
	str	q2, [pCRow3, #32]

	add	pCRow3, pCRow3, #48

.endm

/******************************************************************************/
.macro INIT12x2
	fmov		s8,  wzr
	fmov		s9,  wzr
	fmov		s10, s8
	fmov		s11, s9
	fmov		s12, s8
	fmov		s13, s9
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL12x2_I
	ldp	q0, q1, [pA], #32

	ldp	s4, s5, [pB], #8

	fmul	v8.4s, v0.4s, v4.s[0]
	fmul	v11.4s, v0.4s, v5.s[0]
	ldr	q6, [pA] 
	add	pA, pA, #16


	fmul	v9.4s, v1.4s, v4.s[0]
	fmul	v12.4s, v1.4s, v5.s[0]
	ldp	q0, q1, [pA], #32

	ldp	s2, s3, [pB], #8

	fmul	v10.4s, v6.4s, v4.s[0]
	fmul	v13.4s, v6.4s, v5.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
	ldr	q7, [pA]
	add	pA, pA, #16

.endm

.macro KERNEL12x2_M1
	fmla	v8.4s, v0.4s, v4.s[0]
	fmla	v11.4s, v0.4s, v5.s[0]
	ldr	q7, [pA, #32] 
	ldp	s2, s3,  [pB], #8


	fmla	v9.4s, v1.4s, v4.s[0]
	fmla	v12.4s, v1.4s, v5.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v6.4s, v4.s[0]
	fmla	v13.4s, v6.4s, v5.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL12x2_M2
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v3.s[0]
	ldr	q6, [pA, #32] 
	ldp	s4, s5, [pB], #8


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v3.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v7.4s, v2.s[0]
	fmla	v13.4s, v7.4s, v3.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL12x2_E
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v3.s[0]


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v3.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v10.4s, v7.4s, v2.s[0]
	fmla	v13.4s, v7.4s, v3.s[0]
.endm

.macro KERNEL12x2_SUB
	ldp	q0, q1, [pA], #32
	ldp	s2, s3, [pB], #8

	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v3.s[0]

	ldr	q6, [pA]
	add	pA, pA, #16


	fmla	v9.4s, v1.4s, v2.s[0]
	fmla	v12.4s, v1.4s, v3.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v10.4s, v6.4s, v2.s[0]
	fmla	v13.4s, v6.4s, v3.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE12x2
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow0]
	ldr	q2, [pCRow0, #32]
	fmla	v0.4s, v8.4s, alphaV0
	fmla	v1.4s, v9.4s, alphaV0
	fmla	v2.4s, v10.4s, alphaV0
	stp 	q0, q1, [pCRow0]
	str	q2, [pCRow0, #32]
	add	pCRow0, pCRow0, #48

	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	ldp	q0, q1, [pCRow1]
	ldr	q2, [pCRow1, #32] 
	fmla	v0.4s, v11.4s, alphaV0
	fmla	v1.4s, v12.4s, alphaV0
	fmla	v2.4s, v13.4s, alphaV0
	stp 	q0, q1, [pCRow1]
	str	q2, [pCRow1, #32]
	add	pCRow1, pCRow1, #48

.endm
/******************************************************************************/
.macro INIT12x1
	fmov		s8,  wzr
	fmov		s9,  wzr
	fmov		s10, s8
.endm

/*
v8-v31 store C
v0-v2 A
v4-v5 B
v3,v6-v7 tmp
*/
.macro KERNEL12x1_I
	ldp	q0, q1, [pA], #32

	ldr	s4, [pB]
	add	pB, pB, #4

	fmul	v8.4s, v0.4s, v4.s[0]
	ldr	q6, [pA] 
	add	pA, pA, #16


	fmul	v9.4s, v1.4s, v4.s[0]
	ldp	q0, q1, [pA], #32

	ldr	s2, [pB]
	add	pB, pB, #4

	fmul	v10.4s, v6.4s, v4.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
	ldr	q7, [pA]
	add	pA, pA, #16

.endm

.macro KERNEL12x1_M1
	fmla	v8.4s, v0.4s, v4.s[0]
	ldr	q7, [pA, #32] 
	ldr	s2,  [pB]
	add	pB, pB, #4


	fmla	v9.4s, v1.4s, v4.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v6.4s, v4.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL12x1_M2
	fmla	v8.4s, v0.4s, v2.s[0]
	fmla	v11.4s, v0.4s, v3.s[0]
	ldr	q6, [pA, #32] 
	ldr	s4, [pB]
	add	pB, pB, #4

	fmla	v9.4s, v1.4s, v2.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]

	ldp	q0, q1, [pA], #32
	add	pA, pA, #16

	fmla	v10.4s, v7.4s, v2.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE+64]
.endm

.macro KERNEL12x1_E
	fmla	v8.4s, v0.4s, v2.s[0]


	fmla	v9.4s, v1.4s, v2.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v10.4s, v7.4s, v2.s[0]
.endm

.macro KERNEL12x1_SUB
	ldp	q0, q1, [pA], #32
	ldr	s2, [pB] 

	fmla	v8.4s, v0.4s, v2.s[0]

	ldr	q6, [pA]
	add	pA, pA, #16


	fmla	v9.4s, v1.4s, v2.s[0]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v10.4s, v6.4s, v2.s[0]
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
.endm

.macro SAVE12x1
	fmov	alpha0, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	ldp	q0, q1, [pCRow0]
	ldr	q2, [pCRow0, #32]
	fmla	v0.4s, v8.4s, alphaV0
	fmla	v1.4s, v9.4s, alphaV0
	fmla	v2.4s, v10.4s, alphaV0
	stp 	q0, q1, [pCRow0]
	str	q2, [pCRow0, #32]
	add	pCRow0, pCRow0, #48

.endm

/******************************************************************************/

.macro INIT8x4
	fmov		s16, wzr
	fmov		s17, wzr
	fmov		s20, wzr
	fmov		s21, s16
	fmov		s24, wzr
	fmov		s25, s16
	fmov		s28, wzr
	fmov		s29, s16
.endm

.macro KERNEL8x4_I
	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	q0, [pA], #16
	ldr	q1, [pA], #16

	fmul	v16.4s, v0.4s, v8.s[0]
	fmul	v17.4s, v1.4s, v8.s[0]
	fmul	v20.4s, v0.4s, v9.s[0]
	fmul	v21.4s, v1.4s, v9.s[0]
	fmul	v24.4s, v0.4s, v10.s[0]
	fmul	v25.4s, v1.4s, v10.s[0]
	fmul	v28.4s, v0.4s, v11.s[0]
	fmul	v29.4s, v1.4s, v11.s[0]

	ldp	s12, s13, [pB], #8
	ldp	s14, s15, [pB], #8

	ldr	q4, [pA], #16
	ldr	q5, [pA], #16
.endm

.macro KERNEL8x4_M1
	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v17.4s, v1.4s, v8.s[0]
	fmla	v20.4s, v0.4s, v9.s[0]
	fmla	v21.4s, v1.4s, v9.s[0]
	fmla	v24.4s, v0.4s, v10.s[0]
	fmla	v25.4s, v1.4s, v10.s[0]
	fmla	v28.4s, v0.4s, v11.s[0]
	fmla	v29.4s, v1.4s, v11.s[0]

	ldp	s12, s13, [pB], #8
	ldp	s14, s15, [pB], #8

	ldr	q4, [pA], #16
	ldr	q5, [pA], #16
.endm

.macro KERNEL8x4_M2
	fmla	v16.4s, v4.4s, v12.s[0]
	fmla	v17.4s, v5.4s, v12.s[0]
	fmla	v20.4s, v4.4s, v13.s[0]
	fmla	v21.4s, v5.4s, v13.s[0]
	fmla	v24.4s, v4.4s, v14.s[0]
	fmla	v25.4s, v5.4s, v14.s[0]
	fmla	v28.4s, v4.4s, v15.s[0]
	fmla	v29.4s, v5.4s, v15.s[0]

	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	q0, [pA], #16
	ldr	q1, [pA], #16
.endm

.macro KERNEL8x4_E
	fmla	v16.4s, v4.4s, v12.s[0]
	fmla	v17.4s, v5.4s, v12.s[0]
	fmla	v20.4s, v4.4s, v13.s[0]
	fmla	v21.4s, v5.4s, v13.s[0]
	fmla	v24.4s, v4.4s, v14.s[0]
	fmla	v25.4s, v5.4s, v14.s[0]
	fmla	v28.4s, v4.4s, v15.s[0]
	fmla	v29.4s, v5.4s, v15.s[0]
.endm

.macro KERNEL8x4_SUB
	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	q0, [pA], #16
	ldr	q1, [pA], #16

	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v17.4s, v1.4s, v8.s[0]
	fmla	v20.4s, v0.4s, v9.s[0]
	fmla	v21.4s, v1.4s, v9.s[0]
	fmla	v24.4s, v0.4s, v10.s[0]
	fmla	v25.4s, v1.4s, v10.s[0]
	fmla	v28.4s, v0.4s, v11.s[0]
	fmla	v29.4s, v1.4s, v11.s[0]
.endm

.macro SAVE8x4
	fmov	alpha0, alpha

	ldp	q0, q1, [pCRow0]
	fmla	v0.4s, v16.4s, alphaV0
	fmla	v1.4s, v17.4s, alphaV0
	stp 	q0, q1, [pCRow0]

	add	pCRow0, pCRow0, #32

	ldp	q0, q1, [pCRow1]
	fmla	v0.4s, v20.4s, alphaV0
	fmla	v1.4s, v21.4s, alphaV0
	stp 	q0, q1, [pCRow1]

	add	pCRow1, pCRow1, #32

	ldp	q0, q1, [pCRow2]
	fmla	v0.4s, v24.4s, alphaV0
	fmla	v1.4s, v25.4s, alphaV0
	stp 	q0, q1, [pCRow2]

	add	pCRow2, pCRow2, #32

	ldp	q0, q1, [pCRow3]
	fmla	v0.4s, v28.4s, alphaV0
	fmla	v1.4s, v29.4s, alphaV0
	stp 	q0, q1, [pCRow3]

	add	pCRow3, pCRow3, #32
.endm

/******************************************************************************/

.macro INIT4x4
	fmov		s16, wzr
	fmov		s20, wzr
	fmov		s24, wzr
	fmov		s28, wzr
.endm

.macro KERNEL4x4_I
	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	q0, [pA], #16

	fmul	v16.4s, v0.4s, v8.s[0]
	fmul	v20.4s, v0.4s, v9.s[0]
	fmul	v24.4s, v0.4s, v10.s[0]
	fmul	v28.4s, v0.4s, v11.s[0]

	ldp	s12, s13, [pB], #8
	ldp	s14, s15, [pB], #8

	ldr	q1, [pA], #16
.endm

.macro KERNEL4x4_M1
	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v20.4s, v0.4s, v9.s[0]
	fmla	v24.4s, v0.4s, v10.s[0]
	fmla	v28.4s, v0.4s, v11.s[0]

	ldp	s12, s13, [pB], #8
	ldp	s14, s15, [pB], #8

	ldr	q1, [pA], #16
.endm

.macro KERNEL4x4_M2
	fmla	v16.4s, v1.4s, v12.s[0]
	fmla	v20.4s, v1.4s, v13.s[0]
	fmla	v24.4s, v1.4s, v14.s[0]
	fmla	v28.4s, v1.4s, v15.s[0]

	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	q0, [pA], #16
.endm

.macro KERNEL4x4_E
	fmla	v16.4s, v1.4s, v12.s[0]
	fmla	v20.4s, v1.4s, v13.s[0]
	fmla	v24.4s, v1.4s, v14.s[0]
	fmla	v28.4s, v1.4s, v15.s[0]
.endm

.macro KERNEL4x4_SUB
	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	q0, [pA], #16

	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v20.4s, v0.4s, v9.s[0]
	fmla	v24.4s, v0.4s, v10.s[0]
	fmla	v28.4s, v0.4s, v11.s[0]
.endm

.macro SAVE4x4
	fmov	alpha0, alpha

	ldr 	q0, [pCRow0]
	fmla	v0.4s, v16.4s, alphaV0
	str 	q0, [pCRow0]

	add	pCRow0, pCRow0, #16

	ldr 	q0, [pCRow1]
	fmla	v0.4s, v20.4s, alphaV0
	str 	q0, [pCRow1]

	add	pCRow1, pCRow1, #16

	ldr 	q0, [pCRow2]
	fmla	v0.4s, v24.4s, alphaV0
	str 	q0, [pCRow2]

	add	pCRow2, pCRow2, #16

	ldr 	q0, [pCRow3]
	fmla	v0.4s, v28.4s, alphaV0
	str 	q0, [pCRow3]

	add	pCRow3, pCRow3, #16
.endm

/******************************************************************************/

.macro INIT2x4
	fmov		s16, wzr
	fmov		s20, s16
	fmov		s24, s20
	fmov		s28, s16
.endm

.macro KERNEL2x4_SUB
	ldp	s8, s9, [pB], #8
	ldp	s10, s11, [pB], #8

	ldr	d0, [pA], #8

	fmla	v16.2s, v0.2s, v8.s[0]
	fmla	v20.2s, v0.2s, v9.s[0]
	fmla	v24.2s, v0.2s, v10.s[0]
	fmla	v28.2s, v0.2s, v11.s[0]
.endm

.macro SAVE2x4
	fmov	alpha0, alpha

	ldr	d0, [pCRow0]
	fmla	v0.2s, v16.2s, alphaV0
	str	d0, [pCRow0]

	add	pCRow0, pCRow0, #8

	ldr	d1, [pCRow1]
	fmla	v1.2s, v20.2s, alphaV0
	str	d1, [pCRow1]

	add	pCRow1, pCRow1, #8

	ldr	d0, [pCRow2]
	fmla	v0.2s, v24.2s, alphaV0
	str	d0, [pCRow2]

	add	pCRow2, pCRow2, #8

	ldr	d1, [pCRow3]
	fmla	v1.2s, v28.2s, alphaV0
	str	d1, [pCRow3]

	add	pCRow3, pCRow3, #8
.endm

/******************************************************************************/

.macro INIT1x4
	fmov		s16, wzr
	fmov		s20, s16
.endm

.macro KERNEL1x4_SUB
	ldr	s0, [pA]
	add	pA, pA, #4

	ld1	{v8.2s, v9.2s}, [pB]
	add	pB, pB, #16

	fmla	v16.2s, v8.2s, v0.s[0]
	fmla	v20.2s, v9.2s, v0.s[0]
.endm

.macro SAVE1x4
	fmov	alpha0, alpha

	ld1	{v8.s}[0], [pCRow0]
	ld1	{v8.s}[1], [pCRow1]
	fmla	v8.2s, v16.2s, alphaV0
	st1	{v8.s}[0], [pCRow0]
	st1	{v8.s}[1], [pCRow1]

	add	pCRow0, pCRow0, #4
	add	pCRow1, pCRow1, #4

	ld1	{v12.s}[0], [pCRow2]
	ld1	{v12.s}[1], [pCRow3]
	fmla	v12.2s, v20.2s, alphaV0
	st1	{v12.s}[0], [pCRow2]
	st1	{v12.s}[1], [pCRow3]

	add	pCRow2, pCRow2, #4
	add	pCRow3, pCRow3, #4
.endm

/******************************************************************************/

.macro INIT8x2
	fmov	s16, wzr
	fmov	s17, s16
	fmov	s20, s17
	fmov	s21, s16
.endm

.macro KERNEL8x2_SUB
	ld1	{v8.2s}, [pB]
	add	pB, pB, #8
	ld1	{v0.4s}, [pA]
	add	pA, pA, #16
	ld1	{v1.4s}, [pA]
	add	pA, pA, #16

	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v17.4s, v1.4s, v8.s[0]

	fmla	v20.4s, v0.4s, v8.s[1]
	fmla	v21.4s, v1.4s, v8.s[1]
.endm

.macro SAVE8x2
	fmov	alpha0, alpha

	add	pCRow1, pCRow0, LDC

	ld1	{v0.4s, v1.4s}, [pCRow0]
	fmla	v0.4s, v16.4s, alphaV0
	fmla	v1.4s, v17.4s, alphaV0
	st1 	{v0.4s, v1.4s}, [pCRow0]

	add	pCRow2, pCRow1, LDC

	ld1	{v4.4s, v5.4s}, [pCRow1]
	fmla	v4.4s, v20.4s, alphaV0
	fmla	v5.4s, v21.4s, alphaV0
	st1 	{v4.4s, v5.4s}, [pCRow1]

	add	pCRow0, pCRow0, #32
.endm

/******************************************************************************/

.macro INIT4x2
	fmov	s16, wzr
	fmov	s17, s16
	fmov	s20, s17
	fmov	s21, s16
.endm

.macro KERNEL4x2_SUB
	ld1	{v8.2s}, [pB]
	add	pB, pB, #8
	ld1	{v0.2s, v1.2s}, [pA]
	add	pA, pA, #16

	fmla	v16.2s, v0.2s, v8.s[0]
	fmla	v17.2s, v1.2s, v8.s[0]
	fmla	v20.2s, v0.2s, v8.s[1]
	fmla	v21.2s, v1.2s, v8.s[1]
.endm

.macro SAVE4x2
	fmov	alpha0, alpha

	ld1	{v8.2s, v9.2s}, [pCRow0]
	fmla	v8.2s, v16.2s, alphaV0
	fmla	v9.2s, v17.2s, alphaV0
	st1	{v8.2s, v9.2s}, [pCRow0]

	add	pCRow1, pCRow0, LDC
	ld1	{v12.2s, v13.2s}, [pCRow1]
	fmla	v12.2s, v20.2s, alphaV0
	fmla	v13.2s, v21.2s, alphaV0
	st1	{v12.2s, v13.2s}, [pCRow1]

	add	pCRow0, pCRow0, #16
.endm

/******************************************************************************/

.macro INIT2x2
	fmov		s16, wzr
	fmov		s20, s16
.endm

.macro KERNEL2x2_SUB
	ld1	{v8.2s}, [pB]
	add	pB, pB, #8

	ld1	{v0.2s}, [pA]
	add	pA, pA, #8

	fmla	v16.2s, v0.2s, v8.s[0]
	fmla	v20.2s, v0.2s, v8.s[1]
.endm

.macro SAVE2x2
	fmov	alpha0, alpha

	ld1	{v8.2s}, [pCRow0]
	fmla	v8.2s, v16.2s, alphaV0
	st1	{v8.2s}, [pCRow0]

	add	pCRow1 , pCRow0, LDC
	ld1	{v12.2s}, [pCRow1]
	fmla	v12.2s, v20.2s, alphaV0
	st1	{v12.2s}, [pCRow1]

	add	pCRow0, pCRow0, #8
.endm

/******************************************************************************/

.macro INIT1x2
	fmov		s16, wzr
.endm

.macro KERNEL1x2_SUB
	ld1	{v8.2s} , [pB]
	add	pB , pB, #8

	ldr	s0 , [pA]
	add	pA, pA, #4

	fmla	v16.2s, v8.2s, v0.s[0]
.endm

.macro SAVE1x2
	fmov	alpha0, alpha

	add	pCRow1 , pCRow0, LDC
	ld1	{v8.s}[0], [pCRow0]
	ld1	{v8.s}[1], [pCRow1]
	fmla	v8.2s, v16.2s, alphaV0
	st1	{v8.s}[0], [pCRow0]
	st1	{v8.s}[1], [pCRow1]

	add	pCRow0, pCRow0, #4
.endm

/******************************************************************************/

.macro INIT16x1
	fmov	s16, wzr
	fmov	s17, wzr
	fmov	s18, wzr
	fmov	s19, s16
.endm

.macro KERNEL16x1_SUB
	ldr	s8, [pB]
	add	pB , pB, #4

	ld1	{v0.4s}, [pA]
	add	pA, pA, #16
	ld1	{v1.4s}, [pA]
	add	pA, pA, #16
	ld1	{v2.4s}, [pA]
	add	pA, pA, #16
	ld1	{v3.4s}, [pA]
	add	pA, pA, #16

	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v17.4s, v1.4s, v8.s[0]
	fmla	v18.4s, v2.4s, v8.s[0]
	fmla	v19.4s, v3.4s, v8.s[0]
.endm

.macro SAVE16x1
	fmov	alpha0, alpha

	ld1	{v4.4s, v5.4s, v6.4s, v7.4s}, [pCRow0]
	fmla	v4.4s, v16.4s, alphaV0
	fmla	v5.4s, v17.4s, alphaV0
	fmla	v6.4s, v18.4s, alphaV0
	fmla	v7.4s, v19.4s, alphaV0
	st1 	{v4.4s, v5.4s, v6.4s, v7.4s}, [pCRow0]

	add	pCRow0, pCRow0, #64
.endm

/******************************************************************************/

.macro INIT8x1
	fmov	s16, wzr
	fmov	s17, wzr
.endm

.macro KERNEL8x1_SUB
	ldr	s8, [pB]
	add	pB , pB, #4

	ld1	{v0.4s}, [pA]
	add	pA, pA, #16
	ld1	{v1.4s}, [pA]
	add	pA, pA, #16

	fmla	v16.4s, v0.4s, v8.s[0]
	fmla	v17.4s, v1.4s, v8.s[0]
.endm

.macro SAVE8x1
	fmov	alpha0, alpha

	ld1	{v0.4s, v1.4s}, [pCRow0]
	fmla	v0.4s, v16.4s, alphaV0
	fmla	v1.4s, v17.4s, alphaV0
	st1 	{v0.4s, v1.4s}, [pCRow0]

	add	pCRow0, pCRow0, #32
.endm

/******************************************************************************/

.macro INIT4x1
	fmov	s16, wzr
	fmov	s17, s16
.endm

.macro KERNEL4x1_SUB
	ldr	s8, [pB]
	add	pB , pB, #4

	ld1	{v0.2s, v1.2s}, [pA]
	add	pA , pA, #16

	fmla	v16.2s, v0.2s, v8.s[0]
	fmla	v17.2s, v1.2s, v8.s[0]
.endm

.macro SAVE4x1
	fmov	alpha0, alpha

	ld1	{v8.2s, v9.2s}, [pCRow0]
	fmla	v8.2s, v16.2s, alphaV0
	fmla	v9.2s, v17.2s, alphaV0
	st1	{v8.2s, v9.2s}, [pCRow0]

	add	pCRow0, pCRow0, #16
.endm

/******************************************************************************/

.macro INIT2x1
	fmov		s16, wzr
.endm

.macro KERNEL2x1_SUB
	ldr	s8, [pB]
	add	pB , pB, #4

	ld1	{v0.2s}, [pA]
	add	pA , pA, #8

	fmla	v16.2s, v0.2s, v8.s[0]
.endm

.macro SAVE2x1
	fmov	alpha0, alpha

	ld1	{v8.2s}, [pCRow0]
	fmla	v8.2s, v16.2s, alphaV0
	st1	{v8.2s}, [pCRow0]

	add	pCRow0, pCRow0, #8
.endm

/******************************************************************************/

.macro INIT1x1
	fmov	s16, wzr
.endm

.macro KERNEL1x1_SUB
	ldr	s8, [pB]
	add	pB , pB, #4

	ldr	s0, [pA]
	add	pA , pA, #4

	fmadd 	s16, s0, s8, s16  
.endm

.macro SAVE1x1
	fmov	alpha0, alpha

	ldr 	s8, [pCRow0]
	fmla	s8, s16, alphaV0
	str 	s8, [pCRow0]

	add	pCRow0, pCRow0, #4
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

.Lsgemm_kernel_begin:

	.align 5
	add	sp, sp, #-(11 * 16)
	stp	d8, d9, [sp, #(0 * 16)]
	stp	d10, d11, [sp, #(1 * 16)]
	stp	d12, d13, [sp, #(2 * 16)]
	stp	d14, d15, [sp, #(3 * 16)]
	stp	d16, d17, [sp, #(4 * 16)]
	stp	x18, x19, [sp, #(5 * 16)]
	stp	x20, x21, [sp, #(6 * 16)]
	stp	x22, x23, [sp, #(7 * 16)]
	stp	x24, x25, [sp, #(8 * 16)]
	stp	x26, x27, [sp, #(9 * 16)]
	str	x28, [sp, #(10 * 16)]

	prfm	PLDL1KEEP, [origPB]
	prfm	PLDL1KEEP, [origPA]

	fmov	alpha, s0

	lsl	LDC, LDC, #2			// ldc = ldc * 4

	mov	pB, origPB

	mov	counterJ, origN
	asr 	counterJ, counterJ, #3		// J = J / 8 
	cmp 	counterJ, #0
	ble	.Lsgemm_kernel_L4_BEGIN

/******************************************************************************/

.Lsgemm_kernel_L8_BEGIN:
	mov	pCRow0, pC
	add	pCRow1, pCRow0, LDC
	add	pCRow2, pCRow1, LDC
	add	pCRow3, pCRow2, LDC
	add	pCRow4, pCRow3, LDC
	add	pCRow5, pCRow4, LDC
	add	pCRow6, pCRow5, LDC
	add	pCRow7, pCRow6, LDC

	add	pC, pCRow7, LDC

	mov	pA, origPA			// pA = start of A array

.Lsgemm_kernel_L8_M12_BEGIN:

	mov	counterI, origM
	cmp 	counterI, #11
	ble	.Lsgemm_kernel_L8_M8_BEGIN

	.align 5
.Lsgemm_kernel_L8_M12_20:

	mov	pB, origPB

	asr 	counterL , origK, #3
	cmp	counterL , #2
	blt	.Lsgemm_kernel_L8_M12_32

	KERNEL12x8_I
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2

	subs	counterL, counterL, #2
	ble	.Lsgemm_kernel_L8_M12_22a

	.align 5
.Lsgemm_kernel_L8_M12_22:

	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M12_22

	.align 5
.Lsgemm_kernel_L8_M12_22a:

	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_E

	b	 .Lsgemm_kernel_L8_M12_44

	.align 5
.Lsgemm_kernel_L8_M12_32:

	tst	counterL, #1
	ble	.Lsgemm_kernel_L8_M12_40

	KERNEL12x8_I
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_M2
	KERNEL12x8_M1
	KERNEL12x8_E

	b	.Lsgemm_kernel_L8_M12_44

.Lsgemm_kernel_L8_M12_40:

	INIT12x8

.Lsgemm_kernel_L8_M12_44:

	ands	counterL , origK, #7
	ble	.Lsgemm_kernel_L8_M12_100

	.align 5
.Lsgemm_kernel_L8_M12_46:

	KERNEL12x8_SUB
	subs	counterL, counterL, #1
	bne	.Lsgemm_kernel_L8_M12_46

.Lsgemm_kernel_L8_M12_100:
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]

	SAVE12x8

.Lsgemm_kernel_L8_M12_END:
	subs	counterI, counterI, #12
	cmp	counterI, #11
	bgt	.Lsgemm_kernel_L8_M12_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L8_M8_BEGIN:

	cmp     counterI, #1
	blt	.Lsgemm_kernel_L8_END 

	cmp	counterI, #8
	blt	.Lsgemm_kernel_L8_M4_BEGIN

.Lsgemm_kernel_L8_M8_20:

	mov	pB, origPB

	asr 	counterL , origK, #1		// L = K / 2
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lsgemm_kernel_L8_M8_32

	KERNEL8x8_I				// do one in the K
	KERNEL8x8_M2				// do another in the K

	subs	counterL, counterL, #2
	ble	.Lsgemm_kernel_L8_M8_22a
	.align 5

.Lsgemm_kernel_L8_M8_22:

	KERNEL8x8_M1
	KERNEL8x8_M2

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M8_22

.Lsgemm_kernel_L8_M8_22a:

	KERNEL8x8_M1
	KERNEL8x8_E

	b	 .Lsgemm_kernel_L8_M8_44

.Lsgemm_kernel_L8_M8_32:

	tst	counterL, #1
	ble	.Lsgemm_kernel_L8_M8_40

	KERNEL8x8_I
	KERNEL8x8_E

	b	.Lsgemm_kernel_L8_M8_44

.Lsgemm_kernel_L8_M8_40:

	INIT8x8

.Lsgemm_kernel_L8_M8_44:

	ands	counterL , origK, #1
	ble	.Lsgemm_kernel_L8_M8_100

.Lsgemm_kernel_L8_M8_46:

	KERNEL8x8_SUB

.Lsgemm_kernel_L8_M8_100:

	SAVE8x8

.Lsgemm_kernel_L8_M8_END:

	subs	counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L8_M4_BEGIN:

        cmp	counterI ,#1
	blt	.Lsgemm_kernel_L8_END

	cmp	counterI, #4
	blt	.Lsgemm_kernel_L8_M2_BEGIN

.Lsgemm_kernel_L8_M4_20:

	mov	pB, origPB

	asr 	counterL , origK, #1		// L = K / 2
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lsgemm_kernel_L8_M4_32

	KERNEL4x8_I				// do one in the K
	KERNEL4x8_M2				// do another in the K

	subs	counterL, counterL, #2
	ble	.Lsgemm_kernel_L8_M4_22a
	.align 5

.Lsgemm_kernel_L8_M4_22:

	KERNEL4x8_M1
	KERNEL4x8_M2

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M4_22

.Lsgemm_kernel_L8_M4_22a:

	KERNEL4x8_M1
	KERNEL4x8_E

	b	 .Lsgemm_kernel_L8_M4_44

.Lsgemm_kernel_L8_M4_32:

	tst	counterL, #1
	ble	.Lsgemm_kernel_L8_M4_40

	KERNEL4x8_I
	KERNEL4x8_E

	b	.Lsgemm_kernel_L8_M4_44

.Lsgemm_kernel_L8_M4_40:

	INIT4x8

.Lsgemm_kernel_L8_M4_44:

	ands	counterL , origK, #1
	ble	.Lsgemm_kernel_L8_M4_100

.Lsgemm_kernel_L8_M4_46:

	KERNEL4x8_SUB

.Lsgemm_kernel_L8_M4_100:

	SAVE4x8

.Lsgemm_kernel_L8_M4_END:

	subs	counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L8_M2_BEGIN:

	cmp	counterI, #1
	blt	.Lsgemm_kernel_L8_END
	cmp	counterI, #2			// counterI = counterI / 2
	blt	.Lsgemm_kernel_L8_M1_BEGIN

.Lsgemm_kernel_L8_M2_20:

	INIT2x8

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L8_M2_40

.Lsgemm_kernel_L8_M2_22:

	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB

	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB
	KERNEL2x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M2_22


.Lsgemm_kernel_L8_M2_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L8_M2_100

.Lsgemm_kernel_L8_M2_42:

	KERNEL2x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M2_42

.Lsgemm_kernel_L8_M2_100:

	SAVE2x8

.Lsgemm_kernel_L8_M2_END:

	subs	counterI, counterI, #2


.Lsgemm_kernel_L8_M1_BEGIN:

	cmp	counterI, #1			// counterI = counterI % 2
	blt	.Lsgemm_kernel_L8_END

.Lsgemm_kernel_L8_M1_20:

	INIT1x8

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L8_M1_40

.Lsgemm_kernel_L8_M1_22:
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB

	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB
	KERNEL1x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M1_22


.Lsgemm_kernel_L8_M1_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L8_M1_100

.Lsgemm_kernel_L8_M1_42:

	KERNEL1x8_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L8_M1_42

.Lsgemm_kernel_L8_M1_100:

	SAVE1x8

.Lsgemm_kernel_L8_END:
	add	origPB, origPB, origK, lsl #5	// B = B + K * 4 * 8

	subs	counterJ, counterJ , #1		// j--
	bgt	.Lsgemm_kernel_L8_BEGIN


/******************************************************************************/
/******************************************************************************/

.Lsgemm_kernel_L4_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #7
	ble	.Lsgemm_kernel_L999

	tst	counterJ , #4
	ble	.Lsgemm_kernel_L2_BEGIN

	mov	pCRow0, pC
	add	pCRow1, pCRow0, LDC
	add	pCRow2, pCRow1, LDC
	add	pCRow3, pCRow2, LDC

	add	pC, pCRow3, LDC

	mov	pA, origPA			// pA = start of A array

.Lsgemm_kernel_L4_M12_BEGIN:

	mov	counterI, origM
	cmp 	counterI, #11
	ble	.Lsgemm_kernel_L4_M8_BEGIN

	.align 5
.Lsgemm_kernel_L4_M12_20:

	mov	pB, origPB

	asr 	counterL , origK, #3
	cmp	counterL , #2
	blt	.Lsgemm_kernel_L4_M12_32

	KERNEL12x4_I
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2

	subs	counterL, counterL, #2
	ble	.Lsgemm_kernel_L4_M12_22a

	.align 5
.Lsgemm_kernel_L4_M12_22:

	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M12_22

	.align 5
.Lsgemm_kernel_L4_M12_22a:

	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_E

	b	 .Lsgemm_kernel_L4_M12_44

	.align 5
.Lsgemm_kernel_L4_M12_32:

	tst	counterL, #1
	ble	.Lsgemm_kernel_L4_M12_40

	KERNEL12x4_I
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_M2
	KERNEL12x4_M1
	KERNEL12x4_E

	b	.Lsgemm_kernel_L4_M12_44

.Lsgemm_kernel_L4_M12_40:

	INIT12x4

.Lsgemm_kernel_L4_M12_44:

	ands	counterL , origK, #7
	ble	.Lsgemm_kernel_L4_M12_100

	.align 5
.Lsgemm_kernel_L4_M12_46:

	KERNEL12x4_SUB
	subs	counterL, counterL, #1
	bne	.Lsgemm_kernel_L4_M12_46

.Lsgemm_kernel_L4_M12_100:
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]

	SAVE12x4

.Lsgemm_kernel_L4_M12_END:
	subs	counterI, counterI, #12
	cmp	counterI, #11
	bgt	.Lsgemm_kernel_L4_M12_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L4_M8_BEGIN:

        cmp	counterI, #1
	blt	.Lsgemm_kernel_L4_END

	cmp	counterI, #8
	blt	.Lsgemm_kernel_L4_M4_BEGIN

.Lsgemm_kernel_L4_M8_20:

	mov	pB, origPB

	asr 	counterL , origK, #1		// L = K / 2
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lsgemm_kernel_L4_M8_32

	KERNEL8x4_I				// do one in the K
	KERNEL8x4_M2				// do another in the K

	subs	counterL, counterL, #2
	ble	.Lsgemm_kernel_L4_M8_22a
	.align 5

.Lsgemm_kernel_L4_M8_22:

	KERNEL8x4_M1
	KERNEL8x4_M2

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M8_22

.Lsgemm_kernel_L4_M8_22a:

	KERNEL8x4_M1
	KERNEL8x4_E

	b	 .Lsgemm_kernel_L4_M8_44

.Lsgemm_kernel_L4_M8_32:

	tst	counterL, #1
	ble	.Lsgemm_kernel_L4_M8_40

	KERNEL8x4_I
	KERNEL8x4_E

	b	.Lsgemm_kernel_L4_M8_44

.Lsgemm_kernel_L4_M8_40:

	INIT8x4

.Lsgemm_kernel_L4_M8_44:

	ands	counterL , origK, #1
	ble	.Lsgemm_kernel_L4_M8_100

.Lsgemm_kernel_L4_M8_46:

	KERNEL8x4_SUB

.Lsgemm_kernel_L4_M8_100:

	SAVE8x4

.Lsgemm_kernel_L4_M8_END:

	subs	counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L4_M4_BEGIN:

	cmp	counterI , #1
	blt	.Lsgemm_kernel_L4_END

	cmp	counterI, #4
	blt	.Lsgemm_kernel_L4_M2_BEGIN

.Lsgemm_kernel_L4_M4_20:

	mov	pB, origPB

	asr 	counterL , origK, #1		// L = K / 2
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Lsgemm_kernel_L4_M4_32

	KERNEL4x4_I				// do one in the K
	KERNEL4x4_M2				// do another in the K

	subs	counterL, counterL, #2
	ble	.Lsgemm_kernel_L4_M4_22a
	.align 5

.Lsgemm_kernel_L4_M4_22:

	KERNEL4x4_M1
	KERNEL4x4_M2

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M4_22

.Lsgemm_kernel_L4_M4_22a:

	KERNEL4x4_M1
	KERNEL4x4_E

	b	 .Lsgemm_kernel_L4_M4_44

.Lsgemm_kernel_L4_M4_32:

	tst	counterL, #1
	ble	.Lsgemm_kernel_L4_M4_40

	KERNEL4x4_I
	KERNEL4x4_E

	b	.Lsgemm_kernel_L4_M4_44

.Lsgemm_kernel_L4_M4_40:

	INIT4x4

.Lsgemm_kernel_L4_M4_44:

	ands	counterL , origK, #1
	ble	.Lsgemm_kernel_L4_M4_100

.Lsgemm_kernel_L4_M4_46:

	KERNEL4x4_SUB

.Lsgemm_kernel_L4_M4_100:

	SAVE4x4

.Lsgemm_kernel_L4_M4_END:

	subs	counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L4_M2_BEGIN:

	cmp	counterI , #1
	blt	.Lsgemm_kernel_L4_END

	cmp	counterI, #2			// counterI = counterI / 2
	blt	.Lsgemm_kernel_L4_M1_BEGIN

.Lsgemm_kernel_L4_M2_20:

	INIT2x4

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L4_M2_40

.Lsgemm_kernel_L4_M2_22:

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M2_22


.Lsgemm_kernel_L4_M2_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L4_M2_100

.Lsgemm_kernel_L4_M2_42:

	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M2_42

.Lsgemm_kernel_L4_M2_100:

	SAVE2x4

.Lsgemm_kernel_L4_M2_END:

	subs	counterI, counterI, #2


.Lsgemm_kernel_L4_M1_BEGIN:

 	cmp	counterI, #1			// counterI = counterI % 2
	blt	.Lsgemm_kernel_L4_END

.Lsgemm_kernel_L4_M1_20:

	INIT1x4

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L4_M1_40

.Lsgemm_kernel_L4_M1_22:
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M1_22


.Lsgemm_kernel_L4_M1_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L4_M1_100

.Lsgemm_kernel_L4_M1_42:

	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L4_M1_42

.Lsgemm_kernel_L4_M1_100:

	SAVE1x4

.Lsgemm_kernel_L4_END:
	add	origPB, origPB, origK, lsl #4	// B = B + K * 4 * 4


/******************************************************************************/

.Lsgemm_kernel_L2_BEGIN:   // less than 2 left in N direction

	mov	counterJ , origN
	tst	counterJ , #3
	ble	.Lsgemm_kernel_L999

	tst	counterJ , #2
	ble	.Lsgemm_kernel_L1_BEGIN

	mov	pCRow0, pC			// pCRow0 = pC
        add     pCRow1, pCRow0, LDC
	add	pC,pC,LDC, lsl #1


	mov	pA, origPA			// pA = A


.Lsgemm_kernel_L2_M12_BEGIN:

	mov	counterI, origM

	cmp	counterI, #12
	blt	.Lsgemm_kernel_L2_M8_BEGIN

.Lsgemm_kernel_L2_M12_20:

	INIT12x2

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lsgemm_kernel_L2_M12_40
	.align 5

.Lsgemm_kernel_L2_M12_22:
	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB

	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB
	KERNEL12x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M12_22


.Lsgemm_kernel_L2_M12_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L2_M12_100

.Lsgemm_kernel_L2_M12_42:

	KERNEL12x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M12_42

.Lsgemm_kernel_L2_M12_100:

	SAVE12x2

.Lsgemm_kernel_L2_M12_END:

	subs	counterI, counterI, #12

	cmp	counterI, #11
	bgt	.Lsgemm_kernel_L2_M12_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L2_M8_BEGIN:

        cmp     counterI, #0
	ble	.Lsgemm_kernel_L2_END

        cmp     counterI, #7
	ble	.Lsgemm_kernel_L2_M4_BEGIN

.Lsgemm_kernel_L2_M8_20:

	INIT8x2

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lsgemm_kernel_L2_M8_40
	.align 5

.Lsgemm_kernel_L2_M8_22:
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB

	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB
	KERNEL8x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M8_22


.Lsgemm_kernel_L2_M8_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L2_M8_100

.Lsgemm_kernel_L2_M8_42:

	KERNEL8x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M8_42

.Lsgemm_kernel_L2_M8_100:

	SAVE8x2

.Lsgemm_kernel_L2_M8_END:

	subs	counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L2_M4_BEGIN:

        cmp	counterI, #1
	blt	.Lsgemm_kernel_L2_END

	cmp	counterI, #4
	blt	.Lsgemm_kernel_L2_M2_BEGIN

.Lsgemm_kernel_L2_M4_20:

	INIT4x2

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL,#0
	ble	.Lsgemm_kernel_L2_M4_40
	.align 5

.Lsgemm_kernel_L2_M4_22:
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M4_22


.Lsgemm_kernel_L2_M4_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L2_M4_100

.Lsgemm_kernel_L2_M4_42:

	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M4_42

.Lsgemm_kernel_L2_M4_100:

	SAVE4x2

.Lsgemm_kernel_L2_M4_END:

	subs	counterI, counterI, #4

//------------------------------------------------------------------------------


.Lsgemm_kernel_L2_M2_BEGIN:

        cmp	counterI , #1
	blt	.Lsgemm_kernel_L2_END

	cmp	counterI, #2			// counterI = counterI / 2
	blt	.Lsgemm_kernel_L2_M1_BEGIN

.Lsgemm_kernel_L2_M2_20:

	INIT2x2

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
        cmp	counterL,#0
	ble	.Lsgemm_kernel_L2_M2_40

.Lsgemm_kernel_L2_M2_22:

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M2_22


.Lsgemm_kernel_L2_M2_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L2_M2_100

.Lsgemm_kernel_L2_M2_42:

	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M2_42

.Lsgemm_kernel_L2_M2_100:

	SAVE2x2

.Lsgemm_kernel_L2_M2_END:

	subs	counterI, counterI, #2


.Lsgemm_kernel_L2_M1_BEGIN:

	cmp	counterI, #1			// counterI = counterI % 2
	blt	.Lsgemm_kernel_L2_END

.Lsgemm_kernel_L2_M1_20:

	INIT1x2

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
        cmp     counterL, #0
	ble	.Lsgemm_kernel_L2_M1_40

.Lsgemm_kernel_L2_M1_22:
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M1_22


.Lsgemm_kernel_L2_M1_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L2_M1_100

.Lsgemm_kernel_L2_M1_42:

	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L2_M1_42

.Lsgemm_kernel_L2_M1_100:

	SAVE1x2

.Lsgemm_kernel_L2_END:
	add	origPB, origPB, origK, lsl #3	// B = B + K * 2 * 4

/******************************************************************************/

.Lsgemm_kernel_L1_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #1
	ble	.Lsgemm_kernel_L999 // done


	mov	pCRow0, pC			// pCRow0 = C
	add	pC , pC , LDC			// Update pC to point to next

	mov	pA, origPA			// pA = A

.Lsgemm_kernel_L1_M12_BEGIN:

	mov	counterI, origM

	cmp	counterI, #12
	blt	.Lsgemm_kernel_L1_M8_BEGIN

.Lsgemm_kernel_L1_M12_20:

	INIT12x1

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L1_M12_40
	.align 5

.Lsgemm_kernel_L1_M12_22:
	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB

	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB
	KERNEL12x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M12_22


.Lsgemm_kernel_L1_M12_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L1_M12_100

.Lsgemm_kernel_L1_M12_42:

	KERNEL12x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M12_42

.Lsgemm_kernel_L1_M12_100:

	SAVE12x1

.Lsgemm_kernel_L1_M12_END:

	subs	counterI, counterI, #12

	cmp	counterI, #11
	bgt	.Lsgemm_kernel_L1_M12_20

//------------------------------------------------------------------------------

.Lsgemm_kernel_L1_M8_BEGIN:

        cmp	counterI, #1
	blt	.Lsgemm_kernel_L1_END

	cmp	counterI, #8
	blt	.Lsgemm_kernel_L1_M4_BEGIN

.Lsgemm_kernel_L1_M8_20:

	INIT8x1

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L1_M8_40
	.align 5

.Lsgemm_kernel_L1_M8_22:
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB

	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB
	KERNEL8x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M8_22


.Lsgemm_kernel_L1_M8_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L1_M8_100

.Lsgemm_kernel_L1_M8_42:

	KERNEL8x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M8_42

.Lsgemm_kernel_L1_M8_100:

	SAVE8x1

.Lsgemm_kernel_L1_M8_END:

	subs	counterI, counterI, #8

//------------------------------------------------------------------------------

.Lsgemm_kernel_L1_M4_BEGIN:

        cmp	counterI, #1
	blt	.Lsgemm_kernel_L1_END

	cmp	counterI, #4
	blt	.Lsgemm_kernel_L1_M2_BEGIN

.Lsgemm_kernel_L1_M4_20:

	INIT4x1

	mov	pB, origPB

	asr	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L1_M4_40
	.align 5

.Lsgemm_kernel_L1_M4_22:
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M4_22


.Lsgemm_kernel_L1_M4_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L1_M4_100

.Lsgemm_kernel_L1_M4_42:

	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M4_42

.Lsgemm_kernel_L1_M4_100:

	SAVE4x1

.Lsgemm_kernel_L1_M4_END:

	subs	counterI, counterI, #4

//------------------------------------------------------------------------------

.Lsgemm_kernel_L1_M2_BEGIN:

	cmp	counterI, #1
	blt	.Lsgemm_kernel_L1_END

        cmp	counterI, #2			// counterI = counterI / 2
	blt	.Lsgemm_kernel_L1_M1_BEGIN

.Lsgemm_kernel_L1_M2_20:

	INIT2x1

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L1_M2_40

.Lsgemm_kernel_L1_M2_22:

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M2_22


.Lsgemm_kernel_L1_M2_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L1_M2_100

.Lsgemm_kernel_L1_M2_42:

	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M2_42

.Lsgemm_kernel_L1_M2_100:

	SAVE2x1

.Lsgemm_kernel_L1_M2_END:

	subs	counterI, counterI, #2


.Lsgemm_kernel_L1_M1_BEGIN:

	cmp	counterI, #1			// counterI = counterI % 2
	blt	.Lsgemm_kernel_L1_END

.Lsgemm_kernel_L1_M1_20:

	INIT1x1

	mov	pB, origPB

	asr 	counterL , origK, #3		// counterL = counterL / 8
	cmp	counterL , #0
	ble	.Lsgemm_kernel_L1_M1_40

.Lsgemm_kernel_L1_M1_22:
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M1_22


.Lsgemm_kernel_L1_M1_40:

	ands	counterL , origK, #7		// counterL = counterL % 8
	ble	.Lsgemm_kernel_L1_M1_100

.Lsgemm_kernel_L1_M1_42:

	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Lsgemm_kernel_L1_M1_42

.Lsgemm_kernel_L1_M1_100:

	SAVE1x1

.Lsgemm_kernel_L1_END:

.Lsgemm_kernel_L999:
	mov	x0, #0				// set return value
	ldp	d8, d9, [sp, #(0 * 16)]
	ldp	d10, d11, [sp, #(1 * 16)]
	ldp	d12, d13, [sp, #(2 * 16)]
	ldp	d14, d15, [sp, #(3 * 16)]
	ldp	d16, d17, [sp, #(4 * 16)]
	ldp	x18, x19, [sp, #(5 * 16)]
	ldp	x20, x21, [sp, #(6 * 16)]
	ldp	x22, x23, [sp, #(7 * 16)]
	ldp	x24, x25, [sp, #(8 * 16)]
	ldp	x26, x27, [sp, #(9 * 16)]
	ldr	x28, [sp, #(10 * 16)]
	add	sp, sp, #(11*16)
	ret

	EPILOGUE

