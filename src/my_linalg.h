// ���� my_linalg.h :
// ����������� ���������� ��������� 
// ���������� �������� �������.
// ������ v0.01 ���� 2011
// ������������� ������� ���� ������� �������� �������.
// 1. (��������) ������� ����.
// 2. (��������� ��������) �������� ���������� �� � ��.
// ------------------------------------------------------
// ��������� � ������� ����.
// ���������������: 
// a) ������� � ����������� �������.
// b) ������������ � ��������������.
// ��� ����������� �������:
// ������� �������������� ������ ��������������. (������������ ������������ ������������).
// ��� ����� ������������� ������, � ������� ������� � �����
// ���������� ���������� �������� ����������� �������������������.
// ����������� ��������� ��������� �������� �������� ����������� ������ 
// � ������� �������������� ����� ��������� (������ �������������� ������������)
// ��� �������������� (� ������ ����� ��������� ����� ��������� �������� ����������� 
// ������ ) ������������ ������� �������� �������  ����� � ���������� ��������, ���
// ����������� �� ������, �� ���������� �� ��������������.

#pragma once
#ifndef MY_LINALG_H
#define MY_LINALG_H 1

#include "sparse_gauss.c" // ����� ������ ��� ����������� �������
// ��� �������
// calculateSPARSEgaussArray - ������ ����������� ���� ������� ���������� ������.
// initIMatrix - ��������� ����������� ������.
// setValueIMatrix, addValueIMatrix - ��������� � ���������� �������� � ����������� �������.

#define Real double

// ��� ����������� �������� ����������� �������:
// ��������� ������� ������� ����:
typedef struct tagNONZEROELEM {
	Real aij;
	int key; // ����� �������
    struct tagNONZEROELEM* next; // ��������� �� ��������� ������� � ������
} NONZEROELEM; 

// ��� ��������� �������� �������
typedef struct tagSIMPLESPARSE {
	// ������ ����������
   NONZEROELEM** root; // ������ ��������� ��������� ������ ������
   int n; // ����� ��������� ���������
   //int POOL_SIZE; // ������ ������� ��������� ���������:
   // n - ��������, POOL_SIZE - ������ �� ������.
   //int incCLUSTER_SIZE; // ��� ��������� ������
} SIMPLESPARSE;

// ���� ������ ������� ����
typedef struct Tequation {
	Real ap, ae, an, aw, as, b;
	int iP, iE, iN, iW, iS;
} equation;

// ������ ���� ������� ���������� ������
// ��� ������ �������� �������� � ��� 
// ����� ������������� �������.
// �������� ��� �������������� ������� A nodesxnodes.
void eqsolve_simple_gauss(Real **A, int nodes, Real *b, Real *x);

// ������ ���� ������� ���������� ����������
// ������� ����������� ������ ���� ������������.
// ��� ����� ������������� �������. ��������������
// ������� ����� ������� � �� ��������� ��� �����.
// �� �������������� � ������ � ��� ���� ����������� 
// ������ ���������� ������.
// ��. ������ ���������� � �������: ��� ��� ��������������
// � ��������� ����������. 1986 ���.
void eqsolv_simple_holesskii(Real **A, int nodes, Real *b, Real *x);

// ������� �������� ������� �������
// ���������� ������. ���������� � ������������ ���� 
// �������� ������ ���� ���.
// ������� A �������� � � �� ������������ �������� �������.
void inverse_matrix_simple(Real **A, int nodes, bool flag);

// ������� ������������ ���� ����������
// ������ nodesxnodes C=A*B. ��������� 
// ������������ � ������� B.
void multiply_matrix_simple(Real **A, Real **B, int nodes);

// ��������� ��������� ������� ������������ ��� ��������������� ��� �������
// ������ �������� ����������� ��������:

// 1. ��������� ���������� ������ ������� nxn:
//              t=m*p.
// �� ��������� ������ ��������� �������� � ������� t. 
void multi_m(Real **m, Real **p, Real **t, int n);

// 2. ���������������� ���������� ������� m ��������
// nxn. �� ��������� ������ ��������� ����������������
// �������� � m.
void tr_m(Real **m, int n);

// 3. ���������� ������������ ���������������
// ������� ��� ������������ ������� A �������� 
// nxn. ������� ������������� �������� A[f][g].
Real max_el(Real **A, int n, int& f, int& g);

// 4. �������� ������ ������� � ������: A=B.
// ������� ���������� �������� nxn
void matr_copy(Real **A, Real **B, int n);

// 5. ������� ��������� ���� ���������� ������ ������������ 
// ���� ������� nxn (����� ���������):
//                A=hi*A.
// ����� hi - �������������� ������� ��������:
// hi[f][f] = cosfi;
// hi[g][g] = cosfi;
// hi[f][g] = +sinfi;
// hi[g][f] = -sinfi;
// ����� f � g ������� ��������� ���������.
// ��������� ���������� � ����.
// �� ��������� ������ ��������� �������� � �������� ������� A.
// ������ ������� hi� ��������� ������ ��� ���� ������ ������ ��������
// (cosfi � sinfi), ��� ��������� ����������� ��������� ������ � ��������������.
// rab - ������� ������ ����������� 2xn. �� ������������ ��� �������� ���������.
void multi_m_left(Real **A, Real **rab, int n, int f, int g, Real cosfi, Real sinfi);

// 6. ������� ��������� ���� ���������� ������ ������������ 
// ���� ������� nxn (������ ���������):
//                A=A*hi.
// ����� hi - �������������� ������� ��������:
// hi[f][f] = cosfi;
// hi[g][g] = cosfi;
// hi[f][g] = -sinfi;
// hi[g][f] = +sinfi;
// ����� f � g ������� ��������� ���������.
// ��������� ���������� � ����.
// �� ��������� ������ ��������� �������� � �������� ������� A.
// ������ ������� hi ��������� ������ ��� ���� ������ ������ ��������
// ��� ��������� ����������� ��������� ������ � ��������������.
// rab - ������� ������ ����������� 2xn. �� ������������ ��� �������� ���������.
void multi_m_right(Real **A, Real **rab, int n, int f, int g, Real cosfi, Real sinfi); 
   
// ������� �������� ������ ������ �������� �� � ����
//      A-lambda_scal*E=0
// ������� ���������� �������� � �� �������� ������������,
// ��� �������� ��������������� ��������� epsilon.
void jacobi_matrix_simple(Real **A, Real **U, Real *lambda, int nodes, Real epsilon);

// ����������� ����������.
void BubbleSortGSEP1(Real *a, int *mask, int n);

// ������ ���������� ������������ �������� ����������� ��������
//   GSEP1: 
void GSEP1(Real **A, Real **B, Real **U, Real *lambda, int *mask, int nodes, Real epsilon);

// ����� ������ ��� ��������� ������� A ��������
//              nodes x 2*icolx+1, ���
//   2*icolx+1 - ������ �����. ��� ��� ��� �������
//  A ��������� ���������� �� ��� ��������� ��������
//  ������� ���������� ������ ������ �����.
//  b - ������ ������ ����� ����, x - ������ �������.
//  ��������� ��������� ���������� � ����.
//  ��� ������������ ����������� �������� ��������������
//  ������ �, ������� �������� ����� ������.
//  ����� ���� ������� 1777-1855.
//  � ���������� ������ ������� � ��������.
//  ��� ����������� ��������� �� ��������� � ���
//  ������ ����� ����� ���� ������ ��� nodes 
//  ����� � ��� ���� ������� ������������� ��������� ���������� 
//  � ���� ������ ������������ (����������), ���� �� ��������� 
//  ����������� �������� ���������� ������ �����.
void eqsolve_lenta_gauss(Real **A, int nodes, int icolx, Real *b, Real *x);

// ����� (�����) ������-������� ����������-������� SOR
// ��� ������� ���� � ������� �������� � nxn
// �������� ��������������, �� � ������������ 
// �������������. ������� � ��������������
// ��������� ���������� (�������������).
// b - ������ �����, x - ���������� �������, 
// eps - �������� ����������� �������.
// omega - ����������� ����������� �������� ����������.
void Seidel(Real **A, Real *b, Real *x, int n, Real eps, Real omega);

// ���������� ������������ �� ���� 
// ������������ �����.
//Real fmax(Real fA, Real fB);

// ����������� ��� ��������� �������� �������� 
// � ������ ����� �� ���� ������� ����� ������� �������.
void SOR(equation* &sl, Real* &x, int n);

// ����� ���������� ����������
// ��� ����� ������������� ������� ����.

// ��������� ������� �� ������
Real* MatrixByVector(Real** H,Real* V,int n);

// ��������� ����� �������
Real NormaV(double *V, int n);

// ��������� ������������ ���� ��������
Real Scal(Real *v1, Real *v2, int n);

//----------����� ����������� ����������---------------
/* ������� ���������:
*  A - ������������� ������� ����,
*  dV - ������ ������ �����, 
*  x - ��������� ����������� � ������� ��� NULL.
*  n - ����������� ���� Anxn.
*  ������� A ���������� ������������ ����������� � 
*  ������������ (������������ ������������ ������������).
*  ���������� �������� ���������� 1000, �.�. ��������������,
*  ��� ���� ������� �� ������� �� 1000 �������� �� ��� � �� �������.
*  �������� ������ �� ������� ������� � ���������� ���������:
*  dterminatedTResudual.
*/
Real *SoprGrad(Real **A, Real *dV, Real *x, int n);

// ����������� ������������ ������������ ����������� ������� 
// � CRS ������� ��������:

// ��������� ������� �� ������
// ������������ ������ �������� CRS
// ����������� ������� A (val, col_ind, row_ptr) ���������� �������� nxn.
// ����� ��������� ����� ����� ����������� � ����� n.
void MatrixCRSByVector(Real* val, int* col_ind, int* row_ptr, Real* V, Real* &tmp, int n);

// ��������� ����������������� ������� �� ������
// (������������, ��������, � ������ BiCG - ������������ ����������)
// ��� �������� (�� ����������������� �������) ������������ ������ �������� CRS
// ����������� ������� A (val, col_ind, row_ptr) ���������� �������� nxn.
// ����� ��������� ����� ����� ����������� � ����� n.
Real* MatrixTransposeCRSByVector(Real* val, int* col_ind, int* row_ptr, Real* V, int n);


/* ������� ���������:
*  val, col_ind, row_ptr - ����������� ������� ���� � ������� CRS,
*  dV - ������ ������ �����, 
*  x - ��������� ����������� � ������� ��� NULL.
*  n - ����������� ���� Anxn.
*  ����������� ������� A (val, col_ind, row_ptr) ���������� �������� nxn.
*  ����� ��������� ����� ����� ����������� � ����� n.
*  ������� A ���������� ������������ ����������� � 
*  ������������ (������������ ������������ ������������).
*  ���������� �������� ���������� 1000, �.�. ��������������,
*  ��� ���� ������� �� ������� �� 1000 �������� �� ��� � �� �������.
*  �������� ������ �� ������� ������� � ���������� ���������:
*  dterminatedTResudual.
*/
Real *SoprGradCRS(Real *val, int* col_ind, int* row_ptr, Real *dV, Real *x, int n);

// ����� ������������ ����������
// ��� �������� �������������� ������� � (val, col_ind, row_ptr).
// ����������������� �� ������ ��������, ������ : "������
// ������� ���� ������� �����������".
// dV - ������ ����� ����,
// x - ��������� ����������� � ������� ��� NULL.
// n - ����������� � nxn.
// ���������� �������� ���������� 2000.
// �������� ������ �� ������� ������� � ���������� ���������:
//  dterminatedTResudual.
Real *BiSoprGradCRS(Real *val, int* col_ind, int* row_ptr, Real *dV, Real *x, int n, int maxit);

// ������ ��� �� ����������� ���������������� ������� L.
// ������������ ������������ ����������� �������
// ���� A ������������ �������� ����������� ��������� 
// A~=L*transpose(L); L - ������ ����������� �������.
// L - �������� � ��������� ����:
// 1. ldiag - ������������ �������� L.
// 2. lltr - ��������������� �������� � �������� �������,
// �.�. �������� ����������. 
// 3. jptr - �������������� ������ �������� ��� lltr, 
// 4. iptr - ���������� � ������ ��������� ������ ��� lltr.
// f - ������ ������ ����� �������� nodes.
// ���������� ������ z=inverse(L)*f;
// ������ f ��������.
// ������ (CSIR - ������):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* inverseL(Real* f, Real* ldiag, Real* lltr, int* jptr, int* iptr, int n);

// ������ ��� �� ����������� ���������������� ������� L.
// ������������ ������������ ����������� �������
// ���� A ������������ �������� ����������� ��������� 
// A~=L*transpose(L); L - ������ ����������� �������.
// L - �������� � ��������� ����:
// 1. val - ������������ � ��������������� �������� L.
// � ���������� �������. 
// 3. indx - �������������� ������ ����� ��� val, 
// 4. pntr - ���������� � ������ ���������� �������.
// f - ������ ������ ����� �������� nodes.
// ���������� ������ z=inverse(L)*f;
// ������ f ��������.
// ������ (CSIR - ������):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseL_ITL(Real* f, Real* val, int* indx, int* pntr, Real* &z, int n);

// �������� ��� �� ����������� ����������������� ������� U.
// ������������ ������������ ����������� �������
// ���� A ������������ �������� ����������� ��������� 
// A~=L*transpose(L); L - ������ ����������� �������.
// U=transpose(L);
// U - �������� � ��������� ����:
// 1. udiag - ������������ �������� U.
// 2. uutr - ��������������� �������� � ���������� �������,
// �.�. �������� ������������. 
// ��� ������� �����������, ��:
// 3. jptr - �������������� ������ �������� ��� lltr, 
// 4. iptr - ���������� � ������ ��������� ������ ��� lltr.
// f - ������ ������ ����� �������� nodes.
// ���������� ������ z=inverse(U)*f;
// ������ f ��������.
// ������ (CSIR - ������):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0   
// ------------------------------------------
// udiag==ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// uutr==lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* inverseU(Real* f, Real* udiag, Real* uutr, int* jptr, int* iptr, int n);

// �������� ��� �� ����������� ����������������� ������� U.
// ������������ ������������ ����������� �������
// ���� A ������������ �������� ����������� ��������� 
// A~=L*transpose(L); L - ������ ����������� �������.
// U=transpose(L); - ������� ����������� �������.
// U - �������� � ��������� ����:
// 1. val - ������������ � ��������������� �������� U (� ��������� �������).
// 2. indx - �������������� ������ ��������, 
// 3. pntr - ���������� � ������ ��������� ������ ��� val.
// f - ������ ������ ����� �������� nodes.
// ���������� ������ z=inverse(U)*f;
// ������ f ��������.
// ������ (CSIR_ITL - ������):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseU_ITL(Real* f, Real* val, int* indx, int* pntr, Real* &z, int n);

// ��������� �� ������� CSIR � ������ CSIR_ITL
// �������:
// CSIR: ldiag, lltr, jptr, iptr
// CSIR_ITL: val, indx, pntr
// ������:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// ������ CSIR:
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
//��������� ����������� ������ CSIR_ITL
//val : 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
//indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
//pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
void convertCSIRtoCSIR_ITL(Real *ldiag, Real *lltr, int *jptr, int *iptr, int n, int nz, Real* &val, int* &indx, int* &pntr, int nnz);

// �������� ���������� ���������
// ��� ������������ ����������� ������������
// ������� � �������� nxn.
// n - ����������� ������� ����
// ������� val ���������� � � ��� ������������
// �������� ���������� ��������� IC(0).
// ������:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
//������ CSIR_ITL (������� ����������� �������� ���������).
// val : 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
// pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
// ��������� ������������ ��� ����������:
// ���������� ������ val (indx � pntr �������� ��� ���������):
// val (factorization)= 
// 3.0
// 1.0
// 0.3333333333333333
// 0.3333333333333333
// 3.3166247903554
// 0.6030226891555273
// 0.30151134457776363
// 0.6030226891555273
// 3.1622776601683795
// 0.6324555320336759
// 2.932575659723036
// 0.34099716973523675
// 3.4472773213410837
// 0.2578524458667825
// 2.8284271247461903
// 2.7310738989293286
//-------------------------------------------
void IC0Factor_ITL(Real* val, int* indx, int* pntr, int n);

// ���������������� �������� ���������� ���������.
// (���������� ������� IC0Factor_ITL).
void IC0FactorModify_ITL(Real* val, int* indx, int* pntr, int n);

// ��������� �� ������� CSIR_ITL � ������ CSIR (�������� ��������������)
// ������ ��� ��� ������� �������������� ���������� �������!!!
// �������:
// CSIR_ITL: val, indx, pntr
// CSIR: ldiag, lltr, jptr, iptr
// ������:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
//��������� ����������� ������ CSIR_ITL
//val : 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
//indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
//pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
// ������ CSIR:
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void convertCSIR_ITLtoCSIR(Real* ldiag, Real* lltr, int* jptr, int* iptr, int n, int nz, Real* val, int* indx, int* pntr, int nnz);

// �������� ���������� ��������� IC(0).
// ������� ������ ������ ����������� ������������ ������� � ������� CSIR.
// ������ ��������� ���� �������������� � ������� CSIR_ITL ���������� �������� ITL.
void ICFactor0(Real* ldiag, Real* lltr, int* jptr, int* iptr, int n, int nz);

// ��������� ������������ ������������ �����������  ������� �� ������ 
// ������������ ������ �������� CSIR. � ���� ��������� �������� ������ ��������������� �������� altr. 
// ����������� SPD ������� A (adiag, altr, jptr, iptr) ���������� �������� nxn.
// ����� ��������� ����� ����� ����������� � ����� n.
// ������:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// ������ CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void SPDMatrixCSIRByVector(Real* adiag, Real* altr, int* jptr, int* iptr, Real* V, Real* &tmp, int n);

// ��������� �������������� ������������ �����������  ������� �� ������ 
// ������������ ������ �������� CSIR.  
// ����������� ������� A (adiag, altr, autr, jptr, iptr) ���������� �������� nxn.
// ����� ��������� ����� ����� ����������� � ����� n.
// ��������� adiag �������� ��������. ������ ����������� altr �������� ���������.
// ������� ����������� �������� �� �������� autr. ������� ������� (������� ��������� 
// ��������� ) �������������� ������������. ������ jptr - ������ �������� ��� ������� 
// ������������, ������ iptr - ���������� ��� ���������� ����� ������ ��� ������� ������������.
// ������:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   1.0   10.0   2.0   0.0   0.0   0.0    
// 2.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 2.0   2.0   0.0   0.0   3.0   0.0   8.0 
// ------------------------------------------
// ������ CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 1.0  2.0 1.0 2.0  1.0 1.0  2.0 2.0 3.0
// autr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* MatrixCSIRByVector(Real* adiag, Real* altr, Real* autr, int* jptr, int* iptr, Real* V, int n);

// ��������� ����������������� �������������� ������������ �����������  ������� �� ������ 
// ������������ ������ �������� CSIR.  
// ����������� ������� A (adiag, altr, autr, jptr, iptr) ���������� �������� nxn. �������� 
// ������ �������� �������, � ���������� � ����������������� �������.
// ����� ��������� ����� ����� ����������� � ����� n.
// ��������� adiag �������� ��������. ������ ����������� altr �������� ���������.
// ������� ����������� �������� �� �������� autr. ������� ������� (������� ��������� 
// ��������� ) �������������� ������������. ������ jptr - ������ �������� ��� ������� 
// ������������, ������ iptr - ���������� ��� ���������� ����� ������ ��� ������� ������������.
// ������:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   1.0   10.0   2.0   0.0   0.0   0.0    
// 2.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 2.0   2.0   0.0   0.0   3.0   0.0   8.0 
// ------------------------------------------
// ������ CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 1.0  2.0 1.0 2.0  1.0 1.0  2.0 2.0 3.0
// autr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* MatrixTransposeCSIRByVector(Real* adiag, Real* altr, Real* autr, int* jptr, int* iptr, Real* V, int n);


/* ����� ���������� ���������� �������� � ������� [1952]
*  ������� ���������:
*  adiag, altr, jptr, iptr - ����������� ������� ���� � ������� CSIR,
*  dV - ������ ������ �����, 
*  x - ��������� ����������� � ������� ��� NULL.
*  n - ����������� ���� Anxn.
*  nz - ����������� �������� altr, jptr.
*  ����������� ������� A (adiag, altr, jptr, iptr) ���������� �������� nxn.
*  ����� ��������� ����� ����� ����������� � ����� n.
*  ������� A ���������� ������������ ����������� � 
*  ������������ (������������ ������������ ������������).
*  �������� ������ ������ ����������� � ���������� altr � adiag.
*  ���������� �������� ���������� 1000, �.�. ��������������,
*  ��� ���� ������� �� ������� �� 1000 �������� �� ��� � �� �������.
*  �������� ������ �� ������� ������� � ���������� ���������:
*  dterminatedTResudual.
*  � �������� ������������������� �������� �������� ���������� ���������:
*  M^(-1)==transpose(L)^(-1)*L^(-1); // ���������� �������������������.
*/
Real *SoprGradCSIR(Real* adiag, Real* altr, int* jptr, int* iptr, Real *dV, Real *x, int n, int nz0);

// ������� ���������� ���� ������������� ������� ���� �.
// ������� ���� � ������� � CSIR ������� : adiag, altr, jptr, iptr.
// �������� ���������� ��������� ��� � ������������ � ���������� � ����:
// A = L*transpose(L); � ������� �����������. ������� jptr �  iptr �������� ���� ��.
// ����� ������� : A~=inverse(L)*A*inverse(transpose(L)) ���� ����������� � ������������ ����������.
// ������ ����� ��������������� ������� ����� ���: dV~=inverse(L)*dV.
// ������� ���� ����� ����� A~*x~=dV~; => x~=transpose(L)*x; => x=inverse(transpose(L))*x~;
// ������������������ �������� ������������ ��������� ��������� ���������� �������� ��� ������� ����,
// �������� ������������ �������������� ������� ����.
Real *SoprGradCSIR2(Real* adiag, Real* altr, int* jptr, int* iptr, Real *dV, Real *x, int n, int nz0);

/* ����� ���������� ���������� �������� � ������� [1952]
*  ������� ���������:
*  M - ����������� ������� ���� � ������� SIMPLESPARSE,
*  dV - ������ ������ �����, 
*  x - ��������� ����������� � ������� ��� NULL.
*  n - ����������� ���� Anxn.
*
*  ����������� ������� M ���������� �������� nxn.
*  ����� ��������� ����� ����� ����������� � ����� n.
*  ������� M �������������� ������������ ����������� � 
*  ������������ (������������ ������������ ������������).
*  �������� ������ ��������� ��������. 
*  ���������� �������� ���������� 1000, �.�. ��������������,
*  ��� ���� ������� �� ������� �� 1000 �������� �� ��� � �� �������.
*  �������� ������ �� ������� ������� � ���������� ���������:
*  dterminatedTResudual.
*  � �������� ������������������� �������� �������� ���������� ���������:
*  K^(-1)==transpose(L)^(-1)*L^(-1); // ���������� �������������������.
*  ������� ���������� (�������� ��� ������� ������).
*/
void ICCG(SIMPLESPARSE &M, Real *dV, Real* &x, int n);

// �������� �.�. ����������� [1993]
// ��� �������� �������������� ������.
// ���������������� �� ����������
// "��������� ������ ������� ������ ���������" [2004]
// �������������� ���������������� ������������ ������������.
Real* SoloveichikAlgCSIR_SPD(int isize, // ������ ���������� �������
						Real* adiag, Real* altr, int* jptr, int* iptr, // ������� ����
                         Real *dV,  // ������ ������ �����
                         const Real *dX0, // ������ ���������� �����������
                         bool bconsole_message); // �������� �� �������� ������� �� ������� ?

// �������� �.�. ����������� [1993]
// ��� �������� �������������� ������.
// ����� ������������ ��� ������������ � ������������ ����������� �������.
// � ������������������� �������� ����������� ���������.
// ���������������� �� ����������
// "��������� ������ ������� ������ ���������" [2004]
// �������������� ���������������� ������������ ������������.
Real* SoloveichikAlgCSIR_SPDgood(int isize, int nz0,// ������ ���������� �������
						Real* adiag, Real* altr, int* jptr, int* iptr, // ������� ����
                         Real *dV,  // ������ ������ �����
                         const Real *dX0, // ������ ���������� �����������
                         bool bconsole_message); // �������� �� �������� ������� �� ������� ?

// �������� �.�. ����������� [1993]
// ��� �������� �������������� ������.
// ���������������� �� ����������
// "��������� ������ ������� ������ ���������" [2004]
// �������������� ���������������� ������������ ������������.
void SoloveichikAlgCRS(int isize, // ������ ���������� �������
						Real *val, int* col_ind, int* row_ptr, // ������� ����
                         Real *dV,  // ������ ������ �����
                         const Real* &dX0, // ������ ���������� �����������
                         bool bconsole_message, int maxit); // �������� �� �������� ������� �� ������� ?

// �������������� ����������� �������
void initsimplesparse(SIMPLESPARSE &M, int nodes);

// ��������� ��������� ������� �
// ���������� ����������� ������� M
void addelmsimplesparse(SIMPLESPARSE &M, Real aij, int i, int j, bool bset);

// ������������ ������ ��� ������� SIMPLESPARSE
void simplesparsefree(SIMPLESPARSE &M);


// ����������� ���������� ������ �������� ����������� �������
// � ������ CRS. ����� nodes - ���������.
void simplesparsetoCRS(SIMPLESPARSE &M, Real* &val, int* &col_ind, int* &row_ptr, int nodes);

// ���������� �� ������� ������.
// ����������� ���������� ������ �������� ����������� �������
// � ������ CSIR. ����� nodes - ���������.
// ��� �������� ������ ��� SPD ������.
// ������������ ������������ ����������� ������,
// �������� ������ ������ �����������.
void simplesparsetoCSIR(SIMPLESPARSE &M, Real* &adiag, Real* &altr, int* &jptr, int* &iptr, int nodes);

// ������ ������� � �������
void printM_and_CSIR(SIMPLESPARSE &sparseM, int  n);

// ���������� �� ������� ������.
// ����������� ���������� ������ �������� ����������� �������
// � ������ CSIR_ITL. ����� nodes - ���������.
// ��� �������� ������ ��� SPD ������.
// ������������ ������������ ����������� ������,
// �������� ������ ������� �����������.
// ������ ���������� ������ ������.
void simplesparsetoCSIR_ITLSPD(SIMPLESPARSE &M, Real* &val, int* &indx, int* &pntr, int nodes);

/* �������� LU ���������� ��� �������������� ������
*  ������ � nxn=
*    9.0 0.0 0.0 3.0 1.0 0.0 1.0
*    0.0 11.0 2.0 1.0 0.0 0.0 2.0 
*    0.0 1.0 10.0 2.0 0.0 0.0 0.0 
*    2.0 1.0 2.0 9.0 1.0 0.0 0.0 
*    1.0 0.0 0.0 1.0 12.0 0.0 1.0 
*    0.0 0.0 0.0 0.0 0.0 8.0 0.0
*    2.0  2.0 0.0 0.0 3.0 0.0 8.0
*-----------------------------------------
*  ������������� (� ���� ���� ������ ��������� �� ���� ���������):
*  ������� ����������� ������� �������� ���������, � ������ ������
*  �������� ������������� �� �������� ������� ��������.
*  U_val :   1.0, 1.0, 3.0, 9.0,   2.0, 1.0, 2.0, 11.0,   2.0, 10.0, 1.0, 9.0, 1.0,12.0, 8.0, 8.0
*  U_ind :   6, 4, 3, 0,  6, 3, 2, 1,  3,2, 4,3, 6,4, 5, 6
*  U_ptr :   0, 4, 8, 10, 12, 14, 15, 16
*  ������ ����������� ������� �������� �����������, � ������ �������
*  �������� ������������� �� �������� ������� �����.
*  L_val :  2.0, 1.0, 2.0, 9.0,    2.0, 1.0, 1.0, 11.0,  2.0, 10.0, 1.0, 9.0,  3.0, 12.0, 8.0, 8.0
*  L_ind :  6, 4, 3, 0,  6, 3, 2, 1,   3, 2,  4,3,  6, 4, 5, 6
*  L_ptr :  0, 4, 8, 10, 12, 14, 15, 16
*----------------------------------------------
* ��������� ILU ����������:
* U_val : 1.0, 1.0, 3.0, 9.0, 2.0, 1.0, 2.0, 11.0, 2.0, 10.0, 1.0, 9.0, 1.0, 12.0, 8.0, 8.0.
* L_val : 0.222, 0.111, 0.222, 1.0, -1.273, 0.091, 0.091, 1.0, 0.2, 1.0, 0.111, 1.0, -0.417, 1.0, 1.0, 1.0.
*/
void ILU0_Decomp_ITL(Real* &U_val, int* &U_ind, int* &U_ptr, Real* &L_val, int* &L_ind, int* &L_ptr, int n);

/* ����� ������������ ����������
* ��� �������� �������������� ������� � (val, col_ind, row_ptr).
* ����������������� �� ������ ��������, ������ : "������
* ������� ���� ������� �����������".
* dV - ������ ����� ����,
* x - ��������� ����������� � ������� ��� NULL.
* n - ����������� � nxn.
* ���������� �������� ���������� 2000.
* �������� ������ �� ������� ������� � ���������� ���������:
*  dterminatedTResudual.
* ������ ����� ����������. ���� ������� ������ ������ r_tilda, �� 
* ������� ����� ����� ����������. ����������� �� ����� ������� r_tilda:
* ������� ����� ��������� ������������ Scal(r,r_tilda,n) != 0.0.
*/
Real *BiSoprGrad(IMatrix *xO, SIMPLESPARSE &M,  Real *dV, Real *x, int n);

// �������� �.�. ����������� [1993]
// ��� �������� �������������� ������.
// ���������������� �� ����������
// "��������� ������ ������� ������ ���������" [2004]
// �������������� ���������������� ������������ ������������.
// �������� ILU �������������������.
void SoloveichikAlg( IMatrix *xO, SIMPLESPARSE &M,// ����������� ������� ����
                         Real *dV,  // ������ ������ �����
                         Real* &dX0, // ������ ���������� �����������
                         bool bconsole_message, // �������� �� �������� ������� �� ������� ?
                         int imaxiter); // ����������� ���������� ���-�� ��������

// ����� ��� ��� ������ Bi-CGStab
Real  *Bi_CGStab(int n, Real *val, int* col_ind, int* row_ptr, Real *dV, Real *dX0, int maxit);

#endif