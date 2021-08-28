// ALICEFLOW0_02.cpp :
// 2D �������� ��� ��������� ����� ������ (���������� � ������������ ����������). 
// ���������� ��������� SIMPLE 1972 ���� �� ����������� ������.
// ������ ���-��� 1983 ���.
// ��������� ��������� ��� ������� �����:
// 1. meshin.txt
// 2. solver.txt
// ����� meshin.txt � solver.txt ���������������� � ������� � ��������� DavisTest Delphi
//  ���������� �� ����� ������� � ����� Embarcadero Delphi.
// 
// �������������� ������ ����� Upwind (������ ������ ��� ���������).
// ��� �������� �������� ������������ �������� ICCG, ��� �������� ����� �������.
// ����� ��������� ������� ������������� � ������������ ������ ������������� �������.
// ������������ ����������������� ��������� �����.
// 
// ������������ � ������� 2021 (14-16 ������� 2021). 
// 28 ������ 2011 �������� � �������������� �������.
// 
// ����� �� 17 ������� 2019: ���������� �� SPARSEM � ������ ell �������. ������ �������� ������ ���� ��� ��� ell.
// ����� �������� ������ � ell �������� �� 3D �������. ������� �������� �� ��� ����� � ����������� � �������������� ICCG � SEIDEL
// �� ��� �����, ����� ��������� �� 3D ����.
// ������� ����� ��������� ��������� �������� 8� 610��.
// 
// ����� �� 18 ������� 2019: ������� ������� � AliceFlow0_02(AliceFlow2D) �� AliceFlow_0_59(AliceFlow3D). �������������� AliceFlow0_02.
// 
// ���������� ������������ ICCG ������� ���� �� ��������� ��� �������� ap / 0.9999; �� �������. 24.08.2021. ��������� amg1r5 ������ ��� ��������.
// ����� ICCG ���������� ������ ������������ ���� ��������� ap / 0.9999 ������� � amg ����� �����������.
//
// ����� ����������� �������� ���������:
// A. ����������� ������������ �� OpenGL ������ �� ����� �������.
// B. ���������� ����������� ��� ������� ������� ��������������. ���������� � ������� �� � �������.
// C. ����� ������ - ��������� ��������� - ����������� � ������������ ���������� (���� ����� �����).
// D. ����� SMARTER ��� ���������.
// E. ��������� meshgen2 ������ � AliceFlow.v.0.02 2D �� QUAD Pave ��������� �����.
// F. ����������������, dynamic mesh, VOF �����.
// �� ������� ����� ��������� ����� �������, C/C++ Visual Studio, Intel Parallel Studio. ���������� ��������������. �������� ������� �������.
// �.�. ��� ���� ���� ����� ��������� ��������� ������ DavisTest, � ��� ���� � ����������� �. F �� �������� ����� ���������. �� �� ��������.
// �� ��������� �� �������������� Visual Studio.



#include "stdafx.h"
#include <iostream> // ������� cout, cin
#include <stdlib.h> // ��� ������� exit, atoi, atof
#include <string.h> // strcat, strstr
#include <math.h> // �������������� �������
#include "windows.h" // ��� ������� WinExec
#include <time.h> // ��� �������� ����������

#define Real double // float
Real Cp_active = 1.0;
Real Pe_max = 0.0;
Real Re_max = 0.0;
Real glRCh = 0.05;

#include <algorithm>

#pragma once
#ifndef ALICEFLOW0_02_CPP
#define ALICEFLOW0_02_CPP 1

int ilim = 0;

int inumcore = 2; // ����� ���� ����������

// true - �������� ���������� ���������� � �������.
const bool DEBUG = false;
//#define SOLOVEICHIC_SPEED 1



typedef struct T_PAIR_OMP
{
	int s, e;
} PAIR_OMP;

PAIR_OMP s_par[10] = { {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1} };

const int MAX_STRING_LENGTH_ELL_THERMAL = 7; // 27 ����������.
//const int MAX_STRING_LENGTH_ELL_MECHANICAL = 81; // 81 ���������.
int MAX_STRING_LENGTH_ELL = MAX_STRING_LENGTH_ELL_THERMAL;

Real** data_ell = nullptr;
int** coll_ell = nullptr;
int nnz_ell = 0;

#include "my_meshin_v0_02.c" // ���������� �������� ��������� �����

// ����������� ������� ��� ���������� ��������� ���������-�������� �� ����������� �����
#include "pamendment2.c" 

// ��� ��������� � my_elmatr_quad_f.c
//#include "my_linalg.c" // ���������� ������� �������� �������
// ��� �������: 
// eqsolve_simple_gauss - ������ ���� ������� ���������� ������
// eqsolv_simple_holesskii - ������ ���� ������� ���������� ����������

// ������� ���������� � ��������� tecplot360
#include "my_export_tecplot2.c"



// ������������ ����� �������� ��� ������ ���� �� �������,
// ����� ����� �� �������, 
// ����� ��������� ���������� ��� ����,
// ����� ����� ��� ������������� ���������-��������.
int MAXIT, MAXTS, iSOLVER, iSHEME;



// false - ������������
// true - ��������������.
bool BTIMEDEP;

Real TAU; // ��� �� �������

int maxnod; // ������������ ����� ���� (����������� �������)
int maxelm; // ����������� ���������� ����� ��������� (��������� ������ ��� �����) (����� ����������� �������)
// �� ��������� �������� ������������. ��� �������� ������������ ��� ���������� ����� � ������.
int nve=4; // ����� ������� ���������� �������� (��������� �������� ������ ��� ����������).

int **nvtx=nullptr; // ������ ����� ��� ������� �������� (������������ ������)
int **sosed = nullptr; // �������� ����������� ������ ��� ������� ������������ ������
Real *x = nullptr, *y = nullptr; // ������� ����������
// ��� ����������� ������������
int **nvtxcell = nullptr; // ����� ��� ����������� �������.
int ncell; // ���������� ������ ��� ����������� �������.

bool **boundary = nullptr; // ������ ��� ��������� ����������� ������� (���. ����������).
bool **neiman = nullptr; // ������ ���� �� �������� ����� ������� �������.  (���. ����������).
int **norm = nullptr; // ������� � ������� ��������� ������� (-1) - ���� ���������� ����.
Real **potent = nullptr; // ������ ������� ����������� (������� �������)
Real **prop = nullptr; // �������� ����������
Real dgx=0.0, dgy=0.0; // ��������� ���������� �������
Real temp_ref=0.0; // ������� �������� ����������� � ����������� ����������.


Real** sumanb = nullptr;
Real** tau = nullptr; // �����������.
Real** Flux_gran = nullptr, ** Flux_gran_relx = nullptr;



// ��������� ������ ����������:
Real *alpha = nullptr;


equation **slau = nullptr; // ������������ ������� ����


// ���������  ������ ��� ������� �������������
void my_free2(equation**& slauparam, int maxvalparam, int maxelmparam, Real*& alphaloc) {
	
	for (int i = 0; i < maxvalparam; ++i) {
		switch (i) {
		case Temp: delete[] slauparam[Temp]; slauparam[Temp] = nullptr; break;
		case Vx: delete[] slauparam[Vx]; slauparam[Vx] = nullptr; break;
		case Vy: delete[] slauparam[Vy]; slauparam[Vy] = nullptr;  break;
		case PAm: delete[] slauparam[PAm]; slauparam[PAm] = nullptr; break;
		default: slauparam[i] = nullptr; break;
		}
	}
	delete[] alphaloc;
	delete[] slauparam;

} // my_free2

// ���������  ������ ��� ������� �������������
void my_malloc2(equation** &slauparam, int maxvalparam, int maxelmparam, Real* &alphaloc) {
	slauparam = new equation*[maxvalparam];
	for (int i=0; i<maxvalparam; ++i) {
		switch (i) {
			case Temp : slauparam[Temp] = new equation[maxelmparam]; break;
			case Vx : slauparam[Vx] = new equation[maxelmparam]; break;
			case Vy : slauparam[Vy] = new equation[maxelmparam]; break;
			case PAm : slauparam[PAm] = new equation[maxelmparam]; break;
			default : slauparam[i] = nullptr; break;
		}
	}
	alphaloc = new Real[MAXVAL];
	

} // my_malloc2

void my_init() {
	alpha[Temp] = 1.0;
	alpha[Vx] = 0.7;
	alpha[Vy] = 0.7;
	alpha[Press] = 0.3;
	alpha[PAm] = 0.3;

	for (int iP = 0; iP < maxelm; iP++) {
		Flux_gran_relx[E][iP] = 0.0;
		Flux_gran_relx[W][iP] = 0.0;
		Flux_gran_relx[N][iP] = 0.0;
		Flux_gran_relx[S][iP] = 0.0;

		Flux_gran[E][iP] = 0.0;
		Flux_gran[W][iP] = 0.0;
		Flux_gran[N][iP] = 0.0;
		Flux_gran[S][iP] = 0.0;

		potent[Vxcor][iP] = 0.0;
		potent[Vycor][iP] = 0.0;
		potent[Press][iP] = 0.0;
    }
} // my_init

void marena() {
   for (int iP=0; iP<maxelm; iP++) {
	   Real xc, yc;
	   Real Pi=3.141;
	   xc=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
	   yc=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
	   potent[Vx][iP]=2*Pi*sin(Pi*xc)*sin(Pi*xc)*sin(Pi*yc)*cos(Pi*yc);
	   potent[Vy][iP]=-2*Pi*sin(Pi*xc)*sin(Pi*yc)*sin(Pi*yc)*cos(Pi*xc);
   }
} //marena



// ������ ���� ���������
// �������� ��������� ����������������.
// res - ������������ �������.
void solve(int iVar, Real& res, Real***& B) {
	ilim = 0;
	if (iVar == PAm) {
		Pe_max = 0.0;
		Re_max = 0.0;
	}
	Real alphaP = 0.3;

	// ������ ������� ����
	Real res1 = 0.0;

	Real avgx = 0.0; 
	Real avgy = 0.0;
	switch (iVar) {
	case PAm:

#pragma omp parallel for
		for (int i_1 = 0; i_1 < maxelm; ++i_1) {
			int tid = omp_get_thread_num();

			calc_tau(i_1, slau, nvtx, boundary, potent, sumanb, Flux_gran, Flux_gran_relx, tau, x, y, prop, sosed, neiman, norm, nve, alpha, B[tid]);//�����������.
		}
		if (0) {
			// ��� ������� �������� ������� �� � ���� ������ �� ��������� �����������.
			for (int i_1 = 0; i_1 < maxelm; ++i_1) {
				avgx += tau[E][i_1] + tau[W][i_1];
				avgy += tau[N][i_1] + tau[S][i_1];
			}
			avgx = avgx / (2.0 * maxelm);
			avgy = avgy / (2.0 * maxelm);
			for (int i_1 = 0; i_1 < maxelm; ++i_1) {
				// ������� �������������
				tau[E][i_1] = 2.0 * tau[E][i_1] * avgx / (tau[E][i_1] + avgx);
				tau[W][i_1] = 2.0 * tau[W][i_1] * avgx / (tau[W][i_1] + avgx);

				tau[N][i_1] = 2.0 * tau[N][i_1] * avgy / (tau[N][i_1] + avgy);
				tau[S][i_1] = 2.0 * tau[S][i_1] * avgy / (tau[S][i_1] + avgy);
			}
		}
#pragma omp parallel for
		for (int i_1 = 0; i_1 < maxelm; ++i_1) {
			int tid = omp_get_thread_num();

			my_elmatr_quad_PAm(i_1, slau, nvtx, boundary, potent, sumanb, Flux_gran, Flux_gran_relx, tau, x, y, prop, sosed, neiman, norm, nve, alpha, B[tid]); // �������� ��������
		}

		if (0) {
			// ��������� ����� ���� ��� ������������.

			// ������������ ������������������.
			for (int i = 0; i < maxelm; ++i) {
				slau[PAm][i].ae /= slau[PAm][i].ap;
				slau[PAm][i].aw /= slau[PAm][i].ap;
				slau[PAm][i].an /= slau[PAm][i].ap;
				slau[PAm][i].as /= slau[PAm][i].ap;
				slau[PAm][i].b /= slau[PAm][i].ap;
				slau[PAm][i].ap = 1.0;
			}
		}

#pragma omp parallel for reduction(+: res1)
		for (int i_1 = 0; i_1 < maxelm; ++i_1) {
			if (!((boundary[PAm][i_1]) && (!neiman[PAm][i_1]))) {
				//res=fmax(res,fabs(slau[PAm][i_1].b));
				res1 += fabs(slau[PAm][i_1].b);
			}
		}
		res = (Real)(res1 / maxelm);

		// ������������� ����
		// ��� ICCG ������� ��������� ������� �������
		for (int i = 0; i < maxelm; ++i) if ((boundary[PAm][i]) && (!neiman[PAm][i])) {
			// ������ ��� �������� �������� ��������
#pragma omp parallel for
			for (int j = 0; j < maxelm; ++j) {
				if (j != i) {
					if (slau[PAm][j].iE == i) {
						slau[PAm][j].b -= slau[PAm][j].ae * slau[PAm][i].b / slau[PAm][i].ap;
						slau[PAm][j].ae = 0.0;
						slau[PAm][j].iE = -1;
					}
					if (slau[PAm][j].iW == i) {
						slau[PAm][j].b -= slau[PAm][j].aw * slau[PAm][i].b / slau[PAm][i].ap;
						slau[PAm][j].aw = 0.0;
						slau[PAm][j].iW = -1;
					}
					if (slau[PAm][j].iN == i) {
						slau[PAm][j].b -= slau[PAm][j].an * slau[PAm][i].b / slau[PAm][i].ap;
						slau[PAm][j].an = 0.0;
						slau[PAm][j].iN = -1;
					}
					if (slau[PAm][j].iS == i) {
						slau[PAm][j].b -= slau[PAm][j].as * slau[PAm][i].b / slau[PAm][i].ap;
						slau[PAm][j].as = 0.0;
						slau[PAm][j].iS = -1;
					}
				}
			}
		}

#pragma omp parallel for
		for (int i_1 = 0; i_1 < maxelm; ++i_1) {
			if (!((boundary[PAm][i_1]) /* && (!neiman[PAm][i_1])*/)) {
				// ��� ��� ���������� ICCG �������� ��� ������� !!!
				// ���� ���� slau[PAm][i_1].ap / 0.999999; �� ��� ������� �����.
				slau[PAm][i_1].ap = slau[PAm][i_1].ap / 0.9999;// 999;//0.99; 0.9999
			}
		}

		break;
	default:
		if ((iVar == Vx) || (iVar == Vy) || (iVar == Temp)) {
#pragma omp parallel for
			for (int i_1 = 0; i_1 < maxelm; ++i_1) {
				int tid = omp_get_thread_num();

				my_elmatr_quad_F(i_1, slau[iVar], iVar, BTIMEDEP, TAU, iSHEME, nvtx, boundary, potent, sumanb, Flux_gran, x, y, prop, sosed, neiman, norm, nve, alpha, dgx, dgy, temp_ref, B[tid]); // ������������

			}
		}
		res = 0.0;
		break;
	}


	/*
	for (int i=0; i<maxelm; ++i) {
		printf("%.3f %.3f %.3f %.3f %.3f = %.3f \n", slau[iVar][i].ap, slau[iVar][i].ae, slau[iVar][i].an, slau[iVar][i].as, slau[iVar][i].aw, slau[iVar][i].b);
	}
	printf("\n");
	system("pause");
	//*/

	Real* rthdsd; // ������ ����� ������� ���������
	rthdsd = new Real[maxelm];


	SIMPLESPARSE sparseM; // ����������� �������
	IMatrix sparseS;

	

    
    // ����������� ������� � ������� CSIR
    //Real *adiag, *altr;
    //int *jptr, *iptr;

    //system("pause");
	if (iVar == PAm) {
		// ��������� ������ � ������������� ��� 
		// ���������� ����������� �������.
		//initsimplesparse(sparseM, maxelm);
		//initIMatrix(&sparseS, maxelm);
	}
	else {

		if (iVar == Temp) {
			// ��������� ������ � ������������� ��� 
			// ���������� ����������� �������.
			//initsimplesparse(sparseM, maxelm);
			//initIMatrix(&sparseS, maxelm);
		}
		else {

#ifdef SOLOVEICHIC_SPEED
			// ��������� ������ � ������������� ��� 
			// ���������� ����������� �������.
			initsimplesparse(sparseM, maxelm);
			initIMatrix(&sparseS, maxelm);
#endif
		}
	}
	
#pragma omp parallel for
	for (int i=0; i<maxelm; ++i) {
		switch (iVar) {
			case Vx : 
#ifdef SOLOVEICHIC_SPEED
				addelmsimplesparse(sparseM, slau[iVar][i].ap/alpha[iVar], slau[iVar][i].iP, slau[iVar][i].iP, true);
				setValueIMatrix(&sparseS,slau[iVar][i].iP, slau[iVar][i].iP, slau[iVar][i].ap/alpha[iVar]);
#endif
					  if (!(boundary[Vx][i])) {
						  rthdsd[slau[iVar][i].iP] = slau[iVar][i].b + (1 - alpha[iVar]) * slau[iVar][i].ap * potent[Vxcor][slau[iVar][i].iP] / alpha[iVar];
						  slau[iVar][i].ap /= alpha[iVar];
					  }
					  else {
						  rthdsd[slau[iVar][i].iP] = slau[iVar][i].b;
					  }
					  break;
			case Vy : 
#ifdef SOLOVEICHIC_SPEED
				addelmsimplesparse(sparseM, slau[iVar][i].ap/alpha[iVar], slau[iVar][i].iP, slau[iVar][i].iP, true);
                setValueIMatrix(&sparseS,slau[iVar][i].iP, slau[iVar][i].iP,slau[iVar][i].ap/alpha[iVar]);
#endif
					  if (!(boundary[Vy][i])) {
						  rthdsd[slau[iVar][i].iP] = slau[iVar][i].b + (1 - alpha[iVar]) * slau[iVar][i].ap * potent[Vycor][slau[iVar][i].iP] / alpha[iVar];
						  slau[iVar][i].ap /= alpha[iVar];
					  }
					  else {
						  rthdsd[slau[iVar][i].iP] = slau[iVar][i].b;
					  }
					  break;
			case Temp :

				//addelmsimplesparse(sparseM, slau[iVar][i].ap, slau[iVar][i].iP, slau[iVar][i].iP, true);
				//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iP, slau[iVar][i].ap);

				//rthdsd[slau[iVar][i].iP] = slau[iVar][i].b;

				

				if (!(boundary[Temp][i])) {
					//addelmsimplesparse(sparseM, slau[iVar][i].ap / alpha[iVar], slau[iVar][i].iP, slau[iVar][i].iP, true);
					//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iP, slau[iVar][i].ap / alpha[iVar]);

					rthdsd[slau[iVar][i].iP] = slau[iVar][i].b + (1 - alpha[iVar]) * slau[iVar][i].ap * potent[Temp][slau[iVar][i].iP] / alpha[iVar];
					slau[iVar][i].ap /= alpha[iVar];
				}
				else {
					//addelmsimplesparse(sparseM, slau[iVar][i].ap, slau[iVar][i].iP, slau[iVar][i].iP, true);
					//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iP, slau[iVar][i].ap);

					rthdsd[slau[iVar][i].iP] = slau[iVar][i].b;
				}

				break;
			default : 
				//alphaP = 0.9999;
				//slau[iVar][i].b+= (1 - alphaP) * slau[iVar][i].ap * potent[PAm][slau[iVar][i].iP] / alphaP;
				//slau[iVar][i].ap /= alphaP;

				rthdsd[slau[iVar][i].iP] = slau[iVar][i].b;

				     addelmsimplesparse_Stress_ell(slau[iVar][i].ap, slau[iVar][i].iP, slau[iVar][i].iP, true,true);
                     // setValueIMatrix(&sparseS,slau[iVar][i].iP, slau[iVar][i].iP,slau[iVar][i].ap);
				      
				      break;
		}
		
		Real nonzeroEPS=1e-37; // ��� ��������� ������������� ����

		if (iVar == PAm) {
			//nonzeroEPS = -1;

			if ((slau[iVar][i].iE > -1) && ((slau[iVar][i].ae) > nonzeroEPS)) {
				addelmsimplesparse_Stress_ell(-slau[iVar][i].ae, slau[iVar][i].iP, slau[iVar][i].iE, true,true);
				//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iE, -slau[iVar][i].ae);
			}
			if ((slau[iVar][i].iN > -1) && ((slau[iVar][i].an) > nonzeroEPS)) {
				addelmsimplesparse_Stress_ell( -slau[iVar][i].an, slau[iVar][i].iP, slau[iVar][i].iN, true,true);
				//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iN, -slau[iVar][i].an);
			}
			if ((slau[iVar][i].iS > -1) && ((slau[iVar][i].as) > nonzeroEPS)) {
				addelmsimplesparse_Stress_ell( -slau[iVar][i].as, slau[iVar][i].iP, slau[iVar][i].iS, true,true);
				//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iS, -slau[iVar][i].as);
			}
			if ((slau[iVar][i].iW > -1) && ((slau[iVar][i].aw) > nonzeroEPS)) {
				addelmsimplesparse_Stress_ell( -slau[iVar][i].aw, slau[iVar][i].iP, slau[iVar][i].iW, true, true);
				//setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iW, -slau[iVar][i].aw);
			}
		}
		else {

			if (iVar == Temp) {
				/*
				* // �� �������� � ������������� ������.
				if ((slau[iVar][i].iE > -1) && (fabs(slau[iVar][i].ae) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].ae, slau[iVar][i].iP, slau[iVar][i].iE, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iE, -slau[iVar][i].ae);
				}
				if ((slau[iVar][i].iN > -1) && (fabs(slau[iVar][i].an) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].an, slau[iVar][i].iP, slau[iVar][i].iN, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iN, -slau[iVar][i].an);
				}
				if ((slau[iVar][i].iS > -1) && (fabs(slau[iVar][i].as) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].as, slau[iVar][i].iP, slau[iVar][i].iS, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iS, -slau[iVar][i].as);
				}
				if ((slau[iVar][i].iW > -1) && (fabs(slau[iVar][i].aw) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].aw, slau[iVar][i].iP, slau[iVar][i].iW, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iW, -slau[iVar][i].aw);
				}*/
			}
			else {

#ifdef SOLOVEICHIC_SPEED
				if ((slau[iVar][i].iE > -1) && (fabs(slau[iVar][i].ae) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].ae, slau[iVar][i].iP, slau[iVar][i].iE, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iE, -slau[iVar][i].ae);
				}
				if ((slau[iVar][i].iN > -1) && (fabs(slau[iVar][i].an) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].an, slau[iVar][i].iP, slau[iVar][i].iN, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iN, -slau[iVar][i].an);
				}
				if ((slau[iVar][i].iS > -1) && (fabs(slau[iVar][i].as) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].as, slau[iVar][i].iP, slau[iVar][i].iS, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iS, -slau[iVar][i].as);
				}
				if ((slau[iVar][i].iW > -1) && (fabs(slau[iVar][i].aw) > nonzeroEPS)) {
					addelmsimplesparse(sparseM, -slau[iVar][i].aw, slau[iVar][i].iP, slau[iVar][i].iW, true);
					setValueIMatrix(&sparseS, slau[iVar][i].iP, slau[iVar][i].iW, -slau[iVar][i].aw);
				}
#endif
			}
		}

		
		
	}

	if (iVar == PAm) {
		patch_Ell(maxelm);
	}

	/*if (iVar==PAm) {
		printM_and_CSIR(sparseM, maxelm); 
	    system("pause");
	    printf("\n");
	}// */ // debug 


    // ���� ������� ����:

    //1. ������ ������ ������� ����:

    // 1.1. ��� ������� ������.

	// ������ ���� ������� ���������� ����������
	// ��� ������������ ������������ ����������� ������� s.
	// ��� ����� ������������� �������. � ��� ���� ����� 
	// ���������� ��� ����� ������� ���������� ������.
	//eqsolv_simple_holesskii(s, nodes, rthdsd, potent);
	// ������ ���� ������� ���������� ������
	// ��� ������ �������� �������� � ��� 
	// ����� ������������� �������.
	//eqsolve_simple_gauss(s, nodes, rthdsd, potent[Temp]);

    // 1.2. ��� ����������� ������.

	// ����� �.�. ������ ��� ����������� �������� �������������� �������
    //calculateSPARSEgaussArray(&sparseS, potent, rthdsd);


	//2. ������������ ������ ������� ����:

	// 2.1. ��� ������� ������.

	// ����� ������-�������-����������-������� SOR
	// ��� ������������� ������� ����.
	// ����� ��������������� ��������� �����������.
	//Seidel(s, rthdsd, potent, nodes, 1e-5, 1.855);
	// ������ ������������� ���� �������
	// ���������� ����������.
	// ��� ������������ ������������ ����������� ������� s.
	// ����� �������� ����������, �.�. ���� �������� ������
	// ���������� ����������� �� ��� ����� ��������� ������ NULL.
	//potent=SoprGrad(s, rthdsd, NULL, nodes);
		
    // 2.2. ��� ����������� ������.

	// 2.2.1 ��� ��� CRS (CSIR, CSIR_ITL) ������� ������������ ������.
    ///*
    //simplesparsetoCRS(sparseM, val, col_ind, row_ptr, nodes);
	//printM_and_CRS(sparseM,val,col_ind,row_ptr,nodes);
    //potent=SoprGradCRS(val, col_ind, row_ptr, rthdsd, NULL, nodes);//*/

	// ��� (Congruate Gradients) ��� CSIR ������� ��������
	// ������������ ������������ ����������� ������.
    //simplesparsetoCSIR(sparseM, adiag, altr, jptr, iptr, nodes);
	//printM_and_CSIR(sparseM, adiag, altr, jptr, iptr,  nodes);
	//int inz=(int)((sparseM.n-nodes)/2.0);
	//potent=SoprGradCSIR(adiag, altr, jptr, iptr, rthdsd, NULL, nodes, inz);
	//potent=SoloveichikAlgCSIR_SPD(nodes, adiag, altr, jptr, iptr, rthdsd, NULL, true);
	//potent=SoloveichikAlgCSIR_SPDgood(nodes, inz, adiag, altr, jptr, iptr, rthdsd, NULL, true);
	if (iVar==PAm) {



		 for (int i=0; i<maxelm; ++i) potent[PAm][i]=0.0;
		 if (0) {
			 // ������� �������������������, ������� ���������.
			 SOR(slau[PAm], potent[PAm], rthdsd, maxelm, PAm);

		 }
		 else {

			// SOR(slau[PAm], potent[PAm], maxelm);
			 ICCG(sparseM, rthdsd, potent[PAm], maxelm);
		 }
	}
	//*/

	// 2.2.2 ������ ��� CRS (CSIR, CSIR_ITL) ������� �������� �������������� ������.
	if (iVar != PAm) {
		// ����������� ������� � ������� CRS
		Real* val = nullptr;
		int* col_ind = nullptr, * row_ptr = nullptr;

		if (iVar == Temp) {
			//simplesparsetoCRS(sparseM, val, col_ind, row_ptr, maxelm);
		}
		else {

#ifdef SOLOVEICHIC_SPEED
			simplesparsetoCRS(sparseM, val, col_ind, row_ptr, maxelm);
#endif
		}
		int maxiter = MAXIT;
		//system("pause");

		//if (iVar==PAm) maxiter=2000;
		switch (iSOLVER) {
		case BICGCRS: potent[iVar] = BiSoprGradCRS(val, col_ind, row_ptr, rthdsd, potent[iVar], maxelm, maxiter); break;
		case SOLOVALGCRS:

			if (iVar == Temp) {
				//SoloveichikAlg(&sparseS, sparseM, rthdsd, potent[iVar], true, maxiter);
				SOR(slau[iVar], potent[iVar], rthdsd, maxelm, Temp);
			}
			else {
#ifdef SOLOVEICHIC_SPEED
				//SoloveichikAlgCRS(maxelm, val, col_ind, row_ptr, rthdsd, potent[iVar], true, maxiter); 
				  //system("pause");
				SoloveichikAlg(&sparseS, sparseM, rthdsd, potent[iVar], true, maxiter);
				//system("pause"); system("pause"); system("pause");
			   //system("pause");
#else
				if (1 || (fabs(dgx * dgy) < 1.0e-20)) {
					// ��� ��������� �������� ���������� ���� ����������� ������.

					SOR(slau[iVar], potent[iVar], rthdsd, maxelm, iVar);
				}
#endif
			}
			break;
		case BICGSTABCRS: potent[iVar] = Bi_CGStab(maxelm, val, col_ind, row_ptr, rthdsd, potent[iVar], maxiter); break;
		default: potent[iVar] = Bi_CGStab(maxelm, val, col_ind, row_ptr, rthdsd, potent[iVar], maxiter); break;
		}

		if (iVar == Temp) {
			delete[] val; delete[] col_ind; delete[] row_ptr;
	}
		else {

#ifdef SOLOVEICHIC_SPEED
			delete[] val; delete[] col_ind; delete[] row_ptr;
#endif
		}
	}
	
	// ������ � ��������������� ILU ��������������������
	// �� ������������ ������� ����������� ������ ��������
	// �������� ���������� � ���� ������������.
	//potent=BiSoprGrad(&sparseS, sparseM,  rthdsd, NULL, nodes);
	// ��� ��������� �.�. ����������� ������������������ ��������� ���������� ?
	//potent[Temp]=SoloveichikAlg( &sparseS, sparseM, rthdsd, NULL, true);

	// ������������ ������
	delete[] rthdsd;

	if (iVar == PAm) {
		//simplesparsefree(sparseM, maxelm);
		//system("pause");
		//freeIMatrix(&sparseS);
		//system("pause");
	}
	else {
		if (iVar == Temp) {
			//simplesparsefree(sparseM, maxelm);
			//system("pause");
			//freeIMatrix(&sparseS);
			//system("pause");
		}
		else {

#ifdef SOLOVEICHIC_SPEED
			simplesparsefree(sparseM, maxelm);
			//system("pause");
			freeIMatrix(&sparseS);
			//system("pause");
#endif
	}
	}

} // solve

void my_version_SIMPLE_Algorithm(Real& continity, Real*** &B) {
	// ��������� ���� �������� �������������� ��������.
	
	Real res = 0.0;
	// ������� ������������� ���� ���������:
	
	// ���� ���������� ������������
	if (DEBUG) printf("Vx \n");
	solve(Vx, res, B);
	// ������� ���������� ������� � ��������� tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();
	if (DEBUG) printf("Vy \n");
	solve(Vy, res, B);
	// ������� ���������� ������� � ��������� tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();

#pragma omp parallel for
	for (int i = 0; i < maxelm; ++i) {
		// ���������� � ���� �������� ���������������� ��������� �������������. 
		// High Order Relaxation Factor 0.25
		Real fHORF = 0.25;
		if (!(boundary[Vx][i] && (!neiman[Vx][i]))) {
			potent[Vx][i] = fHORF * potent[Vx][i] + (1.0 - fHORF) * potent[Vxcor][i];
		}
		if (!(boundary[Vy][i] && (!neiman[Vy][i]))) {
			potent[Vy][i] = fHORF * potent[Vy][i] + (1.0 - fHORF) * potent[Vycor][i];
		}
	}

	//exporttecplotxy360( nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent, rhie_chow);
	//system("pause");



	// ������ ��������� ��� �������� ��������:
	// ����� �������� ������������� �������� �������� ����.
	if (DEBUG) printf("PAm\n");
    solve(PAm,continity, B);

	// ������� ���������� ������� � ��������� tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();

	bool bfreePressure = true;
	for (int i = 0; i < maxelm; i++) if (((boundary[PAm][i]) && (!neiman[PAm][i]))) {
		bfreePressure = false;
		//std::cout << i << " incomming";
		//system("pause");
	}
	if (fabs(dgx * dgy) > 1.0e-20) {
		// ���� ���� ���� ������.
		bfreePressure = false;
	}
	if (bfreePressure) {
		Real sum = 0.0;
		Real vol = 0.0;
		for (int i = 0; i < maxelm; i++) {
			// ���������� �������� �������� ������������ ������:
			Real dx = 0.0, dy = 0.0; // ������� ������������ ������
			volume(i, nve, nvtx, x, y, dx, dy);

			sum += potent[PAm][i] * dx * dy;
			vol += dx * dy;
		}

		for (int i = 0; i < maxelm; i++) {
			potent[PAm][i] -= sum / vol;
		}
	}
	// ��������� ��������:
#pragma omp parallel for
	for (int i = 0; i < maxelm; ++i) {
		if (!(boundary[PAm][i] && (!neiman[PAm][i]))) {

			//����� P+deltaP
			// P = (1-alphaP)*P+alphaP*(P+deltaP);
			potent[Press][i] += alpha[PAm] * potent[PAm][i];
		}
	}
	// ��������� ��������:

#pragma omp parallel for
	for (int i=0; i<maxelm; ++i) {
		int tid = omp_get_thread_num();

		correct(i, slau, Vx, nvtx, boundary, potent, sumanb, tau, x, y, sosed, neiman, norm, prop, nve, alpha, B[tid], Flux_gran_relx, Flux_gran);
        correct(i, slau, Vy, nvtx, boundary, potent, sumanb, tau, x, y, sosed, neiman, norm, prop, nve, alpha, B[tid], Flux_gran_relx, Flux_gran);
	}
	// ������� ���������� ������� � ��������� tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();
#pragma omp parallel for
    for (int i=0; i<maxelm; ++i) {
		// ���������� ���� �������� 
		// ��������������� ���������
		// �������������, �����, �����
		// ����� ����������� � ���� ������ 
		// ����������.
		potent[Vxcor][i] = potent[Vx][i];
        potent[Vycor][i] = potent[Vy][i];
	}
    if (DEBUG) printf("Temp\n");
	solve(Temp,res,B);

	if (DEBUG) getchar();

	//exporttecplotxy360( nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent, rhie_chow);
	//system("pause");

} // my_version_SIMPLE_Algorithm


int _tmain(int argc, _TCHAR* argv[])
{
	// ����� �������.
	unsigned int calculation_main_start_time_global_Depend = 0; // ������ ����� ��.
	unsigned int calculation_main_end_time = 0; // ��������� ����� ��.
	unsigned int calculation_main_seach_time = 0; // ����� ���������� ������� ���� � ��.

	calculation_main_start_time_global_Depend = clock(); // ������ ������ �����.


	// ��������� ��������� � ������� Windows
	// ������ ����������� ����� � ������� �����.
	//setlocale(LC_ALL, "");
	system("mode con cols=166 lines=12000");

	std::cout << "AliceFlow.v.0.02 2D CFD solver 2011 (revised 2021)\n";
	std::cout << "convection scheme Upwind\n";
	std::cout << "full calculation on 2 thread (auto decomposition)\n";
	std::cout << " my version SIMPLE algorithm (SIMPLE 1972 S.Patankar, B.Spalding)...\n";
	std::cout << "Rhie - Chow 1983.\n";
	std::cout << "SEIDEL projection method for Speed\n";
	std::cout << "ICCG Incompleate Cholesky Conjuate Gradient method for PRESSURE\n";
	std::cout << "Steady State CFD calculation\n";
	std::cout << "DavisTest Delphi (Pascal), AliceFlow.v.0.02 (C/C++) copyright kirill7785@mail.ru\n";
	std::cout << "for support send message on kirill7785@mail.ru Kirill Ivanov (kirill7785)\n";

	meshin("meshin.txt", nve, maxsos, MAXVAL, maxnod, maxelm, x, y, nvtx, boundary, potent, sumanb, Flux_gran, Flux_gran_relx, tau, prop, sosed, neiman, norm, nvtxcell, ncell);
    my_malloc2(slau, MAXVAL, maxelm, alpha); // ��������� ������ ��� ������������ ����
	my_init();
	loadcas("solver.txt", MAXIT, MAXTS, iSOLVER, iSHEME, BTIMEDEP, TAU, dgx, dgy, temp_ref);

	
	omp_set_num_threads(inumcore);
	nested_dissection(nvtx, x, y, nve, maxnod, maxelm, prop, sosed,  nvtxcell, ncell, MAXVAL, neiman, boundary, norm,  potent);

	Real*** B = new Real**[inumcore];
	for (int tid = 0; tid < inumcore; ++tid) {
		B[tid] = new Real * [3];
		for (int l = 0; l < 3; l++) B[tid][l] = new Real[3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[tid][i][j] = 0.0;
	}

	data_ell = new Real * [maxelm];
	for (int i = 0; i < maxelm; ++i) {
		data_ell[i] = new Real[MAX_STRING_LENGTH_ELL];
	}
	coll_ell = new int* [maxelm];
	for (int i = 0; i < maxelm; ++i) {
		coll_ell[i] = new int[MAX_STRING_LENGTH_ELL];
#pragma omp parallel for
		for (int j = 0; j < MAX_STRING_LENGTH_ELL; ++j) {
			coll_ell[i][j] = -1; // �������������!!!
		}
	}

	init_qass(maxelm, sosed, nvtx, x, y, nve);

	//marena();
	//int i;
	//for (int i=0; i<100; i++) solve(Temp);

    FILE *fpcont; // ���� � ������� ����� ������������ �������
	errno_t err;
	// �������� ����� ��� ������ �������� �����������
	if ((err = fopen_s( &fpcont, "continity.txt", "w")) !=0) {
	   printf("Create File continity.txt Error\n");
       system("pause");
       exit(0);
	}
	else {
         fprintf(fpcont, " Evalution residual\n");
         fprintf(fpcont, " iter \t\t continity\n");
         
		 // 16.08.2021 ���������� 300 ��������.
		 std::cout<<" n=" << maxelm << std::endl;

		 Real continity=1.0; // �������������
		 bool bfirst0 = true;
		 Real continity0 = 1.0;
		 Real continity0ld = 1.0;
	     for (int i=0; i<2500; i++) {
			 // getchar();
			 nnz_ell = 0;

#pragma omp parallel for
			 for (int k = 0; k < maxelm; ++k) {
				 for (int j = 0; j < MAX_STRING_LENGTH_ELL; ++j) {
					 if (coll_ell[k][j] == -1) break; // �������� ��� ���� �������.
					 //data_ell[k][j] = 0.0;
					 coll_ell[k][j] = -1;
				 }
			 }

			 continity0ld = continity;
		     my_version_SIMPLE_Algorithm(continity, B);
			 if (continity > continity0ld) {
				 glRCh = glRCh / 1.001;
			 }
			 if (bfirst0) {
				 continity0 = continity;
				 if (fabs(continity) > 1.0e-17) {
					 bfirst0 = false;
				 }
			 }
			 if (!bfirst0) {
				 fprintf(fpcont, "%d %e\n", i + 1, continity/continity0); // ������ � ���
				 //if (i % 50 == 0)
				 Real alldistit = (continity / continity0) / 1.0e-3;
				 if (i % 10 == 0) {
					 std::cout << i + 1 << " " << continity / continity0 << "  " << (int)(alldistit / (continity0ld / continity)) << " Pe="<< Pe_max<< "  Re="<< Re_max << std::endl;
				 }
				 if (continity / continity0 < 1.0e-3) break; // ������� ��������.
			 }
			 else {
				 fprintf(fpcont, "%d %e\n", i + 1, 1.0); // ������ � ���
				 //if (i % 50 == 0)
					 std::cout << i + 1 << " " << 1.0 << std::endl;
			 }
			 //exporttecplotxy360( nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent, rhie_chow);
	     }
		 //std::cout << 2000 << " " << continity << std::endl;

		 fclose(fpcont); // �������� ����� ��� ������ �������. 
	
		 for (int i = 0; i < maxelm; ++i) {
			 delete[] data_ell[i];
			 delete[] coll_ell[i];
		 }
		 delete[] data_ell;
		 delete[] coll_ell;
					
		

		 // ������� ���������� ������� � ��������� tecplot360
	     exporttecplotxy360( nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	

		 my_free2(slau, MAXVAL, maxelm, alpha);
		 for (int i = 0; i < 8; ++i) {

			 if (i < 3) {
				 delete[] sumanb[i];
			 }
			 if (i < 4) {
				 delete[] tau[i];
				 delete[] nvtx[i];
				 delete[] nvtxcell[i];
				 delete[] Flux_gran[i];
				 delete[] Flux_gran_relx[i];
			 }
			 delete[] sosed[i];
			 if (i < 5) {
				 delete[] prop[i];
			 }
			 if (i < 7) {
				 delete[] potent[i];
			 }
		 }
		 delete[] norm[0];
		 delete[] norm[1];
		 delete[] norm;
		 norm = nullptr;
		 delete[] Flux_gran;
		 delete[] Flux_gran_relx;
		 delete[] sumanb;		 
		 Flux_gran=nullptr;
		 Flux_gran_relx = nullptr;
		 sumanb = nullptr;
		 delete[] tau;
		 tau = nullptr;
		 delete[] nvtxcell;
		 delete[] nvtx;
		 delete[] sosed;
		 delete[] prop;
		 delete[] potent;
		 delete[] x;
		 delete[] y;
		 nvtxcell = nullptr;
		 x = nullptr;
		 y = nullptr;
		 prop = nullptr;
		 nvtx = nullptr;
		 sosed = nullptr;
		 potent = nullptr;
		 for (int j = 0; j < MAXVAL; j++) {
			 delete[] boundary[j];
			 delete[] neiman[j];
		 }
		 delete[] boundary;
		 delete[] neiman;
		 delete[] qass;

		 calculation_main_end_time = clock();
		 calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;

		 // ����� ����� ����������.
		 int im = 0, is = 0, ims = 0;
		 im = (int)(calculation_main_seach_time / 60000); // ������
		 is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // �������
		 ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

		 if (im == 0) {
			 printf("time calculation is:  %d second %d millisecond\n", is, 10 * ims);
		 }
		 else {
			 printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);
		 }

		 
	}
	
	for (int tid = 0; tid < inumcore; ++tid) {
		for (int l = 0; l < 3; l++) delete[] B[tid][l];
		delete[] B[tid];
	}
	delete[] B;

	std::cout << "calculation compleate..." << std::endl;
	system("pause");

	return 0;
}

#endif