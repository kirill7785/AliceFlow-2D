// ALICEFLOW0_02.cpp :
// 2D Решатель для уравнений Навье Стокса (ламинарных в стационарной постановке). 
// Реализация алгоритма SIMPLE 1972 года на совмещённых сетках.
// подход Рхи-Чоу 1983 год.
// Программа считывает два входных файла:
// 1. meshin.txt
// 2. solver.txt
// Файлы meshin.txt и solver.txt подготавливаются и пишутся в программе DavisTest Delphi
//  написанной на языке Паскаль в среде Embarcadero Delphi.
// 
// Поддерживается только схема Upwind (против потока для конвекции).
// Для поправки давления используется решатель ICCG, для скорости метод Зейделя.
// Форма расчётной области прямоугольная с произвольным числом прямоугольных вырезов.
// Используется структурированная расчётная сетка.
// 
// пересмотрено в августе 2021 (14-16 августа 2021). 
// 28 апреля 2011 Движение в изотермической каверне.
// 
// Планы на 17 августа 2019: избавиться от SPARSEM в пользу ell формата. Память выделять только один раз для ell.
// Взять функкции работы с ell форматом из 3D солвера. Сделать рабиение на две части и разделитель и распараллелить ICCG и SEIDEL
// на две части, взять заготовки из 3D кода.
// Текущее время бенчмарка обтекания квадрата 8с 610мс.
// 
// Планы на 18 августа 2019: сделать экспорт в AliceFlow0_02(AliceFlow2D) из AliceFlow_0_59(AliceFlow3D). Протестировать AliceFlow0_02.
// 
// Обнаружена расходимость ICCG солвера если не применять для давления ap / 0.9999; На каверне. 24.08.2021. Попробуем amg1r5 солвер для давления.
// метод ICCG достаточно надёжно отрабатывает если применять ap / 0.9999 поэтому с amg можно повременить.
//
// Планы дальнейшего развития программы:
// A. Графическая визуализация на OpenGL онлайн во время расчёта.
// B. Нахождение производных тех величин которые рассчитываются. Нахождение и экспорт их в техплот.
// C. Новая физика - уравнения Симоненко - Зеньковской в стационарной постановке (если будет время).
// D. Схема SMARTER для конвекции.
// E. Генератор meshgen2 расчёт в AliceFlow.v.0.02 2D на QUAD Pave расчётной сетке.
// F. Нестационарность, dynamic mesh, VOF метод.
// Но главное чтобы программа бысро считала, C/C++ Visual Studio, Intel Parallel Studio. Постоянное профилирование. Скорость расчёта главное.
// Т.к. уже есть один общий медленный двумерный солвер DavisTest, в нем хоть и реализованы п. F но скорость счета черепашья. Он не подходит.
// Не отступать от профилировщика Visual Studio.



#include "stdafx.h"
#include <iostream> // функции cout, cin
#include <stdlib.h> // Для функции exit, atoi, atof
#include <string.h> // strcat, strstr
#include <math.h> // математические функции
#include <windows.h> // для функции WinExec
#include <time.h> // для скорости выполнения

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

int inumcore = 2; // число ядер процессора

// true - печатает отладочную информацию в консоль.
const bool DEBUG = false;
//#define SOLOVEICHIC_SPEED 1



typedef struct T_PAIR_OMP
{
	int s, e;
} PAIR_OMP;

PAIR_OMP s_par[10] = { {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1} };

const int MAX_STRING_LENGTH_ELL_THERMAL = 7; // 27 диагоналей.
//const int MAX_STRING_LENGTH_ELL_MECHANICAL = 81; // 81 диагональ.
int MAX_STRING_LENGTH_ELL = MAX_STRING_LENGTH_ELL_THERMAL;

Real** data_ell = nullptr;
int** coll_ell = nullptr;
int nnz_ell = 0;

#include "my_meshin_v0_02.cpp" // считывание входного сеточного файла

// составление матрицы для обощённого уравнения конвекции-диффузии на совмещённой сетке
#include "pamendment2.cpp" 

// уже подключён в my_elmatr_quad_f.c
//#include "my_linalg.c" // самописные функции линейной алгебры
// Для функций: 
// eqsolve_simple_gauss - решает СЛАУ методом исключения Гаусса
// eqsolv_simple_holesskii - решает СЛАУ методом разложения Холесского

// экспорт результата в программу tecplot360
#include "my_export_tecplot2.cpp"



// максимальное число итераций для одного шага по времени,
// число шагов по времени, 
// выбор решающего устройства для СЛАУ,
// выбор схемы для аппроксимации конвекции-диффузии.
int MAXIT, MAXTS, iSOLVER, iSHEME;



// false - стационарный
// true - нестационарный.
bool BTIMEDEP;

Real TAU; // шаг по времени

int maxnod; // максимальный номер узла (размерность массива)
int maxelm; // максимально допустимое число элементов (элементов меньше чем узлов) (число контрольных объёмов)
// по умолчанию элементы треугольники. Тип элемента определяется при считывании файла с сеткой.
int nve=4; // число узловых переменных элемента (программа работает только для квадратных).

int **nvtx=nullptr; // список узлов для каждого элемента (контрольного объёма)
int **sosed = nullptr; // соседние контрольные объёмы для каждого контрольного объёма
Real *x = nullptr, *y = nullptr; // узловые координаты
// для графической визуализации
int **nvtxcell = nullptr; // связи для контрольных объёмов.
int ncell; // количество связей для контрольных объёмов.

bool **boundary = nullptr; // истина для граничных контрольных объёмов (лог. переменная).
bool **neiman = nullptr; // истина если на грагнице стоит условие Неймана.  (лог. переменная).
int **norm = nullptr; // нормаль к границе расчётной области (-1) - если внутренний узел.
Real **potent = nullptr; // массив узловых потенциалов (искомых функций)
Real **prop = nullptr; // свойства материалов
Real dgx=0.0, dgy=0.0; // ускорение свободного падения
Real temp_ref=0.0; // опорное значение температуры в приближении Буссинеска.


Real** sumanb = nullptr;
Real** tau = nullptr; // псевдовремя.
Real** Flux_gran = nullptr, ** Flux_gran_relx = nullptr;



// параметры нижней релаксации:
Real *alpha = nullptr;


equation **slau = nullptr; // коэффициенты матрицы СЛАУ


// выделение  памяти под матрицу коэффициентов
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

// выделение  памяти под матрицу коэффициентов
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



// решает одно уравнение
// например уравнение теплопроводности.
// res - возвращаемая невязка.
void solve(int iVar, Real& res, Real***& B) {
	ilim = 0;
	if (iVar == PAm) {
		Pe_max = 0.0;
		Re_max = 0.0;
	}
	Real alphaP = 0.3;

	// сборка матрицы СЛАУ
	Real res1 = 0.0;

	Real avgx = 0.0; 
	Real avgy = 0.0;
	switch (iVar) {
	case PAm:

#pragma omp parallel for
		for (int i_1 = 0; i_1 < maxelm; ++i_1) {
			int tid = omp_get_thread_num();

			calc_tau(i_1, slau, nvtx, boundary, potent, sumanb, Flux_gran, Flux_gran_relx, tau, x, y, prop, sosed, neiman, norm, nve, alpha, B[tid]);//псевдовремя.
		}
		if (0) {
			// При высокой скорости течения ни в коем случае не осреднять псевдовремя.
			for (int i_1 = 0; i_1 < maxelm; ++i_1) {
				avgx += tau[E][i_1] + tau[W][i_1];
				avgy += tau[N][i_1] + tau[S][i_1];
			}
			avgx = avgx / (2.0 * maxelm);
			avgy = avgy / (2.0 * maxelm);
			for (int i_1 = 0; i_1 < maxelm; ++i_1) {
				// Среднее гармоническое
				tau[E][i_1] = 2.0 * tau[E][i_1] * avgx / (tau[E][i_1] + avgx);
				tau[W][i_1] = 2.0 * tau[W][i_1] * avgx / (tau[W][i_1] + avgx);

				tau[N][i_1] = 2.0 * tau[N][i_1] * avgy / (tau[N][i_1] + avgy);
				tau[S][i_1] = 2.0 * tau[S][i_1] * avgy / (tau[S][i_1] + avgy);
			}
		}
#pragma omp parallel for
		for (int i_1 = 0; i_1 < maxelm; ++i_1) {
			int tid = omp_get_thread_num();

			my_elmatr_quad_PAm(i_1, slau, nvtx, boundary, potent, sumanb, Flux_gran, Flux_gran_relx, tau, x, y, prop, sosed, neiman, norm, nve, alpha, B[tid]); // поправка давления
		}

		if (0) {
			// Включение этого тоже даёт расходимость.

			// Диагональное предобуславливание.
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

		// симметризация СЛАУ
		// для ICCG солвера обработка условий Дирихле
		for (int i = 0; i < maxelm; ++i) if ((boundary[PAm][i]) && (!neiman[PAm][i])) {
			// только для нулевого значения поправки
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
				// Это даёт сходимость ICCG решателя для каверны !!!
				// если нету slau[PAm][i_1].ap / 0.999999; то для каверны АВОСТ.
				slau[PAm][i_1].ap = slau[PAm][i_1].ap / 0.9999;// 999;//0.99; 0.9999
			}
		}

		break;
	default:
		if ((iVar == Vx) || (iVar == Vy) || (iVar == Temp)) {
#pragma omp parallel for
			for (int i_1 = 0; i_1 < maxelm; ++i_1) {
				int tid = omp_get_thread_num();

				my_elmatr_quad_F(i_1, slau[iVar], iVar, BTIMEDEP, TAU, iSHEME, nvtx, boundary, potent, sumanb, Flux_gran, x, y, prop, sosed, neiman, norm, nve, alpha, dgx, dgy, temp_ref, B[tid]); // стационарный

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

	Real* rthdsd; // правая часть системы уравнений
	rthdsd = new Real[maxelm];


	SIMPLESPARSE sparseM; // разреженная матрица
#ifdef SOLOVEICHIC_SPEED
	IMatrix sparseS;
#endif
	

    
    // разреженная матрица в формате CSIR
    //Real *adiag, *altr;
    //int *jptr, *iptr;

    //system("pause");
	if (iVar == PAm) {
		// выделение памяти и инициализация для 
		// простейшей разреженной матрицы.
		//initsimplesparse(sparseM, maxelm);
		//initIMatrix(&sparseS, maxelm);
	}
	else {

		if (iVar == Temp) {
			// выделение памяти и инициализация для 
			// простейшей разреженной матрицы.
			//initsimplesparse(sparseM, maxelm);
			//initIMatrix(&sparseS, maxelm);
		}
		else {

#ifdef SOLOVEICHIC_SPEED
			// выделение памяти и инициализация для 
			// простейшей разреженной матрицы.
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
		
		Real nonzeroEPS=1e-37; // для отделения вещественного нуля

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
				* // Не работает в многопоточном режиме.
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


    // Блок решения СЛАУ:

    //1. Прямые методы решения СЛАУ:

    // 1.1. Для плотных матриц.

	// решает СЛАУ методом разложения Холесского
	// для симметричной положительно определённой матрицы s.
	// без учёта разреженности матрицы. В два раза более 
	// эффективна чем метод прямого исключения Гаусса.
	//eqsolv_simple_holesskii(s, nodes, rthdsd, potent);
	// решает СЛАУ методом исключения Гаусса
	// без выбора главного элемента и без 
	// учёта разреженности матрицы.
	//eqsolve_simple_gauss(s, nodes, rthdsd, potent[Temp]);

    // 1.2. Для разреженных матриц.

	// Метод К.Ф. Гаусса для разреженной возможно несимметричной матрицы
    //calculateSPARSEgaussArray(&sparseS, potent, rthdsd);


	//2. Итерационные методы решения СЛАУ:

	// 2.1. Для плотных матриц.

	// Метод Гаусса-Зейделя-Ричардсона-Либмана SOR
	// для неразреженной матрицы СЛАУ.
	// метод характеризуется медленной сходимостью.
	//Seidel(s, rthdsd, potent, nodes, 1e-5, 1.855);
	// Решает неразреженную СЛАУ методом
	// сопряжённых градиентов.
	// для симметричной положительно определённой матрицы s.
	// Метод является уточняющим, т.е. если известен вектор
	// начального приближения то его нужно поставить вместо NULL.
	//potent=SoprGrad(s, rthdsd, NULL, nodes);
		
    // 2.2. Для разреженных матриц.

	// 2.2.1 МСГ для CRS (CSIR, CSIR_ITL) формата симметричный случай.
    ///*
    //simplesparsetoCRS(sparseM, val, col_ind, row_ptr, nodes);
	//printM_and_CRS(sparseM,val,col_ind,row_ptr,nodes);
    //potent=SoprGradCRS(val, col_ind, row_ptr, rthdsd, NULL, nodes);//*/

	// МСГ (Congruate Gradients) для CSIR формата хранения
	// симметричный положительно определённый случай.
    //simplesparsetoCSIR(sparseM, adiag, altr, jptr, iptr, nodes);
	//printM_and_CSIR(sparseM, adiag, altr, jptr, iptr,  nodes);
	//int inz=(int)((sparseM.n-nodes)/2.0);
	//potent=SoprGradCSIR(adiag, altr, jptr, iptr, rthdsd, NULL, nodes, inz);
	//potent=SoloveichikAlgCSIR_SPD(nodes, adiag, altr, jptr, iptr, rthdsd, NULL, true);
	//potent=SoloveichikAlgCSIR_SPDgood(nodes, inz, adiag, altr, jptr, iptr, rthdsd, NULL, true);
	if (iVar==PAm) {



		 for (int i=0; i<maxelm; ++i) potent[PAm][i]=0.0;
		 if (0) {
			 // Сильная недоитерированность, решение нефизично.
			 SOR(slau[PAm], potent[PAm], rthdsd, maxelm, PAm);

		 }
		 else {

			// SOR(slau[PAm], potent[PAm], maxelm);
			 ICCG(sparseM, rthdsd, potent[PAm], maxelm);
		 }
	}
	//*/

	// 2.2.2 Методы для CRS (CSIR, CSIR_ITL) формата возможно несимметричный случай.
	if (iVar != PAm) {
		// разреженная матрица в формате CRS
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
					// Для компонент скорости показывает себя единственно хорошо.

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
	
	// методы с интегрированным ILU предобуславливателем
	// На симметричной матрице наблюдается резкое ухудшене
	// скорости сходимости и даже расходимость.
	//potent=BiSoprGrad(&sparseS, sparseM,  rthdsd, NULL, nodes);
	// Для алгоритма Ю.Г. Соловейчика предобуславливание замедляет сходимость ?
	//potent[Temp]=SoloveichikAlg( &sparseS, sparseM, rthdsd, NULL, true);

	// освобождение памяти
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
	// начальное поле давления предполагается заданным.
	
	Real res = 0.0;
	// находим промежуточное поле скоростей:
	
	// цикл устранения нелинейности
	if (DEBUG) printf("Vx \n");
	solve(Vx, res, B);
	// экспорт результата расчёта в программу tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();
	if (DEBUG) printf("Vy \n");
	solve(Vy, res, B);
	// экспорт результата расчёта в программу tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();

#pragma omp parallel for
	for (int i = 0; i < maxelm; ++i) {
		// Релаксация к полю скорости удовлетворяющему уравнению неразрывности. 
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



	// решаем уравнение для поправки давления:
	// можно добавить инициализацию поправки давления нулём.
	if (DEBUG) printf("PAm\n");
    solve(PAm,continity, B);

	// экспорт результата расчёта в программу tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();

	bool bfreePressure = true;
	for (int i = 0; i < maxelm; i++) if (((boundary[PAm][i]) && (!neiman[PAm][i]))) {
		bfreePressure = false;
		//std::cout << i << " incomming";
		//system("pause");
	}
	if (fabs(dgx * dgy) > 1.0e-20) {
		// Типо тест Валь Девиса.
		bfreePressure = false;
	}
	if (bfreePressure) {
		Real sum = 0.0;
		Real vol = 0.0;
		for (int i = 0; i < maxelm; i++) {
			// вычисление размеров текущего контрольного объёма:
			Real dx = 0.0, dy = 0.0; // размеры контрольного объёма
			volume(i, nve, nvtx, x, y, dx, dy);

			sum += potent[PAm][i] * dx * dy;
			vol += dx * dy;
		}

		for (int i = 0; i < maxelm; i++) {
			potent[PAm][i] -= sum / vol;
		}
	}
	// коррекция давления:
#pragma omp parallel for
	for (int i = 0; i < maxelm; ++i) {
		if (!(boundary[PAm][i] && (!neiman[PAm][i]))) {

			//новое P+deltaP
			// P = (1-alphaP)*P+alphaP*(P+deltaP);
			potent[Press][i] += alpha[PAm] * potent[PAm][i];
		}
	}
	// коррекция скорости:

#pragma omp parallel for
	for (int i=0; i<maxelm; ++i) {
		int tid = omp_get_thread_num();

		correct(i, slau, Vx, nvtx, boundary, potent, sumanb, tau, x, y, sosed, neiman, norm, prop, nve, alpha, B[tid], Flux_gran_relx, Flux_gran);
        correct(i, slau, Vy, nvtx, boundary, potent, sumanb, tau, x, y, sosed, neiman, norm, prop, nve, alpha, B[tid], Flux_gran_relx, Flux_gran);
	}
	// экспорт результата расчёта в программу tecplot360
	//exporttecplotxy360(nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent);
	//getchar();
#pragma omp parallel for
    for (int i=0; i<maxelm; ++i) {
		// запоминаем поле скорости 
		// удовлетворяющее уравнению
		// неразрывности, затем, чтобы
		// потом осуществить к нему нижнюю 
		// релаксацию.
		potent[Vxcor][i] = potent[Vx][i];
        potent[Vycor][i] = potent[Vy][i];
	}
    if (DEBUG) printf("Temp\n");
	solve(Temp,res,B);

	if (DEBUG) getchar();

	//exporttecplotxy360( nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent, rhie_chow);
	//system("pause");

} // my_version_SIMPLE_Algorithm


int main()
{
	// Замер времени.
	unsigned int calculation_main_start_time_global_Depend = 0; // начало счёта мс.
	unsigned int calculation_main_end_time = 0; // окончание счёта мс.
	unsigned int calculation_main_seach_time = 0; // время выполнения участка кода в мс.

	calculation_main_start_time_global_Depend = clock(); // момент начала счёта.


	// Поддержка кириллицы в консоли Windows
	// Меняет разделитель целой и дробной части.
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

	char name_mesh_file0[11] = "meshin.txt";
	char* name_mesh_file=new char[11];
	for (int i = 0; i < 10; ++i) name_mesh_file[i] = name_mesh_file0[i];
	name_mesh_file[10] = '\0';


	meshin(name_mesh_file, nve, maxsos, MAXVAL, maxnod, maxelm, x, y, nvtx, boundary, potent, sumanb, Flux_gran, Flux_gran_relx, tau, prop, sosed, neiman, norm, nvtxcell, ncell);
	delete[] name_mesh_file;
	my_malloc2(slau, MAXVAL, maxelm, alpha); // выделение памяти под коэффициенты СЛАУ
	my_init();
	char name_solver_file0[11] = "solver.txt";
	char* name_solver_file = new char[11];
	for (int i = 0; i < 10; ++i) name_solver_file[i] = name_solver_file0[i];
	name_solver_file[10] = '\0';

	loadcas(name_solver_file, MAXIT, MAXTS, iSOLVER, iSHEME, BTIMEDEP, TAU, dgx, dgy, temp_ref);
	delete[] name_solver_file;
	
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
			coll_ell[i][j] = -1; // инициализация!!!
		}
	}

	init_qass(maxelm, sosed, nvtx, x, y, nve);

	//marena();
	//int i;
	//for (int i=0; i<100; i++) solve(Temp);

    FILE *fpcont; // файл в который будут записываться невязки
	errno_t err;
	// создание файла для записи значений температуры
	if ((err = fopen_s( &fpcont, "continity.txt", "w")) !=0) {
	   printf("Create File continity.txt Error\n");
       system("pause");
       exit(0);
	}
	else {
         fprintf(fpcont, " Evalution residual\n");
         fprintf(fpcont, " iter \t\t continity\n");
         
		 // 16.08.2021 достаточно 300 итераций.
		 std::cout<<" n=" << maxelm << std::endl;

		 Real continity=1.0; // инициализация
		 bool bfirst0 = true;
		 Real continity0 = 1.0;
		 Real continity0ld = 1.0;
	     for (int i=0; i<2500; i++) {
			 // getchar();
			 nnz_ell = 0;

#pragma omp parallel for
			 for (int k = 0; k < maxelm; ++k) {
				 for (int j = 0; j < MAX_STRING_LENGTH_ELL; ++j) {
					 if (coll_ell[k][j] == -1) break; // ускоряет для гран условий.
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
				 fprintf(fpcont, "%d %e\n", i + 1, continity/continity0); // печать в мкс
				 //if (i % 50 == 0)
				 Real alldistit = (continity / continity0) / 1.0e-3;
				 if (i % 10 == 0) {
					 std::cout << i + 1 << " " << continity / continity0 << "  " << (int)(alldistit / (continity0ld / continity)) << " Pe="<< Pe_max<< "  Re="<< Re_max << std::endl;
				 }
				 if (continity / continity0 < 1.0e-3) break; // Решение получено.
			 }
			 else {
				 fprintf(fpcont, "%d %e\n", i + 1, 1.0); // печать в мкс
				 //if (i % 50 == 0)
					 std::cout << i + 1 << " " << 1.0 << std::endl;
			 }
			 //exporttecplotxy360( nve, maxelm, ncell, nvtx, nvtxcell, x, y, potent, rhie_chow);
	     }
		 //std::cout << 2000 << " " << continity << std::endl;

		 fclose(fpcont); // закрытие файла для записи невязки. 
	
		 for (int i = 0; i < maxelm; ++i) {
			 delete[] data_ell[i];
			 delete[] coll_ell[i];
		 }
		 delete[] data_ell;
		 delete[] coll_ell;
					
		

		 // экспорт результата расчёта в программу tecplot360
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

		 // Общее время вычисления.
		 int im = 0, is = 0, ims = 0;
		 im = (int)(calculation_main_seach_time / 60000); // минуты
		 is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // секунды
		 ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // миллисекунды делённые на 10

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