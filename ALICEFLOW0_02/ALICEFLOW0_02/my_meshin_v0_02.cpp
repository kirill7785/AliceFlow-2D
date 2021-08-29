// Файл my_meshin_v0_02.cpp
// Считывание сетки и граничных условий
// для системы уравнений Навье-Стокса.

#pragma once
#ifndef MY_MESHIN_V0_02_CPP
#define MY_MESHIN_V0_02_CPP 1

//#include <stdlib.h> // Для функции exit, atoi, atof
//#include <string.h> // strcat, strstr



// Искомые величины
const int Temp = 0; // температура
const int Vx = 1; // горизонтальная скорость
const int Vy = 2; // вертикальная скорость
const int Press = 3; // давление
const int PAm = 4; // поправка давления
const int MAXVAL = 5; // количество искомых величин
const int  Vxcor = 5; // скорости с предыдущей итерации удовлетворяющие 
const int Vycor = 6; // уравнению неразрывности


// схемы для аппроксимации конвекции-диффузии (старые)
const int CR = 1; // Центрально-разностная
const int UDS = 2; // Противопоточная первого порядка
const int COMB = 3; // Комбинированная 
const int POLY = 4; // Полиномиальная C. Патанкара
const int EXP = 5; // экспоненциальная схема
const int BULG = 6; // схема В.К. Булгакова (23) из статьи
const int POW = 7; // показательная

// схемы для аппроксимации конвекции-диффузии (новые)
const int CR2 = 100; // Центральные разности
const int UDS2 = 101; // Противопоточная первого порядка
const int EXP2 = 102; // экспоненциальная схема (точная)
const int KUD = 103; // схема предложенная в диссертации Кудинова Павла Ивановича

const int distsheme = 100; // константа перехода от новой схемы к новой

const int Rho = 0; // плотность
const int Cp = 1; // теплоёмкость
const int Lam = 2; // теплопроводность
const int Mu = 3;  // динамическая вязкость
const int BETA_T = 4; // коэффициент линейного температурного расширения

const int E = 0; // (east) восток
const int N = 1; // север
const int W = 2; // запад
const int S = 3; // юг
const int EE = 4;
const int NN = 5;
const int WW = 6;
const int SS = 7;
const int maxsos = 8; // максимальное число соседей


const int BICGCRS = 1; // бисопряжённые градиенты
const int SOLOVALGCRS = 2; // алгоритм Ю.Г. Соловейчика
const int BICGSTABCRS = 3; // алгоритм Ван Дер Ворста BiCGStab

// Для генерации матрицы СЛАУ требуется переупорядочивание элементов
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
template <typename TVAL>
void swap(TVAL*& v, int i, int j);

// Вот алгоритм PivotList
template <typename TVAL>
int PivotList(TVAL*& list, int first, int last);

// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
template <typename TVAL>
void QuickSort(TVAL*& list, int first, int last);


// вычисляет размеры контрольного объёма
void volume(int iP, int nve, int**& nvtx, Real*& x, Real*& y, Real& dx, Real& dy);

// 17.08.2021
void nested_dissection(int**& nvtx, Real*& x, Real*& y, int nve, int maxnod, int maxelm,
	Real**& prop, int**& sosed, int**& nvtxcell, int& ncell, int maxval, 
	bool**& neiman, bool**& boundary, int**& norm,
	Real**& potent) {
	
	if (inumcore == 2) {

		Real avg = 0.0;
		if (0) {
#pragma omp parallel for reduction(+: avg)
			for (int i = 0; i < maxelm; ++i) {

				Real centx = 0.25 * (x[nvtx[0][i] - 1] + x[nvtx[1][i] - 1] + x[nvtx[2][i] - 1] + x[nvtx[3][i] - 1]);
				avg += centx;
			}
			avg /= (Real)((1.0 * maxelm));
		}
		else {
			// медиана.
			// Разбивает гораздо лучше поровну.
			// Даже для областей усложнённой формы с вырезами.
			Real* list = new Real[maxnod];
#pragma omp parallel for
			for (int i = 0; i < maxnod; ++i) {
				list[i] = x[i];
			}
			//QuickSort(list, 0, maxnod - 1);
			std::sort(list, list+ maxnod - 1);
			avg = list[(maxnod - 1) / 2];
			delete[] list;
		}

		Real fx = 1.0e36;
		int ipos = -1;
		for (int i = 0; i < maxelm; ++i) {
			Real centx = 0.25 * (x[nvtx[0][i]-1] + x[nvtx[1][i]-1] + x[nvtx[2][i]-1] + x[nvtx[3][i]-1]);
			if (fabs(centx - avg) < fx) {
				fx = fabs(centx - avg);
				ipos = i;
			}
		}
		Real dx = 0.0, dy = 0.0;
		volume(ipos, nve, nvtx, x, y, dx, dy);
		dx *= 0.5;
		Real centx0 = 0.25 * (x[nvtx[0][ipos]-1] + x[nvtx[1][ipos]-1] + x[nvtx[2][ipos]-1] + x[nvtx[3][ipos]-1]);

		int* inumerate = new int[maxelm];
#pragma omp parallel for
		for (int i = 0; i < maxelm; ++i) {
			inumerate[i] = -1;
		}

		int iL = 0;
		int iR = 0;
		int iS = 0;


		for (int i = 0; i < maxelm; ++i) {
			Real centx = 0.25 * (x[nvtx[0][i]-1] + x[nvtx[1][i]-1] + x[nvtx[2][i]-1] + x[nvtx[3][i]-1]);

			if (centx < centx0 - dx) {
				inumerate[i] = 1; // Левая часть.
				iL++;
			}
			else if (centx > centx0 + dx) {
				inumerate[i] = 2; // Правая часть.
				iR++;
			}
			else {
				// Разделитель.
				inumerate[i] = 3;
				iS++;
			}
		}



		std::cout << "partition: iL=" << iL << " iR=" << iR << " iS=" << iS << std::endl;
		//getchar();

		s_par[1].s = 0;
		s_par[1].e = iL;
		s_par[2].s = iL;
		s_par[2].e = iL+iR;
		s_par[3].s = iL+iR;
		s_par[3].e = iL+iR+iS;

		for (int i = 0; i < maxelm; ++i) {
			if (inumerate[i] == -1) std::cout<<"error i="<<i<<std::endl;
		}


		int** nvtx_copy = new int* [4];
		for (int i = 0; i < 4; ++i) nvtx_copy[i] = new int[maxelm];

		int* new_numeric = new int[maxelm];
#pragma omp parallel for
		for (int i = 0; i < maxelm; ++i) new_numeric[i] = -1;

		Real** prop_copy = new Real * [5];
		for (int i = 0; i < 5; ++i) prop_copy[i] = new Real[maxelm];

		int** norm_copy = new int* [2];
		norm_copy[0] = new int[maxelm];
		norm_copy[1] = new int[maxelm];


		bool** boundary_copy = new bool* [maxval];
		for (int j = 0; j < maxval; j++) {
			boundary_copy[j] = new bool[maxelm];
		}
		bool** neiman_copy = new bool* [maxval];
		for (int j = 0; j < maxval; j++) {
			neiman_copy[j] = new bool[maxelm];
		}
		Real** potent_copy = new Real* [maxval];
		for (int j = 0; j < maxval; j++) {
			potent_copy[j] = new Real[maxelm];
		}

		int** sosed_copy = new int* [8];
		for (int i_2 = 0; i_2 < 8; ++i_2) sosed_copy[i_2] = new int[maxelm];

		{
			int j_scan = 0;
			for (int i = 0; i < maxelm; ++i) {
				if (inumerate[i] == 1) {
					nvtx_copy[0][j_scan] = nvtx[0][i];
					nvtx_copy[1][j_scan] = nvtx[1][i];
					nvtx_copy[2][j_scan] = nvtx[2][i];
					nvtx_copy[3][j_scan] = nvtx[3][i];

					prop_copy[Rho][j_scan] = prop[Rho][i];
					prop_copy[Cp][j_scan] = prop[Cp][i];
					prop_copy[Lam][j_scan] = prop[Lam][i];
					prop_copy[Mu][j_scan] = prop[Mu][i];
					prop_copy[BETA_T][j_scan] = prop[BETA_T][i];

					for (int j_1 = 0; j_1 < maxval; j_1++) {
						boundary_copy[j_1][j_scan] = boundary[j_1][i];
						neiman_copy[j_1][j_scan] = neiman[j_1][i];
						potent_copy[j_1][j_scan] = potent[j_1][i];
					}

					for (int i_2 = 0; i_2 < 8; ++i_2) sosed_copy[i_2][j_scan] = sosed[i_2][i];

					norm_copy[0][j_scan] = norm[0][i];
					norm_copy[1][j_scan] = norm[1][i];

					new_numeric[i] = j_scan;
					++j_scan;
				}
			}

			for (int i = 0; i < maxelm; ++i) {
				if (inumerate[i] == 2) {
					nvtx_copy[0][j_scan] = nvtx[0][i];
					nvtx_copy[1][j_scan] = nvtx[1][i];
					nvtx_copy[2][j_scan] = nvtx[2][i];
					nvtx_copy[3][j_scan] = nvtx[3][i];

					prop_copy[Rho][j_scan] = prop[Rho][i];
					prop_copy[Cp][j_scan] = prop[Cp][i];
					prop_copy[Lam][j_scan] = prop[Lam][i];
					prop_copy[Mu][j_scan] = prop[Mu][i];
					prop_copy[BETA_T][j_scan] = prop[BETA_T][i];

					for (int j_1 = 0; j_1 < maxval; j_1++) {
						boundary_copy[j_1][j_scan] = boundary[j_1][i];
						neiman_copy[j_1][j_scan] = neiman[j_1][i];
						potent_copy[j_1][j_scan] = potent[j_1][i];
					}

					for (int i_2 = 0; i_2 < 8; ++i_2) sosed_copy[i_2][j_scan] = sosed[i_2][i];

					norm_copy[0][j_scan] = norm[0][i];
					norm_copy[1][j_scan] = norm[1][i];

					new_numeric[i] = j_scan;
					++j_scan;
				}
			}

			for (int i = 0; i < maxelm; ++i) {
				if (inumerate[i] == 3) {
					nvtx_copy[0][j_scan] = nvtx[0][i];
					nvtx_copy[1][j_scan] = nvtx[1][i];
					nvtx_copy[2][j_scan] = nvtx[2][i];
					nvtx_copy[3][j_scan] = nvtx[3][i];

					prop_copy[Rho][j_scan] = prop[Rho][i];
					prop_copy[Cp][j_scan] = prop[Cp][i];
					prop_copy[Lam][j_scan] = prop[Lam][i];
					prop_copy[Mu][j_scan] = prop[Mu][i];
					prop_copy[BETA_T][j_scan] = prop[BETA_T][i];

					for (int j_1 = 0; j_1 < maxval; j_1++) {
						boundary_copy[j_1][j_scan] = boundary[j_1][i];
						neiman_copy[j_1][j_scan] = neiman[j_1][i];
						potent_copy[j_1][j_scan] = potent[j_1][i];
					}

					for (int i_2 = 0; i_2 < 8; ++i_2) sosed_copy[i_2][j_scan] = sosed[i_2][i];

					norm_copy[0][j_scan] = norm[0][i];
					norm_copy[1][j_scan] = norm[1][i];

					new_numeric[i] = j_scan;
					++j_scan;
				}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < maxelm; ++i) {
			nvtx[0][i] = nvtx_copy[0][i];
			nvtx[1][i] = nvtx_copy[1][i];
			nvtx[2][i] = nvtx_copy[2][i];
			nvtx[3][i] = nvtx_copy[3][i];

			prop[Rho][i]=prop_copy[Rho][i];
			prop[Cp][i]=prop_copy[Cp][i];
			prop[Lam][i]=prop_copy[Lam][i];
			prop[Mu][i]=prop_copy[Mu][i];
			prop[BETA_T][i]=prop_copy[BETA_T][i];

			for (int j_1 = 0; j_1 < maxval; j_1++) {
				boundary[j_1][i] = boundary_copy[j_1][i];
				neiman[j_1][i] = neiman_copy[j_1][i];
				potent[j_1][i] = potent_copy[j_1][i];
			}

			for (int i_2 = 0; i_2 < 8; ++i_2) sosed[i_2][i] = sosed_copy[i_2][i];

			norm[0][i] = norm_copy[0][i];
			norm[1][i] = norm_copy[1][i];
		}

#pragma omp parallel for
		for (int i = 0; i < maxelm; ++i) {
			if (sosed[E][i] > 0) sosed[E][i] = 1+ new_numeric[sosed[E][i]-1];
			if (sosed[W][i] > 0) sosed[W][i] = 1+new_numeric[sosed[W][i]-1];
			if (sosed[N][i] > 0) sosed[N][i] = 1+new_numeric[sosed[N][i]-1];
			if (sosed[S][i] > 0) sosed[S][i] = 1+new_numeric[sosed[S][i]-1];
			if (sosed[EE][i] > 0) sosed[EE][i] = 1+new_numeric[sosed[EE][i]-1];
			if (sosed[WW][i] > 0)  sosed[WW][i] = 1+new_numeric[sosed[WW][i]-1];
			if (sosed[NN][i] > 0) sosed[NN][i] = 1+new_numeric[sosed[NN][i]-1];
			if (sosed[SS][i] > 0) sosed[SS][i] = 1+new_numeric[sosed[SS][i]-1];
		}

#pragma omp parallel for
		for (int i = 0; i < ncell; i++) {
			nvtxcell[0][i] = 1+new_numeric[nvtxcell[0][i]-1];
			nvtxcell[1][i] = 1+new_numeric[nvtxcell[1][i]-1];
			nvtxcell[2][i] = 1+new_numeric[nvtxcell[2][i]-1];
			nvtxcell[3][i] = 1+new_numeric[nvtxcell[3][i]-1];
		}

		for (int i = 0; i < 4; ++i) delete[] nvtx_copy[i];
		delete[] nvtx_copy;

		delete[] inumerate;

		delete[] new_numeric;

		for (int i = 0; i < 5; ++i) delete[] prop_copy[i];
		delete[] prop_copy;

		for (int j = 0; j < maxval; j++) {
		delete[] 	boundary_copy[j];
		}
		delete[] 	boundary_copy;

		for (int j = 0; j < maxval; j++) {
			delete[] 	neiman_copy[j];
		}
		delete[] 	neiman_copy;

		for (int j = 0; j < maxval; j++) {
			delete[] 	potent_copy[j];
		}
		delete[] 	potent_copy;

		delete[] norm_copy[0];
		delete[] norm_copy[1];

		delete[] norm_copy;

		for (int i_2 = 0; i_2 < 8; ++i_2) delete[] sosed_copy[i_2];
		delete[] sosed_copy;
	}

}

// Для ускорения сборки матриц.
typedef struct TQASS {
	Real dx, dy;
	int iE, iW, iS, iN;
	Real dxe, dxw, dyn, dys;
	Real dy_dxe, dy_dxw, dx_dyn, dx_dys;
} QASS;

QASS* qass = nullptr;



void init_qass(int maxelm, int** &sosed, int** &nvtx, Real* &x, Real* &y, int nve) {

	for (int iP = 0; iP < maxelm; ++iP) {
		qass[iP].iE = sosed[E][iP] - 1; qass[iP].iN = sosed[N][iP] - 1; qass[iP].iW = sosed[W][iP] - 1; qass[iP].iS = sosed[S][iP] - 1;
		volume(iP, nve, nvtx, x, y, qass[iP].dx, qass[iP].dy);
		

		qass[iP].dxe = qass[iP].dx, qass[iP].dxw = qass[iP].dx, qass[iP].dyn = qass[iP].dy, qass[iP].dys = qass[iP].dy;
		if ((qass[iP].iE > -1)) qass[iP].dxe = 0.25 * (x[nvtx[0][qass[iP].iE] - 1] + x[nvtx[1][qass[iP].iE] - 1] + x[nvtx[2][qass[iP].iE] - 1] + x[nvtx[3][qass[iP].iE] - 1]);
		if ((qass[iP].iE > -1)) qass[iP].dxe -= 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
		if ((qass[iP].iW > -1)) qass[iP].dxw = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
		if ((qass[iP].iW > -1)) qass[iP].dxw -= 0.25 * (x[nvtx[0][qass[iP].iW] - 1] + x[nvtx[1][qass[iP].iW] - 1] + x[nvtx[2][qass[iP].iW] - 1] + x[nvtx[3][qass[iP].iW] - 1]);
		if ((qass[iP].iN > -1)) qass[iP].dyn = 0.25 * (y[nvtx[0][qass[iP].iN] - 1] + y[nvtx[1][qass[iP].iN] - 1] + y[nvtx[2][qass[iP].iN] - 1] + y[nvtx[3][qass[iP].iN] - 1]);
		if ((qass[iP].iN > -1)) qass[iP].dyn -= 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
		if ((qass[iP].iS > -1)) qass[iP].dys = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
		if ((qass[iP].iS > -1)) qass[iP].dys -= 0.25 * (y[nvtx[0][qass[iP].iS] - 1] + y[nvtx[1][qass[iP].iS] - 1] + y[nvtx[2][qass[iP].iS] - 1] + y[nvtx[3][qass[iP].iS] - 1]);
	
		qass[iP].dxe = fabs(qass[iP].dxe);
		qass[iP].dxw = fabs(qass[iP].dxw);
		qass[iP].dyn = fabs(qass[iP].dyn);
		qass[iP].dys = fabs(qass[iP].dys);


		qass[iP].dy_dxe = qass[iP].dy / qass[iP].dxe;
		qass[iP].dy_dxw = qass[iP].dy / qass[iP].dxw;
		qass[iP].dx_dyn = qass[iP].dx / qass[iP].dyn;
		qass[iP].dx_dys = qass[iP].dx / qass[iP].dys;


		if (qass[iP].dx <= 0) {
			std::cout << " dx=" << qass[iP].dx << std::endl;
			system("PAUSE");
		}
		if (qass[iP].dy <= 0) {
			std::cout << " dy=" << qass[iP].dy << std::endl;
			system("PAUSE");
		}

		if (qass[iP].dxe < 1.0e-20) {
			printf("zero dxe=%e iE=%d iP=%d\n", qass[iP].dxe, qass[iP].iE, iP);
			system("PAUSE");
		}

		if (qass[iP].dx != qass[iP].dx) std::cout << "dy=" << qass[iP].dy << " dxe=" << qass[iP].dxe << std::endl;
		if (qass[iP].dy != qass[iP].dy) std::cout << "dy=" << qass[iP].dy << " dxw=" << qass[iP].dxw << std::endl;


		if (qass[iP].dxe!= qass[iP].dxe) std::cout << "dy=" << qass[iP].dy << " dxe=" << qass[iP].dxe << std::endl;
		if (qass[iP].dxw != qass[iP].dxw) std::cout << "dy=" << qass[iP].dy << " dxw=" << qass[iP].dxw << std::endl;

		if (qass[iP].dy_dxe != qass[iP].dy_dxe) std::cout<< "dy="<< qass[iP].dy << " dxe="<< qass[iP].dxe <<std::endl;
		if (qass[iP].dy_dxw != qass[iP].dy_dxw) std::cout << "dy=" << qass[iP].dy << " dxw=" << qass[iP].dxw << std::endl;
		if (qass[iP].dx_dyn != qass[iP].dx_dyn) std::cout << "dx=" << qass[iP].dx << " dxe=" << qass[iP].dyn << std::endl;
		if (qass[iP].dx_dys != qass[iP].dx_dys) std::cout << "dx=" << qass[iP].dx << " dxw=" << qass[iP].dys << std::endl;
		
		
	}
}

// Динамическое выделение памяти
void my_malloc(int nve, int maxsos, int maxval, 
			   int maxnod, int maxelm, int** &nvtx, 
			   bool** &boundary, Real** &potent, Real**& sumanb, Real** & Flux_gran, Real**& Flux_gran_relx, Real** &tau, Real* &x,
			   Real* &y, Real** &prop, int** &sosed,
			   bool** &neiman, int** &norm, int** &nvtxcell, int &ncell) {
	
	qass = new QASS[maxelm]; // для ускорения сборки матрицы.
	int i;
	nvtx=new int*[nve];
    for (i=0; i<nve; i++) nvtx[i]=new int[maxelm];
	nvtxcell=new int*[nve];
    for (i=0; i<nve; i++) nvtxcell[i]=new int[ncell];
	sosed=new int*[maxsos];
    for (i=0; i<maxsos; i++) sosed[i]=new int[maxelm];

    // Vx; Vy; Press; PAm; Temp;
	boundary=new bool*[maxval];
	for (i=0; i<maxval; i++) boundary[i]=new bool[maxelm];
	potent=new Real*[maxval+2]; // 2 для Vxcor, Vycor
    for (i=0; i<(maxval+2); i++) potent[i]=new Real[maxelm];
	sumanb = new Real * [3]; // 2 для Vx, Vy
	for (i = 0; i < (3); i++) sumanb[i] = new Real[maxelm];
	Flux_gran = new Real * [4]; // 2 для Fe + 0.1*FeRhie_Chow
	for (i = 0; i < (4); i++) Flux_gran[i] = new Real[maxelm];
	Flux_gran_relx = new Real * [4]; // 2 для Fe + 0.1*FeRhie_Chow
	for (i = 0; i < (4); i++) Flux_gran_relx[i] = new Real[maxelm];
	tau = new Real * [4]; // псевдовремя.
	for (i = 0; i < (4); i++) tau[i] = new Real[maxelm];
	neiman=new bool*[maxval];
	for (i=0; i<maxval; i++) neiman[i]=new bool[maxelm];

	norm = new int*[2];
    for (i=0; i<2; i++) norm[i]=new int[maxelm];
    
    x=new Real[maxnod];
    y=new Real[maxnod];

    prop=new Real*[5];
    for (i=0; i<5; i++) prop[i]=new Real[maxelm];

} // my_malloc

// печатает на консоль постановку задачи
void print_task(Real* x, Real* y, int maxnod, int maxelm, 
				int** nvtx, Real** &prop, int nve, bool** &boundary,
				Real** &potent, int maxsos, int** &sosed, int **nvtxcell, int &ncell ) {
	int i,j;
	for (i=0; i<maxnod; i++) printf("%d %.2f %.2f\n",i+1,x[i],y[i]);
	printf("\n");
	for (j=0; j<maxelm; j++) { 
	    for (i=0; i<nve; i++) printf("%d ",nvtx[i][j]);
		printf("%.2f %.2f %.2f %.2f %.2f\n",prop[Rho][j],prop[Cp][j],prop[Lam][j],prop[Mu][j],prop[BETA_T][j]);
	}
	printf("\n");
	for (j=0; j<maxelm; j++) { 
	    for (i=0; i<maxsos; i++) printf("%d ",sosed[i][j]);
        printf("\n"); 
	}
	printf("\n");
    for (i=0; i<maxelm; i++) if (boundary[Temp][i]) 
		printf("%d %.2f %.2f %.2f %.2f\n", i+1, potent[Vx][i], potent[Vy][i], potent[Press][i], potent[Temp][i]);
    printf("\n");
    for (j=0; j<ncell; j++) { 
	    for (i=0; i<nve; i++) printf("%d ",nvtxcell[i][j]);
        printf("\n");
	}
	printf("\n");
	//system("pause"); system("pause"); system("pause"); 
} // print_task


// Файл передаётся из программы DavisTest Delphi написанной на Паскале.
// Файл имеет жёсткую структуру:
// секция координат узлов 
// номер узла начиная с единицы, координата x, координата y.
// пустая строка
// секция : ячеек контрольных объёмов.
// номер контрольного объёма начиная с единицы, четыре значения nvtx, 
// восемь значений соседей E,N,W,S,EE,NN,WW,SS ноль если сосед отсутствует,
// плотность, удельная теплоёмкость при постоянном давлении, теплопроводность,
// динамическая вязкость, коэффициент линейного температурного расширения.
// пустая строка
// секция : граничные условия.
// номер узла начиная с единицы, тип условия для Vx: однородный Нейман или Дирихле, 0 если однородный
// Нейман для Vx, значение Дирихле для Vx, тип условия для Vy: однородный Нейман или Дирихле, 0 если однородный
// Нейман для Vy, значение Дирихле для Vy, тип условия для Давления: однородный Нейман или Дирихле, 0 если однородный
// Нейман для Давления, значение Дирихле для Давления, тип условия для Температуры: однородный Нейман или Дирихле, 0 если однородный
// Нейман для Температуры, значение Дирихле для Температуры, два значения для двух внешних номалей E,W,N, или S.
// пустая строка
// секция : четыре целочисленных значения для nvtx_cell в каждой строке - расчётная сетка для Tecplot.
// переход к пустой строке и конец файла
void meshin(char *fname, int &nve, int maxsos, int maxval, int &maxnod, int &maxelm,
			Real* &x, Real* &y, int** &nvtx, bool** &boundary,
			Real** &potent, Real**& sumanb, Real** &Flux_gran, Real** & Flux_gran_relx, Real** &tau, Real** &prop, int** &sosed,
			bool** &neiman, int** &norm, int** &nvtxcell, int &ncell) {
	//char *fname;
	//fname="source.txt";
    FILE *fp;
	errno_t err1;
	if ((err1=fopen_s(&fp,fname,"r"))!=0) {
		printf("No File\n");
	}
	else
	{
		 bool nelf=false;
         int iNc=0; // количество непустых строк в файле
		 int in[4] = {0,0,0,0}; //длина секций
		 int i=0; // счётчик

         // часть 1: подсчёт длин секций.
         do{
			 nelf=false;
             while(fgetc(fp)!='\n' && !feof(fp)) nelf=true; // строка не пустая
             if(nelf) iNc++;
			 else {
				 if ((iNc!=0) && (i<4)) in[i++]=iNc;
				 iNc=0;
			 }
         } while(!feof(fp));

		 nve=4; // прямоугольные ячейки
		 maxsos=8;
         maxnod=in[0]; // число узлов 
		 maxelm=in[1]; // число контрольных объёмов 
         ncell=in[3]; // число связей контрольных объёмов для графической визуализации 
         my_malloc(nve,maxsos,maxval,maxnod,maxelm,nvtx,boundary,potent, sumanb, Flux_gran, Flux_gran_relx, tau, x,y,prop,sosed,neiman,norm,nvtxcell,ncell); // выделение оперативной памяти

         rewind(fp); // переход к началу файла.
         //fclose(fp);
         //fopen_s(&fp,fname,"r");
         
		 //printf("%d %d %d\n",in[0],in[1],in[2]);
		 //system("pause"); system("pause");
		 // часть 2: считывание секций.

		 // часть 2.1: координаты узлов.
         float fin=0.0;
		 for (i=0; i<maxnod; i++) {
             fscanf_s(fp, "%d", &in[0]); // номер узла
			 fscanf_s(fp, "%f", &fin); // x координата
			 x[in[0]-1]=fin;
             fscanf_s(fp, "%f", &fin); // y координата 
			 y[in[0]-1]=fin;
		 }
		 //for (i = 0; i < maxnod; i++) {
			// x[i] += 2.5;// zmeevik
			// y[i] += 2.5;
		// }
		 // часть 2.2: прямоугольники и свойства материала
		 for (i=0; i<maxelm; i++) {
			 int din=0;
             fscanf_s(fp, "%d", &din); // номер контрольного объёма
             fscanf_s(fp, "%d", &din);
             nvtx[0][i]=din;
             fscanf_s(fp, "%d", &din);
             nvtx[1][i]=din;
			 fscanf_s(fp, "%d", &din);
             nvtx[2][i]=din;
             fscanf_s(fp, "%d", &din);
			 nvtx[3][i]=din;

             fscanf_s(fp, "%d", &din);
		     sosed[E][i]=din;
             fscanf_s(fp, "%d", &din);
		     sosed[N][i]=din;
			 fscanf_s(fp, "%d", &din);
		     sosed[W][i]=din;
             fscanf_s(fp, "%d", &din);
			 sosed[S][i]=din;
             fscanf_s(fp, "%d", &din);
			 sosed[EE][i]=din;
             fscanf_s(fp, "%d", &din);
			 sosed[NN][i]=din;
			 fscanf_s(fp, "%d", &din);
			 sosed[WW][i]=din;
             fscanf_s(fp, "%d", &din);
			 sosed[SS][i]=din;

			 fscanf_s(fp, "%f", &fin);
			 prop[Rho][i]=fin;
             fscanf_s(fp, "%f", &fin);
             prop[Cp][i]=fin;
			 Cp_active = fin;
			 fscanf_s(fp, "%f", &fin);
             prop[Lam][i]=fin;
			 fscanf_s(fp, "%f", &fin);
             prop[Mu][i]=fin;
			 fscanf_s(fp, "%f", &fin);
             prop[BETA_T][i]=fin;
			 //std::cout << prop[Rho][i] << " " << prop[Cp][i] << " " << prop[Lam][i] << " " << prop[Mu][i] << " " << prop[BETA_T][i] << std::endl;
			 //getchar();

		 }
		 
		 // часть 2.3: граничные условия
	     for (i=0; i<maxelm; i++) {
			 int j=0;
			 for (j=0; j<maxval; j++) {
				 boundary[j][i]=false; // инициализация
		         potent[j][i]=0.0; // инициализация
			     neiman[j][i]=false; // по умолчанию все условия Дирихле.
			 }
			 for (j=0; j<2; j++) norm[j][i]=-1; // инициализация внутренний узел
             potent[Vxcor][i]=0.0; // инициализация
             potent[Vycor][i]=0.0;
	     } 
		 for (i=0; i<in[2]; i++) {
             fscanf_s(fp, "%d", &in[0]); // номер узла
             int j=0;
			 for (j=0; j<maxval; j++) boundary[j][in[0]-1]=true; // этот контрольный объём принадлежит границе области
			 int din;
			 // Vx
			 fscanf_s(fp, "%d", &din);
			 if (din==0) neiman[Vx][in[0]-1]=true;
             fscanf_s(fp, "%f", &fin); // условие Дирихле
			 if (din==1) {
				 potent[Vx][in[0]-1]=fin;
				 // нулевое поле скорости удовлетворяет уравнению неразрывности
				 // но не удовлетворяет краевым условиям.
                 potent[Vxcor][in[0]-1]=fin; 
			 }
             // Vy
			 fscanf_s(fp, "%d", &din);
			 if (din==0) neiman[Vy][in[0]-1]=true;
             fscanf_s(fp, "%f", &fin); // условие Дирихле
			 if (din==1) {
				 potent[Vy][in[0]-1]=fin;
				 // нулевое поле скорости удовлетворяет уравнению неразрывности
				 // но не удовлетворяет краевым условиям.
                 potent[Vycor][in[0]-1]=fin; 
			 }
			 // Press, PAm
             fscanf_s(fp, "%d", &din);
			 if (din==0) { 
				 neiman[Press][in[0]-1]=true;
				 neiman[PAm][in[0]-1]=true;
			 }
             fscanf_s(fp, "%f", &fin); // условие Дирихле
			 if (din==1) {
				 potent[Press][in[0]-1]=fin;
                 potent[PAm][in[0]-1]=0.0;
			 }
			 // Temp
             fscanf_s(fp, "%d", &din);
			 if (din==0) neiman[Temp][in[0]-1]=true;
             fscanf_s(fp, "%f", &fin); // условие Дирихле
             if (din==1) potent[Temp][in[0]-1]=fin;
			 // внешние нормали
             fscanf_s(fp, "%d", &din);
			 norm[0][in[0]-1]=din-1;
			 fscanf_s(fp, "%d", &din);
			 norm[1][in[0]-1]=din-1;
		 }

		 // часть 2.4 : определение сетки для техплота.
		 for (i=0; i<ncell; i++) {
			 int din=0;
             fscanf_s(fp, "%d", &din);
             nvtxcell[0][i]=din;
             fscanf_s(fp, "%d", &din);
             nvtxcell[1][i]=din;
			 fscanf_s(fp, "%d", &din);
             nvtxcell[2][i]=din;
             fscanf_s(fp, "%d", &din);
			 nvtxcell[3][i]=din;
		 }

         fclose(fp); // закрытие файла
		 if (maxelm<1000) print_task(x, y, maxnod,maxelm, nvtx, prop, nve, boundary, potent,maxsos, sosed, nvtxcell,ncell);
	}
} // meshin

// Считываает настройки солвера.
void loadcas(char *fname, int &maxit, int &maxts, int &isolver, int &isheme, 
			 bool &btimedep, Real &tau, Real &dgx, Real &dgy, Real &temp_ref) {
    //fname="solver.txt";
    FILE *fp;
	errno_t err1;
	if ((err1=fopen_s(&fp,fname,"r"))!=0) {
		printf("No File\n");
	}
	else
	{
		// default:
		btimedep=false; tau=0.1;
		maxit=60; maxts=20; isolver=BICGSTABCRS; isheme=UDS;

        int din=0;
		float fin=(float)1e-1;
        fscanf_s(fp, "%f", &fin);
		dgx=(Real)fin; // ускорение свободного падения по оси Ox
        fscanf_s(fp, "%f", &fin);
		dgy=(Real)fin; // ускорение свободного падения по оси Oy
		fscanf_s(fp, "%f", &fin);
		temp_ref=(Real)fin; // опорное значение темпратуры
		//printf("gx=%f, gy=%f, ref_temp=%f\n",dgx,dgy,temp_ref);
		//system("pause"); // debug
        fscanf_s(fp, "%d", &din);
		//if (din==1) btimedep=true;
		fin=(float)0.1;
        fscanf_s(fp, "%f", &fin);
		tau=(Real)fin;
        fscanf_s(fp, "%d", &din);
		maxit = 2000;// din;
		fscanf_s(fp, "%d", &din);
		maxts=din;
		fscanf_s(fp, "%d", &din);
		isolver= SOLOVALGCRS;
		fscanf_s(fp, "%d", &din);
		isheme = UDS;// din;

		//printf("%d %d %d %d /n",maxit, maxts, isolver, isheme);
		//system("pause");

	}
} // loadcas



#endif