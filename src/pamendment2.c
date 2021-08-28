// pamendment2.c 
// текущий рабочий вариант 
// файла в котором составляется дискретный
// аналог для уравнения поправки давления.

#pragma once
#ifndef PAMENDMENTv_0_02_C
#define PAMENDMENTv_0_02_C 1



// аппроксимация обобщённого уравнения конвекции-диффузии
// на совмещённой сетке
#include "my_elmatr_quad_f.c"


// коррекция скорости
void correct(int iP, equation** &sl, int iVar, 
			 int** &nvtx, bool** &boundary, Real** &potent, Real** &sumanb, Real** &tau,
			 Real* &x, Real* &y, int** &sosed, bool** &neiman,
			 int** &norm, Real** &prop, int nve, Real* &alpha, Real** &B, 
	Real**  & Flux_gran_relx, Real**& Flux_gran) {

    // iP - номер центрального контрольного объёма
	int iE, iN, iW, iS; // номера соседних контрольных объёмов
	iE = qass[iP].iE; iN = qass[iP].iN; iW = qass[iP].iW; iS = qass[iP].iS;

    if (!(boundary[iVar][iP]) && (!neiman[iVar][iP])) {   

        // условие Неймана 
		// в случае если переменная равна true
		bool bE=false, bN=false, bW=false, bS=false;
        if ((boundary[iVar][iP]) && (neiman[iVar][iP])) {

			// внешние нормали
			switch (norm[0][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // первая нормаль

			// внешние нормали
            switch (norm[1][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // вторая нормаль
		}    

		
	    // вычисление размеров текущего контрольного объёма:
	
	    Real dx= qass[iP].dx, dy= qass[iP].dy; // размеры контрольного объёма
		//volume(iP, nve, nvtx, x, y, dx, dy);
		//printf("%.2f %.2f\n",dx,dy); // debug GOOD
		//system("pause");

        Real dxe=1.0, dxw=1.0, dyn=1.0, dys=1.0;
	    //if ((!bE)&&(iE>-1)) dxe=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
		if ((!bE) && (iE > -1)) {
			//dxe-=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
			dxe = qass[iP].dxe;
		}
		else { iE=iP; dxe=dx; }
	    //if ((!bW) && (iW > -1)) dxw=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
		if ((!bW) && (iW > -1)) {
			//dxw-=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
			dxw = qass[iP].dxw;
		}
		else { iW=iP; dxw=dx; }
	    //if ((!bN) && (iN > -1)) dyn=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
		if ((!bN) && (iN > -1)) {
			//dyn-=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
			dyn = qass[iP].dyn;
		}
		else { iN=iP; dyn=dy; }
	    //if ((!bS) && (iS > -1)) dys=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
		if ((!bS) && (iS > -1)) {
			//dys-=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
			dys = qass[iP].dys;
		}
		else { iS=iP; dys=dy; }

	    Real feplus, fwplus, fnplus, fsplus;
	    feplus=0.5*dx/dxe;
	    fwplus=0.5*dx/dxw;
	    fnplus=0.5*dy/dyn;
	    fsplus=0.5*dy/dys;

		// плотность аппроксимируется средним гармоническим
		Real rhoe, rhow, rhon, rhos;
		Real rP, rE, rN, rW, rS;

		rP = prop[Rho][iP];
		if ((!bE) && (iE > -1)) rE = prop[Rho][iE]; else rE = rP;
		if ((!bN) && (iN > -1)) rN = prop[Rho][iN]; else rN = rP;
		if ((!bW) && (iW > -1)) rW = prop[Rho][iW]; else rW = rP;
		if ((!bS) && (iS > -1)) rS = prop[Rho][iS]; else rS = rP;

		rhoe = rE * rP / ((1 - feplus) * rP + feplus * rE);
		rhow = rW * rP / ((1 - fwplus) * rP + fwplus * rW);
		rhon = rN * rP / ((1 - fnplus) * rP + fnplus * rN);
		rhos = rS * rP / ((1 - fsplus) * rP + fsplus * rS);



		Real dl=dx;
		if (iVar==Vx) dl=dy;
		Real dlr = dy;
		if (iVar == Vx) dlr = dx;


		//int i; // счётчик цикла for
		
		//Real** B = new Real * [3];
		//for (int l = 0; l < 3; l++) B[l] = new Real[3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = 0.0;

		Real PAmP, PAmW, PAmE, PAmS, PAmN;
		PAmP=potent[PAm][iP];
        if ((!bW) && (iW > -1)) PAmW=potent[PAm][iW];
		else { // квадратичная экстраполяция
               // узла iW нету.
				// надо восстановить давление PAmW в узле iW с помощью
				// квадратичной экстраполяции
                Real xP, xE, xEE, xP2, xE2, xEE2, xW;
				int iEE=sosed[EE][iP]-1;
				// при условии отсутствия узла iW узел iEE всегда существует.
                xEE=0.25*(x[nvtx[0][iEE]-1]+x[nvtx[1][iEE]-1]+x[nvtx[2][iEE]-1]+x[nvtx[3][iEE]-1]);
				xE=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xP2=xP*xP; xE2=xE*xE;  xEE2=xEE*xEE;
				
				B[0][0]=1.0; B[0][1]=xP; B[0][2]=xP2;
				B[1][0]=1.0; B[1][1]=xE; B[1][2]=xE2;
				B[2][0]=1.0; B[2][1]=xEE; B[2][2]=xEE2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				xW=xP-dxw;
				Real  PAmPloc, PAmEloc, PAmEEloc;
				 PAmPloc=potent[PAm][iP]; PAmEloc=potent[PAm][iE]; PAmEEloc=potent[PAm][iEE];
				PAmW=1.0*(B[0][0]*PAmPloc+B[0][1]*PAmEloc+B[0][2]*PAmEEloc);
				PAmW+=xW*(B[1][0]*PAmPloc+B[1][1]*PAmEloc+B[1][2]*PAmEEloc);
				PAmW+=xW*xW*(B[2][0]*PAmPloc+B[2][1]*PAmEloc+B[2][2]*PAmEEloc); 
		}
        if ((!bE) && (iE > -1)) PAmE=potent[PAm][iE];
		else { // квадратичная экстраполяция
                // узла iE нету.
				// надо восстановить давление PAmE в узле iE с помощью
				// квадратичной экстраполяции
                Real xP, xW, xWW, xP2, xW2, xWW2, xE;
                // при условии отсутствия узла iE узел iWW всегда существует.
				int iWW = sosed[WW][iP]-1;
                xWW=0.25*(x[nvtx[0][iWW]-1]+x[nvtx[1][iWW]-1]+x[nvtx[2][iWW]-1]+x[nvtx[3][iWW]-1]);
				xW=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xP2=xP*xP; xW2=xW*xW;  xWW2=xWW*xWW;
				
				B[0][0]=1.0; B[0][1]=xWW; B[0][2]=xWW2;
				B[1][0]=1.0; B[1][1]=xW; B[1][2]=xW2;
				B[2][0]=1.0; B[2][1]=xP; B[2][2]=xP2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				xE=xP+dxe;
				Real PAmWWloc, PAmWloc, PAmPloc;
			    PAmWWloc=potent[PAm][iWW]; PAmWloc=potent[PAm][iW]; PAmPloc=potent[PAm][iP];
				PAmE=1.0*(B[0][0]*PAmWWloc+B[0][1]*PAmWloc+B[0][2]*PAmPloc);
				PAmE+=xE*(B[1][0]*PAmWWloc+B[1][1]*PAmWloc+B[1][2]*PAmPloc);
				PAmE+=xE*xE*(B[2][0]*PAmWWloc+B[2][1]*PAmWloc+B[2][2]*PAmPloc);
		}
		if ((!bS) && (iS > -1)) PAmS=potent[PAm][iS];
		else { // квадратичная экстраполяция
                // узла iS нету.
				// надо восстановить поправку давления PAmS в узле iS с помощью
				// квадратичной экстраполяции
                Real yP, yN, yNN, yP2, yN2, yNN2, yS;
				int iNN=sosed[NN][iP]-1;
				// при условии отсутствия узла iS узел iNN всегда существует.
                yNN=0.25*(y[nvtx[0][iNN]-1]+y[nvtx[1][iNN]-1]+y[nvtx[2][iNN]-1]+y[nvtx[3][iNN]-1]);
				yN=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				yP2=yP*yP; yN2=yN*yN;  yNN2=yNN*yNN;
				
				B[0][0]=1.0; B[0][1]=yP; B[0][2]=yP2;
				B[1][0]=1.0; B[1][1]=yN; B[1][2]=yN2;
				B[2][0]=1.0; B[2][1]=yNN; B[2][2]=yNN2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				yS=yP-dys;
				Real  PAmPloc, PAmNloc, PAmNNloc;
				PAmPloc=potent[PAm][iP]; PAmNloc=potent[PAm][iN]; PAmNNloc=potent[PAm][iNN];
				PAmS=1.0*(B[0][0]*PAmPloc+B[0][1]*PAmNloc+B[0][2]*PAmNNloc);
				PAmS+=yS*(B[1][0]*PAmPloc+B[1][1]*PAmNloc+B[1][2]*PAmNNloc);
				PAmS+=yS*yS*(B[2][0]*PAmPloc+B[2][1]*PAmNloc+B[2][2]*PAmNNloc);
		}
		if ((!bN) && (iN > -1)) PAmN=potent[PAm][iN];
		else { // квадратичная экстраполяция
                // узла iN нету.
				// надо восстановить давление PN в узле iN с помощью
				// квадратичной экстраполяции
                Real yP, yS, ySS, yP2, yS2, ySS2, yN;
                // при условии отсутствия узла iN узел iSS всегда существует.
				int iSS=sosed[SS][iP]-1;
                ySS=0.25*(y[nvtx[0][iSS]-1]+y[nvtx[1][iSS]-1]+y[nvtx[2][iSS]-1]+y[nvtx[3][iSS]-1]);
				yS=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				yP2=yP*yP; yS2=yS*yS;  ySS2=ySS*ySS;
				
				B[0][0]=1.0; B[0][1]=ySS; B[0][2]=ySS2;
				B[1][0]=1.0; B[1][1]=yS; B[1][2]=yS2;
				B[2][0]=1.0; B[2][1]=yP; B[2][2]=yP2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				yN=yP+dyn;
				Real PAmSSloc, PAmSloc, PAmPloc;
			    PAmSSloc=potent[PAm][iSS]; PAmSloc=potent[PAm][iS]; PAmPloc=potent[PAm][iP];
				PAmN=1.0*(B[0][0]*PAmSSloc+B[0][1]*PAmSloc+B[0][2]*PAmPloc);
				PAmN+=yN*(B[1][0]*PAmSSloc+B[1][1]*PAmSloc+B[1][2]*PAmPloc);
				PAmN+=yN*yN*(B[2][0]*PAmSSloc+B[2][1]*PAmSloc+B[2][2]*PAmPloc); 
		}
		
		//for (int l = 0; l < 3; l++) delete B[l];
		//delete[] B;

		Real deltaP=0.0;
		Real deltaP2 = 0.0;
	
		switch (iVar) {
			case Vx : 
				if ((!bW) && (!bE) && (iW > -1) && (iE > -1)) {
					deltaP2 = ((1.0 - fwplus) * PAmP + fwplus * PAmW);
					deltaP2 -= (feplus * PAmE + (1 - feplus) * PAmP);
				}
					  break;
			case Vy : 
				if ((!bS) && (!bN) && (iS > -1) && (iN > -1)) {
					deltaP2 = ((1.0 - fsplus) * PAmP + (fsplus)*PAmS);
					deltaP2 -= (fnplus * PAmN + (1 - fnplus) * PAmP);
				}
				      break;
		}
		/*switch (iVar) {
		case Vx:
			if ((!bW) && (!bE) && (iW > -1) && (iE > -1)) {
				deltaP = tau[W][iP]*((1.0 - fwplus) * PAmP + fwplus * PAmW)/rhow;
				deltaP -= tau[E][iP]*(feplus * PAmE + (1 - feplus) * PAmP)/rhoe;
			}
			break;
		case Vy:
			if ((!bS) && (!bN) && (iS > -1) && (iN > -1)) {
				deltaP = tau[S][iP] * ((1.0 - fsplus) * PAmP + (fsplus)*PAmS)/rhos;
				deltaP -= tau[N][iP] * (fnplus * PAmN + (1 - fnplus) * PAmP)/rhon;
			}
			break;
		}*/

		//de = alpha[iVar] * dy / apue;
		//Real taue = rhoe * de * dy;
		// taue*deltaP/(dx*rhoe);
		//alpha[iVar] * dy * dy*deltaP/(dx*apue);

		// коррекция скорости не должна подвергаться нижней релаксации.
		potent[iVar][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];// sl[iVar][iP].ap;//alpha[iVar]*

		/*
		// Портит конфузор !!!
		if (iVar == Vx) {
			if ((!bW) && (!bE) && (iW > -1) && (iE > -1)) {
				// Скорректируем также и поток жидкости на грани контрольного объёма.
				Flux_gran_relx[E][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
				Flux_gran_relx[W][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
				Flux_gran[E][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
				Flux_gran[W][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
			}
		}
		if (iVar == Vy) {
			if ((!bS) && (!bN) && (iS > -1) && (iN > -1)) {
				// Скорректируем также и поток жидкости на грани контрольного объёма.
				Flux_gran_relx[N][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
				Flux_gran_relx[S][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
				Flux_gran[N][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
				Flux_gran[S][iP] += alpha[iVar] * dl * (deltaP2) / sumanb[iVar][iP];
			}
		}*/
	
		/*switch (iVar) {
		case Vx: 
			if ((fabs(tau[E][iP]) > 1.0e-30) && (fabs(tau[W][iP]) > 1.0e-30)) {
				potent[iVar][iP] += deltaP / dx;
			}
			else {
				potent[iVar][iP] += alpha[iVar] * dl * (deltaP) / sumanb[iVar][iP];
			}
			break;
		case Vy: 
			if ((fabs(tau[N][iP]) > 1.0e-30) && (fabs(tau[S][iP]) > 1.0e-30)) {
				potent[iVar][iP] += deltaP / dy;
			}
			else {
				potent[iVar][iP] += alpha[iVar] * dl * (deltaP) / sumanb[iVar][iP];
			}
			break;
		}*/

		// коррекция скорости не должна подвергаться нижней релаксации.
		//--->potent[iVar][iP] += alpha[iVar] * dl * dl *  (deltaP) / (dlr * sumanb[iVar][iP]);

		//potent[iVar][iP] += (deltaP) / (dlr);
}
    
} // correct


void calc_tau(int iP, equation**& sl, int**& nvtx, bool**& boundary,
	Real**& potent, Real**& sumanb, Real**& Flux_gran, Real**& Flux_gran_relx,
	Real**& tau, Real*& x, Real*& y, Real**& prop,
	int**& sosed, bool**& neiman, int**& norm, int nve,
	Real*& alpha, Real**& B)
{
	if (!((boundary[PAm][iP]) && (!neiman[PAm][iP]))) {

		Real eps = 1e-37; // для отделения вещественного нуля.

	// iP - номер центрального контрольного объёма
		int iE, iN, iW, iS; // номера соседних контрольных объёмов
		iE = qass[iP].iE; iN = qass[iP].iN; iW = qass[iP].iW; iS = qass[iP].iS;
		sl[PAm][iP].iE = iE; sl[PAm][iP].iN = iN; sl[PAm][iP].iP = iP; sl[PAm][iP].iS = iS; sl[PAm][iP].iW = iW;


		// либо внутренний узел, либо условие Неймана.

		// условие Неймана 
		// в случае если переменная равна true
		bool bE = false, bN = false, bW = false, bS = false;
		if ((boundary[PAm][iP]) && (neiman[PAm][iP])) {
			// внешние нормали
			switch (norm[0][iP]) {
			case E: bE = true; break;
			case N: bN = true; break;
			case W: bW = true; break;
			case S: bS = true; break;
			} // первая нормаль

			// внешние нормали
			switch (norm[1][iP]) {
			case E: bE = true; break;
			case N: bN = true; break;
			case W: bW = true; break;
			case S: bS = true; break;
			} // вторая нормаль
		}


		// вычисление размеров текущего контрольного объёма:
		Real dx = qass[iP].dx, dy = qass[iP].dy; // размеры контрольного объёма
		//volume(iP, nve, nvtx, x, y, dx, dy);


		Real dxe = qass[iP].dxe, dxw = qass[iP].dxw, dyn = qass[iP].dyn, dys = qass[iP].dys;
		//if ((!bE) && (iE > -1)) dxe = 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
		//if ((!bE) && (iE > -1)) dxe -= 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]); else dxe = dx;
		//if ((!bW) && (iW > -1)) dxw = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
		//if ((!bW) && (iW > -1)) dxw -= 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]); else dxw = dx;
		//if ((!bN) && (iN > -1)) dyn = 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
		//if ((!bN) && (iN > -1)) dyn -= 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]); else dyn = dy;
		//if ((!bS)&&(iS>-1)) dys = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
		//if ((!bS) && (iS > -1)) dys -= 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]); else dys = dy;

		Real feplus, fwplus, fnplus, fsplus;
		feplus = 0.5 * dx / dxe;
		fwplus = 0.5 * dx / dxw;
		fnplus = 0.5 * dy / dyn;
		fsplus = 0.5 * dy / dys;

		//printf("%.2f %.2f %.2f %.2f\n",feplus, fwplus, fnplus, fsplus);
		//system("pause");

		Real apue = 1.0, apuw = 1.0, apvn = 1.0, apvs = 1.0;

		// Подлежит удалению.
		//if (!bE) apue = sl[Vx][iE].ap * sl[Vx][iP].ap / ((1 - feplus) * sl[Vx][iP].ap + feplus * sl[Vx][iE].ap); else apue = sl[Vx][iP].ap;
		//if (!bW) apuw = sl[Vx][iW].ap * sl[Vx][iP].ap / ((1 - fwplus) * sl[Vx][iP].ap + fwplus * sl[Vx][iW].ap); else apuw = sl[Vx][iP].ap;
		//if (!bN) apvn = sl[Vy][iN].ap * sl[Vy][iP].ap / ((1 - fnplus) * sl[Vy][iP].ap + fnplus * sl[Vy][iN].ap); else apvn = sl[Vy][iP].ap;
		//if (!bS) apvs = sl[Vy][iS].ap * sl[Vy][iP].ap / ((1 - fsplus) * sl[Vy][iP].ap + fsplus * sl[Vy][iS].ap); else apvs = sl[Vy][iP].ap;

		if ((!bE) && (iE > -1)) apue = sumanb[Vx][iE] * sumanb[Vx][iP] / ((1 - feplus) * sumanb[Vx][iP] + feplus * sumanb[Vx][iE]); else apue = sumanb[Vx][iP];
		if ((!bW) && (iW > -1)) apuw = sumanb[Vx][iW] * sumanb[Vx][iP] / ((1 - fwplus) * sumanb[Vx][iP] + fwplus * sumanb[Vx][iW]); else apuw = sumanb[Vx][iP];
		if ((!bN) && (iN > -1)) apvn = sumanb[Vy][iN] * sumanb[Vy][iP] / ((1 - fnplus) * sumanb[Vy][iP] + fnplus * sumanb[Vy][iN]); else apvn = sumanb[Vy][iP];
		if ((!bS) && (iS > -1)) apvs = sumanb[Vy][iS] * sumanb[Vy][iP] / ((1 - fsplus) * sumanb[Vy][iP] + fsplus * sumanb[Vy][iS]); else apvs = sumanb[Vy][iP];




		Real de, dw, dn, ds;
		de = alpha[Vx] * dy / apue; dw = alpha[Vx] * dy / apuw;
		dn = alpha[Vy] * dx / apvn; ds = alpha[Vy] * dx / apvs;

		// плотность аппроксимируется средним гармоническим
		Real rhoe, rhow, rhon, rhos;
		Real rP, rE, rN, rW, rS;

		rP = prop[Rho][iP];
		if ((!bE) && (iE > -1)) rE = prop[Rho][iE]; else rE = rP;
		if ((!bN) && (iN > -1)) rN = prop[Rho][iN]; else rN = rP;
		if ((!bW) && (iW > -1)) rW = prop[Rho][iW]; else rW = rP;
		if ((!bS) && (iS > -1)) rS = prop[Rho][iS]; else rS = rP;

		rhoe = rE * rP / ((1 - feplus) * rP + feplus * rE);
		rhow = rW * rP / ((1 - fwplus) * rP + fwplus * rW);
		rhon = rN * rP / ((1 - fnplus) * rP + fnplus * rN);
		rhos = rS * rP / ((1 - fsplus) * rP + fsplus * rS);

		tau[E][iP] = (rhoe * de * dy);
		tau[W][iP] = (rhow * dw * dy);
		tau[N][iP] = (rhon * dn * dx);
		tau[S][iP] = (rhos * ds * dx);
	}
	else {
		tau[E][iP] = 0.0;
		tau[W][iP] = 0.0;
		tau[N][iP] = 0.0;
		tau[S][iP] = 0.0;
	}
}

// Составляет матрицу для уравнения 
// поправки давления
void my_elmatr_quad_PAm(int iP, equation**& sl, int** &nvtx, bool** &boundary,
	Real** &potent, Real** &sumanb, Real** &Flux_gran, Real** & Flux_gran_relx, 
	Real** &tau, Real* &x, Real* &y, Real** &prop,
	int** &sosed, bool** &neiman, int** &norm, int nve,
	Real* &alpha, Real** &B) {

	Real eps = 1e-37; // для отделения вещественного нуля.

	// iP - номер центрального контрольного объёма
	int iE, iN, iW, iS; // номера соседних контрольных объёмов
	iE = qass[iP].iE; iN = qass[iP].iN; iW = qass[iP].iW; iS = qass[iP].iS;
	sl[PAm][iP].iE = iE; sl[PAm][iP].iN = iN; sl[PAm][iP].iP = iP; sl[PAm][iP].iS = iS; sl[PAm][iP].iW = iW;

	
		// либо внутренний узел, либо условие Неймана.

		// условие Неймана 
		// в случае если переменная равна true
		bool bE = false, bN = false, bW = false, bS = false;
		if ((boundary[PAm][iP]) && (neiman[PAm][iP])) {
			// внешние нормали
			switch (norm[0][iP]) {
			case E: bE = true; break;
			case N: bN = true; break;
			case W: bW = true; break;
			case S: bS = true; break;
			} // первая нормаль

			// внешние нормали
			switch (norm[1][iP]) {
			case E: bE = true; break;
			case N: bN = true; break;
			case W: bW = true; break;
			case S: bS = true; break;
			} // вторая нормаль
		}


		// вычисление размеров текущего контрольного объёма:
		Real dx = qass[iP].dx, dy = qass[iP].dy; // размеры контрольного объёма
		//volume(iP, nve, nvtx, x, y, dx, dy);


		Real dxe = qass[iP].dxe, dxw = qass[iP].dxw, dyn = qass[iP].dyn, dys = qass[iP].dys;
		//if ((!bE) && (iE > -1)) dxe = 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
		//if ((!bE) && (iE > -1)) dxe -= 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]); else dxe = dx;
		//if ((!bW) && (iW > -1)) dxw = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
		//if ((!bW) && (iW > -1)) dxw -= 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]); else dxw = dx;
		//if ((!bN) && (iN > -1)) dyn = 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
		//if ((!bN) && (iN > -1)) dyn -= 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]); else dyn = dy;
		//if ((!bS)&&(iS>-1)) dys = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
		//if ((!bS) && (iS > -1)) dys -= 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]); else dys = dy;

		Real feplus, fwplus, fnplus, fsplus;
		feplus = 0.5 * dx / dxe;
		fwplus = 0.5 * dx / dxw;
		fnplus = 0.5 * dy / dyn;
		fsplus = 0.5 * dy / dys;

		//printf("%.2f %.2f %.2f %.2f\n",feplus, fwplus, fnplus, fsplus);
		//system("pause");

		Real apue = 1.0, apuw = 1.0, apvn = 1.0, apvs = 1.0;

		// Подлежит удалению.
		//if (!bE) apue = sl[Vx][iE].ap * sl[Vx][iP].ap / ((1 - feplus) * sl[Vx][iP].ap + feplus * sl[Vx][iE].ap); else apue = sl[Vx][iP].ap;
		//if (!bW) apuw = sl[Vx][iW].ap * sl[Vx][iP].ap / ((1 - fwplus) * sl[Vx][iP].ap + fwplus * sl[Vx][iW].ap); else apuw = sl[Vx][iP].ap;
		//if (!bN) apvn = sl[Vy][iN].ap * sl[Vy][iP].ap / ((1 - fnplus) * sl[Vy][iP].ap + fnplus * sl[Vy][iN].ap); else apvn = sl[Vy][iP].ap;
		//if (!bS) apvs = sl[Vy][iS].ap * sl[Vy][iP].ap / ((1 - fsplus) * sl[Vy][iP].ap + fsplus * sl[Vy][iS].ap); else apvs = sl[Vy][iP].ap;

		if ((!bE) && (iE > -1)) apue = sumanb[Vx][iE] * sumanb[Vx][iP] / ((1 - feplus) * sumanb[Vx][iP] + feplus * sumanb[Vx][iE]); else apue = sumanb[Vx][iP];
		if ((!bW) && (iW > -1)) apuw = sumanb[Vx][iW] * sumanb[Vx][iP] / ((1 - fwplus) * sumanb[Vx][iP] + fwplus * sumanb[Vx][iW]); else apuw = sumanb[Vx][iP];
		if ((!bN) && (iN > -1)) apvn = sumanb[Vy][iN] * sumanb[Vy][iP] / ((1 - fnplus) * sumanb[Vy][iP] + fnplus * sumanb[Vy][iN]); else apvn = sumanb[Vy][iP];
		if ((!bS) && (iS > -1)) apvs = sumanb[Vy][iS] * sumanb[Vy][iP] / ((1 - fsplus) * sumanb[Vy][iP] + fsplus * sumanb[Vy][iS]); else apvs = sumanb[Vy][iP];


		Real de, dw, dn, ds;
		de = alpha[Vx] * dy / apue; dw = alpha[Vx] * dy / apuw;
		dn = alpha[Vy] * dx / apvn; ds = alpha[Vy] * dx / apvs;

		// плотность аппроксимируется средним гармоническим
		Real rhoe, rhow, rhon, rhos;
		Real rP, rE, rN, rW, rS;

		rP = prop[Rho][iP];
		if ((!bE) && (iE > -1)) rE = prop[Rho][iE]; else rE = rP;
		if ((!bN) && (iN > -1)) rN = prop[Rho][iN]; else rN = rP;
		if ((!bW) && (iW > -1)) rW = prop[Rho][iW]; else rW = rP;
		if ((!bS) && (iS > -1)) rS = prop[Rho][iS]; else rS = rP;

		rhoe = rE * rP / ((1 - feplus) * rP + feplus * rE);
		rhow = rW * rP / ((1 - fwplus) * rP + fwplus * rW);
		rhon = rN * rP / ((1 - fnplus) * rP + fnplus * rN);
		rhos = rS * rP / ((1 - fsplus) * rP + fsplus * rS);

		
		Real mnull = 0.0;
		Real mnull2 = 0.0;

		Real Fw = 0.0, Fe = 0.0, Fs = 0.0, Fn = 0.0;
		Real Fw_cor = 0.0, Fe_cor = 0.0, Fs_cor = 0.0, Fn_cor = 0.0;
		Real FwRhie_Chow = 0.0, FeRhie_Chow = 0.0, FsRhie_Chow = 0.0, FnRhie_Chow = 0.0;

		//int l; // счётчик цикла for
		
		//Real** B = new Real * [3];
		//for (int l = 0; l < 3; l++) B[l] = new Real[3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = 0.0;

		if ((!bE) && (iE > -1)) {
			Fe = rhoe * dy * (feplus * potent[Vx][iE] + (1.0 - feplus) * potent[Vx][iP]);
			Fe_cor = rhoe * dy * (feplus * potent[Vxcor][iE] + (1.0 - feplus) * potent[Vxcor][iP]);
			

			Real koef = rhoe * dy * dy * alpha[Vx];
			int iEE = sosed[EE][iP] - 1;
			Real dxee = dxe;
			Real PEE = 0.0; // давление в узле iEE
			if (iEE > -1) {
				// если узел существует
				dxee = 0.25 * (x[nvtx[0][iEE] - 1] + x[nvtx[1][iEE] - 1] + x[nvtx[2][iEE] - 1] + x[nvtx[3][iEE] - 1]);
				dxee -= 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
				PEE = potent[Press][iEE];
			}
			else {
				// если узел несуществует 
				// недостающее Давление PEE
				// будет определено с помощью 
				// квадратичной экстраполяции.
				Real xE, xP, xW, xE2, xP2, xW2, xEE;
				xE = 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
				xP = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
				// при условии отсутствия узла iEE узел iW всегда существует.
				xW = 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]);
				xE2 = xE * xE; xP2 = xP * xP; xW2 = xW * xW;

				B[0][0] = 1.0; B[0][1] = xW; B[0][2] = xW2;
				B[1][0] = 1.0; B[1][1] = xP; B[1][2] = xP2;
				B[2][0] = 1.0; B[2][1] = xE; B[2][2] = xE2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				xEE = xE + dxe;
				Real PE, PP, PW;
				PE = potent[Press][iE]; PP = potent[Press][iP]; PW = potent[Press][iW];
				PEE = 1.0 * (B[0][0] * PW + B[0][1] * PP + B[0][2] * PE);
				PEE += xEE * (B[1][0] * PW + B[1][1] * PP + B[1][2] * PE);
				PEE += xEE * xEE * (B[2][0] * PW + B[2][1] * PP + B[2][2] * PE);
			}

			Real dxE = qass[iE].dx, dyE = qass[iE].dy;
			//volume(iE, nve, nvtx, x, y, dxE, dyE);
			Real feeplus = 0.5 * dxE / dxee, fpplus = 0.5 * dxE / dxe;

			if ((!bW) && (iW > -1)) {
				// iE
				FeRhie_Chow += koef * (feplus) * (feeplus * PEE + (1.0 - feeplus) * potent[Press][iE] 
					- (1.0 - fpplus) * potent[Press][iE] - fpplus * potent[Press][iP]) / (sumanb[Vx][iE]);
				// iP
				FeRhie_Chow += koef * (1.0 - feplus) * (feplus * potent[Press][iE] + (1.0 - feplus) * potent[Press][iP] 
					- (1.0 - fwplus) * potent[Press][iP] - fwplus * potent[Press][iW]) / (sumanb[Vx][iP]);
				// e
				FeRhie_Chow -= koef * (potent[Press][iE] - potent[Press][iP]) / apue;
			}
			else {
				// узла iW нету.
				// надо восстановить давление PW в узле iW с помощью
				// квадратичной экстраполяции
				Real xP, xE, xEE, xP2, xE2, xEE2, xW;
				xEE = 0.25 * (x[nvtx[0][iEE] - 1] + x[nvtx[1][iEE] - 1] + x[nvtx[2][iEE] - 1] + x[nvtx[3][iEE] - 1]);
				// при условии отсутствия узла iWW узел iE всегда существует.
				xE = 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
				xP = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
				xP2 = xP * xP; xE2 = xE * xE;  xEE2 = xEE * xEE;

				B[0][0] = 1.0; B[0][1] = xP; B[0][2] = xP2;
				B[1][0] = 1.0; B[1][1] = xE; B[1][2] = xE2;
				B[2][0] = 1.0; B[2][1] = xEE; B[2][2] = xEE2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				xW = xP - dxw;
				Real PE, PP, PW;
				PE = potent[Press][iE]; PP = potent[Press][iP];
				PW = 1.0 * (B[0][0] * PP + B[0][1] * PE + B[0][2] * PEE);
				PW += xW * (B[1][0] * PP + B[1][1] * PE + B[1][2] * PEE);
				PW += xW * xW * (B[2][0] * PP + B[2][1] * PE + B[2][2] * PEE);


				// неравномерная сетка
				// E
				FeRhie_Chow += koef * (feplus) * (feeplus * PEE + (1.0 - feeplus) * potent[Press][iE] 
					- fpplus * potent[Press][iE] - (1.0 - fpplus) * potent[Press][iP]) / (sumanb[Vx][iE]);
				// P
				FeRhie_Chow += koef * (1.0 - feplus) * (feplus * potent[Press][iE] + (1.0 - feplus) * potent[Press][iP] 
					- (1.0 - fwplus)  * potent[Press][iP] - fwplus * PW) / (sumanb[Vx][iP]);
				// e
				FeRhie_Chow -= koef * (potent[Press][iE] - potent[Press][iP]) / apue;
				
				FeRhie_Chow *= mnull;
			}

		}
		else {
			Fe = rhoe * dy * potent[Vx][iP];
			Fe_cor = rhoe * dy * potent[Vxcor][iP];

			// Если жидкость не течёт через стенку
			// то и нормального к стенке перепада давления там нет.
			// Вблизи твёрдой стенки предположим что
			// давление гладкая функция.
			// В дальнейшем это требуется проверить экспериментально. TODO
			//if (fabs(Fe) < eps) FeRhie_Chow = 0.0;
			//else
			{

				// узлов iE и iEE нету
				// давление в них определяется
				// с помощью квадратичной экстраполяции
				Real PEE = 0.0; // давление в узле iEE
				Real PE = 0.0; // давление в узле iE
				int iWW = sosed[WW][iP] - 1;
				Real xWW, xW, xP, xWW2, xW2, xP2, xE, xEE;
				// при условии отсутствия узла iE и iEE узел iWW всегда существует.
				xWW = 0.25 * (x[nvtx[0][iWW] - 1] + x[nvtx[1][iWW] - 1] + x[nvtx[2][iWW] - 1] + x[nvtx[3][iWW] - 1]);
				xW = 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]);
				xP = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
				xWW2 = xWW * xWW; xP2 = xP * xP; xW2 = xW * xW;

				B[0][0] = 1.0; B[0][1] = xWW; B[0][2] = xWW2;
				B[1][0] = 1.0; B[1][1] = xW; B[1][2] = xW2;
				B[2][0] = 1.0; B[2][1] = xP; B[2][2] = xP2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				Real dxee = dxe;
				xE = xP + dxe;
				xEE = xE + dxe;
				Real PWW, PW, PP;
				PWW = potent[Press][iWW]; PW = potent[Press][iW]; PP = potent[Press][iP];
				PE = 1.0 * (B[0][0] * PWW + B[0][1] * PW + B[0][2] * PP);
				PE += xE * (B[1][0] * PWW + B[1][1] * PW + B[1][2] * PP);
				PE += xE * xE * (B[2][0] * PWW + B[2][1] * PW + B[2][2] * PP);

				PEE = 1.0 * (B[0][0] * PWW + B[0][1] * PW + B[0][2] * PP);
				PEE += xEE * (B[1][0] * PWW + B[1][1] * PW + B[1][2] * PP);
				PEE += xEE * xEE * (B[2][0] * PWW + B[2][1] * PW + B[2][2] * PP);

				Real koef = rhoe * dy * dy * alpha[Vx];

				Real feeplus = 0.5, fpplus = 0.5;

				// неравномерная сетка
				FeRhie_Chow += koef * (feplus) * ((1.0 - feeplus) * PEE +  feeplus *PE 
				- fpplus * PE - (1.0 - fpplus) * potent[Press][iP]) / (sumanb[Vx][iP]);
				FeRhie_Chow += koef * (1.0 - feplus) * (feplus * PE + (1.0 - feplus) * potent[Press][iP] 
				- (1.0 - fwplus) * potent[Press][iP] - fwplus * PW) / (sumanb[Vx][iP]);
				FeRhie_Chow -= koef * (PE - potent[Press][iP]) / apue;
			}

			// Поток задан пользователем, какой смысл его коректировать.
			FeRhie_Chow *= mnull2;;
		}

		if ((!bW) && (iW > -1)) {
			Fw = rhow * dy * ((1.0 - fwplus) * potent[Vx][iP] + fwplus * potent[Vx][iW]);
			Fw_cor = rhow * dy * ((1.0 - fwplus) * potent[Vxcor][iP] + fwplus * potent[Vxcor][iW]);
			// Схема по потоку
			//if (Fw > 0.0) Fw = rhow * dy * potent[Vx][iP]; else Fw = rhow * dy * potent[Vx][iW];

			Real koef = rhow * dy * dy * alpha[Vx];
			int iWW = sosed[WW][iP] - 1;
			Real dxww = dxw;
			Real PWW = 0.0;
			if (iWW > -1) {
				// если узел существует
				dxww = 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]);
				dxww -= 0.25 * (x[nvtx[0][iWW] - 1] + x[nvtx[1][iWW] - 1] + x[nvtx[2][iWW] - 1] + x[nvtx[3][iWW] - 1]);
				PWW = potent[Press][iWW];
			}
			else {
				// если узел iWW не существует,
				// то недостающее давление PWW
				// будет определено с помощью 
				// квадратичной интерполляции
				Real xW, xP, xE, xW2, xP2, xE2, xWW;
				xW = 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]);
				// при условии отсутствия узла iWW узел iE всегда существует.
				xE = 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
				xP = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
				xE2 = xE * xE; xP2 = xP * xP; xW2 = xW * xW;

				B[0][0] = 1.0; B[0][1] = xW; B[0][2] = xW2;
				B[1][0] = 1.0; B[1][1] = xP; B[1][2] = xP2;
				B[2][0] = 1.0; B[2][1] = xE; B[2][2] = xE2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				xWW = xW - dxw;
				Real PE, PP, PW;
				PE = potent[Press][iE]; PP = potent[Press][iP]; PW = potent[Press][iW];
				PWW = 1.0 * (B[0][0] * PW + B[0][1] * PP + B[0][2] * PE);
				PWW += xWW * (B[1][0] * PW + B[1][1] * PP + B[1][2] * PE);
				PWW += xWW * xWW * (B[2][0] * PW + B[2][1] * PP + B[2][2] * PE);
			}

			Real dxW = qass[iW].dx, dyW = qass[iW].dy;
			//volume(iW, nve, nvtx, x, y, dxW, dyW);
			Real fpplus = 0.5 * dxW / dxw, fwwplus = 0.5 * dxW / dxww;

			if ((!bE) && (iE > -1)) {
				// P
				FwRhie_Chow += koef * (1.0 - fwplus)* (feplus * potent[Press][iE] + (1.0 - feplus) * potent[Press][iP] 
					- (1.0 - fwplus) * potent[Press][iP] - fwplus * potent[Press][iW]) / (sumanb[Vx][iP]);
				// W
				FwRhie_Chow += koef * fwplus * (fpplus*potent[Press][iP] + (1.0 - fpplus) *  potent[Press][iW]
					- (1.0 - fwwplus) * potent[Press][iW] - fwwplus * PWW) / (sumanb[Vx][iW]);
				// w
				FwRhie_Chow -= koef * (potent[Press][iP] - potent[Press][iW]) / apuw;
			}
			else {
				// узла iE нету.
				// надо восстановить давление PE в узле iE с помощью
				// квадратичной экстраполяции
				Real xP, xW, xWW, xP2, xW2, xWW2, xE;
				// при условии отсутствия узла iE узел iWW всегда существует.
				xWW = 0.25 * (x[nvtx[0][iWW] - 1] + x[nvtx[1][iWW] - 1] + x[nvtx[2][iWW] - 1] + x[nvtx[3][iWW] - 1]);
				xW = 0.25 * (x[nvtx[0][iW] - 1] + x[nvtx[1][iW] - 1] + x[nvtx[2][iW] - 1] + x[nvtx[3][iW] - 1]);
				xP = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
				xP2 = xP * xP; xW2 = xW * xW;  xWW2 = xWW * xWW;

				B[0][0] = 1.0; B[0][1] = xWW; B[0][2] = xWW2;
				B[1][0] = 1.0; B[1][1] = xW; B[1][2] = xW2;
				B[2][0] = 1.0; B[2][1] = xP; B[2][2] = xP2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				xE = xP + dxe;
				Real PE, PW, PP;
				PW = potent[Press][iW]; PP = potent[Press][iP];
				PE = 1.0 * (B[0][0] * PWW + B[0][1] * PW + B[0][2] * PP);
				PE += xE * (B[1][0] * PWW + B[1][1] * PW + B[1][2] * PP);
				PE += xE * xE * (B[2][0] * PWW + B[2][1] * PW + B[2][2] * PP);


				FwRhie_Chow += koef * (1.0 - fwplus) * (feplus * PE + (1.0 - feplus) * potent[Press][iP] 
					- (1.0 - fwplus) * potent[Press][iP] - fwplus * potent[Press][iW]) / (sumanb[Vx][iP]);
				FwRhie_Chow += koef * fwplus * (fpplus *  potent[Press][iP] + (1.0 - fpplus) *  potent[Press][iW]
					- (1.0 - fwwplus) * potent[Press][iW] - fwwplus * PWW) / (sumanb[Vx][iW]);
				FwRhie_Chow -= koef * (potent[Press][iP] - potent[Press][iW]) / apuw;
			
				FwRhie_Chow *= mnull;
			}

		}
		else {
			Fw = rhow * dy * potent[Vx][iP];
			Fw_cor = rhow * dy * potent[Vxcor][iP];

			//if (fabs(Fw) < eps) FwRhie_Chow = 0.0;
			//else
			{
				// узлов iW и iWW нету
				// давление в них определяется 
				// с помощью квадратичной экстраполяции
				Real PWW = 0.0; // давление в узле iWW
				Real PW = 0.0; // давление в узле iW
				int iEE = sosed[EE][iP] - 1;
				Real xP, xE, xEE, xP2, xE2, xEE2, xW, xWW;
				// при условии отсутствия узлов iW и iWW узел iEE всегда существует
				xEE = 0.25 * (x[nvtx[0][iEE] - 1] + x[nvtx[1][iEE] - 1] + x[nvtx[2][iEE] - 1] + x[nvtx[3][iEE] - 1]);
				xE = 0.25 * (x[nvtx[0][iE] - 1] + x[nvtx[1][iE] - 1] + x[nvtx[2][iE] - 1] + x[nvtx[3][iE] - 1]);
				xP = 0.25 * (x[nvtx[0][iP] - 1] + x[nvtx[1][iP] - 1] + x[nvtx[2][iP] - 1] + x[nvtx[3][iP] - 1]);
				xEE2 = xEE * xEE; xP2 = xP * xP; xE2 = xE * xE;

				B[0][0] = 1.0; B[0][1] = xP; B[0][2] = xP2;
				B[1][0] = 1.0; B[1][1] = xE; B[1][2] = xE2;
				B[2][0] = 1.0; B[2][1] = xEE; B[2][2] = xEE2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				Real dxww = dxw;
				xW = xP - dxw;
				xWW = xW - dxw;
				Real PEE, PE, PP;
				PEE = potent[Press][iEE]; PE = potent[Press][iE]; PP = potent[Press][iP];
				PW = 1.0 * (B[0][0] * PP + B[0][1] * PE + B[0][2] * PEE);
				PW += xW * (B[1][0] * PP + B[1][1] * PE + B[1][2] * PEE);
				PW += xW * xW * (B[2][0] * PP + B[2][1] * PE + B[2][2] * PEE);

				PWW = 1.0 * (B[0][0] * PP + B[0][1] * PE + B[0][2] * PEE);
				PWW += xWW * (B[1][0] * PP + B[1][1] * PE + B[1][2] * PEE);
				PWW += xWW * xWW * (B[2][0] * PP + B[2][1] * PE + B[2][2] * PEE);

				Real koef = rhow * dy * dy * alpha[Vx];

				Real fpplus = 0.5, fwwplus = 0.5;

				FwRhie_Chow += koef * (1.0 - fwplus)  * (feplus * PE + (1.0 - feplus) * potent[Press][iP] 
				- (1.0 - fwplus) * potent[Press][iP] - fwplus * PW) / (sumanb[Vx][iP]);
				FwRhie_Chow += koef * fwplus  * ( fpplus* potent[Press][iP] +  (1.0 - fpplus)* PW 
				- (1.0 - fwwplus) * PW - fwwplus * PWW) / (sumanb[Vx][iP]);
				FwRhie_Chow -= koef * (potent[Press][iP] - PW) / apuw;
			}
		
			// Поток задан пользователем, какой смысл его коректировать.
			FwRhie_Chow  *= mnull2;
         }

		if ((!bN) && (iN > -1)) {
			Fn = rhon * dx * (fnplus * potent[Vy][iN] + (1.0 - fnplus) * potent[Vy][iP]);
			Fn_cor = rhon * dx * (fnplus * potent[Vycor][iN] + (1.0 - fnplus) * potent[Vycor][iP]);
			// Схема по потоку
			//if (Fn > 0.0) Fn = rhon * dx * potent[Vy][iN]; else Fn = rhon * dx * potent[Vy][iP];

			Real koef = rhon * dx * dx * alpha[Vy];
			int iNN = sosed[NN][iP] - 1;
			Real dynn = dyn;
			Real PNN = 0.0; // давление в узле iNN
			if (iNN > -1) {
				// если узел существует
				dynn = 0.25 * (y[nvtx[0][iNN] - 1] + y[nvtx[1][iNN] - 1] + y[nvtx[2][iNN] - 1] + y[nvtx[3][iNN] - 1]);
				dynn -= 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
				PNN = potent[Press][iNN];
			}
			else {
				// если узел несуществует 
				// недостающее Давление PNN
				// будет определено с помощью 
				// квадратичной экстраполяции.
				Real yN, yP, yS, yN2, yP2, yS2, yNN;
				yN = 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
				yP = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
				// при условии отсутствия узла iNN узел iS всегда существует.
				yS = 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]);
				yN2 = yN * yN; yP2 = yP * yP; yS2 = yS * yS;

				B[0][0] = 1.0; B[0][1] = yS; B[0][2] = yS2;
				B[1][0] = 1.0; B[1][1] = yP; B[1][2] = yP2;
				B[2][0] = 1.0; B[2][1] = yN; B[2][2] = yN2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				yNN = yN + dyn;
				Real PN, PP, PS;
				PN = potent[Press][iN]; PP = potent[Press][iP]; PS = potent[Press][iS];
				PNN = 1.0 * (B[0][0] * PS + B[0][1] * PP + B[0][2] * PN);
				PNN += yNN * (B[1][0] * PS + B[1][1] * PP + B[1][2] * PN);
				PNN += yNN * yNN * (B[2][0] * PS + B[2][1] * PP + B[2][2] * PN);
			}

			Real dxN = qass[iN].dx, dyN = qass[iN].dy;
			//volume(iN, nve, nvtx, x, y, dxN, dyN);
			Real fnnplus = 0.5 * dyN / dynn, fpplus = 0.5 * dyN / dyn;

			if ((!bS) && (iS > -1)) {
				// неравномерная сетка
				FnRhie_Chow += koef * (fnplus) * (fnnplus * PNN + (1.0 - fnnplus) * potent[Press][iN] 
					- (1.0 - fpplus) * potent[Press][iN] - fpplus * potent[Press][iP]) / (sumanb[Vy][iN]);
				FnRhie_Chow += koef * (1.0 - fnplus) * (fnplus * potent[Press][iN] + (1.0 - fnplus) * potent[Press][iP]
					- (1.0 - fsplus) * potent[Press][iP] - fsplus * potent[Press][iS]) / (sumanb[Vy][iP]);
				FnRhie_Chow -= koef * (potent[Press][iN] - potent[Press][iP]) / apvn;
			}
			else {
				// узла iS нету.
				// надо восстановить давление PS в узле iS с помощью
				// квадратичной экстраполяции
				Real yP, yN, yNN, yP2, yN2, yNN2, yS;
				// при условии отсутствия узла iS узел iNN всегда существует.
				yNN = 0.25 * (y[nvtx[0][iNN] - 1] + y[nvtx[1][iNN] - 1] + y[nvtx[2][iNN] - 1] + y[nvtx[3][iNN] - 1]);
				yN = 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
				yP = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
				yP2 = yP * yP; yN2 = yN * yN;  yNN2 = yNN * yNN;

				B[0][0] = 1.0; B[0][1] = yP; B[0][2] = yP2;
				B[1][0] = 1.0; B[1][1] = yN; B[1][2] = yN2;
				B[2][0] = 1.0; B[2][1] = yNN; B[2][2] = yNN2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				yS = yP - dys;
				Real PS, PP, PN;
				PN = potent[Press][iN]; PP = potent[Press][iP];
				PS = 1.0 * (B[0][0] * PP + B[0][1] * PN + B[0][2] * PNN);
				PS += yS * (B[1][0] * PP + B[1][1] * PN + B[1][2] * PNN);
				PS += yS * yS * (B[2][0] * PP + B[2][1] * PN + B[2][2] * PNN);

				// неравномерная сетка
				FnRhie_Chow += koef * (fnplus) * (fnnplus * PNN + (1.0 - fnnplus) * potent[Press][iN] 
					- (1.0 - fpplus) * potent[Press][iN] - fpplus * potent[Press][iP]) / (sumanb[Vy][iN]);
				FnRhie_Chow += koef * (1.0 - fnplus) * (fnplus * potent[Press][iN] + (1.0 - fnplus) * potent[Press][iP]
					- (1.0 - fsplus) * potent[Press][iP] - fsplus * PS) / (sumanb[Vy][iP]);
				FnRhie_Chow -= koef * (potent[Press][iN] - potent[Press][iP]) / apvn;
				
				FnRhie_Chow *= mnull;
			}

		}
		else {
			Fn = rhon * dx * potent[Vy][iP];
			Fn_cor = rhon * dx * potent[Vycor][iP];

			//if (fabs(Fn) < eps)  FnRhie_Chow = 0.0;
			//else
			{

				// узлов iN и iNN нету
				// давление в них определяется
				// с помощью квадратичной экстраполяции
				Real PNN = 0.0; // давление в узле iNN
				Real PN = 0.0; // давление в узле iN
				int iSS = sosed[SS][iP] - 1;
				Real ySS, yS, yP, ySS2, yS2, yP2, yN, yNN;
				// при условии отсутствия узла iN и iNN узел iSS всегда существует.
				ySS = 0.25 * (y[nvtx[0][iSS] - 1] + y[nvtx[1][iSS] - 1] + y[nvtx[2][iSS] - 1] + y[nvtx[3][iSS] - 1]);
				yS = 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]);
				yP = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
				ySS2 = ySS * ySS; yP2 = yP * yP; yS2 = yS * yS;

				B[0][0] = 1.0; B[0][1] = ySS; B[0][2] = ySS2;
				B[1][0] = 1.0; B[1][1] = yS; B[1][2] = yS2;
				B[2][0] = 1.0; B[2][1] = yP; B[2][2] = yP2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				Real dynn = dyn;
				yN = yP + dyn;
				yNN = yN + dyn;
				Real PSS, PS, PP;
				PSS = potent[Press][iSS]; PS = potent[Press][iS]; PP = potent[Press][iP];
				PN = 1.0 * (B[0][0] * PSS + B[0][1] * PS + B[0][2] * PP);
				PN += yN * (B[1][0] * PSS + B[1][1] * PS + B[1][2] * PP);
				PN += yN * yN * (B[2][0] * PSS + B[2][1] * PS + B[2][2] * PP);

				PNN = 1.0 * (B[0][0] * PSS + B[0][1] * PS + B[0][2] * PP);
				PNN += yNN * (B[1][0] * PSS + B[1][1] * PS + B[1][2] * PP);
				PNN += yNN * yNN * (B[2][0] * PSS + B[2][1] * PS + B[2][2] * PP);

				Real koef = rhon * dx * dx * alpha[Vy];

				Real fnnplus = 0.5, fpplus = 0.5;

				// неравномерная сетка
				FnRhie_Chow += koef * (fnplus) * (fnnplus * PNN + (1.0 - fnnplus) * PN 
				- (1.0 - fpplus) * PN -  fpplus* potent[Press][iP]) / (sumanb[Vy][iP]);
				FnRhie_Chow += koef * (1.0 - fnplus) * (fnplus * PN + (1.0 - fnplus) * potent[Press][iP] 
				- (1.0 - fsplus) * potent[Press][iP] - fsplus * potent[Press][iS]) / (sumanb[Vy][iP]);
				FnRhie_Chow -= koef * (PN - potent[Press][iP]) / sumanb[Vy][iP];
			}
		
			// Поток задан пользователем, какой смысл его коректировать.
			FnRhie_Chow  *= mnull2;
        }

		if ((!bS) && (iS > -1)) {
			Fs = rhos * dx * ((1.0 - fsplus) * potent[Vy][iP] + fsplus * potent[Vy][iS]);
			Fs_cor = rhos * dx * ((1.0 - fsplus) * potent[Vycor][iP] + fsplus * potent[Vycor][iS]);
			// Схема по потоку
			//if (Fs > 0.0) Fs = rhos * dx * potent[Vy][iP]; else Fs = rhos * dx * potent[Vy][iS];

			Real koef = rhos * dx * dx * alpha[Vy];
			int iSS = sosed[SS][iP] - 1;
			Real dyss = dys;
			Real PSS = 0.0; // давление в узле iSS
			if (iSS > -1) {
				dyss = 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]);
				dyss -= 0.25 * (y[nvtx[0][iSS] - 1] + y[nvtx[1][iSS] - 1] + y[nvtx[2][iSS] - 1] + y[nvtx[3][iSS] - 1]);
				PSS = potent[Press][iSS];
			}
			else {
				// если узел несуществует 
				 // недостающее Давление PSS
				 // будет определено с помощью 
				 // квадратичной экстраполяции.
				Real yN, yP, yS, yN2, yP2, yS2, ySS;
				// при условии отсутствия узла iSS узел iN всегда существует.
				yN = 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
				yP = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
				yS = 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]);
				yN2 = yN * yN; yP2 = yP * yP; yS2 = yS * yS;

				B[0][0] = 1.0; B[0][1] = yS; B[0][2] = yS2;
				B[1][0] = 1.0; B[1][1] = yP; B[1][2] = yP2;
				B[2][0] = 1.0; B[2][1] = yN; B[2][2] = yN2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				ySS = yS - dys;
				Real PN, PP, PS;
				PN = potent[Press][iN]; PP = potent[Press][iP]; PS = potent[Press][iS];
				PSS = 1.0 * (B[0][0] * PS + B[0][1] * PP + B[0][2] * PN);
				PSS += ySS * (B[1][0] * PS + B[1][1] * PP + B[1][2] * PN);
				PSS += ySS * ySS * (B[2][0] * PS + B[2][1] * PP + B[2][2] * PN);
			}

			Real dxS = qass[iS].dx, dyS = qass[iS].dy;
			//volume(iS, nve, nvtx, x, y, dxS, dyS);
			Real fpplus = 0.5 * dyS / dys, fssplus = 0.5 * dyS / dyss;

			if ((!bN) && (iN>-1)) {
				FsRhie_Chow += koef * (1.0 - fsplus) * (fnplus * potent[Press][iN] + (1.0 - fnplus) * potent[Press][iP] 
					- (1.0 - fsplus) * potent[Press][iP] - fsplus * potent[Press][iS]) / (sumanb[Vy][iP]);
				FsRhie_Chow += koef * fsplus * (fpplus * potent[Press][iP] +  (1.0 - fpplus) * potent[Press][iS]
					- (1.0 - fssplus) * potent[Press][iS] - fssplus * PSS) / (sumanb[Vy][iS]);
				FsRhie_Chow -= koef * (potent[Press][iP] - potent[Press][iS]) / apvs;
			}
			else {
				// узла iN нету.
				// надо восстановить давление PN в узле iN с помощью
				// квадратичной экстраполяции
				Real yP, yS, ySS, yP2, yS2, ySS2, yN;
				// при условии отсутствия узла iE узел iWW всегда существует.
				ySS = 0.25 * (y[nvtx[0][iSS] - 1] + y[nvtx[1][iSS] - 1] + y[nvtx[2][iSS] - 1] + y[nvtx[3][iSS] - 1]);
				yS = 0.25 * (y[nvtx[0][iS] - 1] + y[nvtx[1][iS] - 1] + y[nvtx[2][iS] - 1] + y[nvtx[3][iS] - 1]);
				yP = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
				yP2 = yP * yP; yS2 = yS * yS;  ySS2 = ySS * ySS;

				B[0][0] = 1.0; B[0][1] = ySS; B[0][2] = ySS2;
				B[1][0] = 1.0; B[1][1] = yS; B[1][2] = yS2;
				B[2][0] = 1.0; B[2][1] = yP; B[2][2] = yP2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				yN = yP + dyn;
				Real PN, PS, PP;
				PS = potent[Press][iS]; PP = potent[Press][iP];
				PN = 1.0 * (B[0][0] * PSS + B[0][1] * PS + B[0][2] * PP);
				PN += yN * (B[1][0] * PSS + B[1][1] * PS + B[1][2] * PP);
				PN += yN * yN * (B[2][0] * PSS + B[2][1] * PS + B[2][2] * PP);

				FsRhie_Chow += koef * (1.0 - fsplus) * (fnplus * PN + (1.0 - fnplus) * potent[Press][iP] 
					- (1.0 - fsplus) * potent[Press][iP] - fsplus * PS) / (sumanb[Vy][iP]);
				FsRhie_Chow += koef * fsplus * (fpplus*potent[Press][iP] + (1.0 - fpplus) *  PS
					- (1.0 - fssplus) * PS - fssplus * PSS) / (sumanb[Vy][iS]);
				FsRhie_Chow -= koef * (potent[Press][iP] - PS) / apvs;
			
				FsRhie_Chow *= mnull;
			}

		}
		else {
			Fs = rhos * dx * potent[Vy][iP];
			Fs_cor = rhos * dx * potent[Vycor][iP];

			//if (fabs(Fs) < eps) FsRhie_Chow = 0.0;
			//else
			{

				// узлов iS и iSS нету
				// давление в них определяется 
				// с помощью квадратичной экстраполяции
				Real PSS = 0.0; // давление в узле iSS
				Real PS = 0.0; // давление в узле iS
				int iNN = sosed[NN][iP] - 1;
				Real yP, yN, yNN, yP2, yN2, yNN2, yS, ySS;
				// при условии отсутствия узлов iS и iSS узел iNN всегда существует
				yNN = 0.25 * (y[nvtx[0][iNN] - 1] + y[nvtx[1][iNN] - 1] + y[nvtx[2][iNN] - 1] + y[nvtx[3][iNN] - 1]);
				yN = 0.25 * (y[nvtx[0][iN] - 1] + y[nvtx[1][iN] - 1] + y[nvtx[2][iN] - 1] + y[nvtx[3][iN] - 1]);
				yP = 0.25 * (y[nvtx[0][iP] - 1] + y[nvtx[1][iP] - 1] + y[nvtx[2][iP] - 1] + y[nvtx[3][iP] - 1]);
				yNN2 = yNN * yNN; yP2 = yP * yP; yN2 = yN * yN;

				B[0][0] = 1.0; B[0][1] = yP; B[0][2] = yP2;
				B[1][0] = 1.0; B[1][1] = yN; B[1][2] = yN2;
				B[2][0] = 1.0; B[2][1] = yNN; B[2][2] = yNN2;
				inverse_matrix_simple(B, 3, false); // обращает матрицу
				Real dyss = dys;
				yS = yP - dys;
				ySS = yS - dys;
				Real PNN, PN, PP;
				PNN = potent[Press][iNN]; PN = potent[Press][iN]; PP = potent[Press][iP];
				PS = 1.0 * (B[0][0] * PP + B[0][1] * PN + B[0][2] * PNN);
				PS += yS * (B[1][0] * PP + B[1][1] * PN + B[1][2] * PNN);
				PS += yS * yS * (B[2][0] * PP + B[2][1] * PN + B[2][2] * PNN);

				PSS = 1.0 * (B[0][0] * PP + B[0][1] * PN + B[0][2] * PNN);
				PSS += ySS * (B[1][0] * PP + B[1][1] * PN + B[1][2] * PNN);
				PSS += ySS * ySS * (B[2][0] * PP + B[2][1] * PN + B[2][2] * PNN);

				Real koef = rhos * dx * dx * alpha[Vy];

				Real fpplus = 0.5, fssplus = 0.5;

				FsRhie_Chow += koef * (1.0 - fsplus) * (fnplus * PN + (1.0 - fnplus) * potent[Press][iP] 
				- (1.0 - fsplus) * potent[Press][iP] - fsplus * PS) / (sumanb[Vy][iP]);
				FsRhie_Chow += koef * fsplus * (  fpplus *potent[Press][iP] + (1.0 - fpplus) *PS
				- (1.0 - fssplus) * PS - fssplus * PSS) / (sumanb[Vy][iP]);
				FsRhie_Chow -= koef * (potent[Press][iP] - PS) / apvs;
			}
		
			// Поток задан пользователем, какой смысл его коректировать.
			FsRhie_Chow *= mnull2;
        }

		//for (int l = 0; l < 3; l++) delete B[l];
		//delete[] B;


		 // коэффициенты диффузии:
		 Real GP1, GE1, GW1, GN1, GS1;
		 Real Ge1, Gw1, Gn1, Gs1;
		 int iVarDiffusion = Mu;
		 //if (iVar == Temp) iVarDiffusion = Lam;

		 GP1 = prop[iVarDiffusion][iP];
		 if (!bE) GE1 = prop[iVarDiffusion][iE]; else GE1 = GP1;
		 if (!bN) GN1 = prop[iVarDiffusion][iN]; else GN1 = GP1;
		 if (!bW) GW1 = prop[iVarDiffusion][iW]; else GW1 = GP1;
		 if (!bS) GS1 = prop[iVarDiffusion][iS]; else GS1 = GP1;

		 Ge1 = GE1 * GP1 / ((1 - feplus) * GP1 + feplus * GE1);
		 Gw1 = GW1 * GP1 / ((1 - fwplus) * GP1 + fwplus * GW1);
		 Gn1 = GN1 * GP1 / ((1 - fnplus) * GP1 + fnplus * GN1);
		 Gs1 = GS1 * GP1 / ((1 - fsplus) * GP1 + fsplus * GS1);

		 // Диффузионная составляющая потока:
		 Real De1 = 1.0, Dw1 = 1.0, Dn1 = 1.0, Ds1 = 1.0;
		 if ((!bE)&&(iE>-1)) De1 = Ge1 * qass[iP].dy_dxe;
		 if ((!bW) && (iW > -1)) Dw1 = Gw1 * qass[iP].dy_dxw;
		 if ((!bN) && (iN > -1)) Dn1 = Gn1 * qass[iP].dx_dyn;
		 if ((!bS) && (iS > -1)) Ds1 = Gs1 * qass[iP].dx_dys;

		 // Числа Пекле:
		 Real Pe1=0.0, Pw1=0.0, Pn1=0.0, Ps1=0.0;
		 if ((!bE) && (iE > -1))  Pe1 = fabs(Fe / De1);
		 if ((!bW) && (iW > -1)) Pw1 = fabs(-Fw / Dw1);
		 if ((!bN) && (iN > -1)) Pn1 = fabs(Fn / Dn1);
		 if ((!bS) && (iS > -1)) Ps1 = fabs(-Fs / Ds1);


		 Real lim_arg = 0.25 * (fabs(Pe1) + fabs(Pw1) + fabs(Ps1) + fabs(Pn1));

		 if (lim_arg > 10000000) {
			 std::wcout << De1 << " " << Dw1 << " " << Dn1 << " " << Ds1 << " " << std::endl;
			 std::wcout << Ge1 << " " << Gw1 << " " << Gn1 << " " << Gs1 << " " << std::endl;
			 std::wcout << Fe << " " << Fw << " " << Fn << " " << Fs << " " << std::endl;
			 getchar();
		 }

		 Real lim_Re = sqrt(dx * dy) * sqrt(potent[Vxcor][iP] * potent[Vxcor][iP] + potent[Vycor][iP] * potent[Vycor][iP])*prop[Rho][iP] / prop[Mu][iP];
		 if (lim_arg > Pe_max) Pe_max = lim_arg;
		 if (lim_Re > Re_max) Re_max = lim_Re;

		 Real RCh = lim(lim_arg);// 0.01;//0.05;  0.1;
		Real m1 = 1.0;



		Flux_gran[E][iP] =  (Fe + RCh * FeRhie_Chow) + m1 * RCh * (1.0 - alpha[Vx]) * (Flux_gran_relx[E][iP] - Fe_cor);
		Flux_gran[W][iP] =  (Fw + RCh * FwRhie_Chow) + m1 * RCh * (1.0 - alpha[Vx]) * (Flux_gran_relx[W][iP] - Fw_cor);
		Flux_gran[N][iP] =  (Fn + RCh * FnRhie_Chow) + m1 * RCh * (1.0 - alpha[Vy]) * (Flux_gran_relx[N][iP] - Fn_cor);
		Flux_gran[S][iP] =  (Fs + RCh * FsRhie_Chow) + m1 * RCh * (1.0 - alpha[Vy]) * (Flux_gran_relx[S][iP] - Fs_cor);


		Flux_gran_relx[E][iP] = Flux_gran[E][iP];
		Flux_gran_relx[W][iP] = Flux_gran[W][iP];
		Flux_gran_relx[N][iP] = Flux_gran[N][iP];
		Flux_gran_relx[S][iP] = Flux_gran[S][iP];


		if ((boundary[PAm][iP]) && (!neiman[PAm][iP])) {

			// граничное условие Дирихле.
			sl[PAm][iP].ae = 0.0;
			sl[PAm][iP].aw = 0.0;
			sl[PAm][iP].an = 0.0;
			sl[PAm][iP].as = 0.0;
			sl[PAm][iP].ap = 1.0;
			sl[PAm][iP].b = potent[PAm][iP];			

			tau[E][iP] = 0.0;
			tau[W][iP] = 0.0;
			tau[N][iP] = 0.0;
			tau[S][iP] = 0.0;
		}
		/*else  if ((boundary[PAm][iP]) && (neiman[PAm][iP])) {

			switch (norm[0][iP]) {
			case E: //E узла нет
				// граничное условие Неймана однородное.
				sl[PAm][iP].ae = 0.0;
				sl[PAm][iP].aw = -1.0;
				sl[PAm][iP].an = 0.0;
				sl[PAm][iP].as = 0.0;
				sl[PAm][iP].ap = 1.0;
				sl[PAm][iP].b = 0.0;
				break;
			case N: //N узла нет
				// граничное условие Неймана однородное.
				sl[PAm][iP].ae = 0.0;
				sl[PAm][iP].aw = 0.0;
				sl[PAm][iP].an = 0.0;
				sl[PAm][iP].as = -1.0;
				sl[PAm][iP].ap = 1.0;
				sl[PAm][iP].b = 0.0;

				break;
			case W: //W узла нет
				// граничное условие Неймана однородное.
				sl[PAm][iP].ae = -1.0;
				sl[PAm][iP].aw = 0.0;
				sl[PAm][iP].an = 0.0;
				sl[PAm][iP].as = 0.0;
				sl[PAm][iP].ap = 1.0;
				sl[PAm][iP].b = 0.0;

				break;
			case S: //S узла нет
				// граничное условие Неймана однородное.
				sl[PAm][iP].ae = 0.0;
				sl[PAm][iP].aw = 0.0;
				sl[PAm][iP].an = -1.0;
				sl[PAm][iP].as = 0.0;
				sl[PAm][iP].ap = 1.0;
				sl[PAm][iP].b = 0.0;

				break;
			} // первая нормаль


		}*/
		else
		{
			

			if ((!bE) && (iE>-1)) sl[PAm][iP].ae = (tau[E][iP])*dy/dxe; else sl[PAm][iP].ae = 0.0;
			if ((!bW) && (iW > -1)) sl[PAm][iP].aw = (tau[W][iP])*dy/dxw; else sl[PAm][iP].aw = 0.0;
			if ((!bN) && (iN > -1)) sl[PAm][iP].an = (tau[N][iP])*dx/dyn; else sl[PAm][iP].an = 0.0;
			if ((!bS) && (iS > -1)) sl[PAm][iP].as = (tau[S][iP])*dx/dys; else sl[PAm][iP].as = 0.0;
			sl[PAm][iP].ap = sl[PAm][iP].ae + sl[PAm][iP].aw + sl[PAm][iP].an + sl[PAm][iP].as;

			

		//Real deltat = 1e-6 * dx * dy * alpha[Press] / ((sl[PAm][iP].ae + sl[PAm][iP].aw + sl[PAm][iP].an + sl[PAm][iP].as) * (1 - alpha[Press]));
		//sl[PAm][iP].b = (Fw - Fe + Fs - Fn);// +0.1 * (FwRhie_Chow - FeRhie_Chow + FsRhie_Chow - FnRhie_Chow);
		
		// Должно быть -130, +530 Re=400.
		// без *dx*dy делю на dxe -271000 397520

		sl[PAm][iP].b = 0.0;
		if (((!bW) && (!bE) && (iW > -1) && (iE > -1))&&((!bS) && (!bN) && (iS > -1) && (iN > -1))) {//*dx*dy
				
			    sl[PAm][iP].b += (Flux_gran[W][iP] - Flux_gran[E][iP]);
				sl[PAm][iP].b += (Flux_gran[S][iP] - Flux_gran[N][iP]);
		}
		else {
			Flux_gran[W][iP] = 0.0;
			Flux_gran[E][iP] = 0.0;
			Flux_gran[S][iP] = 0.0;
			Flux_gran[N][iP] = 0.0;

			Flux_gran_relx[W][iP] = 0.0;
			Flux_gran_relx[E][iP] = 0.0;
			Flux_gran_relx[S][iP] = 0.0;
			Flux_gran_relx[N][iP] = 0.0;
		}
		// Это имеет смысл.
		//sl[PAm][iP].b = ((bE || bW ? 0.0 : Fw - Fe) + (bN || bS ? 0.0 : Fs - Fn));
		//sl[PAm][iP].b += ((bE || bW ? 0.0 : RCh *(FwRhie_Chow - FeRhie_Chow)) + (bN || bS ? 0.0 : RCh *(FsRhie_Chow - FnRhie_Chow)));
		

        

	}


} // my_elmatr_quad_PAm

#endif