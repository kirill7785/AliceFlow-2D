// ���� my_elmatr_quad_f.c 
// ������ ������� ��� �����������
// ��������� ���������-�������� �� 
// ����������� �����.

#pragma once
#ifndef MY_ELMATR_QUAD_F_C
#define MY_ELMATR_QUAD_F_C 1

#include "amg1r5.c" // �������������� ������������� ����� ����� ���� � ������ ������� 1984 ���. 
#include "my_linalg.c" // ���������� ������� �������� �������
// ��� �������: 
// eqsolve_simple_gauss - ������ ���� ������� ���������� ������
// eqsolv_simple_holesskii - ������ ���� ������� ���������� ����������

// ������������� ������������� ����� A(|P|)
#include "my_approx_convective.c" // �� ��������� ����� BULG

#include "my_approx_convective2.c" // ������������� ������������� � ������������� �����

#define Real double // float

// ������� ��������
#define Temp 0 // �����������
#define Vx 1 // �������������� ��������
#define Vy 2 // ������������ ��������
#define Press 3 // ��������
#define PAm  4 // �������� �������� 

#define distsheme 100 // ��������� �������� �� ����� ����� � �����

#define Rho 0 // ���������
#define Cp 1 // �����������
#define Lam 2 // ����������������
#define Mu 3 // ������������ ��������
#define BETA_T 4 // ����������� ��������� �������������� ����������

#define E 0 // (east) ������
#define N 1 // �����
#define W 2 // �����
#define S 3 // ��
#define EE 4
#define NN 5 
#define WW 6
#define SS 7



// ��������� ������� ������������ ������
void volume(int iP, int nve, int** &nvtx, Real* &x, Real * &y, Real &dx, Real &dy) {
	// ���������� �������� �������� ������������ ������:
        int i;
	    int *ii=new int[4];// ������ ������ � ����������������� � ��������� ����������.
	    for (i=0; i<4; i++) ii[i]=-1;

	    Real fxmin=1e100, fxmax=-1e100;
	    Real fymin=1e100, fymax=-1e100;

	    // �������������
	    for (i=0; i<4; i++) {
		    if (x[nvtx[i][iP]-1]<fxmin) fxmin=x[nvtx[i][iP]-1];
		    if (y[nvtx[i][iP]-1]<fymin) fymin=y[nvtx[i][iP]-1];

		    if (x[nvtx[i][iP]-1]>fxmax) fxmax=x[nvtx[i][iP]-1];
		    if (y[nvtx[i][iP]-1]>fymax) fymax=y[nvtx[i][iP]-1];
	    }
        for (i=0; i<4; i++) {
		    if ((fabs(x[nvtx[i][iP]-1]-fxmin)<1e-20) && (fabs(y[nvtx[i][iP]-1]-fymin)<1e-20)) ii[0]=nvtx[i][iP]-1;
            if ((fabs(x[nvtx[i][iP]-1]-fxmin)<1e-20) && (fabs(y[nvtx[i][iP]-1]-fymax)<1e-20)) ii[2]=nvtx[i][iP]-1;
		    if ((fabs(x[nvtx[i][iP]-1]-fxmax)<1e-20) && (fabs(y[nvtx[i][iP]-1]-fymin)<1e-20)) ii[1]=nvtx[i][iP]-1;
		    if ((fabs(x[nvtx[i][iP]-1]-fxmax)<1e-20) && (fabs(y[nvtx[i][iP]-1]-fymax)<1e-20)) ii[3]=nvtx[i][iP]-1;
	    }


	
	    dx=0.0, dy=0.0; // ������� ������������ ������
	    dx=x[ii[1]]-x[ii[0]]; dy=y[ii[2]]-y[ii[0]];
		//printf("%.2f %.2f\n",dx,dy); // debug GOOD
		//getchar();
		delete[] ii; // ������������ ����������� ������
}




// �������� ���� ��������� ������� ���� ��� ����������� ��������� 
// ��������� - ��������, ��� ������������ ������������ ������.
// ��� ������������� �����.
void my_elmatr_quad_F(int iP, equation*& sl, int iVar, bool btimedep, Real tauparam, int ishconvection,
	int**& nvtx, bool**& boundary, Real**& potent, Real**& sumanb, Real**& Flux_gran, Real*& x, Real*& y, Real**& prop,
	int**& sosed, bool**& neiman, int**& norm, int nve, Real*& alpha, Real dgx, Real dgy, Real temp_ref, Real** &B) {
    // btimedep==false - ������������, ����� (true) ��������������
	// tauparam - ��� �� �������.
    const bool imitation_time=false; // false - ������ ����� �� ������������, true - ������������

    // iP - ����� ������������ ������������ ������
	int iE, iN, iW, iS; // ������ �������� ����������� �������
	iE=qass[iP].iE; iN= qass[iP].iN; iW= qass[iP].iW; iS= qass[iP].iS;
	sl[iP].iE=iE; sl[iP].iN=iN; sl[iP].iP=iP; sl[iP].iS=iS; sl[iP].iW=iW;

	if ((boundary[iVar][iP]) && (!neiman[iVar][iP])) {
		// ��������� ������� �������.
        sl[iP].ae=0.0;
	    sl[iP].aw=0.0;
	    sl[iP].an=0.0;
	    sl[iP].as=0.0;
		sl[iP].ap=1.0;
		if (iVar == Vx) {
			sumanb[Vx][iP] = 1.0;
		}
		if (iVar == Vy) {
			sumanb[Vy][iP] = 1.0;
		}
		sl[iP].b=potent[iVar][iP];
	}
	else if ((boundary[iVar][iP]) && (neiman[iVar][iP])) {

		// ������� �������
		switch (norm[0][iP]) {
		case E: 
			// ���������� ������� �������.
			sl[iP].ae = 0.0;
			sl[iP].aw = 1.0;
			sl[iP].an = 0.0;
			sl[iP].as = 0.0;
			sl[iP].ap = 1.0;
			if (iVar == Vx) {
				sumanb[Vx][iP] = 1.0;
			}
			if (iVar == Vy) {
				sumanb[Vy][iP] = 1.0;
			}
			sl[iP].b = 0.0;
			break;
		case N: 
			// ���������� ������� �������.
			sl[iP].ae = 0.0;
			sl[iP].aw = 0.0;
			sl[iP].an = 0.0;
			sl[iP].as = 1.0;
			sl[iP].ap = 1.0;
			if (iVar == Vx) {
				sumanb[Vx][iP] = 1.0;
			}
			if (iVar == Vy) {
				sumanb[Vy][iP] = 1.0;
			}
			sl[iP].b = 0.0;
			break;
		case W:  
			// ���������� ������� �������.
			sl[iP].ae = 1.0;
			sl[iP].aw = 0.0;
			sl[iP].an = 0.0;
			sl[iP].as = 0.0;
			sl[iP].ap = 1.0;
			if (iVar == Vx) {
				sumanb[Vx][iP] = 1.0;
			}
			if (iVar == Vy) {
				sumanb[Vy][iP] = 1.0;
			}
			sl[iP].b = 0.0;
			break;
		case S: 
			// ���������� ������� �������.
			sl[iP].ae = 0.0;
			sl[iP].aw = 0.0;
			sl[iP].an = 1.0;
			sl[iP].as = 0.0;
			sl[iP].ap = 1.0;
			if (iVar == Vx) {
				sumanb[Vx][iP] = 1.0;
			}
			if (iVar == Vy) {
				sumanb[Vy][iP] = 1.0;
			}
			sl[iP].b = 0.0;
			
			break;
		} // ������ �������
	}
	else
	{
		// ���� ���������� ����, ���� ������� �������.

		// ������� ������� 
		// � ������ ���� ���������� ����� true
		bool bE=false, bN=false, bW=false, bS=false;
        if ((boundary[iVar][iP]) && (neiman[iVar][iP])) {

			// ������� �������
			switch (norm[0][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // ������ �������

			// ������� �������
            switch (norm[1][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // ������ �������
		}

		
	    // ���������� �������� �������� ������������ ������:
		Real dx= qass[iP].dx, dy= qass[iP].dy;// ����� �������� ������������� ������
	    //volume(iP, nve, nvtx, x, y, dx, dy);
	   
		//printf("%.2f %.2f\n",dx,dy); // debug GOOD
		//getchar();

	    Real dxe= qass[iP].dxe, dxw= qass[iP].dxw, dyn= qass[iP].dyn, dys= qass[iP].dys;
	    //if ((!bE) && (iE>-1)) dxe=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
	    //if ((!bE) && (iE > -1)) dxe-=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
	    //if ((!bW) && (iW > -1)) dxw=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
	    //if ((!bW) && (iW > -1)) dxw-=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
	    //if ((!bN) && (iN > -1)) dyn=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
	    //if ((!bN) && (iN > -1)) dyn-=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
	    //if ((!bS) && (iS > -1)) dys=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
	    //if ((!bS) && (iS > -1)) dys-=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);

	    Real feplus, fwplus, fnplus, fsplus;
	    feplus=0.5*dx/dxe;
	    fwplus=0.5*dx/dxw;
	    fnplus=0.5*dy/dyn;
	    fsplus=0.5*dy/dys;

		//printf("%.2f %.2f %.2f %.2f\n",feplus, fwplus, fnplus, fsplus);
		//getchar();

	    // ��������� ���������������� ������� �������������
	    Real rhoe, rhow, rhon, rhos;
	    Real rP, rE, rN, rW, rS;
	    switch (iVar) {
		    case Temp : rP=prop[Rho][iP]*prop[Cp][iP];
			            if (!bE) rE=prop[Rho][iE]*prop[Cp][iE]; else rE=rP;
                        if (!bN) rN=prop[Rho][iN]*prop[Cp][iN]; else rN=rP;
			            if (!bW) rW=prop[Rho][iW]*prop[Cp][iW]; else rW=rP;
                        if (!bS) rS=prop[Rho][iS]*prop[Cp][iS]; else rS=rP;
			        break;
		    default : rP=prop[Rho][iP];
			          if (!bE) rE=prop[Rho][iE]; else rE=rP;
                      if (!bN) rN=prop[Rho][iN]; else rN=rP;
			          if (!bW) rW=prop[Rho][iW]; else rW=rP;
                      if (!bS) rS=prop[Rho][iS]; else rS=rP;
			        break;
	    }
	    rhoe=rE*rP/(feplus*rP+(1-feplus)*rE);
	    rhow=rW*rP/(fwplus*rP+(1-fwplus)*rW);
	    rhon=rN*rP/(fnplus*rP+(1-fnplus)*rN);
	    rhos=rS*rP/(fsplus*rP+(1-fsplus)*rS);

		
	    // �������� ������������ �������� �� ����� ��.
        Real ue, uw, vn, vs;
	    if (!bE) ue=feplus*potent[Vxcor][iE]+(1-feplus)*potent[Vxcor][iP]; else ue=potent[Vxcor][iP];
        if (!bW) uw=fwplus*potent[Vxcor][iP]+(1-fwplus)*potent[Vxcor][iW]; else uw=potent[Vxcor][iP];
	    if (!bN) vn=fnplus*potent[Vycor][iN]+(1-fnplus)*potent[Vycor][iP]; else vn=potent[Vycor][iP];
        if (!bS) vs=fsplus*potent[Vycor][iP]+(1-fsplus)*potent[Vycor][iS]; else vs=potent[Vycor][iP];
    

	    // ������
	    Real Fe=0.0, Fw=0.0, Fn=0.0, Fs=0.0;
	   /* if (!bE) Fe = rhoe * ue * dy;
	    if (!bW) Fw=rhow*uw*dy;
	    if (!bN) Fn=rhon*vn*dx;
	    if (!bS) Fs=rhos*vs*dx;
		*/
		
		Fw = Flux_gran[W][iP];
		Fe = Flux_gran[E][iP];
		Fs = Flux_gran[S][iP];
		Fn = Flux_gran[N][iP];

	    // ������������ ��������:
	    Real GP, GE, GW, GN, GS;
	    Real Ge, Gw, Gn, Gs;
		int iVarDiffusion = Mu;
		if (iVar == Temp) iVarDiffusion = Lam;

		GP = prop[iVarDiffusion][iP];
		if (!bE) GE = prop[iVarDiffusion][iE]; else GE = GP;
		if (!bN) GN = prop[iVarDiffusion][iN]; else GN = GP;
		if (!bW) GW = prop[iVarDiffusion][iW]; else GW = GP;
		if (!bS) GS = prop[iVarDiffusion][iS]; else GS = GP;

        Ge=GE*GP/((1 - feplus)*GP+ feplus *GE);
	    Gw=GW*GP/((1 - fwplus)*GP+ fwplus *GW);
	    Gn=GN*GP/((1 - fnplus)*GP+ fnplus *GN);
	    Gs=GS*GP/((1 - fsplus)*GP+ fsplus *GS);

	    // ������������ ������������ ������:
	    Real De=1.0, Dw=1.0, Dn=1.0, Ds=1.0;
	    if (!bE) De=Ge* qass[iP].dy_dxe;
	    if (!bW) Dw=Gw* qass[iP].dy_dxw;
	    if (!bN) Dn=Gn* qass[iP].dx_dyn;
	    if (!bS) Ds=Gs* qass[iP].dx_dys;

		Real rev_sign = 0.0;

		if (ishconvection == UDS) {

			// ���������� ������������� ����������� �������:
			if ((bE) || (iE == -1)) sl[iP].ae = 0.0; else {
				sl[iP].ae = De  + fmax(-Fe, 0);
				rev_sign += De + fmax(Fe, 0);
			}
			if ((bW) || (iW == -1)) sl[iP].aw = 0.0; else {
				sl[iP].aw = Dw  + fmax(Fw, 0);
				rev_sign += Dw + fmax(-Fw, 0);
			}
			if ((bN) || (iN == -1)) sl[iP].an = 0.0; else {
				sl[iP].an = Dn  + fmax(-Fn, 0);
				rev_sign += Dn + fmax(Fn, 0);
			}
			if ((bS) || (iS == -1)) sl[iP].as = 0.0; else {
				sl[iP].as = Ds  + fmax(Fs, 0);
				rev_sign += Ds + fmax(-Fs, 0);
			}

		}
		else {

			// ����� �����:
			Real Pe, Pw, Pn, Ps;
			Pe = Fe / De;
			Pw = -Fw / Dw;
			Pn = Fn / Dn;
			Ps = -Fs / Ds;

			

			if (ishconvection < distsheme) {
				// ���������� ������������� ����������� �������:
				if ((bE) || (iE == -1)) sl[iP].ae = 0.0; else {
					sl[iP].ae = De * ApproxConvective(fabs(Pe), ishconvection) + fmax(-Fe, 0);
					//rev_sign += De * ApproxConvective(fabs(Pe), ishconvection) + fmax(Fe, 0);;
				}
				if ((bW) || (iW == -1)) sl[iP].aw = 0.0; else {
					sl[iP].aw = Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(Fw, 0);
					//rev_sign += Dw * ApproxConvective(fabs(Pw), ishconvection) + fmax(-Fw, 0);;
				}
				if ((bN) || (iN == -1)) sl[iP].an = 0.0; else {
					sl[iP].an = Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(-Fn, 0);
					//rev_sign += Dn * ApproxConvective(fabs(Pn), ishconvection) + fmax(Fn, 0);;
				}
				if ((bS) || (iS == -1)) sl[iP].as = 0.0; else {
					sl[iP].as = Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(Fs, 0);
					//rev_sign += Ds * ApproxConvective(fabs(Ps), ishconvection) + fmax(-Fs, 0);;
				}
			}
			else
			{
				// ���������� ������������� ����������� �������:
				if (bE) sl[iP].ae = 0.0; else sl[iP].ae = -Fe * fC(Pe, ishconvection) + De * fD(Pe, ishconvection);
				if (bW) sl[iP].aw = 0.0; else sl[iP].aw = Fw * fC(Pw, ishconvection) + Dw * fD(Pw, ishconvection);
				if (bN) sl[iP].an = 0.0; else sl[iP].an = -Fn * fC(Pn, ishconvection) + Dn * fD(Pn, ishconvection);
				if (bS) sl[iP].as = 0.0; else sl[iP].as = Fs * fC(Ps, ishconvection) + Ds * fD(Ps, ishconvection);
			}

		}

        Real tau, apzero;
	    Real Fold=0.0; // �������� ������� � ����������� ���������� ����
	    if (btimedep) {
		   // ��������������
		   tau=tauparam;
		   Fold=potent[iVar][iP];
	    }
	    else {
	    	// ������������

		    // ����� �����������:
		    tau=rP*dx*dy*alpha[iVar]/((sl[iP].ae+sl[iP].aw+sl[iP].an+sl[iP].as)*(1.0-alpha[iVar]));
			switch (iVar) {
				// ����� ������������� � ���������
				// ��������������� ��������� �������������.
				case Vx : Fold=potent[Vxcor][iP]; break;
				case Vy : Fold=potent[Vycor][iP]; break;
				default : Fold=potent[iVar][iP]; break;
			}
	    	
	    }

		if ((!btimedep) && (!imitation_time)) {
			// ������������ � �������� ����� �� ������� �� ������������
			apzero=0.0;
		}
		else apzero=rP*dx*dy/tau;

		apzero = 0.0;
  
	    // ������������ ����
	    Real dSc=0.0, dSp=0.0;

		// ��� Vx � Vy ��������� ������� �������� ����������� �� ����������� �����
        Real zP, zW, zE, zS, zN;
        //int i; // ������� ����� for
		

		//Real** B = new Real * [3];
		//for (int l = 0; l < 3; l++) B[l] = new Real[3];
		//for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = 0.0;

        zP=potent[Press][iP]; 
		if ((!bW)  && (iW > -1)) zW = potent[Press][iW];
		else {
			    // ������������ �������������
                // ���� iW ����.
				// ���� ������������ �������� zW � ���� iW � �������
				// ������������ �������������
                Real xP, xE, xEE, xP2, xE2, xEE2, xW;
				int iEE=sosed[EE][iP]-1;
				// ��� ������� ���������� ���� iW ���� iEE ������ ����������.
                xEE=0.25*(x[nvtx[0][iEE]-1]+x[nvtx[1][iEE]-1]+x[nvtx[2][iEE]-1]+x[nvtx[3][iEE]-1]);
				xE=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xP2=xP*xP; xE2=xE*xE;  xEE2=xEE*xEE;
				
				B[0][0]=1.0; B[0][1]=xP; B[0][2]=xP2;
				B[1][0]=1.0; B[1][1]=xE; B[1][2]=xE2;
				B[2][0]=1.0; B[2][1]=xEE; B[2][2]=xEE2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				xW=xP-dxw;
				Real  PPloc, PEloc, PEEloc;
				 PPloc=potent[Press][iP]; PEloc=potent[Press][iE]; PEEloc=potent[Press][iEE];
				zW=1.0*(B[0][0]*PPloc+B[0][1]*PEloc+B[0][2]*PEEloc);
				zW+=xW*(B[1][0]*PPloc+B[1][1]*PEloc+B[1][2]*PEEloc);
				zW+=xW*xW*(B[2][0]*PPloc+B[2][1]*PEloc+B[2][2]*PEEloc); 
		}
		if ((!bE)  && (iE > -1)) zE = potent[Press][iE];
		else { // �������� � ���� iE ����� ����������
			   // � ������� ������������ ������������
			
                // ���� iE ����.
				// ���� ������������ �������� zE � ���� iE � �������
				// ������������ �������������
                Real xP, xW, xWW, xP2, xW2, xWW2, xE;
                // ��� ������� ���������� ���� iE ���� iWW ������ ����������.
				int iWW = sosed[WW][iP]-1;
                xWW=0.25*(x[nvtx[0][iWW]-1]+x[nvtx[1][iWW]-1]+x[nvtx[2][iWW]-1]+x[nvtx[3][iWW]-1]);
				xW=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xP2=xP*xP; xW2=xW*xW;  xWW2=xWW*xWW;
				
				B[0][0]=1.0; B[0][1]=xWW; B[0][2]=xWW2;
				B[1][0]=1.0; B[1][1]=xW; B[1][2]=xW2;
				B[2][0]=1.0; B[2][1]=xP; B[2][2]=xP2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				xE=xP+dxe;
				Real PWWloc, PWloc, PPloc;
			    PWWloc=potent[Press][iWW]; PWloc=potent[Press][iW]; PPloc=potent[Press][iP];
				zE=1.0*(B[0][0]*PWWloc+B[0][1]*PWloc+B[0][2]*PPloc);
				zE+=xE*(B[1][0]*PWWloc+B[1][1]*PWloc+B[1][2]*PPloc);
				zE+=xE*xE*(B[2][0]*PWWloc+B[2][1]*PWloc+B[2][2]*PPloc);
		}
        if ((!bS) && (iS > -1)) zS = potent[Press][iS];
		else {// ������������ �������������
                // ���� iS ����.
				// ���� ������������ �������� �������� zS � ���� iS � �������
				// ������������ �������������
                Real yP, yN, yNN, yP2, yN2, yNN2, yS;
				int iNN=sosed[NN][iP]-1;
				// ��� ������� ���������� ���� iS ���� iNN ������ ����������.
                yNN=0.25*(y[nvtx[0][iNN]-1]+y[nvtx[1][iNN]-1]+y[nvtx[2][iNN]-1]+y[nvtx[3][iNN]-1]);
				yN=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				yP2=yP*yP; yN2=yN*yN;  yNN2=yNN*yNN;
				
				B[0][0]=1.0; B[0][1]=yP; B[0][2]=yP2;
				B[1][0]=1.0; B[1][1]=yN; B[1][2]=yN2;
				B[2][0]=1.0; B[2][1]=yNN; B[2][2]=yNN2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				yS=yP-dys;
				Real  PPloc, PNloc, PNNloc;
				PPloc=potent[Press][iP]; PNloc=potent[Press][iN]; PNNloc=potent[Press][iNN];
				zS=1.0*(B[0][0]*PPloc+B[0][1]*PNloc+B[0][2]*PNNloc);
				zS+=yS*(B[1][0]*PPloc+B[1][1]*PNloc+B[1][2]*PNNloc);
				zS+=yS*yS*(B[2][0]*PPloc+B[2][1]*PNloc+B[2][2]*PNNloc);
		}
		if ((!bN) && (iN > -1)) zN = potent[Press][iN];
		else {   // ������������ �������������
                // ���� iN ����.
				// ���� ������������ �������� zN � ���� iN � �������
				// ������������ �������������
                Real yP, yS, ySS, yP2, yS2, ySS2, yN;
                // ��� ������� ���������� ���� iN ���� iSS ������ ����������.
				int iSS=sosed[SS][iP]-1;
                ySS=0.25*(y[nvtx[0][iSS]-1]+y[nvtx[1][iSS]-1]+y[nvtx[2][iSS]-1]+y[nvtx[3][iSS]-1]);
				yS=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				yP2=yP*yP; yS2=yS*yS;  ySS2=ySS*ySS;
				
				B[0][0]=1.0; B[0][1]=ySS; B[0][2]=ySS2;
				B[1][0]=1.0; B[1][1]=yS; B[1][2]=yS2;
				B[2][0]=1.0; B[2][1]=yP; B[2][2]=yP2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				yN=yP+dyn;
				Real PSSloc, PSloc, PPloc;
			    PSSloc=potent[Press][iSS]; PSloc=potent[Press][iS]; PPloc=potent[Press][iP];
				zN=1.0*(B[0][0]*PSSloc+B[0][1]*PSloc+B[0][2]*PPloc);
				zN+=yN*(B[1][0]*PSSloc+B[1][1]*PSloc+B[1][2]*PPloc);
				zN+=yN*yN*(B[2][0]*PSSloc+B[2][1]*PSloc+B[2][2]*PPloc); 
		}
        

       
		Real body_force=0.0; // �������� ���� ��������� � ����� ������� ��� �������� (�����������)
		Real tp1, tp2;

	    switch (iVar) {
	    	case Vx : // ����� ������ ������
					  //if (Fe>0.0) tp1=zP; else tp1=zE;
					  //if (Fw>0.0) tp2=zW; else tp2=zP;
					  tp1 = (1.0 - feplus) * zP + feplus * zE;
					  tp2 = (1.0 - fwplus) * zP + fwplus * zW;
					  dSc=(tp2-tp1)*dy;
					  if (fabs(prop[BETA_T][iP]) > 1e-37) {
						  // ���������� ����
                         body_force=-prop[Rho][iP]*prop[BETA_T][iP]*dgx*(potent[Temp][iP]-temp_ref);
					  }
					  else 
					  {
						  // ���� �������
                         body_force=prop[Rho][iP]*dgx;
					  }
                      dSc+=body_force*dx*dy;				       
					  break;
            case Vy : // ����� ������ ������
					  //if (Fn>0.0) tp1=zP; else tp1=zN;
					  //if (Fs>0.0) tp2=zS; else tp2=zP;
					  tp1 = (1.0 - fnplus) * zP + fnplus * zN;
					  tp2 = (1.0 - fsplus) * zP + fsplus * zS;
					  dSc=(tp2-tp1)*dx;
				      if (fabs(prop[BETA_T][iP]) > 1e-37) {
						  // ���������� ����
                         body_force=-prop[Rho][iP]*prop[BETA_T][iP]*dgy*(potent[Temp][iP]-temp_ref);
					  }
					  else 
					  {
						  // ���� �������
                         body_force=prop[Rho][iP]*dgy;
					  }
                      dSc+=body_force*dx*dy; 
					  break;
	    	default : dSc=0.0; break;
	    }
	
		//for (int l = 0; l < 3; l++) delete B[l];
		//delete[] B;

		// ���������� ���������� ����� ������� ������� ������������� � ������� ������ ���������� ���������.
		Real attrs = 0.0;
		if (0&&((!bW) && (iW > -1)) && ((!bE) && (iE > -1)) && ((!bS) && (iS > -1)) && ((!bN) && (iN > -1))) {
			attrs += -fmax(Fe, 0.0) * (ue - potent[Vxcor][iP]) + fmax(-Fe, 0.0) * (ue - potent[Vxcor][iE]);
			attrs += -fmax(-Fw, 0.0) * (uw - potent[Vxcor][iP]) + fmax(Fw, 0.0) * (uw - potent[Vxcor][iW]);
			attrs += -fmax(Fn, 0.0) * (vn - potent[Vycor][iP]) + fmax(-Fn, 0.0) * (vn - potent[Vycor][iN]);
			attrs += -fmax(-Fs, 0.0) * (vs - potent[Vycor][iP]) + fmax(Fs, 0.0) * (vs - potent[Vycor][iS]);
		}

		// ����� �����:
		Real Pe, Pw, Pn, Ps;
		Pe = Fe / De;
		Pw = -Fw / Dw;
		Pn = Fn / Dn;
		Ps = -Fs / Ds;

		Real lim_arg = 0.25 * (fabs(Pe) + fabs(Pw) + fabs(Ps) + fabs(Pn));

		//if (ilim < 20) {
			//std::cout << "Pe= "<< lim_arg << "   limiter= " << lim(lim_arg) << std::endl;
			//getchar();
			//ilim++;
		//}

		Real Kap = 0.0*lim(lim_arg);
		Kap = 0.0;
		// ��� ������ �� ��������. � ������� ��� ���������� ���������.

		if (iVar == Vx) {
			sumanb[Vx][iP] = sl[iP].ae + sl[iP].aw + sl[iP].an + sl[iP].as + Kap*(Fe - Fw + Fn - Fs);
			//sumanb[Vx][iP] = rev_sign + Kap * fabs(Fe - Fw + Fn - Fs);
		}
		if (iVar == Vy) {
			sumanb[Vy][iP] = sl[iP].ae + sl[iP].aw + sl[iP].an + sl[iP].as + Kap*(Fe - Fw + Fn - Fs);
			//sumanb[Vy][iP] = rev_sign + Kap * fabs(Fe - Fw + Fn - Fs);
		}

        sl[iP].ap=sl[iP].ae+sl[iP].aw+sl[iP].an+sl[iP].as+apzero-dSp*dx*dy; // ������������ ������� �������
		//sl[iP].ap =  rev_sign + apzero - dSp * dx * dy; // ������������ ������� �������
		sl[iP].ap += Kap*(Fe-Fw+Fn-Fs);
	    sl[iP].b=dSc+apzero*Fold+attrs;// ������ ����� //-Fold*(Fe-Fw+Fn-Fs); // ���� ���� �������� ���������� ������� 

	}

} // my_elmatr_quad_F



#endif