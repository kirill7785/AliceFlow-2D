// ���� pamendment.c ����������
// ������ ������� ��� ���������
// �������� ��������.
// ���������� ������� ���-���.
// ���������� 4 ��� 2011.
// ��� ������ �������� �������� ������
// � ����������� ������ ��� ����������� �����.
// ���� �������, �������� ��������. ����������� ���� pamendment2.c

#pragma once
#ifndef PAMENDMENTv_0_01_C
#define PAMENDMENTv_0_01_C 1

#include "my_linalg.c" // ���������� ������� �������� �������
// ��� �������: 
// eqsolve_simple_gauss - ������ ���� ������� ���������� ������
// eqsolv_simple_holesskii - ������ ���� ������� ���������� ����������

// ������������� ����������� ��������� ���������-��������
// �� ����������� �����
#include "my_elmatr_quad_f.c"


// ��������� ��������
void correct(int iP, equation** &sl, int iVar, 
			 int** nvtx, bool** boundary, Real** potent,
			 Real* x, Real* y, int** sosed, bool** neiman,
			 int** norm, int nve) {
    // iP - ����� ������������ ������������ ������
	int iE, iN, iW, iS; // ������ �������� ����������� �������
	iE=sosed[E][iP]-1; iN=sosed[N][iP]-1; iW=sosed[W][iP]-1; iS=sosed[S][iP]-1;

    if (!(boundary[iVar][iP]) && (!neiman[iVar][iP])) {   

        // ������� ������� 
		// � ������ ���� ���������� ����� true
		bool bE=false, bN=false, bW=false, bS=false;
        if ((boundary[iVar][iP]) && (neiman[iVar][iP])) {
			switch (norm[0][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // ������ �������

            switch (norm[1][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // ������ �������
		}    

		
	    // ���������� �������� �������� ������������ ������:
	
	    Real dx=0.0, dy=0.0; // ������� ������������ ������
		volume(iP, nve, nvtx, x, y, dx, dy);
		//printf("%.2f %.2f\n",dx,dy); // debug GOOD
		//getchar();

        Real dxe=1.0, dxw=1.0, dyn=1.0, dys=1.0;
	    if (!bE) dxe=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
	    if (!bE) dxe-=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]); 
		else { iE=iP; dxe=dx; }
	    if (!bW) dxw=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
	    if (!bW) dxw-=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]); 
		else { iW=iP; dxw=dx; }
	    if (!bN) dyn=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
	    if (!bN) dyn-=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]); 
		else { iN=iP; dyn=dy; }
	    if (!bS) dys=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
	    if (!bS) dys-=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]); 
		else { iS=iP; dys=dy; }

	    Real feplus, fwplus, fnplus, fsplus;
	    feplus=0.5*dx/dxe;
	    fwplus=0.5*dx/dxw;
	    fnplus=0.5*dy/dyn;
	    fsplus=0.5*dy/dys;

		Real dl=dx;
		if (iVar==Vx) dl=dy;

		int i; // ������� ����� for
		

		Real PAmP, PAmW, PAmE, PAmS, PAmN;
		PAmP=potent[PAm][iP];
        PAmW=potent[PAm][iW];
        if (!bE) PAmE=potent[PAm][iE];
		else { // ������������ �������������
                // ���� iE ����.
				// ���� ������������ �������� PAmE � ���� iE � �������
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
				Real PAmWW, PAmW, PAmP;
			    PAmWW=potent[PAm][iWW]; PAmW=potent[PAm][iW]; PAmP=potent[PAm][iP];
				PAmE=1.0*(B[0][0]*PAmWW+B[0][1]*PAmW+B[0][2]*PAmP);
				PAmE+=xE*(B[1][0]*PAmWW+B[1][1]*PAmW+B[1][2]*PAmP);
				PAmE+=xE*xE*(B[2][0]*PAmWW+B[2][1]*PAmW+B[2][2]*PAmP);
		}
		PAmS=potent[PAm][iS];
		PAmN=potent[PAm][iN];
		

		Real deltaP=0.0;
		switch (iVar) {
			case Vx : deltaP=(fwplus*PAmP+(1-fwplus)*PAmW);
				      deltaP-=(feplus*PAmE+(1-feplus)*PAmP); 
					  break;
			case Vy : deltaP=(fsplus*PAmP+(1-fsplus)*PAmS);
				      deltaP-=(fnplus*PAmN+(1-fnplus)*PAmP);
				      break;
		}

		// ��������� �������� �� ������ ������������ ������ ����������.
		potent[iVar][iP]+=dl*(deltaP)/sl[iVar][iP].ap;//alpha[iVar]*
	}
    
} // correct

// ���������� ������� ��� ��������� 
// �������� ��������
void my_elmatr_quad_PAm(int iP, equation** &sl, int** nvtx, bool** boundary,
						Real** potent, Real* x, Real* y, Real** prop, 
					    int** sosed, bool** neiman, int** norm, int nve,
						Real* alpha, Real** &rhie_chow) {
    // iP - ����� ������������ ������������ ������
	int iE, iN, iW, iS; // ������ �������� ����������� �������
	iE=sosed[E][iP]-1; iN=sosed[N][iP]-1; iW=sosed[W][iP]-1; iS=sosed[S][iP]-1;
	sl[PAm][iP].iE=iE; sl[PAm][iP].iN=iN; sl[PAm][iP].iP=iP; sl[PAm][iP].iS=iS; sl[PAm][iP].iW=iW;


	
    if ((boundary[PAm][iP]) && (!neiman[PAm][iP])) {
		// ��������� ������� �������.
        sl[PAm][iP].ae=0.0;
	    sl[PAm][iP].aw=0.0;
	    sl[PAm][iP].an=0.0;
	    sl[PAm][iP].as=0.0;
		sl[PAm][iP].ap=1.0;
		sl[PAm][iP].b=potent[PAm][iP];
	}
	else
	{
	
		// ���� ���������� ����, ���� ������� �������.

		// ������� ������� 
		// � ������ ���� ���������� ����� true
		bool bE=false, bN=false, bW=false, bS=false;
        if ((boundary[PAm][iP]) && (neiman[PAm][iP])) {
			switch (norm[0][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // ������ �������

            switch (norm[1][iP]) {
				case E : bE=true; break;
				case N : bN=true; break;
				case W : bW=true; break;
				case S : bS=true; break;
			} // ������ �������
		}    

	     

		
	    // ���������� �������� �������� ������������ ������:
        
	    // ����� �������� ������������� ������
	    Real dx=0.0, dy=0.0; // ������� ������������ ������
	    volume(iP, nve, nvtx, x, y, dx, dy);
		//printf("%.2f %.2f\n",dx,dy); // debug GOOD
		//getchar();

        Real dxe=1.0, dxw=1.0, dyn=1.0, dys=1.0;
	    if (!bE) dxe=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
	    if (!bE) dxe-=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]); else dxe=dxw;
	    if (!bW) dxw=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
	    if (!bW) dxw-=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]); else dxw=dxe;
	    if (!bN) dyn=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
	    if (!bN) dyn-=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]); else dyn=dys;
	    if (!bS) dys=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
	    if (!bS) dys-=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]); else dys=dyn;

	    Real feplus, fwplus, fnplus, fsplus;
	    feplus=0.5*dx/dxe;
	    fwplus=0.5*dx/dxw;
	    fnplus=0.5*dy/dyn;
	    fsplus=0.5*dy/dys;

		//printf("%.2f %.2f %.2f %.2f\n",feplus, fwplus, fnplus, fsplus);
		//getchar();

		Real apue=1.0, apuw=1.0, apvn=1.0, apvs=1.0;
		if (!bE) apue=sl[Vx][iE].ap*sl[Vx][iP].ap/(feplus*sl[Vx][iP].ap+(1-feplus)*sl[Vx][iE].ap); else apue=sl[Vx][iP].ap;
		if (!bW) apuw=sl[Vx][iW].ap*sl[Vx][iP].ap/(fwplus*sl[Vx][iP].ap+(1-fwplus)*sl[Vx][iW].ap); else apuw=sl[Vx][iP].ap;
		if (!bN) apvn=sl[Vy][iN].ap*sl[Vy][iP].ap/(fnplus*sl[Vy][iP].ap+(1-fnplus)*sl[Vy][iN].ap); else apvn=sl[Vy][iP].ap;
		if (!bS) apvs=sl[Vy][iS].ap*sl[Vy][iP].ap/(fsplus*sl[Vy][iP].ap+(1-fsplus)*sl[Vy][iS].ap); else apvs=sl[Vy][iP].ap;

		Real de, dw, dn, ds;
		de=dy/apue; dw=dy/apuw; 
		dn=dx/apvn; ds=dx/apvs;

        // ��������� ���������������� ������� �������������
	    Real rhoe, rhow, rhon, rhos;
	    Real rP, rE, rN, rW, rS;

        rP=prop[Rho][iP];
		if (!bE) rE=prop[Rho][iE]; else rE=rP;
        if (!bN) rN=prop[Rho][iN]; else rN=rP;
		if (!bW) rW=prop[Rho][iW]; else rW=rP;
        if (!bS) rS=prop[Rho][iS]; else rS=rP;

		rhoe=rE*rP/(feplus*rP+(1.0-feplus)*rE);
	    rhow=rW*rP/(fwplus*rP+(1.0-fwplus)*rW);
	    rhon=rN*rP/(fnplus*rP+(1.0-fnplus)*rN);
	    rhos=rS*rP/(fsplus*rP+(1.0-fsplus)*rS);

        if (!bE) sl[PAm][iP].ae=rhoe*de*dy; else sl[PAm][iP].ae=0.0;
	    if (!bW) sl[PAm][iP].aw=rhow*dw*dy; else sl[PAm][iP].aw=0.0;
	    if (!bN) sl[PAm][iP].an=rhon*dn*dx; else sl[PAm][iP].an=0.0;
	    if (!bS) sl[PAm][iP].as=rhos*ds*dx; else sl[PAm][iP].as=0.0;
		sl[PAm][iP].ap=sl[PAm][iP].ae+sl[PAm][iP].aw+sl[PAm][iP].an+sl[PAm][iP].as;


        Real Fw=0.0, Fe=0.0, Fs=0.0, Fn=0.0; 
		Real FwRhie_Chow=0.0, FeRhie_Chow=0.0, FsRhie_Chow=0.0, FnRhie_Chow=0.0;

		if (!bE) {
            Fe=rhoe*dy*(feplus*potent[Vx][iE]+(1.0-feplus)*potent[Vx][iP]);
			///*
			Real koef=rhoe*dy*dy/(2.0*apue);
			int iEE=sosed[EE][iP]-1;
			Real dxee=dxe;
			Real PEE=0.0; // �������� � ���� iEE
			if (iEE>-1) {
				// ���� ���� ����������
				dxee=0.25*(x[nvtx[0][iEE]-1]+x[nvtx[1][iEE]-1]+x[nvtx[2][iEE]-1]+x[nvtx[3][iEE]-1]);
	            dxee-=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
				PEE=potent[Press][iEE];
			}
			else {
				// ���� ���� ������������ 
				// ����������� �������� PEE
				// ����� ���������� � ������� 
				// ������������ �������������.
				Real xE, xP, xW, xE2, xP2, xW2, xEE;
				xE=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				// ��� ������� ���������� ���� iEE ���� iW ������ ����������.
				xW=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
				xE2=xE*xE; xP2=xP*xP; xW2=xW*xW;
				
				int l; // ������� ����� for
				
				B[0][0]=1.0; B[0][1]=xW; B[0][2]=xW2;
				B[1][0]=1.0; B[1][1]=xP; B[1][2]=xP2;
				B[2][0]=1.0; B[2][1]=xE; B[2][2]=xE2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				xEE=xE+dxe;
				Real PE, PP, PW;
				PE=potent[Press][iE]; PP=potent[Press][iP]; PW=potent[Press][iW];
				PEE=1.0*(B[0][0]*PW+B[0][1]*PP+B[0][2]*PE);
				PEE+=xEE*(B[1][0]*PW+B[1][1]*PP+B[1][2]*PE);
				PEE+=xEE*xEE*(B[2][0]*PW+B[2][1]*PP+B[2][2]*PE);
			}
			if (!bW) {
				FeRhie_Chow+=koef*(PEE-potent[Press][iP])/(dxe+dxee);
                FeRhie_Chow+=koef*(potent[Press][iE]-potent[Press][iW])/(dxe+dxw);
			    FeRhie_Chow-=koef*2.0*(potent[Press][iE]-potent[Press][iP])/dxe;
			}
			else {
                // ���� iW ����.
				// ���� ������������ �������� PW � ���� iW � �������
				// ������������ �������������
                Real xP, xE, xEE, xP2, xE2, xEE2, xW;
                xEE=0.25*(x[nvtx[0][iEE]-1]+x[nvtx[1][iEE]-1]+x[nvtx[2][iEE]-1]+x[nvtx[3][iEE]-1]);
                // ��� ������� ���������� ���� iWW ���� iE ������ ����������.
				xE=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xP2=xP*xP; xE2=xE*xE;  xEE2=xEE*xEE;
				
				B[0][0]=1.0; B[0][1]=xP; B[0][2]=xP2;
				B[1][0]=1.0; B[1][1]=xE; B[1][2]=xE2;
				B[2][0]=1.0; B[2][1]=xEE; B[2][2]=xEE2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				xW=xP-dxw;
				Real PE, PP, PW;
				PE=potent[Press][iE]; PP=potent[Press][iP]; 
				PW=1.0*(B[0][0]*PP+B[0][1]*PE+B[0][2]*PEE);
				PW+=xW*(B[1][0]*PP+B[1][1]*PE+B[1][2]*PEE);
				PW+=xW*xW*(B[2][0]*PP+B[2][1]*PE+B[2][2]*PEE); 

				// ���������� �������
                //FeRhie_Chow+=koef*(PEE-potent[Press][iP])/(dxe+dxee);
                //FeRhie_Chow-=koef*1.5*(potent[Press][iE]-potent[Press][iP])/dxe;
				// ����� �������
				FeRhie_Chow+=koef*(PEE-potent[Press][iP])/(dxe+dxee);
                FeRhie_Chow+=koef*(potent[Press][iE]-PW)/(dxe+dxw);
			    FeRhie_Chow-=koef*2.0*(potent[Press][iE]-potent[Press][iP])/dxe;
			}
			//*/
		}
		else {
			Fe=rhoe*dy*potent[Vx][iP];
            
            // ����� iE � iEE ����
			// �������� � ��� ������������
			// � ������� ������������ �������������
            Real PEE=0.0; // �������� � ���� iEE
			Real PE=0.0; // �������� � ���� iE
            int iWW=sosed[WW][iP]-1;
			Real xWW, xW, xP, xWW2, xW2, xP2, xE, xEE;
			// ��� ������� ���������� ���� iE � iEE ���� iWW ������ ����������.
            xWW=0.25*(x[nvtx[0][iWW]-1]+x[nvtx[1][iWW]-1]+x[nvtx[2][iWW]-1]+x[nvtx[3][iWW]-1]);
			xW=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
			xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
			xWW2=xWW*xWW; xP2=xP*xP; xW2=xW*xW;
			
			B[0][0]=1.0; B[0][1]=xWW; B[0][2]=xWW2;
			B[1][0]=1.0; B[1][1]=xW; B[1][2]=xW2;
			B[2][0]=1.0; B[2][1]=xP; B[2][2]=xP2;
			inverse_matrix_simple(B, 3, false); // �������� �������
			Real dxee=dxe;
			xE=xP+dxe;
			xEE=xE+dxe;
			Real PWW, PW, PP;
			PWW=potent[Press][iWW]; PW=potent[Press][iW]; PP=potent[Press][iP]; 
			PE=1.0*(B[0][0]*PWW+B[0][1]*PW+B[0][2]*PP);
			PE+=xE*(B[1][0]*PWW+B[1][1]*PW+B[1][2]*PP);
			PE+=xE*xE*(B[2][0]*PWW+B[2][1]*PW+B[2][2]*PP);

			PEE=1.0*(B[0][0]*PWW+B[0][1]*PW+B[0][2]*PP);
			PEE+=xEE*(B[1][0]*PWW+B[1][1]*PW+B[1][2]*PP);
			PEE+=xEE*xEE*(B[2][0]*PWW+B[2][1]*PW+B[2][2]*PP);

            Real koef=rhoe*dy*dy/(2.0*apue); 

			FeRhie_Chow+=koef*(PEE-potent[Press][iP])/(dxe+dxee);
            FeRhie_Chow+=koef*(PE-potent[Press][iW])/(dxe+dxw);
			FeRhie_Chow-=koef*2.0*(PE-potent[Press][iP])/dxe;
		}

		if (!bW) {
            Fw=rhow*dy*(fwplus*potent[Vx][iP]+(1.0-fwplus)*potent[Vx][iW]);
			///*
			Real koef=rhow*dy*dy/(2.0*apuw);
			int iWW=sosed[WW][iP]-1;
			Real dxww=dxw;
            Real PWW=0.0;
			if (iWW>-1) {
				// ���� ���� ����������
				dxww=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
	            dxww-=0.25*(x[nvtx[0][iWW]-1]+x[nvtx[1][iWW]-1]+x[nvtx[2][iWW]-1]+x[nvtx[3][iWW]-1]);
				PWW=potent[Press][iWW];
			}
			else {
				// ���� ���� iWW �� ����������,
				// �� ����������� �������� PWW
				// ����� ���������� � ������� 
				// ������������ �������������
                Real xW, xP, xE, xW2, xP2, xE2, xWW;
                xW=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
                // ��� ������� ���������� ���� iWW ���� iE ������ ����������.
				xE=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xE2=xE*xE; xP2=xP*xP; xW2=xW*xW;
				
				B[0][0]=1.0; B[0][1]=xW; B[0][2]=xW2;
				B[1][0]=1.0; B[1][1]=xP; B[1][2]=xP2;
				B[2][0]=1.0; B[2][1]=xE; B[2][2]=xE2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				xWW=xW-dxw;
				Real PE, PP, PW;
				PE=potent[Press][iE]; PP=potent[Press][iP]; PW=potent[Press][iW];
				PWW=1.0*(B[0][0]*PW+B[0][1]*PP+B[0][2]*PE);
				PWW+=xWW*(B[1][0]*PW+B[1][1]*PP+B[1][2]*PE);
				PWW+=xWW*xWW*(B[2][0]*PW+B[2][1]*PP+B[2][2]*PE);
			}
			if (!bE) {
				FwRhie_Chow+=koef*(potent[Press][iE]-potent[Press][iW])/(dxe+dxw);
                FwRhie_Chow+=koef*(potent[Press][iP]-PWW)/(dxw+dxww);
			    FwRhie_Chow-=koef*2.0*(potent[Press][iP]-potent[Press][iW])/dxw;
			}
			else {
                // ���� iE ����.
				// ���� ������������ �������� PE � ���� iE � �������
				// ������������ �������������
                Real xP, xW, xWW, xP2, xW2, xWW2, xE;
                // ��� ������� ���������� ���� iE ���� iWW ������ ����������.
                xWW=0.25*(x[nvtx[0][iWW]-1]+x[nvtx[1][iWW]-1]+x[nvtx[2][iWW]-1]+x[nvtx[3][iWW]-1]);
				xW=0.25*(x[nvtx[0][iW]-1]+x[nvtx[1][iW]-1]+x[nvtx[2][iW]-1]+x[nvtx[3][iW]-1]);
				xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
				xP2=xP*xP; xW2=xW*xW;  xWW2=xWW*xWW;
				
				B[0][0]=1.0; B[0][1]=xWW; B[0][2]=xWW2;
				B[1][0]=1.0; B[1][1]=xW; B[1][2]=xW2;
				B[2][0]=1.0; B[2][1]=xP; B[2][2]=xP2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				xE=xP+dxe;
				Real PE, PW, PP;
			    PW=potent[Press][iW]; PP=potent[Press][iP];
				PE=1.0*(B[0][0]*PWW+B[0][1]*PW+B[0][2]*PP);
				PE+=xE*(B[1][0]*PWW+B[1][1]*PW+B[1][2]*PP);
				PE+=xE*xE*(B[2][0]*PWW+B[2][1]*PW+B[2][2]*PP); 


                // ������ ������� 
                //FwRhie_Chow+=koef*(potent[Press][iP]-PWW)/(dxw+dxww);
                //FwRhie_Chow-=koef*1.5*(potent[Press][iP]-potent[Press][iW])/dxw;
                // ����� �������
                FwRhie_Chow+=koef*(PE-potent[Press][iW])/(dxe+dxw);
                FwRhie_Chow+=koef*(potent[Press][iP]-PWW)/(dxw+dxww);
			    FwRhie_Chow-=koef*2.0*(potent[Press][iP]-potent[Press][iW])/dxw;
			}
			//*/
		}
		else {
			Fw=rhow*dy*potent[Vx][iP];

            // ����� iW � iWW ����
			// �������� � ��� ������������ 
			// � ������� ������������ �������������
			Real PWW=0.0; // �������� � ���� iWW
			Real PW=0.0; // �������� � ���� iW
			int iEE=sosed[EE][iP]-1;
			Real xP, xE, xEE, xP2, xE2, xEE2, xW, xWW;
			// ��� ������� ���������� ����� iW � iWW ���� iEE ������ ����������
            xEE=0.25*(x[nvtx[0][iEE]-1]+x[nvtx[1][iEE]-1]+x[nvtx[2][iEE]-1]+x[nvtx[3][iEE]-1]);
			xE=0.25*(x[nvtx[0][iE]-1]+x[nvtx[1][iE]-1]+x[nvtx[2][iE]-1]+x[nvtx[3][iE]-1]);
			xP=0.25*(x[nvtx[0][iP]-1]+x[nvtx[1][iP]-1]+x[nvtx[2][iP]-1]+x[nvtx[3][iP]-1]);
			xEE2=xEE*xEE; xP2=xP*xP; xE2=xE*xE;
			
			B[0][0]=1.0; B[0][1]=xP; B[0][2]=xP2;
			B[1][0]=1.0; B[1][1]=xE; B[1][2]=xE2;
			B[2][0]=1.0; B[2][1]=xEE; B[2][2]=xEE2;
			inverse_matrix_simple(B, 3, false); // �������� �������
			Real dxww=dxw;
			xW=xP-dxw;
			xWW=xW-dxw;
			Real PEE, PE, PP;
			PEE=potent[Press][iEE]; PE=potent[Press][iE]; PP=potent[Press][iP]; 
			PW=1.0*(B[0][0]*PP+B[0][1]*PE+B[0][2]*PEE);
			PW+=xW*(B[1][0]*PP+B[1][1]*PE+B[1][2]*PEE);
			PW+=xW*xW*(B[2][0]*PP+B[2][1]*PE+B[2][2]*PEE);

			PWW=1.0*(B[0][0]*PP+B[0][1]*PE+B[0][2]*PEE);
			PWW+=xWW*(B[1][0]*PP+B[1][1]*PE+B[1][2]*PEE);
			PWW+=xWW*xWW*(B[2][0]*PP+B[2][1]*PE+B[2][2]*PEE);
 
			Real koef=rhow*dy*dy/(2.0*apuw);

            FwRhie_Chow+=koef*(potent[Press][iE]-PW)/(dxe+dxw);
            FwRhie_Chow+=koef*(potent[Press][iP]-PWW)/(dxw+dxww);
			FwRhie_Chow-=koef*2.0*(potent[Press][iP]-PW)/dxw;
		}

		if (!bN) {
            Fn=rhon*dx*(fnplus*potent[Vy][iN]+(1.0-fnplus)*potent[Vy][iP]);
			///*
			Real koef=rhon*dx*dx/(2.0*apvn);
			int iNN=sosed[NN][iP]-1;
			Real dynn=dyn;
			Real PNN=0.0; // �������� � ���� iNN
			if (iNN>-1) {
				// ���� ���� ����������
				dynn=0.25*(y[nvtx[0][iNN]-1]+y[nvtx[1][iNN]-1]+y[nvtx[2][iNN]-1]+y[nvtx[3][iNN]-1]);
	            dynn-=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
				PNN=potent[Press][iNN];
			}
			else {
                // ���� ���� ������������ 
				// ����������� �������� PNN
				// ����� ���������� � ������� 
				// ������������ �������������.
				Real yN, yP, yS, yN2, yP2, yS2, yNN;
				yN=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				// ��� ������� ���������� ���� iNN ���� iS ������ ����������.
				yS=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
				yN2=yN*yN; yP2=yP*yP; yS2=yS*yS;
				
				B[0][0]=1.0; B[0][1]=yS; B[0][2]=yS2;
				B[1][0]=1.0; B[1][1]=yP; B[1][2]=yP2;
				B[2][0]=1.0; B[2][1]=yN; B[2][2]=yN2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				yNN=yN+dyn;
				Real PN, PP, PS;
				PN=potent[Press][iN]; PP=potent[Press][iP]; PS=potent[Press][iS];
				PNN=1.0*(B[0][0]*PS+B[0][1]*PP+B[0][2]*PN);
				PNN+=yNN*(B[1][0]*PS+B[1][1]*PP+B[1][2]*PN);
				PNN+=yNN*yNN*(B[2][0]*PS+B[2][1]*PP+B[2][2]*PN);
			}
			if (!bS) {
				FnRhie_Chow+=koef*(PNN-potent[Press][iP])/(dyn+dynn);
                FnRhie_Chow+=koef*(potent[Press][iN]-potent[Press][iS])/(dyn+dys);
			    FnRhie_Chow-=koef*2.0*(potent[Press][iN]-potent[Press][iP])/dyn;
			}
			else {
				// ���� iS ����.
				// ���� ������������ �������� PS � ���� iS � �������
				// ������������ �������������
                Real yP, yN, yNN, yP2, yN2, yNN2, yS;
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
				Real PS, PP, PN;
				PN=potent[Press][iN]; PP=potent[Press][iP]; 
				PS=1.0*(B[0][0]*PP+B[0][1]*PN+B[0][2]*PNN);
				PS+=yS*(B[1][0]*PP+B[1][1]*PN+B[1][2]*PNN);
				PS+=yS*yS*(B[2][0]*PP+B[2][1]*PN+B[2][2]*PNN);


				// ���������� �������
                //FnRhie_Chow+=koef*(PNN-potent[Press][iP])/(dyn+dynn);
                //FnRhie_Chow-=koef*1.5*(potent[Press][iN]-potent[Press][iP])/dyn;
				// ����� �������:
                FnRhie_Chow+=koef*(PNN-potent[Press][iP])/(dyn+dynn);
                FnRhie_Chow+=koef*(potent[Press][iN]-PS)/(dyn+dys);
			    FnRhie_Chow-=koef*2.0*(potent[Press][iN]-potent[Press][iP])/dyn;
			}
			//*/
		}
		else {
			Fn=rhon*dx*potent[Vy][iP];

            // ����� iN � iNN ����
			// �������� � ��� ������������
			// � ������� ������������ �������������
            Real PNN=0.0; // �������� � ���� iNN
			Real PN=0.0; // �������� � ���� iN
            int iSS=sosed[SS][iP]-1;
			Real ySS, yS, yP, ySS2, yS2, yP2, yN, yNN;
			// ��� ������� ���������� ���� iN � iNN ���� iSS ������ ����������.
            ySS=0.25*(y[nvtx[0][iSS]-1]+y[nvtx[1][iSS]-1]+y[nvtx[2][iSS]-1]+y[nvtx[3][iSS]-1]);
			yS=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
			yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
			ySS2=ySS*ySS; yP2=yP*yP; yS2=yS*yS;
			
			B[0][0]=1.0; B[0][1]=ySS; B[0][2]=ySS2;
			B[1][0]=1.0; B[1][1]=yS; B[1][2]=yS2;
			B[2][0]=1.0; B[2][1]=yP; B[2][2]=yP2;
			inverse_matrix_simple(B, 3, false); // �������� �������
			Real dynn=dyn;
			yN=yP+dyn;
			yNN=yN+dyn;
			Real PSS, PS, PP;
			PSS=potent[Press][iSS]; PS=potent[Press][iS]; PP=potent[Press][iP]; 
			PN=1.0*(B[0][0]*PSS+B[0][1]*PS+B[0][2]*PP);
			PN+=yN*(B[1][0]*PSS+B[1][1]*PS+B[1][2]*PP);
			PN+=yN*yN*(B[2][0]*PSS+B[2][1]*PS+B[2][2]*PP);

			PNN=1.0*(B[0][0]*PSS+B[0][1]*PS+B[0][2]*PP);
			PNN+=yNN*(B[1][0]*PSS+B[1][1]*PS+B[1][2]*PP);
			PNN+=yNN*yNN*(B[2][0]*PSS+B[2][1]*PS+B[2][2]*PP);

			Real koef=rhon*dx*dx/(2.0*apvn);

			FnRhie_Chow+=koef*(PNN-potent[Press][iP])/(dyn+dynn);
            FnRhie_Chow+=koef*(PN-potent[Press][iS])/(dyn+dys);
			FnRhie_Chow-=koef*2.0*(PN-potent[Press][iP])/dyn;
		}

		if (!bS) {
            Fs=rhos*dx*(fsplus*potent[Vy][iP]+(1.0-fsplus)*potent[Vy][iS]);
			///*
			Real koef=rhos*dx*dx/(2*apvs);
			int iSS=sosed[SS][iP]-1;
			Real dyss=dys;
			Real PSS=0.0; // �������� � ���� iSS
			if (iSS>-1) {
				dyss=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
	            dyss-=0.25*(y[nvtx[0][iSS]-1]+y[nvtx[1][iSS]-1]+y[nvtx[2][iSS]-1]+y[nvtx[3][iSS]-1]);
				PSS=potent[Press][iSS];
			}
			else {
               // ���� ���� ������������ 
				// ����������� �������� PSS
				// ����� ���������� � ������� 
				// ������������ �������������.
				Real yN, yP, yS, yN2, yP2, yS2, ySS;
				// ��� ������� ���������� ���� iSS ���� iN ������ ����������.
				yN=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				yS=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
				yN2=yN*yN; yP2=yP*yP; yS2=yS*yS;
				
				B[0][0]=1.0; B[0][1]=yS; B[0][2]=yS2;
				B[1][0]=1.0; B[1][1]=yP; B[1][2]=yP2;
				B[2][0]=1.0; B[2][1]=yN; B[2][2]=yN2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				ySS=yS-dys;
				Real PN, PP, PS;
				PN=potent[Press][iN]; PP=potent[Press][iP]; PS=potent[Press][iS];
				PSS=1.0*(B[0][0]*PS+B[0][1]*PP+B[0][2]*PN);
				PSS+=ySS*(B[1][0]*PS+B[1][1]*PP+B[1][2]*PN);
				PSS+=ySS*ySS*(B[2][0]*PS+B[2][1]*PP+B[2][2]*PN);
			}
			if (!bN) {
				FsRhie_Chow+=koef*(potent[Press][iN]-potent[Press][iS])/(dyn+dys);
                FsRhie_Chow+=koef*(potent[Press][iP]-PSS)/(dys+dyss);
			    FsRhie_Chow-=koef*2.0*(potent[Press][iP]-potent[Press][iS])/dys;
			}
			else {
                // ���� iN ����.
				// ���� ������������ �������� PN � ���� iN � �������
				// ������������ �������������
                Real yP, yS, ySS, yP2, yS2, ySS2, yN;
                // ��� ������� ���������� ���� iE ���� iWW ������ ����������.
                ySS=0.25*(y[nvtx[0][iSS]-1]+y[nvtx[1][iSS]-1]+y[nvtx[2][iSS]-1]+y[nvtx[3][iSS]-1]);
				yS=0.25*(y[nvtx[0][iS]-1]+y[nvtx[1][iS]-1]+y[nvtx[2][iS]-1]+y[nvtx[3][iS]-1]);
				yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
				yP2=yP*yP; yS2=yS*yS;  ySS2=ySS*ySS;
				
				B[0][0]=1.0; B[0][1]=ySS; B[0][2]=ySS2;
				B[1][0]=1.0; B[1][1]=yS; B[1][2]=yS2;
				B[2][0]=1.0; B[2][1]=yP; B[2][2]=yP2;
				inverse_matrix_simple(B, 3, false); // �������� �������
				yN=yP+dyn;
				Real PN, PS, PP;
			    PS=potent[Press][iS]; PP=potent[Press][iP];
				PN=1.0*(B[0][0]*PSS+B[0][1]*PS+B[0][2]*PP);
				PN+=yN*(B[1][0]*PSS+B[1][1]*PS+B[1][2]*PP);
				PN+=yN*yN*(B[2][0]*PSS+B[2][1]*PS+B[2][2]*PP); 

                // ���������� �������
                //FsRhie_Chow+=koef*(potent[Press][iP]-PSS)/(dys+dyss);
                //FsRhie_Chow-=koef*1.5*(potent[Press][iP]-potent[Press][iS])/dys;
				// ����� �������
                FsRhie_Chow+=koef*(PN-potent[Press][iS])/(dyn+dys);
                FsRhie_Chow+=koef*(potent[Press][iP]-PSS)/(dys+dyss);
			    FsRhie_Chow-=koef*2.0*(potent[Press][iP]-potent[Press][iS])/dys;
			}
			//*/
		}
		else {
			Fs=rhos*dx*potent[Vy][iP];

			// ����� iS � iSS ����
			// �������� � ��� ������������ 
			// � ������� ������������ �������������
			Real PSS=0.0; // �������� � ���� iSS
			Real PS=0.0; // �������� � ���� iS
			int iNN=sosed[NN][iP]-1;
			Real yP, yN, yNN, yP2, yN2, yNN2, yS, ySS;
			// ��� ������� ���������� ����� iS � iSS ���� iNN ������ ����������
            yNN=0.25*(y[nvtx[0][iNN]-1]+y[nvtx[1][iNN]-1]+y[nvtx[2][iNN]-1]+y[nvtx[3][iNN]-1]);
			yN=0.25*(y[nvtx[0][iN]-1]+y[nvtx[1][iN]-1]+y[nvtx[2][iN]-1]+y[nvtx[3][iN]-1]);
			yP=0.25*(y[nvtx[0][iP]-1]+y[nvtx[1][iP]-1]+y[nvtx[2][iP]-1]+y[nvtx[3][iP]-1]);
			yNN2=yNN*yNN; yP2=yP*yP; yN2=yN*yN;
			
			B[0][0]=1.0; B[0][1]=yP; B[0][2]=yP2;
			B[1][0]=1.0; B[1][1]=yN; B[1][2]=yN2;
			B[2][0]=1.0; B[2][1]=yNN; B[2][2]=yNN2;
			inverse_matrix_simple(B, 3, false); // �������� �������
			Real dyss=dys;
			yS=yP-dys;
			ySS=yS-dys;
			Real PNN, PN, PP;
			PNN=potent[Press][iNN]; PN=potent[Press][iN]; PP=potent[Press][iP]; 
			PS=1.0*(B[0][0]*PP+B[0][1]*PN+B[0][2]*PNN);
			PS+=yS*(B[1][0]*PP+B[1][1]*PN+B[1][2]*PNN);
			PS+=yS*yS*(B[2][0]*PP+B[2][1]*PN+B[2][2]*PNN);

			PSS=1.0*(B[0][0]*PP+B[0][1]*PN+B[0][2]*PNN);
			PSS+=ySS*(B[1][0]*PP+B[1][1]*PN+B[1][2]*PNN);
			PSS+=ySS*ySS*(B[2][0]*PP+B[2][1]*PN+B[2][2]*PNN);
 
			Real koef=rhos*dx*dx/(2*apvs);

            FsRhie_Chow+=koef*(potent[Press][iN]-PS)/(dyn+dys);
            FsRhie_Chow+=koef*(potent[Press][iP]-PSS)/(dys+dyss);
			FsRhie_Chow-=koef*2.0*(potent[Press][iP]-PS)/dys;
		}

		Real deltat=1e-6*dx*dy*alpha[Press]/((sl[PAm][iP].ae+sl[PAm][iP].aw+sl[PAm][iP].an+sl[PAm][iP].as)*(1-alpha[Press]));
		sl[PAm][iP].b=Fw-Fe+Fs-Fn+deltat*(FwRhie_Chow-FeRhie_Chow+FsRhie_Chow-FnRhie_Chow);
        rhie_chow[0][iP]=Fw-Fe+Fs-Fn;
		rhie_chow[1][iP]=deltat*(FwRhie_Chow-FeRhie_Chow+FsRhie_Chow-FnRhie_Chow);
        

	}


} // my_elmatr_quad_PAm

#endif