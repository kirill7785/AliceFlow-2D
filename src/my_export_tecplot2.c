// ���� my_export_tecplot2.c �������� �����������
// ������������� � ��������� tecplot360

#pragma once
#ifndef MY_EXPORT_TECPLOT2_C
#define MY_EXPORT_TECPLOT2_C 1

// �������� ���������� �����
// ������� ���������� ������� � ��������� tecplot360
// ������������ �������� � ��� ����������� � ��� ���������� �����.
void exporttecplotxy360(int nve, int maxelm, int ncell, int** nvtx, int** nvtxcell, Real* x, Real* y, Real** potent)
{
	FILE *fp;
	errno_t err;
	// �������� ����� ��� ������.
	if ((err = fopen_s( &fp, "ALICEFLOW0_02.PLT", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"ALICEFLOW0_02\"\n");

		// ������ ��� ����������
		fprintf(fp, "VARIABLES = x, y,  \"UX [m/s]\", \"UY [m/s]\", \"Speed [m/s]\", \"Pressure [N/m2]\", \"Temperature [C]\", \"Nested desection\" \n");

		// ������ ���������� � �����
		if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
        if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);

		int i=0; // �������� 
		int j=0; // ����� for

		// ������ x
	    for (i=0; i<maxelm; i++) {	
				fprintf(fp, "%e ", 0.25*(x[nvtx[0][i]-1]+x[nvtx[1][i]-1]+x[nvtx[2][i]-1]+x[nvtx[3][i]-1]));
				if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// ������ y
		for (i=0;i<maxelm; i++) {
		 	fprintf(fp, "%e ", 0.25*(y[nvtx[0][i]-1]+y[nvtx[1][i]-1]+y[nvtx[2][i]-1]+y[nvtx[3][i]-1]));
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		
		
		// ������ �������������� ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[Vx][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");
        // ������ ������������ ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[Vy][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

		// ������ ������ ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", sqrt(potent[Vx][i]*potent[Vx][i]+potent[Vy][i]*potent[Vy][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");

		// ������ ��������
		for (i=0;i<maxelm; i++) {
			fprintf(fp, "%e ", potent[Press][i]);
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");        
		
		// ������ �����������
		for (i = 0; i < maxelm; i++) {
			fprintf(fp, "%e ", potent[Temp][i]);
			if (i % 10 == 0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");

		// Nested desection
		if (inumcore == 2) {
			for (int i = s_par[1].s; i < s_par[1].e; ++i) {
				fprintf(fp, "%e ", 1.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			for (int i = s_par[2].s; i < s_par[2].e; ++i) {
				fprintf(fp, "%e ", 2.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
			for (int i = s_par[3].s; i < s_par[3].e; ++i) {
				fprintf(fp, "%e ", 3.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
		}
		else {
			for (i = 0; i < maxelm; i++) {
				fprintf(fp, "%e ", 1.0);
				if (i % 10 == 0) fprintf(fp, "\n");
			}
		}
		fprintf(fp, "\n");

		// ������ ���������� � ���������� �����
		for (i=0;i<ncell; i++) {
			if (nve==4) fprintf(fp, "%d %d %d %d\n", nvtxcell[0][i], nvtxcell[1][i], nvtxcell[2][i], nvtxcell[3][i]);
		}

		fclose(fp); // �������� �����
        //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_02.PLT",SW_NORMAL);
		WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_02.PLT", SW_NORMAL);
	}
} // exporttecplotxy360


#endif
