// ���� my_approx_convective.c 
// ������������� ������������� �����

#ifndef MY_APPROX_CONVECTIVE_C
#define MY_APPROX_CONVECTIVE_C 1

#include <math.h>

#define Real double

// ����� ��� ������������� ���������-��������
#define CR 1 // ����������-����������
#define UDS 2 // ��������������� ������� �������
#define COMB 3 // ��������������� 
#define POLY 4 // �������������� C. ���������
#define EXP 5 // ���������������� �����
#define BULG 6 // ����� �.�. ��������� (23) �� ������
#define POW 7 // �������������

// ���������� ������������ �� ���� 
// ������������ �����.
/*
Real fmax(Real fA, Real fB) {
	Real r=fB;
	if (fA > fB) r=fA;
	return r;
} // fmax
*/

Real lim(Real Pe) {
	// 0.2 ��� ����� DavisTest Ra=10^5
	   //return 0.1*1.0 / (1.0 + 0.6712 * sqrt(fabs(Pe)));
	//return  1.0 / (1.0 + 7.0 * sqrt(fabs(Pe)));
	//0.06

	// Re = 94; 1.2; 3200;
	// 0.05 ������
	//return (Pe < 2000.0 ? 0.05: (Pe < 3000.0 ? 0.05 - 0.04*0.001*(Pe-2000.0): 0.005));//
	//return (Pe < 2000.0 ? 0.05 : 0.005);

	// �� �������. ��������� 25.08.2021
	// �� �������� ������������� ������ �� �������� glRCh.
	return glRCh;
	//return 0.1;
}

// ������� A(|P|) 
// ��� ��������� ����:
// fabsPe - ������ ����� �����.
// ishconvection - ����� �����:
// CR==1 - ����������-���������� �����,
// UDS==2 - � ���������� ������ ������,
// COMB==3 - ���������������,
// POLY==4 - �� ��������� �������,
// EXP==5 - ���������������� (������),
// BULG==6 - ����� �.�. ��������� (23) �� ������,
// POW==7 - ������������� �����������.
Real ApproxConvective(Real fabsPe, int ishconvection) {

	if (ishconvection == UDS) {

		// 16.08.2021
		// ����� ������������ ������� �� ������� ���������.
		// ������ �������. ������� �����.
		return 1.0;
	}
	else {

		Real r = 1; // �� ��������� � ��������� ������ ������
		Real rCR, rpoly, rpoly2, rpoly5;
		if (ishconvection < 6) {
			rCR = 1.0 - 0.5 * fabsPe;
			rpoly = 1.0 - 0.1 * fabsPe;
			rpoly2 = rpoly * rpoly;
			rpoly5 = rpoly2 * rpoly2 * rpoly; // �� ��������� �������
		}

		switch (ishconvection) {
		case CR: // ���������� ���������� ����� �������� ������ 
				// ��� ��������� ����� ����� ������� �� ������ 1.0.
			if (fabsPe < 1.0) r = rCR;
			// ������� ��� ������ ����� �� ������ ������� 1.0 ��
			// ������� 10.0 �������������� ������� �� �������������
			// ���������� ����� �� ��������� �������.
			if ((fabsPe >= 1.0) && (fabsPe < 10.0)) r = rpoly5;
			else r = 0.0; // ��� ������ ����� ������� 10.0 ������������� ����.
			break;
		case UDS: // � ���������� ������ ������
				//if (fabsPe<1.0) r=1.0; else r=0.0;
			r = 1.0; // UDS ������ �������, ������� �����.
			break;
		case COMB: // ���������������
			r = fmax(0.0, rCR);
			break;
		case POLY: // �� ��������� �������
			if (fabsPe < 10.0) r = rpoly5; else r = 0.0;
			break;
		case EXP: // ���������������� (������)
			if (fabsPe < 10.0) {
				if (fabsPe < 0.01) r = rCR; else r = fabsPe / (exp(fabsPe) - 1.0);
			}
			else r = 0.0;
			break;
		case BULG: // �.�. ��������, �.�. ��������
				// ���������� ��������������� ����������� �����������
				// ����������� (23) �� ������.
			r = 1.0 / (1.0 + 0.6712 * fabsPe * fabsPe);
			break;
		case POW: // ������������� �����������
			r = pow(0.553, fabsPe);
			break;
		default: r = 1.0 / (1.0 + 0.6712 * fabsPe * fabsPe); break;
		}

		return r;
	}
} // ApproxConvective

#endif