// ���� my_approx_convective2.c 
// �������� ������� ��� ������������� 
// ������������� � ������������� �����.

#ifndef MY_APPROX_CONVECTIVE2_C
#define MY_APPROX_CONVECTIVE2_C 1

#define Real double // float

// ����� ��� ������������� ���������-��������
#define CR2 100 // ����������� ��������
#define UDS2 101 // ��������������� ������� �������
#define EXP2 102 // ���������������� ����� (������)
#define KUD 103 // ����� ������������ � ����������� �������� ����� ���������

Real fC(Real Pe, int isheme) {
	Real r=1.0;
	const Real Pemax=12.0;

	switch (isheme) {
		case CR2 : r=0.5; break;
		case UDS2 : if (Pe<0.0) r=1.0; else r=0.0; break;
		case EXP2 : if (fabs(Pe) < 0.01) r=0.5; else r=(exp(0.5*Pe)-1.0)/(exp(Pe)-1.0); break;
		case KUD : if (fabs(Pe)<=Pemax) r=0.5*(1.0-Pe/Pemax); else r=0.5*(1.0-Pe/fabs(Pe)); break;
		default : if (fabs(Pe)<=Pemax) r=0.5*(1.0-Pe/Pemax); else r=0.5*(1.0-Pe/fabs(Pe)); break;
	}
	return (r);
} // ������������� ���������

Real fD(Real Pe, int isheme) {
	Real r=1.0;
	const Real Pemax=12.0;

	switch (isheme) {
		case CR2 : r=1.0; break;
		case UDS2 : r=1.0; break;
		case EXP2 : if (fabs(Pe) < 0.01) r=1.0; else r=(Pe*exp(0.5*Pe))/(exp(Pe)-1.0); break;
		case KUD : if (fabs(Pe)<=Pemax) r=1.0-(Pe*Pe)/(Pemax*Pemax); else r=0.0; break;
		default : if (fabs(Pe)<=Pemax) r=1.0-(Pe*Pe)/(Pemax*Pemax); else r=0.0; break;
	}
	return (r);
} // ������������� ��������

#endif