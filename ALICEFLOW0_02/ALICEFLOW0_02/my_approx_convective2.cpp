// Файл my_approx_convective2.cpp 
// содержит функции для аппроксимации 
// конвективного и диффузионного члена.

#ifndef MY_APPROX_CONVECTIVE2_CPP
#define MY_APPROX_CONVECTIVE2_CPP 1

#define Real double // float

// схемы для аппроксимации конвекции-диффузии
#define CR2 100 // Центральные разности
#define UDS2 101 // Противопоточная первого порядка
#define EXP2 102 // экспоненциальная схема (точная)
#define KUD 103 // схема предложенная в диссертации Кудинова Павла Ивановича

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
} // аппроксимация конвекции

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
} // аппроксимация диффузии

#endif