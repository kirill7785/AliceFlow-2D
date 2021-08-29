// Файл my_linalg.c
// самостоятенльная реализация некоторых функций линейной алгебры.

#pragma once
#ifndef MY_LINALG_C
#define MY_LINALG_C 1


#include <stdio.h> // для функции getchar
#include <stdlib.h> // Для функции exit, atoi, atof
#include <math.h> // математические функции sqrt, fabs
#include <omp.h> // OpenMP
#include "my_linalg.h" // самописные функции линейной алгебры


#define Real double // модель веществекнного числа


const Real dterminatedTResudual = 1e-15; // для МСГ Congruate Gradients

/*  Решает систему уравнений для квадратной
 *  несимметричной матрицы коэффициентов A
 *        A*x=b;
 *  где A размеров nodesxnodes. Нумерация 
 *  элементов везде начинается с нуля.
 *  Процедура представляет собой метод Гаусса
 *  без выбора главного элемента и без учёта
 *  разреженности матрицы.
 *  A и b не сохраняются. 
*/
void eqsolve_simple_gauss(Real **A, int nodes, Real *b, Real *x) {
   int i=0, j=0, k=0; // счётчики цикла for
   const Real epsilon = 1e-100;
   Real M, sum, akk;

   omp_set_num_threads(inumcore); // установка числа потоков

   // приведение к треугольному виду:
   for(k=0; k<nodes; k++){
	   akk=A[k][k];
       if(fabs(akk)<epsilon){
		  // решение не может быть получено, т.к.
		  // на диагонали находится ноль.
	      printf("\nSolution is not exist! Gauss divizion by zero...\n");
	      system("pause");
	      exit(0);
	   }
       #pragma omp parallel for shared(k, nodes, A, b) private(i,j,M) firstprivate(akk)
       for(i=k+1; i<nodes; i++) {
	      
	      M = A[i][k] / akk;
	      for(j=k; j<nodes; j++){
	          A[i][j] -= M * A[k][j];
	      }
	      b[i] -= M*b[k];
       }
   }
   // процесс обратного исключения
   x[nodes-1]=b[nodes-1]/A[nodes-1][nodes-1];
   for(i=nodes-2; i>=0; i--){
       sum = 0.0;
       #pragma omp parallel for shared(A,x,i,nodes) private(j) reduction (+: sum)
       for(j = i+1; j<nodes; j++){
	       sum+= A[i][j]*x[j];
       }
       x[i] = (b[i] - sum) / A[i][i];
   }
} // eqsolve_simple_gauss

/*  Решает систему уравнений для квадратной
 *  симметричной положительно определённой
 *  (с диагональным преобладанием) матрицы
 *  коэффициентов А:
 *        A*x=b;
 *  где A размеров nodesxnodes. Матрица А
 *  предполагается не разреженной. Нумерация 
 *  элементов везде начинается с нуля.
 *  Процедура представляет собой разложение Холесского:
 *        A=L*transpose(L),
 *  после которого выполняются прямое исключение и 
 *  обратная подстановка. A и b не сохраняются. 
*/
void eqsolv_simple_holesskii(Real **A, int nodes, Real *b, Real *x) {
	// Разложение Холесского: замена A верхним и нижним 
	// треугольными множителями.
	A[0][0]=sqrt(A[0][0]);
	A[1][0]/=A[0][0];
	A[0][1]=A[1][0];
	A[1][1]=sqrt(A[1][1]-A[1][0]*A[1][0]);

	omp_set_num_threads(inumcore); // установка числа потоков

	int irow,irow1;
	int icol, icol1;
	Real sum;
	int k;
	for (irow=2; irow<nodes; irow++) {
		irow1=irow-1;
		A[irow][0]/=A[0][0];
        A[0][irow]=A[irow][0];
        #pragma omp parallel for shared(irow1,A) private(icol, icol1, sum, k)
		for (icol=1; icol<=irow1; icol++) {
			icol1=icol-1;
            sum=0.0;   
            for (k=0; k<=icol1; k++) sum+=A[irow][k]*A[icol][k];
			A[irow][icol]=(A[irow][icol]-sum)/A[icol][icol];
			A[icol][irow]=A[irow][icol];
		}
		sum=0.0;
		#pragma omp parallel for shared(A,irow,irow1) private(k) reduction (+: sum)
		for (k=0; k<=irow1; k++) sum+=A[irow][k]*A[irow][k];
		A[irow][irow]=sqrt(A[irow][irow]-sum);
	}
    
	// Прямое исключение. Происходит разрушение правой части
	b[0]/=A[0][0];

	for (irow=1; irow<nodes; irow++) {
		irow1=irow-1;
		sum=0.0;
		#pragma omp parallel for shared(A,b,irow,irow1) private(icol) reduction (+: sum)
		for (icol=0; icol<=irow1; icol++) sum+=A[irow][icol]*b[icol];
        b[irow]=(b[irow]-sum)/A[irow][irow];
	}

	// Обратная подстановка используется верхний треугольный множитель
	x[nodes-1]=b[nodes-1]/A[nodes-1][nodes-1];
	for (k=1; k<=nodes; k++) {
		irow=nodes+1-k-1;
		irow1=irow+1;
		sum=0.0;
        #pragma omp parallel for shared(A,x,irow,irow1,nodes) private(icol) reduction (+: sum)
		for (icol=irow1; icol<nodes; icol++) sum+=A[irow][icol]*x[icol];
		x[irow]=(b[irow]-sum)/A[irow][irow];
	}

} // eqsolv_simple_holesskii

/* Находит обратную матрицу для 
*  квадратной матрицы A nodes*nodes с 
*  ненулевыми элементами на главной диагонали.
*  Решение производится путём метода исключения
*  Гаусса, а именно решая nodes СЛАУ. 
*          A*inv=e
*  Приведение  к треугольному виду делается
*  только один раз.
* Если flag==true, то матрица уже приведена к верхнетреугольному виду.
*/
void inverse_matrix_simple(Real**  A, int nodes, bool flag) {

    const Real epsilon = 1e-100;

	Real **e; // единичная матрица правых частей.
	Real **inv; // будущая обратная матрица

	int i1=0, j1=0, k1=0;
	e = new Real* [nodes];
    for (i1=0; i1<nodes; i1++) e[i1]=new Real[nodes]; 
	inv = new Real* [nodes];
    for (i1=0; i1<nodes; i1++) inv[i1]=new Real[nodes];
    
	// инициализация
	for (i1=0; i1<nodes; i1++) for (j1=0; j1<nodes; j1++) {
		inv[i1][j1]=0.0; // обратная матрица
		e[i1][j1]=0.0; // правые части
	}
	for (i1=0; i1<nodes; i1++) e[i1][i1]=1.0;


    
	if (!flag) { // если матрица ещё не приведена к верхнетреугольному виду
        Real M;
		// приведение к верхне треугольному виду:
        for(k1=0; k1<nodes; k1++){
           for(i1=k1+1; i1<nodes; i1++){
		       // Если на диагонали ноль:
		       if (fabs(A[k1][k1])<epsilon) {
			      // решение не может быть получено, т.к.
			      // на диагонали находится ноль.
				   std::cout << "diag = "  << A[k1][k1] << std::endl;
	              printf("\nSolution is not exist.\n");
	              system("pause");
		          exit(0);
		       }
	           M = A[i1][k1] / A[k1][k1];
	           for(j1=k1; j1<nodes; j1++){
	              A[i1][j1] -= M * A[k1][j1];
	           }
		       // преобразование правых частей:
              for(j1=0; j1<nodes; j1++) e[i1][j1] -= M*e[k1][j1];
           }
        }
	}
	Real *sum=new Real[nodes];

   // процесс обратного исключения
   for(i1=nodes-1; i1>=0; i1--){
	   // инициализация
       for (k1=0; k1<nodes; k1++) sum[k1] = 0.0;

       for(j1 = i1+1; j1<nodes; j1++){
		   for (k1=0; k1<nodes; k1++) {
	           sum[k1]+= A[i1][j1]*inv[j1][k1];
		   }
       }
       for (k1=0; k1<nodes; k1++) {
	        inv[i1][k1]=(e[i1][k1] - sum[k1])/A[i1][i1];
	   }
   }
   for(i1=nodes-1; i1>=0; i1--) delete[] e[i1];
   delete[] e;

   for(k1=0; k1<nodes; k1++){
       for(i1=0; i1<nodes; i1++){
		   if (inv[k1][i1] != inv[k1][i1]) {
			   std::cout << "bad inverse 3*3\n";
		   }
		   A[k1][i1]=inv[k1][i1];
	   }
   }
   for(i1=nodes-1; i1>=0; i1--) delete[] inv[i1];
   delete[] inv;
   delete[] sum;
} // inverse_matrix_simple
 
/* Находит произведение двух квадратных
* матриц A и B размерами nodesxnodes 
*             C=A*B. 
* Результат  записывается в матрицу B.
*/
void multiply_matrix_simple(Real **A, Real **B, int nodes) {
	int i1=0, j1=0, k1=0; // счётчики цикла for
	
	Real **c;
	c = new Real* [nodes];
    for (i1=0; i1<nodes; i1++) c[i1]=new Real[nodes];

	for (i1=0; i1<nodes; i1++) for (j1=0; j1<nodes; j1++) c[i1][j1]=0.0; // инициализация

	// умножение C=A*B:
    for (i1=0; i1 < nodes; i1++)
        for (k1=0; k1 < nodes; k1++)
            for (j1=0; j1 < nodes; j1++)
                c[i1][k1]+=(A[i1][j1])*(B[j1][k1]);

	// копирование результата в B:
    for (i1=0; i1<nodes; i1++) for (j1=0; j1<nodes; j1++) B[i1][j1]=c[i1][j1];

	delete[] c;
} // multiply_matrix_simple


// Следующие несколько функций (шесть, но теперь
// транспонирование не используется, а умножение заменено более быстрым). 
// используются как вспомогательные для решения
// полной проблемы собственных значений:

/* 1. Умножение квадратных матриц размера nxn:
*                t=m*p.
* Нумерация начинается с нуля.
* По окончании работы результат хранится в матрице t.
*/
void multi_m(Real **m, Real **p, Real **t, int n) {
    for (int i = 0; i < n; i++)
       for (int j = 0; j < n; j++) {
           Real s = 0;
           for (int l = 0; l < n; l++)
               s += m[i][l]*p[l][j];
           t[i][j] = s;
    }
} // multi_m 

/* 2. Транспонирование квадратной матрицы m
*  размером nxn. По окончании работы в матрице 
*  m хранится результат транспонирования.
*/
void tr_m(Real **m, int n) {
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++) {
            Real buf = m[i][j];
            m[i][j] = m[j][i];
            m[j][i] = buf;
        }
} // tr_m

/* 3. Возвращает максимальный внедиагональный
* элемент для симметричной матрицы A размером 
* nxn. Позиция максимального элемента A[f][g].
* Это медленная реализация, т.к. она не использует
* информацию о предыдущих поисках максимального
* элемента в матрице А.
*/
Real max_el(Real **A, int n, int& f, int& g) {
   Real max = A[0][1];
   f=0; g=1; // стартовое значение
   for (int j = 1; j < n; j++)
      for (int i = 0; i < j; i++) {
        if (A[i][j] > max) {
            max = A[i][j];
            f = i; g = j;
        }
    }
    return max;
 } // max_el

/* 4. Копирует вторую матрицу в первую: A=B.
* Матрицы квадратные размером nxn
*/
void matr_copy(Real **A, Real **B, int n) {
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
		  A[i][j]=B[i][j];
}

/* 5. Быстрое умножение двух квадратных матриц специального 
* вида размера nxn (левое умножение):
*                A=hiс*A.
* Здесь hic -  несимметричная транспонированная матрица вращения:
* hic[f][f] = cosfi;
* hic[g][g] = cosfi;
* hic[f][g] = +sinfi;
* hic[g][f] = -sinfi;
* Здесь f и g позиции ненулевых элементов.
* Нумерация начинается с нуля.
* По окончании работы результат хранится в исходной матрице A.
* Теперь матрица hiс передаётся только как свои четыре особых элемента
* что позволяет существенно экономить память и быстродействие.
*/
void multi_m_left(Real **A, Real **rab, int n, int f, int g, Real cosfi, Real sinfi) {
	/* Устаревший неэффективный но понятный код:
    for (int i = 0; i < n; i++)
       for (int j = 0; j < n; j++) {
		   if ((i!=f) && (i!=g)) {
			   t[i][j]=A[i][j];
		   }
		   else if (i==f) {
			   //t[i][j]=hic[f][f]*A[f][j]+hic[f][g]*A[g][j];
               t[i][j]=cosfi*A[f][j]+sinfi*A[g][j];
		   }
		   else if (i==g) {
			   //t[i][j]=hic[g][f]*A[f][j]+hic[g][g]*A[g][j];
			   t[i][j]=-sinfi*A[f][j]+cosfi*A[g][j];
		   }
    }
	*/
    
	// Теперь результат умножения возвращается прямо в матрице А
	// В качестве рабочего используется массив rab размерности
	// 2xn. Трудоёмкость операции всего 4*n умножений.
	for (int j = 0; j < n; j++) {
	   rab[0][j]=cosfi*A[f][j]+sinfi*A[g][j];
	   rab[1][j]=-sinfi*A[f][j]+cosfi*A[g][j];
	}
    for (int j = 0; j < n; j++) {
	   A[f][j]=rab[0][j];
	   A[g][j]=rab[1][j];
	}

} // multi_m_left 

/* 6. Быстрое умножение двух квадратных матриц специального 
* вида размера nxn (правое умножение):
*                A=A*hi.
* Здесь hi - несимметричная матрица вращения:
* hi[f][f] = cosfi;
* hi[g][g] = cosfi;
* hi[f][g] = -sinfi;
* hi[g][f] = +sinfi;
* Здесь f и g позиции ненулевых элементов.
* Нумерация начинается с нуля.
* По окончании работы результат хранится в исходной матрице A.
* Теперь матрица hi передаётся только как свои четыре особых элемента
* что позволяет существенно экономить память и быстродействие.
*/
void multi_m_right(Real **A, Real **rab, int n, int f, int g, Real cosfi, Real sinfi) {
	/* Неэффективное умножение
    for (int i = 0; i < n; i++)
       for (int j = 0; j < n; j++) {
		   if ((j!=f) && (j!=g)) {
			   t[i][j]=A[i][j];
		   }
		   else if (j==f) {
			   //t[i][j]=A[i][f]*hi[f][f]+A[i][g]*hi[g][f];
               t[i][j]=A[i][f]*cosfi+A[i][g]*sinfi;
		   }
		   else if (j==g) {
			   //t[i][j]=A[i][f]*hi[f][g]+A[i][g]*hi[g][g];
			   t[i][j]=-A[i][f]*sinfi+A[i][g]*cosfi;
		   }
    }
	*/

	// Теперь результат умножения возвращается прямо в матрице А
	// В качестве рабочего используется массив rab размерности
	// 2xn. Трудоёмкость операции всего 4*n умножений.
	for (int i = 0; i < n; i++) {
	   rab[0][i]=A[i][f]*cosfi+A[i][g]*sinfi; // f
	   rab[1][i]=-A[i][f]*sinfi+A[i][g]*cosfi; // g
	}
    for (int i = 0; i < n; i++) {
		A[i][f]=rab[0][i];
		A[i][g]=rab[1][i];
	}

} // multi_m_right 


/* Оригинальный алгоритм Якоби от 1846 года. 
* Решает полную проблему собственных значений в форме
*              A-lambda_scal*E=0
*  методом вращений. См. Ревизников лекции книга.
*  Симметричная положительно определённая матрица A 
*  размером nodesxnodes. Матрица A в результате работы  
*  портится ( на диагонали у её испорченного варианта будут СЗ).
*  В квадратной матрице U по столбцам хранятся 
*  собственные векторы. В векторе lambda находится список
*  собственных значений.
*  Процесс нахождения векторов и СЗ является итерационным,
*  Его точность характеризуется значением epsilon.
*  На каждой итерации делается 12xnodes матричных умножений.
*  Дополнительная память равна 2xnodes.
*  EIGEM - метод Якоби.
*/
void jacobi_matrix_simple(Real **A, Real **U, Real *lambda, int nodes, Real epsilon) {

	// значение этой постоянной нужно подобрать импирически.
    const Real eps=1e-10; // точность  с которой элементы проверяются на равенство,
	
	int i,j; // счётчики цикла for
    int im , jm; // позиция максимального элемента
    int p = 1; // номер итерации
	Real maxij; // максимальный элемент
	Real fi; // значение угла
	Real cosfi, sinfi; // значение косинуса и синуса угла fi

	/* При быстром умножении матриц этот код не требуется
	// матрица  вращения
	Real **hi=new Real*[nodes];
    for (i = 0; i < nodes; i++) hi[i]=new Real[nodes];

	// инициализация всей матрицы выполняется один раз:
    for (i = 0; i < nodes; i++)
         for (j = 0; j < nodes; j++) {
            if (i == j)
                hi[i][j] = 1.0;
            else hi[i][j] = 0.0;
    }

	// вспомогательная матрица вращения (копия)
	Real **hic=new Real*[nodes];
    for (i = 0; i < nodes; i++) hic[i]=new Real[nodes];

	// инициализация вспомогательной матрицы один раз:
	matr_copy(hic,hi,nodes); // тоже единичная матрица 
	*/

	// вспомогательная матрица для умножения
    // Устаревший неэффективный по памяти вариант.
    //Real **b=new Real*[nodes];
    //for (i = 0; i < nodes; i++) b[i]=new Real[nodes];

	// Рабочий массив.
	Real **rab=new Real*[2];
    for (i = 0; i < 2; i++) rab[i]=new Real[nodes];

    maxij = max_el(A, nodes,im,jm);
    
	// каждую итерацию делается 12xnodes умножений.
    while (fabs(maxij) > epsilon) {
       
       
	   // Вычисление угла:
	   if (fabs(A[im][im]-A[jm][jm])<eps) {
		   // особый случай значения равны
		   fi=3.141/4.0;
	   }
	   else fi= atan(2*maxij/(A[im][im]-A[jm][jm]))/2;
       
       // Нахождение тригонометрических функций
	   // от угла fi:
       cosfi = cos(fi);
	   sinfi = sin(fi);
 
	   /* При быстром умножении этот закоментированный код не используется.
	   // матрица вращения не является симметричной:
	   // инициализация матрицы вращения
       hi[im][im] = cosfi;
       hi[jm][jm] = cosfi;
       hi[im][jm] = -sinfi;
       hi[jm][im] = sinfi;
	   // транспонированный вариант: 
       hic[im][im] = cosfi;
       hic[jm][jm] = cosfi;
       hic[im][jm] = +sinfi; // транспонирование.
       hic[jm][im] = -sinfi;
	   */
 
       //  инициализация матрицы СВ которые хранятся по столбцам
	   if (p==1) {
		   //matr_copy(U,hi,nodes);
           for (i = 0; i < nodes; i++)
               for (j = 0; j < nodes; j++) {
                   if (i == j)
                      U[i][j] = 1.0;
                   else U[i][j] = 0.0;
               }
		   U[im][im] = cosfi;
           U[jm][jm] = cosfi;
           U[im][jm] = -sinfi;
           U[jm][im] = sinfi;

	   } else {
            //multi_m(U,hi,b, nodes);
            multi_m_right(U, rab, nodes, im, jm, cosfi, sinfi); // Быстрое умножение
			//matr_copy(U,b,nodes); // теперь вместо промежуточной матрицы b используется 
			// экономичная матрица rab. Эффективность операции 4xnodes умножений.
	   }
      
       //multi_m(hic, A, b, nodes); // b=transpose(H)*A
	   multi_m_left(A, rab, nodes, im, jm, cosfi, sinfi); // Быстрое умножение: 4xnodes операций умножения
       //multi_m(b, hi, A, nodes); // A=b*H.
	   multi_m_right(A, rab, nodes, im, jm, cosfi, sinfi); // Быстрое умножение: 4xnodes операций умножения
 
	   /* При быстром умножении этот закоментированный код не используется.
	   // восстановление матриц вращения:
       hi[im][im] = 1.0;
       hi[jm][jm] = 1.0;
       hi[im][jm] = 0.0;
       hi[jm][im] = 0.0;
	   // восстановление копии матрицы вращения:
       hic[im][im] = 1.0;
       hic[jm][jm] = 1.0;
       hic[im][jm] = 0.0;
       hic[jm][im] = 0.0;
	   */

	   maxij = max_el(A, nodes,im,jm); // определение максимума ресурсоёмкая операция.
       p++; // переход к следующей итерации

    } // while

    for (i = 0; i < nodes; i++) lambda[i]=A[i][i]; //  СЗ

} // jacobi_matrix_simple




// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
template <typename MY_IND_TYPE>
void swapCSIR(MY_IND_TYPE*& v, Real*& dr, int i, int j)
{
	MY_IND_TYPE tempi;
	Real tempr;

	// change v[i] <-> v[j]
	tempi = v[i];
	v[i] = v[j];
	v[j] = tempi;
	// change dr[i] <-> dr[j]
	tempr = dr[i];
	dr[i] = dr[j];
	dr[j] = tempr;

} // swap

// Вот алгоритм PivotList
template <typename MY_IND_TYPE>
int PivotListCSIR(MY_IND_TYPE*& jptr, Real*& altr, int first, int last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	int PivotValue = jptr[first];
	int PivotPoint = first;

	for (int index = (first + 1); index <= last; index++) {
		if (jptr[index] < PivotValue) {
			PivotPoint++;
			swapCSIR(jptr, altr, PivotPoint, index);
		}
	}

	swapCSIR(jptr, altr, first, PivotPoint);

	return PivotPoint;
} // PivotListCSIR


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
template <typename MY_IND_TYPE>
void QuickSortCSIR(MY_IND_TYPE*& jptr, Real*& altr, int first, int last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	if (0) {
		// BubbleSort
		int numberOfPairs = last - first + 1;
		bool swappedElements = true;
		while (swappedElements) {
			numberOfPairs--;
			swappedElements = false;
			for (int i = first; i <= first + numberOfPairs - 1; i++) {
				if (jptr[i] > jptr[i + 1]) {
					swapCSIR(jptr, altr, i, i + 1);
					swappedElements = true;
				}
			}
		}
	}
	else
	{
		int pivot;

		if (first < last) {
			pivot = PivotListCSIR(jptr, altr, first, last);
			QuickSortCSIR(jptr, altr, first, pivot - 1);
			QuickSortCSIR(jptr, altr, pivot + 1, last);
		}
	}
} // QuickSortCSIR


// Реализация на связном списке
 // Добавляет ненулевой элемент в
 // простейшую разряженную матрицу M
 // Проверки на равенство добавляемого элемента нулю нет, поэтому
 // может добавить и нулевой элемент.
void addelmsimplesparse_Stress_ell(Real aij, int i, int j, bool bset, bool bsetD) {
	const Real MY_ZERO_TOLERANCE = -1.0; // 1.0e-300; -1.0 Отключена проверка. Так нужно т.к. при условиях Дирихле для симметричности портрета нужны нули.

	//NONZEROELEM* p;
	//p = M.root[i];// Корневой элемент строки i
	// линейный поиск элемента с ключом key
	//while ((p != nullptr) && (p->key != j)) p = p->next;
	int iscan = 0;
	while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1) && (coll_ell[i][iscan] != j)) ++iscan;
	if ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] == j)) {
		// элемент найден
		if (bsetD) {
			// Удалить строку.
			data_ell[i][iscan] = aij;
			if (iscan > 0) {
				std::cout << "Fatal error in addelmsimplesparse_Stress_ell\n";
				system("PAUSE");
			}
		}
		else {
			if (fabs(aij) > MY_ZERO_TOLERANCE) {
				if (bset) data_ell[i][iscan] = aij; // установка
				else {
					data_ell[i][iscan] += aij; // добавление
				}
			}
		}
	}
	else
	{
		// если такого элемента нет в списке
		// то добавление элемента в начало списка.
		if (fabs(aij) > MY_ZERO_TOLERANCE) {
			if ((iscan < MAX_STRING_LENGTH_ELL)) {

				coll_ell[i][iscan] = j;
				data_ell[i][iscan] = aij;

#pragma omp atomic
				nnz_ell++;

			}
			else {
				std::cout << "ERROR!!! UVELICHTE MAX_STRING_LENGTH_ELL\n";
				system("PAUSE");
				exit(1);
			}
		}
	}
} // addelmsimplesparse_ell

// Удаляет положительные внедиагональные элементы.
void patch_Ell_zero(int isize)
{
	for (int j = 0; j < isize; ++j) {
		int i = j;

		Real dsum = 0.0;

		int iscan = 0;
		int ics = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1)) {


			if (j != coll_ell[i][iscan]) {
				if (data_ell[i][iscan] > 0.0) {

					data_ell[i][iscan] = 0.0;

					/*
					for (int iscan3 = iscan; iscan3 < MAX_STRING_LENGTH_ELL - 1; ++iscan3) {
						if (coll_ell[i][iscan3] == -1) break;
						data_ell[i][iscan3] = data_ell[i][iscan3 + 1];
						coll_ell[i][iscan3] = coll_ell[i][iscan3 + 1];
					}
					data_ell[i][MAX_STRING_LENGTH_ELL - 1] = 0.0;
					coll_ell[i][MAX_STRING_LENGTH_ELL - 1] = -1;
					iscan--;
					*/
				}
				else {
					dsum += fabs(data_ell[i][iscan]);
					ics++;
				}
			}

			++iscan;
		}

		if (ics > 0) {
			data_ell[i][0] = dsum; // коррекция диагонали для соблюдения баланса.
		}

		/*iscan = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1)) {
			std::cout << "a[" << i << "," << coll_ell[i][iscan] << "]=" << data_ell[i][iscan] << "  ";
			iscan++;
		}
		std::cout << std::endl;
		getchar();*/

	}
	//getchar();
}


// Диагональ должна стоять первой в строке.
void patch_Ell(int isize)
{
	for (int j = 0; j < isize; ++j) {
		int i = j;

		int iscan = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1) && (coll_ell[i][iscan] != j)) ++iscan;
		if ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] == j)) {
			int itmp = coll_ell[i][0];
			coll_ell[i][0] = coll_ell[i][iscan];
			coll_ell[i][iscan] = itmp;

			Real dtmp = data_ell[i][0];
			data_ell[i][0] = data_ell[i][iscan];
			data_ell[i][iscan] = dtmp;
		}
		else {
			std::cout << "error patch_Ell: diagonal otsutstvuet.\n";
			std::cout << j << "  isize=" << isize << "\n";
			int iscan1 = 0;
			while ((iscan1 < MAX_STRING_LENGTH_ELL) && (coll_ell[j][iscan1] != -1)) {
				std::cout << coll_ell[j][iscan1] << " "<<data_ell[j][iscan1] << std::endl;
				++iscan1;
			}
			system("PAUSE");
		}
		/*iscan = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1)) {
			std::cout << "a[" << i << "," << coll_ell[i][iscan] << "]=" << data_ell[i][iscan] << "  ";
			iscan++;
		}
		std::cout << std::endl;
		getchar();*/

	}
	//patch_Ell_zero(isize);
	//getchar();
}



// При обнаружении сбоя включить stable версия данной функции.
// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
template <typename MY_IND_TYPE>
void ell_to_CRS(Real*& val, MY_IND_TYPE*& col_ind, MY_IND_TYPE*& row_ptr, int nodes) {
	bool flag = true;
	int k; // счётчик
	//for (k = 0; k < nodes; k++) if (M.root[k] == nullptr) {
		//flag = false; break;
	//}

	//if (flag) {
	val = new Real[nnz_ell];
	col_ind = new MY_IND_TYPE[nnz_ell];
	row_ptr = new MY_IND_TYPE[nodes + 1];

	/*bool* bcheck = new bool[M.n];
	for (int i_1 = 0; i_1 < M.n; i_1++) {
		bcheck[i_1] = false;
	}
	NONZEROELEM* p_1 = nullptr;
	for (k = 0; k < nodes; k++) {
		p_1 = M.root[k];
		while (p_1 != nullptr) {
			if (bcheck[p_1->key]) {
				printf("ERROR MATRIX CHECK duplicate ja index string=%lld col_ind=%lld\n", k, p_1->key);
				system("pause");
			}
			bcheck[p_1->key] = true;

			p_1 = p_1->next;
		}
		p_1 = M.root[k];
		// Сброс.
		while (p_1 != nullptr) {
			bcheck[p_1->key] = false;
			p_1 = p_1->next;
		}
	}
	delete[] bcheck;*/

	// инициализация
	for (k = 0; k < (nnz_ell); k++) {
		val[k] = 0.0;
		col_ind[k] = 0;
	}
	for (k = 0; k <= nodes; k++) {
		row_ptr[k] = nnz_ell; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	}

	// Быстрая Сортировка Хоара.
	// упорядочивание по строкам
	//QuickSort(...); не требуется,
	// т.к. сама структура хранения 
	// подразумевает упорядочивание по строкам.

	/*
	// заполнение разреженной матрицы
	for (k=0; k<M.n; k++) {
		val[k]=M.a[k].aij;
		col_ind[k]=M.a[k].j;
		row_ptr[M.a[k].i]=myi_min(k,row_ptr[M.a[k].i]);
	}
	*/


	int ik = 0; // счётчик ненулевых элементов СЛАУ
	//NONZEROELEM* p;

	int nnz = 0;
	for (k = 0; k < nodes; k++) {
		int iscan = 0;
		//p = M.root[k];
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[k][iscan] != -1)) {
			//if (ik < M.n) {
				// Защита от записи нулевого коэффициента.
			//if (fabs(data_ell[k][iscan]) > 1.0e-300) При граничных условиях Дирихле для симметричности портрета нужно некоторое небольшое количество нулей.
			{
				nnz++;// 20.03.2021
				val[ik] = data_ell[k][iscan];
				/*if (p->key == k) {

					printf("diag=%e string %lld \n", val[ik], k);
					getchar();
				}*/
				// p->key начинается с нуля
				col_ind[ik] = coll_ell[k][iscan];
				/*if (p->key < 0) {
					printf("%lld\n", row_ptr[k]);
					system("pause");
				}*/
				row_ptr[k] = myi_min(ik, row_ptr[k]);
				ik++;
			}
			/*else {
				if (col_ind[ik] == k) {

					printf("diag=%e string %lld \n", val[ik], k);
					getchar();
				}
			}*/
			//}
			/*else {
				printf("error: ik>=M.n in simplesparsetoCRS\n");
				printf("see module my_linalg.c\n");
				system("pause");
				exit(1);
			}*/
			//p = p->next;
			iscan++;
		}
	}

	// Исправлена ошибка. Истинное количество ненулевых элементов меньше из за
	// того что берутся лишь ненулевые элементы выше некоторого порога.
	printf("reserve: M.n =%lld,  natural count: nnz=%lld\n", nnz_ell, nnz);
	//getchar(); // 20.03.2021
	row_ptr[nodes] = nnz;

	// в каждой строке элементы отсортированы по номерам столбцов:
	for (k = 0; k < nodes; k++) {


		if (0) {
			// BubbleSort
			BubbleSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
		}
		else {
			QuickSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
		}

	}

	//}
} // ell_to_CRS

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только нижний треугольник.
void  ell_to_CSIR(Real*& adiag, Real*& altr, int*& jptr, int*& iptr, int nodes) {
	bool flag = true;
	int k; // счётчик
	for (k = 0; k < nodes; k++) if (coll_ell[k][0] == -1) {

		std::cout << "error ell_to_CSIR:  k=" << k << "  nodes=" << nodes << std::endl;
		

		std::cout << coll_ell[k][0] << " " << coll_ell[k][1] << " " << coll_ell[k][2] << " " << coll_ell[k][3] << " " << coll_ell[k][4] << " " << coll_ell[k][5] << " \n";
		system("PAUSE");

		flag = false;
		//break;
	}

	if (flag) {
		// поддиагональные элементы в altr хранятся построчно
		int nz = (int)(nnz_ell - nodes) / 2; // число ненулевых элементов
		adiag = new Real[nodes]; // диагональные элементы
		altr = new Real[nz]; // поддиагональные элементы
		jptr = new int[nz]; // номера столцов для нижнего треугольника
		iptr = new int[nodes + 1]; // указатели на следующую строку


		// инициализация
		for (k = 0; k < nodes; k++) adiag[k] = 0.0;
		for (k = 0; k < (nz); k++) {
			altr[k] = 0.0;
			jptr[k] = 0;
		}
		for (k = 0; k <= nodes; k++) {
			iptr[k] = nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; k++) {
			val[k]=M.a[k].aij;
			col_ind[k]=M.a[k].j;
			row_ptr[M.a[k].i]= myi_min(k,row_ptr[M.a[k].i]);
		}
		*/
		/*
		int ik=0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
		for (k=0; k<nodes; k++) {
			p=M.root[k];
			while (p!=nullptr) {
				val[ik]=p->aij;
				col_ind[ik]=p->key;
				row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
				p=p->next;
			}
		}
		*/

		int ik = 0, imin = 1, k1; // счётчик ненулевых поддиагональных элементов СЛАУ
		bool bvisit;
		//NONZEROELEM* p;
		for (k = 0; k < nodes; k++) {
			bvisit = false;
			//p = M.root[k];
			int iscan = 0;
			while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[k][iscan] != -1)) {
				if (coll_ell[k][iscan] == k) {
					adiag[k] = data_ell[k][iscan];
					bvisit = false;
				}
				else if (coll_ell[k][iscan] < k) {
					if (ik < (nz)) {
						altr[ik] = data_ell[k][iscan]; // ненулевое значение
						jptr[ik] = coll_ell[k][iscan]; // номер столбца
					}
					else {
						printf("non simmetric matrix ICCG. ell_to_CSIR\n");
						std::cout << "ik=" << ik << " nz=" << nz << " coll_ell[" << k << "][iscan]=" << coll_ell[k][iscan] << "n==" << nodes << std::endl;
						//getchar();
						system("pause");
					}
					bvisit = true;
				}
				imin = min(ik, iptr[k]);
#if doubleintprecision == 1
				//printf("imin=%lld\n",imin);
#else
				//printf("imin=%d\n",imin);
#endif

				iptr[k] = imin;
				if (imin == 0) for (k1 = 0; k1 < k; k1++) iptr[k1] = 0;
				if (bvisit) {
					ik++;
					bvisit = false;
				}
				//p = p->next;
				++iscan;
			}
			if (iscan >= MAX_STRING_LENGTH_ELL) {
				std::cout << "ERROR!!! UVELICHTE MAX_STRING_LENGTH_ELL\n";
				system("PAUSE");
				exit(1);
			}
		}


		for (k = 0; k < nodes; k++) QuickSortCSIR(jptr, altr, iptr[k], iptr[k + 1] - 1);

	}
} // ell_to_CSIR


// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR_ITL. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только верхний треугольник.
// Память выделяется внутри метода.
void ell_to_CSIR_ITLSPD(Real*& val, int*& indx, int*& pntr, int nodes) {
	bool flag = true;
	int k; // счётчик
	for (k = 0; k < nodes; k++) if (coll_ell[k][0] == -1) {
		flag = false;

		std::cout << "error ell_to_CSIR_ITLSPD:  k=" << k << "  nodes=" << nodes << std::endl;
		system("PAUSE");
		//break;
	}

	if (flag) {

		//printM_and_CSIR(M, nodes); // debug

		// поддиагональные элементы в altr хранятся построчно
		int nz = (int)((nnz_ell - nodes) / 2 + nodes); // число ненулевых элементов
		val = new Real[nz]; // диагональные элементы и наддиагональные элементы
		indx = new int[nz]; // номера столцов для нижнего треугольника
		pntr = new int[nodes + 1]; // указатели на следующую строку


		// инициализация
		for (k = 0; k < (nz); k++) {
			val[k] = 0.0;
			indx[k] = 0;
		}
		for (k = 0; k <= nodes; k++) {
			pntr[k] = nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}



		int ik = 0; // счётчик ненулевых поддиагональных элементов СЛАУ
		//NONZEROELEM* p;
		for (k = 0; k < nodes; k++) {

			//p = M.root[k];
			int iscan = 0;
			while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[k][iscan] != -1)) {

				// k - номер диагонального элемента
				if ((coll_ell[k][iscan] >= k) && (coll_ell[k][iscan] < nodes)) {
					if (ik < (nz)) {
						val[ik] = data_ell[k][iscan]; // ненулевое значение
						indx[ik] = coll_ell[k][iscan]; // номер столбца	
					}
					else {
						std::cout << "\nik=" << ik << " nz=" << nz << " nodes=" << nodes << " k=" << k << std::endl;
						printf(" Error non simmetric matrix ICCG. ell_to_CSIR_ITLSPD\n");
						//getchar();
						system("pause");
					}
					pntr[k] = min(ik, pntr[k]);

					ik++;
				}

				//p = p->next;
				iscan++;
			}

			if (iscan >= MAX_STRING_LENGTH_ELL) {
				std::cout << "ERROR!!! UVELICHTE MAX_STRING_LENGTH_ELL\n";
				system("PAUSE");
				exit(1);
			}

		}

		for (k = 0; k < nodes; k++) QuickSortCSIR(indx, val, pntr[k], pntr[k + 1] - 1);


		/*
		FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
	#if doubleintprecision == 1
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");
		// debug
		for (k=0; k<=nodes; k++) {
			fprintf(fp,"%lld ",pntr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<pntr[nodes]; k++) {
			fprintf(fp, "%e %lld\n",val[k],indx[k]);
		}
		fprintf(fp, "nz==%lld\n", nz);
	#else
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");
		// debug
		for (k=0; k<=nodes; k++) {
			fprintf(fp,"%d ",pntr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<pntr[nodes]; k++) {
			fprintf(fp, "%e %d\n",val[k],indx[k]);
		}
		fprintf(fp, "nz==%d\n", nz);
	#endif


		fclose(fp);
	}
	printf("ready");
	getchar();
	*/
	}
} // ell_to_CSIR_ITLSPD



/* Следующие несколько функций применяются для GSEP
*  в целях упорядочивания по возрастанию набора 
*  собственных значений.
*/

// Пузырьковая сортировка.
void BubbleSortGSEP1(Real *a, int *mask, int n) {
   int i=0, j=0, k=0;
   Real x;

   for (i=1; i<n; i++) {
	   for (j=n-1; j>=i; j--) {
		   if (a[j-1] > a[j]) {
			   // swap
			   x=a[j-1];
			   a[j-1]=a[j];
			   a[j]=x;
               k=mask[j-1];
			   mask[j-1]=mask[j];
			   mask[j]=k;
		   }
	   }
   }
} // BubbleSortGSEP1


/* Первая обобщённая симметричная проблема собственных значений
*   GSEP1:  A*x-lambda_scal*B*x=0;
*   Разложение Холесского: B=L*transpose(L);
*   L - нижняя треугольная, transpose(L) - верхняя треугольная.
*/
void GSEP1(Real **A, Real **B, Real **U, Real *lambda, int *mask, int nodes, Real epsilon) {

	// Разложение Холесского: замена B верхним и нижним 
	// треугольными множителями.
	B[0][0]=sqrt(B[0][0]);
	B[1][0]/=B[0][0];
	B[0][1]=B[1][0];
	B[1][1]=sqrt(B[1][1]-B[1][0]*B[1][0]);

	int irow,irow1;
	int icol, icol1;
	Real sum;
	int k;
	for (irow=2; irow<nodes; irow++) {
		irow1=irow-1;
		B[irow][0]/=B[0][0];
        B[0][irow]=B[irow][0];
		for (icol=1; icol<=irow1; icol++) {
			icol1=icol-1;
            sum=0.0;
            for (k=0; k<=icol1; k++) sum+=B[irow][k]*B[icol][k];
			B[irow][icol]=(B[irow][icol]-sum)/B[icol][icol];
			B[icol][irow]=B[irow][icol];
		}
		sum=0.0;
		for (k=0; k<=irow1; k++) sum+=B[irow][k]*B[irow][k];
		B[irow][irow]=sqrt(B[irow][irow]-sum);
	}

	printf("L*LT 10...\n");
   // TODO: дальше и до конца функции идёт чрезвычайно
	// неэффективный кусок кода:

    int i=0, j=0;

	for (i=0; i<nodes; i++) mask[i]=i;

	// нижняя треугольная матрица
    Real **L=new Real*[nodes];
    for (i = 0; i < nodes; i++) L[i]=new Real[nodes];

	/* Этот закоментированный кусок кода относится к медленной реализации:
	// Если использовать медленную реализацию то это надо раскоментировать.
	// инициализация всей матрицы выполняется один раз:
    for (i = 0; i < nodes; i++)
         for (j = 0; j < nodes; j++) {
            if (j > i)
                L[i][j] = 0.0;
            else L[i][j] = B[i][j];
    }
	*/

    /*
    // верхняя треугольная матрица
    Real **LT=new Real*[nodes];
    for (i = 0; i < nodes; i++) LT[i]=new Real[nodes];

	// инициализация всей матрицы выполняется один раз:
    for (i = 0; i < nodes; i++)
         for (j = 0; j < nodes; j++) {
            if (j < i)
                LT[i][j] = 0.0;
            else LT[i][j] = B[i][j];
    }
	*/

    // вспомогательная матрица для умножения
    Real **b=new Real*[nodes];
    for (i = 0; i < nodes; i++) b[i]=new Real[nodes];

	// Ac копия матрицы А
    Real **Ac=new Real*[nodes];
    for (i = 0; i < nodes; i++) Ac[i]=new Real[nodes];
	matr_copy(Ac,A,nodes); // сохранение А TODO временно потом удалить

	// Медленная реализация
    //inverse_matrix_simple(L,nodes); // нахождение L^(-1)
	//multi_m(L,A,b,nodes); // b=(L^(-1))*A;
	//matr_copy(A,b,nodes); // A=(L^(-1))*A;

	// Более быстрая реализация
    // A=(L^(-1))*A;
	for (i=0; i < nodes; i++) {
		A[0][i]/=B[0][0];

	    for (irow=1; irow<nodes; irow++) {
		    irow1=irow-1;
		    sum=0.0;
		    for (icol=0; icol<=irow1; icol++) sum+=B[irow][icol]*A[icol][i];
            A[irow][i]=(A[irow][i]-sum)/B[irow][irow];
	    }
	}

    printf("(L^(-1))*A 20...\n");
	//matr_copy(L,LT,nodes); // L=transpose(L); т.к. матрица L больше не нужна:
	// дальше везде под именем L используется transpose(L).
    
    // L=LT: L=transpose(L);
	// Теперь L верхняя треугольная матрица
    for (i = 0; i < nodes; i++)
         for (j = 0; j < nodes; j++) {
            if (j < i)
                L[i][j] = 0.0;
            else L[i][j] = B[i][j];
    }

    // Эта матрица уже приведена к верхнетреугольному виду,
    // поэтому второй раз её к такому виду приводить не понадобиться поэтому true.
    inverse_matrix_simple(L,nodes,true); // нахождение (transpose(L))^(-1)
	 
	multi_m(A,L,b,nodes); // b=(L^(-1))*A*(transpose(L))^(-1).
    matr_copy(A,b,nodes); // A=(L^(-1))*A*(transpose(L))^(-1).

	printf("C 30...\n");

	jacobi_matrix_simple(A,U,lambda,nodes,epsilon); // нахождение СВ и СЗ с заданной точностью.

    printf("C 90...\n");

	BubbleSortGSEP1(lambda,mask,nodes); // упорядочивание собственных значений.
	multi_m(L,U,b,nodes); // b=((transpose(L))^(-1))*U
    matr_copy(U,b,nodes); // собственные вектора.

	/* проверка найденных собственных значений.
    multi_m(Ac,U,b,nodes); // b=A*U
    matr_copy(L,U,nodes); // L=U
	tr_m(L,nodes);
	multi_m(L,b,Ac,nodes); // Ac=transpose(U)*A*U

	Real *test=new Real[nodes];
    for (int i=0; i<nodes; i++) test[i]=Ac[i][i];
    BubbleSortGSEP1(test,mask,nodes); 
    for (int i=0; i<8; i++) printf("%.2f ",test[i]/3.141/3.141); // собственные значения
	printf("\n");
	*/

	delete[] L;  delete[] b; //delete LT;
} // GSEP1

/* метод Гаусса для ленточной матрицы A размером
*              nodes x 2*icolx+1, где
*   2*icolx+1 - ширина ленты. Под тем что матрица
*  A ленточная понимается то что ненулевые элементы
*  матрицы содержатся только внутри ленты.
*  b - вектор правой части СЛАУ, x - вектор решение.
*  Нумерация элементов начинается с нуля.
*  Для положительно определённых возможно несимметричных
*  матриц А, которые задаются своей лентой.
*  Гаусс Карл Фридрих 1777-1855.
*  В результате работы матрица А портится.
*/
void eqsolve_lenta_gauss(Real **A, int nodes, int icolx, Real *b, Real *x) {

	const Real eps=1e-300; // для сравнения с нулём
	Real dCik, dSum=0.0;
	int max;

	int *move=new int[nodes]; // массив сдвигов.
	int i=0, j=0, k=0; // счётчики цикла for
	for (i=0; i<nodes; i++) move[i]=icolx-i; // инициализация массива сдвигов

	for (i=0; i<nodes; i++) x[i]=0.0; // инициализация

	// прямой ход метода Гаусса
	// приведение к верхнему треугольному виду:

	// по всем столбцам слева направо
	for (k=0; k<nodes; k++) {
        max=min(k+icolx,nodes-1);
		// цикл по всем строкам ниже строки с номером k
		for (i=k+1; i<=max; i++) {
			// применяется только в том случае
			// если элемент ненулевой
			// это должно несколько ускорить счёт.
			if (fabs(A[i][k+move[i]]) > eps) {
               
                if(fabs(A[k][k+move[k]])<eps){
			          // решение не может быть получено, т.к.
			          // на диагонали находится ноль.
	                  printf("\nSolution is not exist! divizion by zero...\n");
	                  system("pause");
		              exit(0);
	            }

                // обработка фиксированной строки с номером i
				dCik=A[i][k+move[i]]/A[k][k+move[k]];
				// преаобразование матрицы к верхнетреугольному виду:
				for (j=k; j<=max; j++) A[i][j+move[i]] -= dCik*A[k][j+move[k]];
				b[i]-= dCik*b[k]; // преобразование правой части
			}
		}
	}

    // Теперь когда матрица приведена к верхнетреугольному виду
	// можно совершить обратный ход метода Гаусса:
	for (k=nodes-1; k>=0; k--) {
        dSum=0.0; // обнуление сумматора
		max=min(k+icolx,nodes-1);
		for (i=k+1; i<=max; i++) dSum+= A[k][i+move[k]]*x[i];
		x[k]=(b[k]-dSum)/A[k][k+move[k]];
	}

}  // eqsolve_lenta_gauss

// Метод (Якоби) Гаусса-Зейделя
// для решения СЛАУ с матрицей А nxn
// возможно несимметричной, но с диагональным 
// преобладанием. Матрица А предполагается
// полностью заполненой (неразреженной).
// b - правая часть, x - уточняемое решение, 
// eps - точность определения решения.
// omega - импирически подобранный параметр релаксации.
void Seidel(Real **A, Real *b, Real *x, int n, Real eps, Real omega) {
	int i,j;
	Real s1, s2, s, v, m=0.0;
	bool bdiag=true;

	// Исследуем сходимость
	for (i=0; i<n; i++) {
		s=0.0;
		for (j=0; j<n; j++) {
			if (j!=i) s+=fabs(A[i][j]);
		}
		if (s>=fabs(A[i][i])) {
			bdiag=false;
		}
	}
	if (!bdiag) {
		printf("net diagonalnogo preobladaniq...");
		system("pause");
	}

	do {
		m=0.0;
		for (i=0; i<n; i++) {
			// Вычисляем суммы
			s1=s2=0.0; 
			for (j=0; j<=i-1; j++) s1+=A[i][j]*x[j];
			for (j=i+1; j<n; j++) s2+=A[i][j]*x[j];
			// Вычисляем новое приближение и погрешность
			v=x[i];
			x[i]=omega*(b[i]-s1-s2)/A[i][i]+(1-omega)*x[i];

			if (fabs(v-x[i])>m) m=fabs(v-x[i]);
		}

	} while (m > eps);

} // Seidel

// возвращает максимальное из двух 
// вещественных чисел.
/*
Real fmax(Real fA, Real fB) {
	Real r=fB;
	if (fA > fB) r=fA;
	return r;
} // fmax 
*/

// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
void SOR(equation* &sl, Real* &x, int n) {
	Real rURF = 1.0;// 1.855; // параметр верхней релаксации
	// пороговое значение невязки
	Real eps = 1e-3;
	Real ptilda;
	Real sE,sW,sN, sS;
	int i=0,j=0, kend=3000; // счётчик цикла for
	Real dmax=1.0;
	while ((dmax>eps) && (j<kend)) {
		dmax=0.0;
	    for (i=0; i<n; i++) {
            if (sl[i].iE>-1) sE=sl[i].ae*x[sl[i].iE]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*x[sl[i].iW]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*x[sl[i].iN]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*x[sl[i].iS]; else sS=0.0;
		    ptilda=(sE+sW+sN+sS+sl[i].b)/sl[i].ap;
		    //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			dmax+=fabs(sl[i].ap*(ptilda-x[sl[i].iP]));
		    x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
	    }
		dmax/=n;
		if (DEBUG) printf("%e \n",dmax);
		j++;
	}

} // SOR

// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
void SOR(equation*& sl, Real*& x, Real*& rthdsd, int n, int iVar) {
	Real rURF = 1.0;// 1.855; // параметр верхней релаксации
	// пороговое значение невязки
	Real eps = 1e-3;
	
	if (iVar == Temp) eps = 1.0e-3/ Cp_active;
	if (iVar == PAm) {
		eps = 1.0e-5;
		rURF = 1.55;
	}
	
	int  j = 0, kend = 3000; // счётчик цикла for
	Real dmax = 1.0;
	while ((dmax > eps) && (j < kend)) {
		dmax = 0.0;

		if (inumcore == 2) {

			Real dmax1 = 0.0;
			Real dmax2 = 0.0;

#pragma omp parallel sections
			{
#pragma omp section
				{
					for (int i = s_par[1].s; i < s_par[1].e; ++i) {
						Real sE, sW, sN, sS;
						Real ptilda;
						if (sl[i].iE > -1) sE = sl[i].ae * x[sl[i].iE]; else sE = 0.0;
						if (sl[i].iW > -1) sW = sl[i].aw * x[sl[i].iW]; else sW = 0.0;
						if (sl[i].iN > -1) sN = sl[i].an * x[sl[i].iN]; else sN = 0.0;
						if (sl[i].iS > -1) sS = sl[i].as * x[sl[i].iS]; else sS = 0.0;
						ptilda = (sE + sW + sN + sS + rthdsd[i]) / sl[i].ap;
						//dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
						dmax1 += fabs(sl[i].ap * (ptilda - x[sl[i].iP]));
						x[sl[i].iP] = x[sl[i].iP] + rURF * (ptilda - x[sl[i].iP]);
					}
				}
#pragma omp section
				{
					for (int i = s_par[2].s; i < s_par[2].e; ++i) {
						Real sE, sW, sN, sS;
						Real ptilda;
						if (sl[i].iE > -1) sE = sl[i].ae * x[sl[i].iE]; else sE = 0.0;
						if (sl[i].iW > -1) sW = sl[i].aw * x[sl[i].iW]; else sW = 0.0;
						if (sl[i].iN > -1) sN = sl[i].an * x[sl[i].iN]; else sN = 0.0;
						if (sl[i].iS > -1) sS = sl[i].as * x[sl[i].iS]; else sS = 0.0;
						ptilda = (sE + sW + sN + sS + rthdsd[i]) / sl[i].ap;
						//dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
						dmax2 += fabs(sl[i].ap * (ptilda - x[sl[i].iP]));
						x[sl[i].iP] = x[sl[i].iP] + rURF * (ptilda - x[sl[i].iP]);
					}
				}
			}

			dmax += dmax1 + dmax2;
			for (int i = s_par[3].s; i < s_par[3].e; ++i) {
				Real sE, sW, sN, sS;
				Real ptilda;
				if (sl[i].iE > -1) sE = sl[i].ae * x[sl[i].iE]; else sE = 0.0;
				if (sl[i].iW > -1) sW = sl[i].aw * x[sl[i].iW]; else sW = 0.0;
				if (sl[i].iN > -1) sN = sl[i].an * x[sl[i].iN]; else sN = 0.0;
				if (sl[i].iS > -1) sS = sl[i].as * x[sl[i].iS]; else sS = 0.0;
				ptilda = (sE + sW + sN + sS + rthdsd[i]) / sl[i].ap;
				//dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
				dmax += fabs(sl[i].ap * (ptilda - x[sl[i].iP]));
				x[sl[i].iP] = x[sl[i].iP] + rURF * (ptilda - x[sl[i].iP]);
			}
			dmax /= n;
		}
		else {

			for (int i = 0; i < n; ++i) {
				Real sE, sW, sN, sS;
				Real ptilda;
				if (sl[i].iE > -1) sE = sl[i].ae * x[sl[i].iE]; else sE = 0.0;
				if (sl[i].iW > -1) sW = sl[i].aw * x[sl[i].iW]; else sW = 0.0;
				if (sl[i].iN > -1) sN = sl[i].an * x[sl[i].iN]; else sN = 0.0;
				if (sl[i].iS > -1) sS = sl[i].as * x[sl[i].iS]; else sS = 0.0;
				ptilda = (sE + sW + sN + sS + rthdsd[i]) / sl[i].ap;
				//dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
				dmax += fabs(sl[i].ap * (ptilda - x[sl[i].iP]));
				x[sl[i].iP] = x[sl[i].iP] + rURF * (ptilda - x[sl[i].iP]);
			}
			dmax /= n;
		}
		if (DEBUG) {
			printf("%e \n", dmax);
		}
		j++;
	}

	//if (iVar == Temp) {
		//std::cout << " j= " << j << std::endl;
	//}

} // SOR

/* Метод Сопряжённых градиентов
*  без учёта разреженности матрицы СЛАУ.
*/

// умножение матрицы на вектор
Real* MatrixByVector(Real** H,Real* V,int n){
	Real* tmp=new Real[n];
	Real sum=0.0;
	for (int i=0;i<n;++i){
		for (int j=0;j<n;++j)
			sum+=V[j]*H[i][j];
		tmp[i]=sum;
		sum=0.0;}
	return tmp;
} // MatrixByVector

// Евклидова норма вектора
Real NormaV(double *V, int n){
	Real norma;
	Real s=0;
	for(int i=0;i<n;i++)
		s+=V[i]*V[i];
	norma=sqrt(s);
	return norma;
} // NormaV

// Скалярное произведение двух векторов
Real Scal(Real *v1, Real *v2, int n){
	Real s=0.0;
	int i; // счётчик цикла for

   // omp_set_num_threads(inumcore);

  #pragma omp parallel for shared(v1, v2) private(i) reduction (+: s) schedule (guided)
	for ( i=0; i<n; i++)
	{ 
		s+=v1[i]*v2[i];
	}
	return s;
} // Scal 

//----------метод сопряженных градиентов---------------
/* Входные параметры:
*  A - неразреженная матрица СЛАУ,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или NULL.
*  n - размерность СЛАУ Anxn.
*  Матрица A полагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*/
Real *SoprGrad(Real **A, Real *dV, Real *x, int n){
	printf("Reshenie metodom sopryjennyh gradientov:\n");
	int k=0;
	int i; // счётчик
	Real *ap=new Real[n],
		 *z=new Real[n], *p=new Real[n];
	Real a, b, nz;

	// шаг 1.1
	//X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	ap=MatrixByVector(A,x,n);
	for (i=0; i<n; i++) z[i]=dV[i]-ap[i];

	if (Scal(z,z,n)!=0){
		// шаг 1.3
	   for (i=0; i<n; i++)	p[i]=z[i];
	   nz=1000.;
	   while ((nz>e) && (k<1000)) {
		   // шаг 2.1
	 	  ap=MatrixByVector(A,p,n);
		  // шаг 2.2
		  //a=Scal(z,p,n)/Scal(z,ap,n);
		  a=Scal(z,p,n)/Scal(ap,p,n); // шаговый множитель
		  // шаг 2.3 и 2.4
		  for (i=0; i<n; i++) {
		      x[i]+=a*p[i]; // очередное приближение
			  z[i]-=a*ap[i]; // невязка k+1-го приближения
		  }
		  // шаг 2.5
		  nz=NormaV(z,n);
		  if (k%10==0) printf("iter residual\n");
		  printf(" %d %e\n", k, nz);
		  // шаг 3.1
		  b=Scal(z,ap,n)/Scal(p,ap,n);
		  // шаг 3.2
		  for (i=0; i<n; i++) {
		     p[i]=z[i]-b*p[i]; // новое направление минимизации
		  }
          // шаг 3.3 
		  k++;
	   } // while

	   // Освобождение памяти
        delete[] ap;
		delete[] z; delete[] p;

	   return x;
	}
	else {
		// Освобождение памяти
		delete[] ap;
		delete[] z; delete[] p;

		return x;
	}
} // SoprGrad

/* Описание стандарта хранения CRS:
*  1. val - ненулевые значения элементов матрицы отсортированные
*  по номерам строк (нумерация начинается с нуля).
*  2. col_ind - соответствующие элементам из val номера столбцов.
*  3. row_ptr - используется для определения начала следующей строки.
*  Пример:
*
*  9.0   0.0   0.0   3.0   1.0   0.0   1.0    
*  0.0   11.0   2.0   1.0   0.0   0.0   2.0    
*  0.0   2.0   10.0   2.0   0.0   0.0   0.0    
*  3.0   1.0   2.0   9.0   1.0   0.0   0.0    
*  1.0   0.0   0.0   1.0   12.0   0.0   1.0    
*  0.0   0.0   0.0   0.0   0.0   8.0   0.0    
*  1.0   2.0   0.0   0.0   1.0   0.0   8.0    
*
*------------- Разреженная матрица ------------ 
* Формат хранения: CRS  
* val:      9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
* col_ind:  0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
* row_ptr:  0 4 8 10 12 14 15 16 
*------------------------------------------------------
*/


// умножение матрицы на вектор
// используется формат хранения CRS
// Разреженная матрица A (val, col_ind, row_ptr) квадратная размером nxn.
// Число уравнений равно числу неизвестных и равно n.
void MatrixCRSByVector(Real* val, int* col_ind, int* row_ptr, Real* V, Real* &tmp, int n)
{

	omp_set_num_threads(inumcore);

    // вектор tmp индексируется начиная с нуля так же как и вектор V
#pragma omp parallel for
	for (int i=0; i<n; ++i) tmp[i]=0.0;
	if (tmp == NULL)
	{
		printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
		system("pause");
		exit(0);  // завершение программы
	}
    
	

    #pragma omp parallel for shared(row_ptr, val, col_ind, V, tmp) schedule (guided)
	for (int i=0; i<n; ++i) {
	    Real sum = 0.0;
		int rowend=row_ptr[i+1];
		int rowbeg=row_ptr[i];
	    for (int j = rowbeg; j<rowend; ++j)
		{
		    	sum += val[j]*V[col_ind[j]];
		}
		tmp[i] = sum;
	}
	
	//return tmp;
} // MatrixCRSByVector

// умножение транспонированной матрицы на вектор
// (используется, например, в методе BiCG - бисопряжённых градиентов)
// для исходной (не транспонированной матрицы) используется формат хранения CRS
// Разреженная матрица A (val, col_ind, row_ptr) квадратная размером nxn.
// Число уравнений равно числу неизвестных и равно n.
Real* MatrixTransposeCRSByVector(Real* val, int* col_ind, int* row_ptr, Real* V, int n)
{
	
	Real* tmp=new double[n]; // вектор индексируется начиная с нуля так же как и вектор V
	if (tmp == NULL)
	{
		printf("malloc: out of memory for vector tmp in MatrixTransposeCRSByVector\n"); // нехватка памяти
		system("pause");
		exit(0);
		return NULL; // завершение программы
	}
	
	
    int i,j; // Счётчики цикла
	int rowend, rowbeg;
    
	for (i=0; i<n; i++) tmp[i]=0.0;

	for (j=0; j<n; j++) {
		rowend=row_ptr[j+1];
		rowbeg=row_ptr[j];
	    for (i = rowbeg; i<rowend; i++)
		{
		    	tmp[col_ind[i]] += val[i]*V[j];
		}
	}
	
	return tmp;
} // MatrixTransposeCRSByVector


/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  val, col_ind, row_ptr - разреженная матрица СЛАУ в формате CRS,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или NULL.
*  n - размерность СЛАУ Anxn.
*  Разреженная матрица A (val, col_ind, row_ptr) квадратная размером nxn.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица A полагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*/
Real *SoprGradCRS(Real *val, int* col_ind, int* row_ptr, Real *dV, Real *x, int n){
	printf("Conjugate Gradients Method...:\n");
	int k=0;
	int i; // счётчик
	Real *ap=new Real[n],
		 *z=new Real[n], *p=new Real[n];
	Real a, b, nz;

    omp_set_num_threads(inumcore);

	// шаг 1.1
	//X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	MatrixCRSByVector(val,col_ind,row_ptr,x,ap,n);
	
    #pragma omp parallel for shared(z,dV,ap) private(i) schedule (guided)
	for (i=0; i<n; i++) z[i]=dV[i]-ap[i];

	if (Scal(z,z,n)!=0){
		// шаг 1.3
       #pragma omp parallel for shared(p,z) private(i) schedule (guided)
	   for (i=0; i<n; i++)	p[i]=z[i];

	   nz=1000.;
	   while ((nz>e) && (k<2*n)) {
		   // шаг 2.1
		  // чтобы избежать утечки памяти
	 	  MatrixCRSByVector(val,col_ind,row_ptr,p,ap,n);
		  // шаг 2.2
		  //a=Scal(z,p,n)/Scal(z,ap,n);
		  a=Scal(z,p,n)/Scal(ap,p,n); // шаговый множитель
		  // шаг 2.3 и 2.4
		  #pragma omp parallel for shared(x,z,p,ap,a) private(i) schedule (guided)
		  for (i=0; i<n; i++) {
		      x[i]+=a*p[i]; // очередное приближение
			  z[i]-=a*ap[i]; // невязка k+1-го приближения
		  }
		  // шаг 2.5
		  nz=NormaV(z,n);
		  if (k%10==0) printf("iter residual\n");
		  printf(" %d %e\n", k, nz);
		  // шаг 3.1
		  b=Scal(z,ap,n)/Scal(p,ap,n);
		  // шаг 3.2
		  #pragma omp parallel for shared(p,z,b) private(i) schedule (guided)
		  for (i=0; i<n; i++) {
		     p[i]=z[i]-b*p[i]; // новое направление минимизации
		  }
          // шаг 3.3 
		  k++;
	   } // while

	   // Освобождение памяти
        delete[] ap;
		delete[] z; delete[] p;

	   return x;
	}
	else {
		// Освобождение памяти
		delete[] ap;
		delete[] z; delete[] p;

		return x;
	}
} // SoprGradCRS

// Метод бисопряжённых градиентов
// для возможно несимметричной матрицы А (val, col_ind, row_ptr).
// Запрограммировано по книжке Баландин, Шурина : "Методы
// решения СЛАУ большой размерности".
// dV - правая часть СЛАУ,
// x - начальное приближение к решению или NULL.
// n - размерность А nxn.
// Количество итераций ограничено 2000.
// Точность выхода по невязке задаётся в глобальной константе:
//  dterminatedTResudual.
Real *BiSoprGradCRS(Real *val, int* col_ind, int* row_ptr, Real *dV, Real *x, int n, int maxit){
	if (DEBUG) printf("BiConjugate Gradients Method...:\n");

	Real *r=new Real[n], *r_tilda=new Real[n];
	Real *p=new Real[n], *p_tilda=new Real[n];
	Real nz; // невязка
	Real *ap=new Real[n];
	Real a,b,dold, dnew;

	int i; // счётчик цикла for
	int k=0; // номер итерации.

	// Начальное приближение:
    //X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;

	MatrixCRSByVector(val,col_ind,row_ptr,x,ap,n);
	for (i=0; i<n; i++) {
		r[i]=dV[i]-ap[i];
		r_tilda[i]=r[i];
		p[i]=r[i];
		p_tilda[i]=r_tilda[i];
	}

	nz=NormaV(r,n); // начальное значение невязки
	dold=Scal(r,r_tilda,n);

    while ((nz>e) && (k<maxit)) {
		MatrixCRSByVector(val,col_ind,row_ptr,p,ap,n);

		a=dold/Scal(ap,p_tilda,n);
		for (i=0; i<n; i++) {
           x[i]+=a*p[i];
		   r[i]-=a*ap[i];
		}
		delete[] ap;
		ap=MatrixTransposeCRSByVector(val,col_ind,row_ptr,p_tilda,n);
        for (i=0; i<n; i++) {
			r_tilda[i]-=a*ap[i];
		}
		dnew=Scal(r,r_tilda,n);
		b=dnew/dold;
		dold=dnew;
		// вычисление невязки.
        nz=NormaV(r,n);
		if (DEBUG) if (k%10==0) printf("iter residual\n");
		printf(" %d %e\n", k, nz);

		if (fabs(b) < 1e-270) {
			printf("\nBiCG divergence detected...\n");
            system("pause");
			exit(0); // выход из приложения.
			break; // выход из цикла while
		}

        for (i=0; i<n; i++) {
			p[i]=r[i]+b*p[i];
			p_tilda[i]=r_tilda[i]+b*p_tilda[i];
		}

		k++; // переход к следующей итерации.
	}

	// Освобождение памяти
	delete[] r; delete[] r_tilda; 
	delete[] p; delete[] p_tilda;
	delete[] ap;

	return x;

} // BiSoprGradCRS

// Прямой ход по разреженной нижнетреугольной матрице L.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// L - хранится в следующем виде:
// 1. ldiag - диагональные элементы L.
// 2. lltr - поддиагональные элементы в строчном формате,
// т.е. хранение построчное. 
// 3. jptr - соотвествующие номера столбцов для lltr, 
// 4. iptr - информация о начале следующей строки для lltr.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(L)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* inverseL(Real* f, Real* ldiag, Real* lltr, int* jptr, int* iptr, int n) {
	Real *z=new Real[n];

    if (z == nullptr)
	{
		printf("malloc: out of memory for vector z in inverse(L)*f \n"); // нехватка памяти
		system("pause");
		exit(0);
		return nullptr; // завершение программы
	}

	int i,j;
	for (i=0; i<n; i++) {
		for (j=iptr[i]; j<iptr[i+1]; j++) {
			f[i]-=z[jptr[j]]*lltr[j];
		}
		z[i]=f[i]/ldiag[i];
	}
	return z;
}//inverseL

// Прямой ход по разреженной нижнетреугольной матрице L.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// L - хранится в следующем виде:
// 1. val - диагональные и поддиагональные элементы L.
// в столбцовом порядке. 
// 3. indx - соотвествующие номера строк для val, 
// 4. pntr - информация о начале следующего столбца.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(L)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseL_ITL0(Real* f, Real* val, int* indx, int* pntr, Real* &z, int n) {
	

    if (z == NULL)
	{
		Real *z=new Real[n];
		if (z==NULL) {
			printf("malloc: out of memory for vector z in inverse(L)*f \n"); // нехватка памяти
		    system("pause");
		    exit(0); // завершение программы
		}
	}

	int i,j;
	for (i=0; i<n; i++) {
        z[i]=f[i]/val[pntr[i]];
		// оьработка i-го столбца
		for (j=pntr[i]+1; j<pntr[i+1]; j++) {
			f[indx[j]]-=z[i]*val[j];
		}
		
	}

}//inverseL_ITL


// Прямой ход по разреженной нижнетреугольной матрице L.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// L - хранится в следующем виде:
// 1. val - диагональные и поддиагональные элементы L.
// в столбцовом порядке. 
// 3. indx - соответствующие номера строк для val, 
// 4. pntr - информация о начале следующего столбца.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(L)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseL_ITL(Real* f, Real* val, int* indx, int* pntr, Real*& z, int n) {

	// Real **fbuf;
	// набор векторов fbuf нужен только в параллельной версии, в серийной версии можно просто передавать nullptr.
	// количество векторов в fbuf равно количеству потоков.

	if (z == nullptr)
	{
		// Попробуем выделить память. 23.03.2019
		z = new Real[n];
		if (z == nullptr) {
			printf("malloc: out of memory for vector z in inverse(L)*f \n"); // нехватка памяти
		   // getchar();
			system("pause");
			exit(0); // завершение программы
		}
	}


	

		//bool bserial=true;

		//if (bserial) {


		if (0) {

			for (int i = 0; i < n; i++) {
				z[i] = f[i] / val[pntr[i]];

				// обработка i-го столбца
				// эта часть не поддаётся распараллеливанию.
				// из за зависимостей по данным для f.
				for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
					f[indx[j]] -= z[i] * val[j];
				}

			}

		}
		else {

			// однопоточное исполнение.
			if (inumcore == 1) {

				for (int i = 0; i < n; i++) {
					z[i] = f[i] / val[pntr[i]];

					// обработка i-го столбца
					// эта часть не поддаётся распараллеливанию.
					// из за зависимостей по данным для f.
					for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
						f[indx[j]] -= z[i] * val[j];
					}

				}

			}
			if (inumcore == 2) {

#ifdef _OPENMP 
				omp_set_num_threads(inumcore); // установка числа потоков
#endif





#pragma omp parallel sections
				{
#pragma omp section
					{
						for (int i = s_par[1].s; i < s_par[1].e; i++) {

							//if (inumerate[i] == 1)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[2].s; i < s_par[2].e; i++) {

							//if (inumerate[i] == 2)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}


				for (int i = s_par[3].s; i < s_par[3].e; i++) {

					//if (inumerate[i] == 3)
					{

						z[i] = f[i] / val[pntr[i]];

						// обработка i-го столбца
						// эта часть не поддаётся распараллеливанию.
						// из за зависимостей по данным для f.
						for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
							f[indx[j]] -= z[i] * val[j];
						}
					}

				}

			}

			if (inumcore == 3)
			{
#pragma omp parallel sections
				{
#pragma omp section
					{
						for (int i = s_par[1].s; i < s_par[1].e; i++) {

							//if (inumerate[i] == 1)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[2].s; i < s_par[2].e; i++) {

							//if (inumerate[i] == 2)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[3].s; i < s_par[3].e; i++) {

							//if (inumerate[i] == 3)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}

#pragma omp parallel sections
				{
#pragma omp section
					{
						for (int i = s_par[4].s; i < s_par[4].e; i++) {

							//if (inumerate[i] == 4)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[5].s; i < s_par[5].e; i++) {

							//if (inumerate[i] == 5)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}

			}

			if ((inumcore == 4) || (inumcore == 5)) {

#ifdef _OPENMP 
				omp_set_num_threads(inumcore); // установка числа потоков
#endif





#pragma omp parallel sections
				{
#pragma omp section
					{
						for (int i = s_par[4].s; i < s_par[4].e; i++) {

							//if (inumerate[i] == 4)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[6].s; i < s_par[6].e; i++) {

							//if (inumerate[i] == 6)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[7].s; i < s_par[7].e; i++) {

							//if (inumerate[i] == 7)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[9].s; i < s_par[9].e; i++) {

							//if (inumerate[i] == 9)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}


#pragma omp parallel sections
				{
#pragma omp section
					{
						for (int i = s_par[5].s; i < s_par[5].e; i++) {

							//if (inumerate[i] == 5)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (int i = s_par[8].s; i < s_par[8].e; i++) {

							//if (inumerate[i] == 8)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}


				for (int i = s_par[3].s; i < s_par[3].e; i++) {

					//if (inumerate[i] == 3)
					{

						z[i] = f[i] / val[pntr[i]];

						// обработка i-го столбца
						// эта часть не поддаётся распараллеливанию.
						// из за зависимостей по данным для f.
						for (int j = pntr[i] + 1; j < pntr[i + 1]; j++) {
							f[indx[j]] -= z[i] * val[j];
						}
					}

				}

			}


		}
	
	/*

}
else {
	// параллельное исполнение.
	// параллельный код требует правильного разрешения зависимостей по данным.

	// Нам понадобиться
	// n=omp_get_max_threads();
	// дополнительных векторов.

	int nt=0;
#pragma omp parallel shared(nt)
		{
			// число нитей.
			nt=omp_get_max_threads();
		}



		for (int i=0; i<nt; i++) {
			for (int j=0; j<n; j++) {
				fbuf[i][j]=0.0; // инициализация.
			}
		}

#pragma omp for  shared(n, z, val, f, fbuf, pntr, indx, fbuf)
		for (int i=0; i<n; i++) {
		   // Проблема в том что здесь используется f[i], а оно может быть обновлённым, что здесь не учитывается !!!
			z[i]=f[i]/val[pntr[i]];
			// обработка i-го столбца
			// эта часть не поддаётся распараллеливанию.
			// из за зависимостей по данным для f.
			for (int j=pntr[i]+1; j<pntr[i+1]; j++) {
				fbuf[omp_get_thread_num()][indx[j]]-=z[i]*val[j];
			}

		}

	}
	*/
}//inverseL_ITL

// Обратный ход по разреженной верхнетреугольной матрице U.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// U=transpose(L);
// U - хранится в следующем виде:
// 1. udiag - диагональные элементы U.
// 2. uutr - наддиагональные элементы в столбцовом формате,
// т.е. хранение постолбцовое. 
// Так портрет симметричен, то:
// 3. jptr - соотвествующие номера столбцов для lltr, 
// 4. iptr - информация о начале следующей строки для lltr.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(U)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// udiag==ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// uutr==lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* inverseU(Real* f, Real* udiag, Real* uutr, int* jptr, int* iptr, int n) {
	Real *z=new Real[n];

    if (z == NULL)
	{
		printf("malloc: out of memory for vector z in inverse(U)*f \n"); // нехватка памяти
		system("pause");
		exit(0);
		return NULL; // завершение программы
	}

	int i,j;
	for (i=(n-1); i>=0; i--) {
        z[i]=f[i]/udiag[i];
		// Обработка i-го столбца над диагональю:
		for (j=iptr[i]; j<iptr[i+1]; j++) {
			f[jptr[j]]-=z[i]*uutr[j];
		}
		
	}
	return z;
}//inverseU

// Обратный ход по разреженной верхнетреугольной матрице U.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// U=transpose(L); - верхняя треугольная матрица.
// U - хранится в следующем виде:
// 1. val - диагональные и наддиагональные элементы U (в строковом формате).
// 2. indx - соотвествующие номера столбцов, 
// 3. pntr - информация о начале следующей строки для val.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(U)*f;
// Вектор f портится.
// пример (CSIR_ITL - формат):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseU_ITL0(Real* f, Real* val, int* indx, int* pntr, Real* &z, int n) {

    if (z == NULL)
	{
		z = new Real[n];
		if (z==NULL) {
			printf("malloc: out of memory for vector z in inverse(U)*f \n"); // нехватка памяти
		    system("pause");
		    exit(0); // завершение программы
		}
	}

	int i,j;
	for (i=(n-1); i>=0; i--) {
        
		// Обработка i-ой строки:
		for (j=pntr[i]+1; j<pntr[i+1]; j++) {
			f[i]-=z[indx[j]]*val[j];
		}
		// делим на диагональный элемент:
        z[i]=f[i]/val[pntr[i]];
		
	}
	
}//inverseU_ITL

// Обратный ход по разреженной верхнетреугольной матрице U.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// U=transpose(L); - верхняя треугольная матрица.
// U - хранится в следующем виде:
// 1. val - диагональные и наддиагональные элементы U (в строковом формате).
// 2. indx - соответствующие номера столбцов, 
// 3. pntr - информация о начале следующей строки для val.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(U)*f;
// Вектор f портится.
// пример (CSIR_ITL - формат):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseU_ITL(Real* f, Real* val, int* indx, int* pntr, Real*& z, int n) {

	if (z == nullptr)
	{
		z = new Real[n];
		if (z == nullptr) {
			printf("malloc: out of memory for vector z in inverse(U)*f \n"); // нехватка памяти
			//getchar();
			system("pause");
			exit(0); // завершение программы
		}
	}




	// 09.08.2021

	if (inumcore == 1) {

		for (int i = (n - 1); i >= 0; --i) {

			// Обработка i-ой строки:
			// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
			//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
			for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
				f[i] -= z[indx[j]] * val[j];
			}
			// делим на диагональный элемент:
			z[i] = f[i] / val[pntr[i]];

		}

	}

	if (inumcore == 2) {

		for (int i = (s_par[3].e - 1); i >= s_par[3].s; --i) {

			//if (inumerate[i] == 3)
			{

				// Обработка i-ой строки:
				// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
				//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
				for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
					f[i] -= z[indx[j]] * val[j];
				}
				// делим на диагональный элемент:
				z[i] = f[i] / val[pntr[i]];

			}
		}

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (int i = (s_par[1].e - 1); i >= s_par[1].s; --i) {

					//if (inumerate[i] == 1) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (int i = (s_par[2].e - 1); i >= s_par[2].s; --i) {

					//if (inumerate[i] == 2) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

	}


	if (inumcore == 3) {

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (int i = (s_par[4].e - 1); i >= s_par[4].s; --i) {

					//if (inumerate[i] == 4) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (int i = (s_par[5].e - 1); i >= s_par[5].s; --i) {

					//if (inumerate[i] == 5) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}



#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (int i = (s_par[1].e - 1); i >= s_par[1].s; --i) {

					//if (inumerate[i] == 1) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (int i = (s_par[2].e - 1); i >= s_par[2].s; --i) {

					//if (inumerate[i] == 2) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}

#pragma omp section 
			{
				for (int i = (s_par[3].e - 1); i >= s_par[3].s; --i) {

					//if (inumerate[i] == 3)
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

	}


	if ((inumcore == 4) || (inumcore == 5)) {

		for (int i = (s_par[3].e - 1); i >= s_par[3].s; --i) {

			//if (inumerate[i] == 3) 
			{

				// Обработка i-ой строки:
				// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
				//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
				for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
					f[i] -= z[indx[j]] * val[j];
				}
				// делим на диагональный элемент:
				z[i] = f[i] / val[pntr[i]];

			}
		}

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (int i = (s_par[5].e - 1); i >= s_par[5].s; --i) {

					//if (inumerate[i] == 5) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (int i = (s_par[8].e - 1); i >= s_par[8].s; --i) {

					//if (inumerate[i] == 8)
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}




#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (int i = (s_par[4].e - 1); i >= s_par[4].s; --i) {

					//if (inumerate[i] == 4)
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (int i = (s_par[6].e - 1); i >= s_par[6].s; --i) {

					//if (inumerate[i] == 6) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}

#pragma omp section 
			{
				for (int i = (s_par[7].e - 1); i >= s_par[7].s; --i) {

					//if (inumerate[i] == 7) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (int i = (s_par[9].e - 1); i >= s_par[9].s; --i) {

					//if (inumerate[i] == 9) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (int j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

	}


}//inverseU_ITL



// Переводит из формата CSIR в формат CSIR_ITL
// форматы:
// CSIR: ldiag, lltr, jptr, iptr
// CSIR_ITL: val, indx, pntr
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
//Формируем разреженный формат CSIR_ITL
//val : 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
//indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
//pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
void convertCSIRtoCSIR_ITL(Real *ldiag, Real *lltr, int *jptr, int *iptr, int n, int nz, Real* &val, int* &indx, int* &pntr, int nnz) {
	int i,j,k;
	//nnz=n+nz; // размер массивов val и indx
	// выделение оперативной памяти:
	val = new Real[nnz];
	indx = new int[nnz];
	pntr = new int[n+1];
	for (i=0; i<=n; i++) pntr[i]=nnz;

	if ((val == NULL) || (indx == NULL) || (pntr == NULL))
	{
		printf("malloc: out of memory in convertCSIRtoCSIR_ITL \n"); // нехватка памяти
		system("pause");
		exit(0); // завершение программы
	}

	// Алгоритм :
	// По порядку для всех столбцов формата CSIR_ITL
	int ic=0; // счётчик ненулевых элементов
	for (k=0; k<n; k++) {
		// добавление диагонального элемента k - го стобца
		val[ic]=ldiag[k];
		indx[ic]=k;
		pntr[k]=min(ic,pntr[k]);
		ic++;

		// добавление остальных элементов k-го столбца
		// сканирование матрицы в CSIR формате:
		for (i=1; i<n; i++) {
			for (j=iptr[i]; j<iptr[i+1]; j++)
				if (jptr[j] == k) {
					// добавление элемента в k-ый столбец
					val[ic]=lltr[j];
					indx[ic]=i;
                    pntr[k]=min(ic,pntr[k]);
					ic++;
				}
		}

	}

} // convertCSIRtoCSIR_ITL

// Неполное разложение Холецкого
// для положительно определённой симметричной
// матрицы А размером nxn.
// n - размерность матрицы СЛАУ
// Матрица val изменяется и в ней возвращается
// неполное разложение Холецкого IC(0):
// val == U верхняя треугольная матрица
// A = transpose(U)*U=L*transpose(L);
// L=transpose(U);
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
//формат CSIR_ITL (верхний треугольник хранится построчно).
// val : 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
// pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
// Результат факторизации без заполнения:
// изменённый массив val (indx и pntr остались без изменений):
// val (factorization)= 
// 3.0
// 1.0
// 0.3333333333333333
// 0.3333333333333333
// 3.3166247903554
// 0.6030226891555273
// 0.30151134457776363
// 0.6030226891555273
// 3.1622776601683795
// 0.6324555320336759
// 2.932575659723036
// 0.34099716973523675
// 3.4472773213410837
// 0.2578524458667825
// 2.8284271247461903
// 2.7310738989293286
//-------------------------------------------
void IC0Factor_ITL(Real* val, int* indx, int* pntr, int n)
{
  int d, g, h, i, j, k;
  Real z;

  for (k = 0; k < n - 1; k++) {
    d = pntr[k];
    z = val[d] = sqrt(val[d]);

    for (i = d + 1; i < pntr[k+1]; i++)
      val[i] /= z;

    for (i = d + 1; i < pntr[k+1]; i++) {
      z = val[i];
      h = indx[i];
      g = i;

      for (j = pntr[h] ; j < pntr[h+1]; j++)
        for ( ; g < pntr[k+1] && indx[g+1] <= indx[j]; g++)
          if (indx[g] == indx[j])
             val[j] -= z * val[g];
    }
  }
  d = pntr[n-1];
  val[d] = sqrt(val[d]);
} // IC0Factor_ITL

// Модифицированное неполное разложение Холецкого.
void IC0FactorModify_ITL(Real* val, int* indx, int* pntr, int n)
{
  int d, g, h, i, j, k;
  Real z, accumulate_fill_in;

  for (k = 0; k < n - 1; k++) {
    d = pntr[k];
    z = val[d] = sqrt(val[d]);

    for (i = d + 1; i < pntr[k+1]; i++)
      val[i] /= z;

    for (i = d + 1; i < pntr[k+1]; i++) {
      z = val[i];
      h = indx[i];
      g = i;

      accumulate_fill_in = 0.0;

      for (j = pntr[h] ; j < pntr[h+1]; j++)
        for ( ; g < pntr[k+1] && indx[g+1] <= indx[j]; g++)
          if (indx[g] == indx[j]) // номера столбцов равны
             val[j] -= z * val[g];
	  else //index does not match accumulate the fill-in value
		  accumulate_fill_in += z * val[g];

	  val[pntr[h]] -= accumulate_fill_in;

    }
  }
  d = pntr[n-1];
  val[d] = sqrt(val[d]);
} // IC0FactorModify_ITL

// Переводит из формата CSIR_ITL в формат CSIR (обратное преобразование)
// Память под все массивы предполагается выделенной заранее!!!
// форматы:
// CSIR_ITL: val, indx, pntr
// CSIR: ldiag, lltr, jptr, iptr
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
//Формируем разреженный формат CSIR_ITL
//val : 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
//indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
//pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
// формат CSIR:
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void convertCSIR_ITLtoCSIR(Real* ldiag, Real* lltr, int* jptr, int* iptr, int n, int nz, Real* val, int* indx, int* pntr, int nnz) {
	int i,j,k;//,k1;
	int imin=1;
	//nz=nnz-n; // размер массивов lltr и jptr
	// память предполагается выделенной заранее!!!
	// jptr и iptr изменяться не будут
	for (i=0; i<n; i++) ldiag[i]=0.0;
	for (i=0; i<nz; i++) {
		lltr[i]=0.0;
		//jptr[i]=0;
	}
	//for (i=0; i<=n; i++) iptr[i]=nz;


	// Алгоритм :
	// По порядку для всех строк формата CSIR
	int ic=0; // счётчик ненулевых элементов
	for (k=0; k<n; k++) {
		// добавление диагонального элемента k - ой строки
		ldiag[k]=val[pntr[k]];

		// добавление остальных элементов k-ой строки
		// сканирование матрицы в CSIR_ITL формате:
		for (i=0; i<n-1; i++) {
			for (j=pntr[i]+1; j<pntr[i+1]; j++)
				if (indx[j] == k) {
					// добавление элемента в k-ую строку
					lltr[ic]=val[j];
					//jptr[ic]=i;
					//imin=min(ic,iptr[k]);
                    //iptr[k]=imin;
					//if (imin==0) {
					//	for (k1=0; k1<k; k1++) iptr[k1]=0;
					//}
					ic++;
				}
		}

	}

} // convertCSIR_ITLtoCSIR

// неполное разложение Холецкого IC(0).
// входные данные нижний треугольник симметричной матрицы в формате CSIR.
// Внутри программы идут преобразования к формату CSIR_ITL библиотеки шаблонов ITL.
void ICFactor0(Real* ldiag, Real* lltr, int* jptr, int* iptr, int n, int nz) {
    
	Real *val;
	int *indx, *pntr;

	// внутри происходит выделение памяти
	// преобразование (прямое и обратное) ресурсоёмкая операция для больших матриц,
	// поэтому от неё нужно отказаться.
	convertCSIRtoCSIR_ITL(ldiag, lltr, jptr, iptr, n, nz, val, indx, pntr, n+nz);
	printf("Incoplete Cholesky 49.9%%...\n");
	IC0Factor_ITL(val, indx, pntr, n);
	printf("Incoplete Cholesky 50%%...\n");
    convertCSIR_ITLtoCSIR(ldiag, lltr, jptr, iptr, n, nz, val, indx, pntr, n+nz);
	printf("Incoplete Cholesky 100%%...\n");

	// освобождение памяти
	delete val; delete indx; delete pntr;
} // ICFactor0


// умножение симметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR. В силу симметрии хранятся только поддиагональные элементы altr. 
// Разреженная SPD матрица A (adiag, altr, jptr, iptr) квадратная размером nxn.
// Число уравнений равно числу неизвестных и равно n.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void  SPDMatrixCSIRByVector0(Real* adiag, Real* altr, int* jptr, int* iptr, Real* V, Real* &tmp, int n)
{
	
	// вектор tmp индексируется начиная с нуля так же как и вектор V
	if (tmp == NULL)
	{
		printf("in SPDMatrixCSIRByVector tmp==NULL\n");
		system("pause");
		tmp =new Real[n];
		if (tmp==NULL) {
			printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
		    system("pause");
		    exit(0); // завершение программы
		}
	}
	
	
    int i,j; // Счётчики цикла
    

	omp_set_num_threads(inumcore);

    #pragma omp parallel for shared(tmp, V, adiag) private(i) schedule (guided)
	for (i=0; i<n; i++) tmp[i]=V[i]*adiag[i];

    
	for (i=0; i<n; i++) {
	    for (j = iptr[i]; j<iptr[i+1]; j++)
		{
		    tmp[i] += V[jptr[j]]*altr[j];
		    tmp[jptr[j]] += V[i]*altr[j];
		}
	}
	
} // SPDMatrixCSIRByVector

// умножение симметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR. В силу симметрии хранятся только поддиагональные элементы altr. 
// Разреженная SPD матрица A (adiag, altr, jptr, iptr) квадратная размером n*n.
// Число уравнений равно числу неизвестных и равно n.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void  SPDMatrixCSIRByVector(Real* adiag, Real* altr, int* jptr, int* iptr, Real* V, Real*& tmp, int n)
{

	// вектор tmp индексируется начиная с нуля так же как и вектор V
	if (tmp == nullptr)
	{
		printf("in SPDMatrixCSIRByVector tmp==nullptr\n");
		//getchar();
		system("pause");
		tmp = new Real[n];
		if (tmp == nullptr) {
			printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
			//getchar();
			system("pause");
			exit(0); // завершение программы
		}
	}


	// int i,j; // Счётчики цикла

 /*
 #ifdef _OPENMP
	 omp_set_num_threads(inumcore);
 #endif
 */

#pragma omp parallel for shared(tmp, V, adiag)  
	for (int i = 0; i < n; ++i) tmp[i] = V[i] * adiag[i];

	// Последовательная секция
	/*
	for (i=0; i<n; i++) {
		for (j = iptr[i]; j<iptr[i+1]; j++)
		{
			tmp[i] += V[jptr[j]]*altr[j];
			tmp[jptr[j]] += V[i]*altr[j];
		}
	}
	*/

	// Часть первая из двух.
#pragma omp parallel for shared(tmp, V, altr, iptr, jptr,n) 
	for (int i = 0; i < n; ++i) {
		for (int j = iptr[i]; j < iptr[i + 1]; ++j)
		{
			tmp[i] += V[jptr[j]] * altr[j];
		}
	}

	if (inumcore == 1) {

		// Вторая часть не поддаётся распараллеливанию
		for (int i = 0; i < n; ++i) {

			// эта часть не поддаётся распараллеливанию.
			//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
			for (int j = iptr[i]; j < iptr[i + 1]; ++j)
			{
				tmp[jptr[j]] += V[i] * altr[j];
			}
		}
	}

	if (inumcore == 2) {


#pragma omp parallel sections
		{
#pragma omp section
			{
				// Вторая часть не поддаётся распараллеливанию
				for (int i = s_par[1].s; i < s_par[1].e; ++i) {

					// эта часть не поддаётся распараллеливанию.
					//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				// Вторая часть не поддаётся распараллеливанию
				for (int i = s_par[2].s; i < s_par[2].e; ++i) {

					// эта часть не поддаётся распараллеливанию.
					//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}

		// Вторая часть не поддаётся распараллеливанию
		for (int i = s_par[3].s; i < s_par[3].e; ++i) {

			// эта часть не поддаётся распараллеливанию.
			//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
			for (int j = iptr[i]; j < iptr[i + 1]; ++j)
			{
				tmp[jptr[j]] += V[i] * altr[j];
			}
		}
	}


	if (inumcore == 3) {


#pragma omp parallel sections
		{
#pragma omp section
			{
				for (int i = s_par[1].s; i < s_par[1].e; ++i) {

					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[2].s; i < s_par[2].e; ++i) {

					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[3].s; i < s_par[3].e; ++i) {

					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}


#pragma omp parallel sections
		{
#pragma omp section
			{
				for (int i = s_par[4].s; i < s_par[4].e; ++i) {

					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[5].s; i < s_par[5].e; ++i) {

					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}



	}

	if ((inumcore == 4) || (inumcore == 5)) {

#pragma omp parallel sections
		{
#pragma omp section
			{
				for (int i = s_par[4].s; i < s_par[4].e; ++i) {
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[6].s; i < s_par[6].e; ++i) {
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[7].s; i < s_par[7].e; ++i) {
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[9].s; i < s_par[9].e; ++i) {
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}

#pragma omp parallel sections
		{
#pragma omp section
			{
				for (int i = s_par[5].s; i < s_par[5].e; ++i) {
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (int i = s_par[8].s; i < s_par[8].e; ++i) {
					for (int j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}

		for (int i = s_par[3].s; i < s_par[3].e; ++i) {
			for (int j = iptr[i]; j < iptr[i + 1]; ++j)
			{
				tmp[jptr[j]] += V[i] * altr[j];
			}
		}
	}

} // SPDMatrixCSIRByVector

// умножение несимметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR.  
// Разреженная матрица A (adiag, altr, autr, jptr, iptr) квадратная размером nxn.
// Число уравнений равно числу неизвестных и равно n.
// Диагональ adiag хранится отдельно. Нижний треугольник altr хранится построчно.
// Верхний треугольник хранится по столбцам autr. Портрет матрицы (позиции ненулевых 
// элементов ) предполагается симметричным. Массив jptr - номера столбцов для нижнего 
// треугольника, массив iptr - показывает где начинаются новые строки для нижнего треугольника.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   1.0   10.0   2.0   0.0   0.0   0.0    
// 2.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 2.0   2.0   0.0   0.0   3.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 1.0  2.0 1.0 2.0  1.0 1.0  2.0 2.0 3.0
// autr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* MatrixCSIRByVector(Real* adiag, Real* altr, Real* autr, int* jptr, int* iptr, Real* V, int n)
{
	
	Real* tmp=new double[n]; // вектор индексируется начиная с нуля так же как и вектор V
	if (tmp == NULL)
	{
		printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
		system("pause");
		exit(0);
		return NULL; // завершение программы
	}
	
	
    int i,j; // Счётчики цикла

	for (i=0; i<n; i++) tmp[i]=V[i]*adiag[i];

    
	for (i=0; i<n; i++) {
	    for (j = iptr[i]; j<iptr[i+1]; j++)
		{
		    	tmp[i] += V[jptr[j]]*altr[j];
		        tmp[jptr[j]] += V[i]*autr[j];
		}
	}
	
	return tmp;
} // MatrixCSIRByVector

// умножение транспонированной несимметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR.  
// Разреженная матрица A (adiag, altr, autr, jptr, iptr) квадратная размером nxn. Хранится 
// именно исходная матрица, а умножается её транспонированный вариант.
// Число уравнений равно числу неизвестных и равно n.
// Диагональ adiag хранится отдельно. Нижний треугольник altr хранится построчно.
// Верхний треугольник хранится по столбцам autr. Портрет матрицы (позиции ненулевых 
// элементов ) предполагается симметричным. Массив jptr - номера столбцов для нижнего 
// треугольника, массив iptr - показывает где начинаются новые строки для нижнего треугольника.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   1.0   10.0   2.0   0.0   0.0   0.0    
// 2.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 2.0   2.0   0.0   0.0   3.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 1.0  2.0 1.0 2.0  1.0 1.0  2.0 2.0 3.0
// autr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
Real* MatrixTransposeCSIRByVector(Real* adiag, Real* altr, Real* autr, int* jptr, int* iptr, Real* V, int n)
{
	
	Real* tmp=new double[n]; // вектор индексируется начиная с нуля так же как и вектор V
	if (tmp == NULL)
	{
		printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
		system("pause");
		exit(0);
		return NULL; // завершение программы
	}
	
	
    int i,j; // Счётчики цикла

	for (i=0; i<n; i++) tmp[i]=V[i]*adiag[i];

    
	for (i=0; i<n; i++) {
	    for (j = iptr[i]; j<iptr[i+1]; j++)
		{
		    	tmp[i] += V[jptr[j]]*autr[j];
		        tmp[jptr[j]] += V[i]*altr[j];
		}
	}
	
	return tmp;
} // MatrixTransposeCSIRByVector


/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  adiag, altr, jptr, iptr - разреженная матрица СЛАУ в формате CSIR,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или NULL.
*  n - размерность СЛАУ Anxn.
*  nz - размерность массивов altr, jptr.
*  Разреженная матрица A (adiag, altr, jptr, iptr) квадратная размером nxn.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица A полагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Хранится только нижний треугольник с диагональю altr и adiag.
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*  В качестве предобуславливателя работает неполное разложение Холецкого:
*  M^(-1)==transpose(L)^(-1)*L^(-1); // обращённый предобуславливатель.
*  
*/
Real *SoprGradCSIR(Real* adiag, Real* altr, int* jptr, int* iptr, Real *dV, Real *x, int n, int nz){

	printf("Reshenie metodom sopryjennyh gradientov:\n");
	int k=0;
	int i; // счётчик
	Real *ap=new Real[n], *vcopy=new Real[n],
		 *z=new Real[n], *p=new Real[n];
    Real a, b, res;
	
	// для неполного разложения Холецкого:
	Real  *ldiag=new Real[n], *lltr=new Real[nz];
	int *jptrsort=new int[nz];
	Real *f=new Real[n];

	Real dold, dnew;
	

	
	// инициализация
	for (i=0; i<n; i++) ldiag[i]=adiag[i];
	for (i=0; i<nz; i++) lltr[i]=altr[i];
	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	printf("Incoplete Cholesky decomposition beginig...:\n");
    ICFactor0(ldiag, lltr, jptr, iptr, n, nz);
	printf("Incoplete Cholesky decomposition finish...:\n");//*/

    
	for (i=0; i<nz; i++) jptrsort[i]=jptr[i];
	for (i=0; i<n; i++) QuickSort(jptrsort, iptr[i], iptr[i+1]-1);
    //printf("jptrsort...\n");
	//for (i=0; i<nz; i++) printf("%d ",jptrsort[i]); system("pause");



	// шаг 1.1
	//X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, ap, n);
	for (i=0; i<n; i++) z[i]=dV[i]-ap[i];
	for (i=0; i<n; i++) vcopy[i]=z[i];
    f=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
    for (i=0; i<n; i++) vcopy[i]=f[i]; delete[] f; 
	f=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
    dnew=Scal(z,f,n);

	if (fabs(dnew)>1e-100){
		// шаг 1.3
	   for (i=0; i<n; i++)	p[i]=f[i];
	   res=1000.;
	   while ((fabs(res)>e) && (k<1000)) {
		   // шаг 2.1
		  SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, ap, n);

		  // шаг 2.2
		  a=dnew/Scal(p,ap,n);// шаговый множитель
		  // шаг 2.3 и 2.4
		  for (i=0; i<n; i++) {
		      x[i]+=a*p[i]; // очередное приближение 
              z[i]-=a*ap[i];// невязка k+1-го приближения
		  }
          for (i=0; i<n; i++) vcopy[i]=z[i]; delete[] f; 
          f=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
          for (i=0; i<n; i++) vcopy[i]=f[i]; delete[] f; 
	      f=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
		  // шаг 2.5
          dold=dnew;
		  dnew=Scal(z,f,n);

		  
		  res=dnew;
		  if (k%10==0) printf("iter residual\n");
		  printf(" %d %e\n", k, res);
		  // шаг 3.1
		  b=dnew/dold;
		  // шаг 3.2
		  for (i=0; i<n; i++) {
		     p[i]=f[i]+b*p[i]; // новое направление минимизации
		  }
          // шаг 3.3
		  k++;
	   } // while

	   // Освобождение памяти
        delete[] ap; delete[] vcopy;
		delete[] z; delete[] p; delete[] f;

	   return x;
	}
	else {
		// Освобождение памяти
		delete[] ap; delete[] vcopy;
		delete[] z; delete[] p; delete[] f;

		return x;
	}
} // SoprGradCSIR


// простая реализация явно преобразующая матрицу СЛАУ А.
// Матрица СЛАУ А задаётся в CSIR формате : adiag, altr, jptr, iptr.
// Неполное разложение Холецкого для А представляет её приближённо в виде:
// A = L*transpose(L); с нулевым заполнением. Массивы jptr и  iptr остаются теми же.
// Тогда матрица : A~=inverse(L)*A*inverse(transpose(L)) тоже симметрична и положительно определена.
// Правая часть преобразованной системы имеет вид: dV~=inverse(L)*dV.
// Решение СЛАУ тогда равно A~*x~=dV~; => x~=transpose(L)*x; => x=inverse(transpose(L))*x~;
// Предобуславливание неполным разлождением Холецкого уменьшает количество итераций при решении СЛАУ,
// улучшает спектральные характеристики матрицы СЛАУ.
Real *SoprGradCSIR2(Real* adiag, Real* altr, int* jptr, int* iptr, Real *dV, Real *x, int n, int nz0){
	printf("Reshenie metodom sopryjennyh gradientov:\n");
	int k=0;
	int i; // счётчик
	Real *ap, *vcopy=new Real[n],
		 *z=new Real[n], *p=new Real[n];
	Real a, b, nz;

    // для неполного разложения Холецкого:
	Real  *ldiag=new Real[n], *lltr=new Real[nz0];
	int *jptrsort=new int[nz0];


    // инициализация
	for (i=0; i<n; i++) ldiag[i]=adiag[i];
	for (i=0; i<nz0; i++) lltr[i]=altr[i];
	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	printf("Incoplete Cholesky decomposition beginig...:\n");
    ICFactor0(ldiag, lltr, jptr, iptr, n, nz0);
	printf("Incoplete Cholesky decomposition finish...:\n");//*/
   


   /*
	ldiag[0]=1.0; ldiag[1]=1.0;  ldiag[2]=1.838477; ldiag[3]=2.00055;
    ldiag[4]=0.590477; ldiag[5]=1.0;  ldiag[6]=1.0;
	lltr[0]=-1.22383866; lltr[1]=-0.5439282932;  lltr[2]=-1.33247070; //*/
    
    /* // переставлены элементы
	ldiag[0]=1.0; ldiag[1]=1.0;  ldiag[2]=1.838477; ldiag[3]=2.00055;
    ldiag[4]=0.590477; ldiag[5]=1.465913;  ldiag[6]=0.37585673;
	lltr[0]=-1.22383866; lltr[1]=-1.33247070;  lltr[2]=-0.5439282932; lltr[3]=-0.1457305633;
    lltr[4]=-0.4998613742; lltr[5]=-1.401073265;  lltr[6]=-0.06498197865;//*/

	for (i=0; i<nz0; i++) jptrsort[i]=jptr[i];
	for (i=0; i<n; i++) QuickSort(jptrsort, iptr[i], iptr[i+1]-1);

	// шаг 1.1
	//X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	//ap=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, n);
	//for (i=0; i<n; i++) z[i]=dV[i]-ap[i];

	for(i=0;i<n;i++) vcopy[i]=x[i]; 
    ap=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
	for(i=0;i<n;i++) vcopy[i]=ap[i]; 
    SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, ap, n);
    for(i=0;i<n;i++) vcopy[i]=ap[i]; delete ap;
    ap=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
    for(i=0;i<n;i++) vcopy[i]=dV[i]; delete dV;
	dV=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
    
    for (i=0; i<n; i++) z[i]=dV[i]-ap[i];

	if (Scal(z,z,n)!=0){
		// шаг 1.3
	   for (i=0; i<n; i++)	p[i]=z[i];
	   nz=1000.;
	   while ((nz>e) && (k<1000)) {
		   // шаг 2.1
	 	  //ap=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, n);

		   delete ap; // освобождение памяти
           for(i=0;i<n;i++) vcopy[i]=p[i];
          ap=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
           for(i=0;i<n;i++) vcopy[i]=ap[i]; 
          SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, ap, n);
          for(i=0;i<n;i++) vcopy[i]=ap[i]; delete ap;
          ap=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);

		  // шаг 2.2
		  //a=Scal(z,p,n)/Scal(z,ap,n);
		  a=Scal(z,p,n)/Scal(ap,p,n); // шаговый множитель
		  // шаг 2.3 и 2.4
		  for (i=0; i<n; i++) {
		      x[i]+=a*p[i]; // очередное приближение
			  z[i]-=a*ap[i]; // невязка k+1-го приближения
		  }
		  // шаг 2.5
		  nz=NormaV(z,n);
		  if (k%10==0) printf("iter residual\n");
		  printf(" %d %e\n", k, nz);
		  // шаг 3.1
		  b=Scal(z,ap,n)/Scal(p,ap,n);
		  // шаг 3.2
		  for (i=0; i<n; i++) {
		     p[i]=z[i]-b*p[i]; // новое направление минимизации
		  }
          // шаг 3.3 
		  k++;
	   } // while

	   // Освобождение памяти
        delete[] ap; delete[] vcopy;
		delete[] z; delete[] p;

		for(i=0;i<n;i++) vcopy[i]=x[i]; delete[] x;
		x=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
	   return x;
	}
	else {
		// Освобождение памяти
		delete[] ap; delete[] vcopy;
		delete[] z; delete[] p;

		return x;
	}
} // SoprGradCSIR2

/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  M - разреженная матрица СЛАУ в формате SIMPLESPARSE,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или NULL.
*  n - размерность СЛАУ Anxn.
*
*  Разреженная матрица M квадратная размером nxn.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица M предполагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Хранятся только ненулевые элементы. 
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*  В качестве предобуславливателя работает неполное разложение Холецкого:
*  K^(-1)==transpose(L)^(-1)*L^(-1); // обращённый предобуславливатель.
*  
*/
void ICCG(SIMPLESPARSE &M, Real *dV, Real* &x, int n)
{

	if (DEBUG) printf("Reshenie metodom sopryjennyh gradientov:\n");
    // матрица СЛАУ
	// в формате CSIR:
	Real *adiag=NULL, *altr=NULL;
	int *jptr=NULL, *iptr=NULL;

	// предобуславливатель:
	// неполным разложением Холесского в
	// формате CSIR_ITL:
	Real *val=NULL;
	int *indx=NULL, *pntr=NULL;
	
	int k=0;
	int i; // счётчик
	Real *ap=new Real[n], *vcopy=new Real[n], *f=new Real[n],
		 *z=new Real[n], *p=new Real[n];
    Real a, b, res;
	

	Real dold, dnew;
	
	
	// инициализация
	// Память выделяется внутри:
	//simplesparsetoCSIR(M, adiag, altr, jptr, iptr, n);
	//simplesparsetoCSIR_ITLSPD(M, val, indx, pntr, n);
	
	// инициализация
	// Память выделяется внутри:

	ell_to_CSIR(adiag, altr, jptr, iptr, n);
	ell_to_CSIR_ITLSPD(val, indx, pntr, n);


	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	if (DEBUG) printf("Incoplete Cholesky decomposition beginig...:\n");
	//IC0Factor_ITL(val, indx, pntr, n);
	IC0FactorModify_ITL(val, indx, pntr, n);
	if (DEBUG) printf("Incoplete Cholesky decomposition finish...:\n");//*/


	// шаг 1.1
	//X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, ap, n);
	for (i=0; i<n; i++) z[i]=dV[i]-ap[i];
	for (i=0; i<n; i++) vcopy[i]=z[i]; 
    inverseL_ITL(vcopy, val, indx, pntr, f, n);
    for (i=0; i<n; i++) vcopy[i]=f[i];  
	inverseU_ITL(vcopy, val, indx, pntr, f, n);
    dnew=Scal(z,f,n);
	

	if (fabs(dnew)>1e-37){
		// шаг 1.3
	   for (i=0; i<n; i++)	p[i]=f[i];
	   res=1000.;
	   while ((fabs(res)>e) && (k<1000)) {
		   // шаг 2.1
		  SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, ap, n);

		  // шаг 2.2
		  if (fabs(dnew) < 1.0e-24) {
			  if (fabs(Scal(p, ap, n)) < 1.0e-24) {
				  a = 1.0;
			  }
			  else {
				  a = 0.0;
			  }
		  }
		  else {
			  a = dnew / Scal(p, ap, n);// шаговый множитель
		  }
		  // шаг 2.3 и 2.4
		  for (i=0; i<n; i++) {
		      x[i]+=a*p[i]; // очередное приближение 
              z[i]-=a*ap[i];// невязка k+1-го приближения
		  }
          for (i=0; i<n; i++) vcopy[i]=z[i];  
          inverseL_ITL(vcopy, val, indx, pntr, f, n);
          for (i=0; i<n; i++) vcopy[i]=f[i]; 
	      inverseU_ITL(vcopy, val, indx, pntr, f, n);
		  // шаг 2.5
          dold=dnew;
		  dnew=Scal(z,f,n);

		  
		  res=dnew;
		  if (DEBUG) {
			  if (k % 10 == 0) printf("iter residual\n");
			  printf(" %d %e\n", k, res);
		  }
		  // шаг 3.1
		  b=dnew/dold;
		  // шаг 3.2
		  for (i=0; i<n; i++) {
		     p[i]=f[i]+b*p[i]; // новое направление минимизации
		  }
          // шаг 3.3
		  k++;
	   } // while

	   // Освобождение памяти
        delete[] ap; delete[] vcopy;
		delete[] z; delete[] p; delete[] f;  
	}
	else {
		// Освобождение памяти
		delete[] ap; delete[] vcopy;
		delete[] z; delete[] p; delete[] f;		
	}

	// Освобождение памяти
	delete[] adiag; 
	delete[] altr;
	delete[] jptr;
	delete[] iptr;

	delete[] val; 
	delete[] indx;
	delete[] pntr;
	
} // ICCG

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
Real* SoloveichikAlgCSIR_SPD(int isize, // размер квадратной матрицы
						Real* adiag, Real* altr, int* jptr, int* iptr, // матрица СЛАУ
                         Real *dV,  // вектор правой части
                         const Real *dX0, // вектор начального приближения
                         bool bconsole_message) // выводить ли значения невязки на консоль ?
{

     int i,k; // счётчики цикла for
     Real *dx, *dax, *dr, *dz, *dp, *dar1, *dres;
     Real dar, dbr, dnz, dscalp;
	 Real kend=1000; // ограничение на максимальное число итераций
	 Real epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new Real[isize]; dax=new Real[isize]; dr= new Real[isize];
    dz=new Real[isize]; dp=new Real[isize]; dar1=new Real[isize];
	dres=new Real[isize]; // вектор результата
   

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==NULL) {
	   for (i=0; i<isize; i++) dx[i]=0.0;
   }
   else {
	   for (i=0; i<isize; i++) dx[i]=dX0[i];
   }

   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dx, dax, isize); // результат занесён в  dax
   for (i=0; i<isize; i++) dr[i]= dV[i] - dax[i];  // начальная невязка
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   for (i=0; i<isize; i++) dz[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dz, dp, isize); // результат занесён в dp

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (dnz > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         for (i=0; i<isize; i++)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         dnz=dnz-dar*dar/dscalp; // норма невязки
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
            printf("%d %e \n",k,dnz);
		 } 
		 SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dr, dar1, isize);// результат занесён в dar1=A*dr
         dbr=-Scal(dp,dar1,isize)*dscalp;
         for (i=0; i<isize; i++)
		 {
            dz[i]=dr[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }
         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e7) 
		 {
            // восстановление начального приближения
            for (i=0; i<isize; i++) if (dX0==NULL) dx[i]=0.0; else dx[i]=dX0[i];
            printf("\n divergence Soloveichik solver \n");
            bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
      // возвращение результата
      for (i=0; i<isize; i++) dres[i]=dx[i];
   }
   else
   {
	   if (dX0 != NULL) {
		   // возвращает начальное приближение
#pragma omp parallel for
		   for (i = 0; i < isize; i++) dres[i] = dX0[i];
	   }
	   else {
#pragma omp parallel for
		   for (i = 0; i < isize; i++) dres[i] = 0.0;
	   }
   }

   // освобождение памяти выделенной под динамические массивы
   delete[] dx; delete[] dax; delete[] dr;
   delete[] dz; delete[] dp; delete[] dar1;

   return dres; 

} // SoloveichikAlgCSIR_SPD

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
Real* SoloveichikAlgCSIR_SPDgood(int isize, int nz0,// размер квадратной матрицы
						Real* adiag, Real* altr, int* jptr, int* iptr, // матрица СЛАУ
                         Real *dV,  // вектор правой части
                         const Real *dX0, // вектор начального приближения
                         bool bconsole_message) // выводить ли значения невязки на консоль ?
{

     int i,k; // счётчики цикла for
     Real *dx, *dax, *dr, *dz, *dp, *dar1, *dres, *df, *vcopy;
     Real dar, dbr, dnz, dscalp;
	 Real kend=1000; // ограничение на максимальное число итераций
	 Real epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new Real[isize]; dr= new Real[isize];
    dz=new Real[isize]; dp=new Real[isize]; dar1=new Real[isize];
	dres=new Real[isize]; vcopy=new Real[isize]; // вектор результата
	df=new Real[isize];
   


	// для неполного разложения Холецкого:
	Real  *ldiag=new Real[isize], *lltr=new Real[nz0];
	int *jptrsort=new int[nz0];


    // инициализация
	for (i=0; i<isize; i++) ldiag[i]=adiag[i];
	for (i=0; i<nz0; i++) lltr[i]=altr[i];
	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	printf("Incoplete Cholesky decomposition beginig...:\n");
    ICFactor0(ldiag, lltr, jptr, iptr, isize, nz0);
	printf("Incoplete Cholesky decomposition finish...:\n");
    

   /*
	ldiag[0]=1.0; ldiag[1]=1.0;  ldiag[2]=1.838477; ldiag[3]=2.00055;
    ldiag[4]=0.590477; ldiag[5]=1.465913;  ldiag[6]=0.37585673;
	lltr[0]=-1.22383866; lltr[1]=-0.5439282932;  lltr[2]=-1.33247070; lltr[3]=-0.4998613742;
    lltr[4]=-0.1457305633; lltr[5]=-0.06498197865;  lltr[6]=-1.401073265;//*/

    //lltr[0]=-1.22383866; lltr[1]=-1.33247070;  lltr[2]=-0.5439282932; lltr[3]=-0.1457305633;
    //lltr[4]=-0.4998613742; lltr[5]=-1.401073265;  lltr[6]=-0.06498197865;

    for (i=0; i<nz0; i++) jptrsort[i]=jptr[i];
	//for (i=0; i<isize; i++) QuickSort(jptrsort, iptr[i], iptr[i+1]-1);

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==NULL) {
	   for (i=0; i<isize; i++) dx[i]=0.0;
   }
   else {
	   for (i=0; i<isize; i++) dx[i]=dX0[i];
   }

   //dax=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dx, isize); // результат занесён в  dax
   for (i=0; i<isize; i++) vcopy[i]=dx[i];
   dax=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, isize);
   for (i=0; i<isize; i++) vcopy[i]=dax[i];
   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, dax, isize);
   for (i=0; i<isize; i++) vcopy[i]=dax[i]; delete dax;
   dax=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);

   for (i=0; i<isize; i++) vcopy[i]=dV[i]; delete dV;
   dV=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);

   for (i=0; i<isize; i++) dr[i]= dV[i] - dax[i];  // начальная невязка
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   for (i=0; i<isize; i++) dz[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   //dp=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dz, isize); // результат занесён в dp
   for (i=0; i<isize; i++) vcopy[i]=dz[i]; 
   dp=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, isize);
   for (i=0; i<isize; i++) vcopy[i]=dp[i]; 
   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, dp, isize);
   for (i=0; i<isize; i++) vcopy[i]=dp[i]; delete[] dp;
   dp=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (fabs(dnz) > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         for (i=0; i<isize; i++)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         //dnz=dnz-dar*dar/dscalp; // норма невязки
		 dnz=Scal( dr, dr, isize);
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
            printf("%d %e \n",k,dnz);
		 } 
		 //dar1=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dr, isize);// результат занесён в dar1=A*dr
          
         for (i=0; i<isize; i++) vcopy[i]=dr[i];
         dar1=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, isize);
		 for (i=0; i<isize; i++) vcopy[i]=dar1[i]; 
         SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, dar1, isize); 
		 for (i=0; i<isize; i++) vcopy[i]=dar1[i]; delete[] dar1;
         dar1=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);


         dbr=-Scal(dp,dar1,isize)*dscalp;
         for (i=0; i<isize; i++)
		 {
            dz[i]=dr[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }
         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e7) 
		 {
            // восстановление начального приближения
            for (i=0; i<isize; i++) if (dX0==NULL) dx[i]=0.0; else dx[i]=dX0[i];
            printf("\n divergence Soloveichik solver \n");
            bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
      // возвращение результата
      //for (i=0; i<isize; i++) dres[i]=dx[i];
	  dres=inverseU(dx, ldiag, lltr, jptrsort, iptr, isize);
   }
   else
   {
	   if (dX0 != NULL) {
		   // возвращает начальное приближение
#pragma omp parallel for
		   for (i = 0; i < isize; i++) dres[i] = dX0[i];
	   }
	   else {
		   // возвращает начальное приближение
#pragma omp parallel for
		   for (i = 0; i < isize; i++) dres[i] = 0.0;
	   }
   }

   // освобождение памяти выделенной под динамические массивы
   delete[] dx; delete[] dax; delete[] dr;
   delete[] dz; delete[] dp; delete[] dar1;
   delete[] vcopy;

   return dres; 

} // SoloveichikAlgCSIR_SPDgood

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
void SoloveichikAlgCRS(int isize, // размер квадратной матрицы
						 Real *val, int* col_ind, int* row_ptr, // матрица СЛАУ
                         Real *dV,  // вектор правой части
                         Real* &dX0, // вектор начального приближения
                         bool bconsole_message, int maxit) // выводить ли значения невязки на консоль ?
{

     int i,k; // счётчики цикла for
     Real *dx, *dax, *dr, *dz, *dp, *dar1, *dres, *dstart;
     Real dar, dbr, dnz, dscalp;
	 Real kend=maxit; // ограничение на максимальное число итераций
	 Real epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new Real[isize]; dax=new Real[isize]; dr= new Real[isize];
    dz=new Real[isize]; dp=new Real[isize]; dar1=new Real[isize];
	dres=new Real[isize], dstart=new Real[isize]; // вектор результата
   

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==NULL) {
	   for (i=0; i<isize; i++) { 
		   dx[i]=0.0;
		   dstart[i]=0.0;
	   }
	   dX0=new Real[isize];

   }
   else {
	   for (i=0; i<isize; i++) {
		   dx[i]=dX0[i];
           dstart[i]=dX0[i];
	   }
   }

   
   MatrixCRSByVector(val,col_ind,row_ptr,dx, dax,isize); // результат занесён в  dax
   for (i=0; i<isize; i++) dr[i]= dV[i] - dax[i];  // начальная невязка
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   for (i=0; i<isize; i++) dz[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   MatrixCRSByVector(val,col_ind,row_ptr,dz, dp, isize);// результат занесён в dp

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (dnz > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         for (i=0; i<isize; i++)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         dnz=dnz-dar*dar/dscalp; // норма невязки
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
            printf("%d %e \n",k,dnz);
		 } 
		 
		 MatrixCRSByVector(val,col_ind,row_ptr,dr, dar1, isize);// результат занесён в dar1=A*dr
         dbr=-Scal(dp,dar1,isize)*dscalp;
         for (i=0; i<isize; i++)
		 {
            dz[i]=dr[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }
         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e7) 
		 {
            // восстановление начального приближения
            for (i=0; i<isize; i++) if (dX0==NULL) dx[i]=0.0; else dx[i]=dstart[i];
            printf("\n divergence Soloveichik solver \n");
            bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
      // возвращение результата
      for (i=0; i<isize; i++) dres[i]=dx[i];
   }
   else
   {
      // возвращает начальное приближение
	  for (i=0; i<isize; i++) dres[i]=dstart[i];
   }

   // освобождение памяти выделенной под динамические массивы
   delete[] dx; delete[] dax; delete[] dr;
   delete[] dz; delete[] dp; delete[] dar1;

   //return dres;
   for (i=0; i<isize; i++) dX0[i]=dres[i];
   delete[] dres; delete[] dstart;

} // SoloveichikAlgCRS


/* Реализация на диннамическом массиве
// инициализирует разреженную матрицу
void initsimplesparse(SIMPLESPARSE &M) {
	M.a=NULL;
	M.n=0;
	M.incCLUSTER_SIZE=10;
	M.POOL_SIZE=0;
} // initsimplesparse
*/

// Реализация на связном списке
// инициализирует разреженную матрицу
void initsimplesparse(SIMPLESPARSE &M, int nodes) {
	M.n=0; // изначально все элементы нулевые 
	M.root=new NONZEROELEM*[nodes];
	int i; // номер строки, номер уравнения в СЛАУ
	for (i=0; i<nodes; i++) M.root[i]=NULL; 
} // initsimplesparse

/* Реализация на массиве.
// Добавляет ненулевой элемент в
// простейшую разряженную матрицу M
void addelmsimplesparse(SIMPLESPARSE &M, Real aij, int i, int j, bool bset) {
	if (M.n==0) {
		// первый элемент
		M.POOL_SIZE+=M.incCLUSTER_SIZE;
		M.n++;
		M.a=new NONZEROELEM[M.POOL_SIZE];
		M.a[0].aij=aij;
		M.a[0].i=i;
		M.a[0].j=j;
	}
	else if (M.n<M.POOL_SIZE) 
	{
		bool flag=false; // элемент не найден
		int i1; // счётчик
		for (i1=0; i1<M.n; i1++) if ((M.a[i1].i==i) && (M.a[i1].j==j)) {
           flag=true;
           if (bset) M.a[i1].aij=aij;  // установка
		   else M.a[i1].aij+=aij; // добавление
		}
		if (!flag) {
			M.a[M.n].aij=aij;
		    M.a[M.n].i=i;
		    M.a[M.n].j=j;
            M.n++;
		} 
	}
	else // M.n==M.POOL_SIZE
	{
        bool flag=false; // элемент не найден
		int i1; // счётчик
		for (i1=0; i1<M.n; i1++) if ((M.a[i1].i==i) && (M.a[i1].j==j)) {
           flag=true;
           if (bset) M.a[i1].aij=aij;  // установка
		   else M.a[i1].aij+=aij; // добавление
		}
		if (!flag) {
           NONZEROELEM* list=new NONZEROELEM[M.POOL_SIZE];
		   for (i1=0; i1<M.n; i1++) list[i1]=M.a[i1]; // копирование
		   delete M.a;
		   M.POOL_SIZE+=M.incCLUSTER_SIZE;
		   M.a=new NONZEROELEM[M.POOL_SIZE];
           for (i1=0; i1<M.n; i1++) M.a[i1]=list[i1]; // обратное копирование
           M.a[M.n].aij=aij;
		   M.a[M.n].i=i;
		   M.a[M.n].j=j;
		   M.n++;

		}
	}
} // addelmsimplesparse
*/



// Реализация на связном списке
// Добавляет ненулевой элемент в
// простейшую разряженную матрицу M
// Проверки на равенство добавляемого элемента нулю нет, поэтому
// может добавить и нулевой элемент.
void addelmsimplesparse(SIMPLESPARSE &M, Real aij, int i, int j, bool bset) {
    NONZEROELEM* p;
	p=M.root[i];
	// линейный поиск элемента с ключём key
	while ((p!=NULL) && (p->key!=j)) p=p->next;
	if (p!=NULL) {
		// элемент найден
		if (bset) p->aij=aij; // установка
		else p->aij+=aij; // добавление
	}
	else 
	{
		// если такого элемента нет в списке
		// то добавление элемента в начало списка.
        NONZEROELEM* q=new NONZEROELEM;
		q->aij=aij;
		q->key=j;
		q->next=M.root[i];
		M.root[i]=q;
		q=NULL;
		M.n++; // количество ненулевых элементов увеличилось на 1. 
	}
} // addelmsimplesparse

// освобождение памяти для матрицы SIMPLESPARSE
void simplesparsefree(SIMPLESPARSE &M, int nodes) {
	int i; // счётчик цикла for
	for (i=0; i<nodes; i++) {
        NONZEROELEM* p9, *q9;
        p9=M.root[i]; q9=p9;
		M.root[i]=NULL;
		while (p9!=NULL) {
			p9=p9->next;
			q9->next=NULL;
			delete q9;
			q9=p9;
		}
	}
	delete M.root;
} // simplesparsefree 

/*
// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
void swap(NONZEROELEM* &v, int i, int j)
{
        NONZEROELEM temp;

		// change v[i] <-> v[j]
		temp = v[i];
		v[i] = v[j];
		v[j] = temp;
} // swap

// Вот алгоритм PivotList
int PivotList(NONZEROELEM* &list, int first, int last) {
	// list обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	int PivotValue = list[first].key;
	int PivotPoint = first;

	for (int index=(first+1); index<=last; index++) {
		if (list[index].key<PivotValue) {
			PivotPoint++;
			swap(list, PivotPoint, index);
		}
	}

	swap(list, first, PivotPoint);

	return PivotPoint;
} // PivotList


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
void QuickSort(NONZEROELEM* &list, int first, int last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	int pivot;

	if (first < last) {
        pivot = PivotList(list, first, last);
        QuickSort(list, first, pivot-1);
		QuickSort(list, pivot+1, last);
	}
} // QuickSort
*/
// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
template <typename TVAL>
void swap(TVAL* &v, int i, int j)
{
	    TVAL temp;

		// change v[i] <-> v[j]
		temp = v[i];
		v[i] = v[j];
		v[j] = temp;
} // swap

// Вот алгоритм PivotList
template <typename TVAL>
int PivotList(TVAL* &list, int first, int last) {
	// list обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	TVAL PivotValue = list[first];
	int PivotPoint = first;

	for (int index=(first+1); index<=last; index++) {
		if (list[index]<PivotValue) {
			PivotPoint++;
			swap(list, PivotPoint, index);
		}
	}

	swap(list, first, PivotPoint);

	return PivotPoint;
} // PivotList


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
template <typename TVAL>
void QuickSort(TVAL* &list, int first, int last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	int pivot;

	if (first < last) {
        pivot = PivotList(list, first, last);
        QuickSort(list, first, pivot-1);
		QuickSort(list, pivot+1, last);
	}
} // QuickSort

// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
void swapCSIR(int* &v, Real* &dr, int i, int j)
{
        int tempi;
		Real tempr;

		// change v[i] <-> v[j]
		tempi = v[i];
		v[i] = v[j];
		v[j] = tempi;
		// change dr[i] <-> dr[j]
		tempr = dr[i];
		dr[i] = dr[j];
		dr[j] = tempr;

} // swap

// Вот алгоритм PivotList
int PivotListCSIR(int* &jptr, Real* &altr, int first, int last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	int PivotValue = jptr[first];
	int PivotPoint = first;

	for (int index=(first+1); index<=last; index++) {
		if (jptr[index]<PivotValue) {
			PivotPoint++;
			swapCSIR(jptr, altr, PivotPoint, index);
		}
	}

	swapCSIR(jptr, altr, first, PivotPoint);

	return PivotPoint;
} // PivotList


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
void QuickSortCSIR(int* &jptr, Real* &altr, int first, int last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	int pivot;

	if (first < last) {
        pivot = PivotListCSIR(jptr, altr, first, last);
        QuickSortCSIR(jptr, altr, first, pivot-1);
		QuickSortCSIR(jptr, altr, pivot+1, last);
	}
} // QuickSortCSIR

/* Реализация на динамическом массиве.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
void simplesparsetoCRS(SIMPLESPARSE &M, Real* &val, int* &col_ind, int* &row_ptr, int nodes) {
	if (M.n!=0) {
		val = new Real[M.n];
		col_ind = new int[M.n];
		row_ptr = new int[nodes+1];

		int k; // счётчик
		// инициализация
        for (k=0; k<(M.n); k++) {
		   val[k]=0.0;
		   col_ind[k]=0;
	    }
        for (k=0; k<=nodes; k++) {
		    row_ptr[k]=M.n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		QuickSort(M.a, 0, M.n-1);

		// заполнение разреженной матрицы
		for (k=0; k<M.n; k++) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
	}
} // simplesparsetoCRS
*/

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
void simplesparsetoCRS(SIMPLESPARSE &M, Real* &val, int* &col_ind, int* &row_ptr, int nodes) {
	bool flag=true;
    int k; // счётчик
	for (k=0; k<nodes; k++) if (M.root[k]==NULL) {
		flag=false; break;
	}

	if (flag) {
		val = new Real[M.n];
		col_ind = new int[M.n];
		row_ptr = new int[nodes+1];

		
		// инициализация
        for (k=0; k<(M.n); k++) {
		   val[k]=0.0;
		   col_ind[k]=0;
	    }
        for (k=0; k<=nodes; k++) {
		    row_ptr[k]=M.n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; k++) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
		*/
		int ik=0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
        for (k=0; k<nodes; k++) {
			p=M.root[k];
			while (p!=NULL) {
				val[ik]=p->aij;
				col_ind[ik]=p->key;
                row_ptr[k]=min(ik,row_ptr[k]);
				ik++;
				p=p->next;
			}
		}

		// в каждой строке элементы отсортированы по номерам столбцов:
        for (k=0; k<nodes; k++) QuickSortCSIR(col_ind, val, row_ptr[k], row_ptr[k+1]-1); 

	}
} // simplesparsetoCRS

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только нижний треугольник.
void simplesparsetoCSIR(SIMPLESPARSE &M, Real* &adiag, Real* &altr, int* &jptr, int* &iptr, int nodes) {
	bool flag=true;
    int k; // счётчик
	for (k=0; k<nodes; k++) if (M.root[k]==NULL) {
		flag=false; break;
	}

	if (flag) {
		// поддиагональные элементы в altr хранятся построчно
		int nz=(int)(M.n-nodes)/2; // число ненулевых элементов
		adiag = new Real[nodes]; // диагональные элементы
		altr = new Real[nz]; // поддиагональные элементы
		jptr = new int[nz]; // номера столцов для нижнего треугольника
		iptr = new int[nodes+1]; // указатели на следующую строку

		
		// инициализация
		for (k=0; k<nodes; k++) adiag[k]=0.0;
        for (k=0; k<(nz); k++) {
		   altr[k]=0.0;
		   jptr[k]=0;
	    }
        for (k=0; k<=nodes; k++) {
		    iptr[k]=nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; k++) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
		*/
		/*
		int ik=0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
        for (k=0; k<nodes; k++) {
			p=M.root[k];
			while (p!=NULL) {
				val[ik]=p->aij;
				col_ind[ik]=p->key;
                row_ptr[k]=min(ik,row_ptr[k]);
				ik++;
				p=p->next;
			}
		}
		*/

		int ik=0, imin=1,k1; // счётчик ненулевых поддиагональных элементов СЛАУ
		bool bvisit;
		NONZEROELEM* p;
        for (k=0; k<nodes; k++) {
			bvisit=false;
			p=M.root[k];
			while (p!=NULL) {
				if (p->key==k) {
					adiag[k]=p->aij;
				}
				else if (p->key<k) {
					if (ik<(nz)) {
						altr[ik]=p->aij; // ненулевое значение
					    jptr[ik]=p->key; // номер столбца
					}
					else {
						printf("non simmetric matrix ICCG. simplesparsetoCSIR\n");
						system("pause");
					}
					bvisit=true;			   
				}
				imin=min(ik,iptr[k]);
				//printf("imin=%d\n",imin);
                iptr[k]=imin;
                if (imin==0) for (k1=0; k1<k; k1++) iptr[k1]=0;	
				if (bvisit) { 
					ik++;
					bvisit=false;
				}
				p=p->next;
			}
		}


		for (k=0; k<nodes; k++) QuickSortCSIR(jptr, altr, iptr[k], iptr[k+1]-1);

	}
} // simplesparsetoCSIR


// печать матрицы в консоль
void printM_and_CSIR(SIMPLESPARSE &sparseM, int  n) {
	int i;
	// печать простейшей формы разреженной матрицы.
	for (i=0; i<n; i++) {
        NONZEROELEM* pelm=sparseM.root[i];
		while (pelm!=NULL) {
			printf("a[%d][%d]=%e  ",i,pelm->key,pelm->aij);
			pelm=pelm->next;
		}
		printf("\n");
	}//*/
	system("pause");
}

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR_ITL. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только верхний треугольник.
// Память выделяется внутри метода.
void simplesparsetoCSIR_ITLSPD(SIMPLESPARSE &M, Real* &val, int* &indx, int* &pntr, int nodes) {
	bool flag=true;
    int k; // счётчик
	for (k=0; k<nodes; k++) if (M.root[k]==NULL) {
		flag=false; break;
	}

	if (flag) {
 
		//printM_and_CSIR(M, nodes); // debug

		// поддиагональные элементы в altr хранятся построчно
		int nz=(int)((M.n-nodes)/2 + nodes); // число ненулевых элементов
		val = new Real[nz]; // диагональные элементы и наддиагональные элементы
		indx = new int[nz]; // номера столцов для нижнего треугольника
		pntr = new int[nodes+1]; // указатели на следующую строку

		
		// инициализация
        for (k=0; k<(nz); k++) {
		   val[k]=0.0;
		   indx[k]=0;
	    }
        for (k=0; k<=nodes; k++) {
		    pntr[k]=nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        

		int ik=0; // счётчик ненулевых поддиагональных элементов СЛАУ
		NONZEROELEM* p;
        for (k=0; k<nodes; k++) {
			
			p=M.root[k];
			while (p!=NULL) {

				// k - номер диагонального элемента
				if (p->key>=k) {
					if (ik<(nz)) {
						val[ik]=p->aij; // ненулевое значение
					    indx[ik]=p->key; // номер столбца	
					}
					else {
						printf(" Error non simmetric matrix ICCG. simplesparsetoCSIR_ITLSPD\n");
						system("pause");
					}
					pntr[k]=min(ik,pntr[k]);

					ik++;
				}

				p=p->next;
			}

		}

		for (k=0; k<nodes; k++) QuickSortCSIR(indx, val, pntr[k], pntr[k+1]-1);

		/* // debug
		for (k=0; k<=nodes; k++) {
		    printf("%d ",pntr[k]);
		}
        printf("\n");
		printf("nz==%d\n", nz);
		system("pause");
		//*/

	}
} // simplesparsetoCSIR_ITLSPD

/* Неполное LU разложение для несимметричных матриц
*  Пример А nxn=
*    9.0 0.0 0.0 3.0 1.0 0.0 1.0
*    0.0 11.0 2.0 1.0 0.0 0.0 2.0 
*    0.0 1.0 10.0 2.0 0.0 0.0 0.0 
*    2.0 1.0 2.0 9.0 1.0 0.0 0.0 
*    1.0 0.0 0.0 1.0 12.0 0.0 1.0 
*    0.0 0.0 0.0 0.0 0.0 8.0 0.0
*    2.0  2.0 0.0 0.0 3.0 0.0 8.0
*-----------------------------------------
*  инициализация (в этом виде данные поступают на вход процедуре)
*  память предполагается выделенной заранее :
*  верхняя треугольная матрица хранится построчно, в каждой строке
*  элементы отсортированы по убыванию номеров столбцов.
*  U_val :   1.0, 1.0, 3.0, 9.0,   2.0, 1.0, 2.0, 11.0,   2.0, 10.0, 1.0, 9.0, 1.0,12.0, 8.0, 8.0
*  U_ind :   6, 4, 3, 0,  6, 3, 2, 1,  3,2, 4,3, 6,4, 5, 6
*  U_ptr :   0, 4, 8, 10, 12, 14, 15, 16
*  нижняя треугольная матрица хранится постолбцово, в каждом столбце
*  элементы отсортированы по убыванию номеров строк.
*  L_val :  2.0, 1.0, 2.0, 9.0,    2.0, 1.0, 1.0, 11.0,  2.0, 10.0, 1.0, 9.0,  3.0, 12.0, 8.0, 8.0
*  L_ind :  6, 4, 3, 0,  6, 3, 2, 1,   3, 2,  4,3,  6, 4, 5, 6
*  L_ptr :  0, 4, 8, 10, 12, 14, 15, 16
*----------------------------------------------
*  Результат ILU разложения:
*  U_val : 1.0, 1.0, 3.0, 9.0, 2.0, 1.0, 2.0, 11.0, 2.0, 10.0, 1.0, 9.0, 1.0, 12.0, 8.0, 8.0.
*  L_val : 0.222, 0.111, 0.222, 1.0, -1.273, 0.091, 0.091, 1.0, 0.2, 1.0, 0.111, 1.0, -0.417, 1.0, 1.0, 1.0.
*/
void ILU0_Decomp_ITL(Real* &U_val, int* &U_ind, int* &U_ptr, Real* &L_val, int* &L_ind, int* &L_ptr, int n)
{
	/*
	// выделение памяти
	int n=7;
	//Real U_val[16] = { 3.0, 1.0, 1.0, 9.0,  2.0, 1.0, 2.0, 11.0, 2.0, 10.0, 1.0, 9.0, 1.0,12.0, 8.0, 8.0};
	//int U_ind[16] = { 3, 4, 6, 0,  2, 3, 6, 1,  3,2, 4,3, 6,4, 5, 6};
	//int U_ptr[8] = {0, 4, 8, 10, 12, 14, 15, 16};

	// Отсортированы в порядке убывания по столбцам.

	// verno
	Real U_val[16] = {  1.0, 1.0, 3.0, 9.0,     2.0, 1.0, 2.0, 11.0,   2.0, 10.0, 1.0, 9.0, 1.0,12.0, 8.0, 8.0};
	int U_ind[16] = { 6, 4, 3, 0,  6, 3, 2, 1,  3,2, 4,3, 6,4, 5, 6};
	int U_ptr[8] = {0, 4, 8, 10, 12, 14, 15, 16};

	//Real U_val[16] = {  9.0, 11.0, 2.0, 10.0, 3.0, 1.0, 2.0, 9.0, 1.0, 1.0, 12.0, 8.0, 1.0, 2.0, 1.0, 8.0 };
	//int U_ind[16] = { 0, 1, 1, 2, 0, 1, 2, 3, 0, 3, 4, 5, 0, 1, 4, 6};
	//int U_ptr[8] = {0, 1, 2, 4, 8, 11, 15, 16};

	//Real L_val[16] = {2.0, 1.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0,  3.0};
	//int L_ind[16] = { 3, 4, 6, 2, 3, 6, 3, 4, 6};
	//int L_ptr[8] = {0, 3, 6, 7, 8, 8, 8, 8};

	// verno
	Real L_val[16] = {2.0, 1.0, 2.0, 9.0,    2.0, 1.0, 1.0, 11.0,  2.0, 10.0, 1.0, 9.0,  3.0, 12.0, 8.0, 8.0};
	int L_ind[16] = { 6, 4, 3, 0,  6, 3, 2, 1,   3, 2,  4,3,  6, 4, 5, 6};
	int L_ptr[8] = {0, 4, 8, 10, 12, 14, 15, 16};

     //Real L_val[16] = {9.0, 11.0, 1.0, 10.0, 2.0, 1.0, 2.0, 9.0, 1.0, 1.0, 12.0, 8.0, 2.0, 2.0, 3.0, 8.0};
	//int L_ind[16] = { 0, 1, 1, 2, 0, 1, 2, 3, 0, 3, 4, 5, 0, 1, 4, 6};
	//int L_ptr[8] = {0, 1, 2, 4, 8, 11, 12, 16};
	*/

	// решение
	int i, j, qn, pn, rn; 
      for (i = 0; i < n - 1; i++) {
	     Real multiplier = U_val[U_ptr[i+1]-1];
    
	     for (j = L_ptr[i]; j < L_ptr[i+1]; j++)
	          L_val[j] /= multiplier;
    
	     for (j = U_ptr[i+1]; j < U_ptr[i+2]-1; j++) {
	         multiplier = U_val[j];
	         qn = j + 1;
	         rn = L_ptr[i+1];
	         for (pn = L_ptr[U_ind[j]]; L_ind[pn] <= i + 1 && pn < L_ptr[U_ind[j]+1]; pn++) {
	              while (U_ind[qn] < L_ind[pn] && qn < U_ptr[i+2]) qn++;

	              if (L_ind[pn] == U_ind[qn] && qn < U_ptr[i+2])
	                     U_val[qn] -= multiplier * L_val[pn];
	         }
	         for (; pn < L_ptr[U_ind[j]+1]; pn++) {
	             while (L_ind[rn] < L_ind[pn] && rn < L_ptr[i+2])  rn++;

	             if (L_ind[pn] == L_ind[rn] && rn < L_ptr[i+2])
	                    L_val[rn] -= multiplier * L_val[pn];
	         }
	      }
      }
	  L_val[L_ptr[n-1]]=1.0;

	  // сортировка по возрастанию
	  for (i = 0; i < n; i++) {
          QuickSortCSIR(U_ind, U_val, U_ptr[i], U_ptr[i+1]-1);
          QuickSortCSIR(L_ind, L_val, L_ptr[i], L_ptr[i+1]-1);
	  }

	/*
	printf("Uval : ");
	for (i=0; i<16; i++) printf("%.3f, ",U_val[i]);
	printf("\n\n Lval: ");
	for (i=0; i<16; i++) printf("%.3f, ",L_val[i]);
	system("pause");
	exit(0);
	*/
} // ILU0_Decomp_ITL

/* Метод бисопряжённых градиентов
* для возможно несимметричной матрицы А (val, col_ind, row_ptr).
* Запрограммировано по книжке Баландин, Шурина : "Методы
* решения СЛАУ большой размерности".
* dV - правая часть СЛАУ,
* x - начальное приближение к решению или NULL.
* n - размерность А nxn.
* Количество итераций ограничено 2000.
* Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
* Иногда метод расходится. Если выбрать другой вектор r_tilda, то 
* процесс может стать сходящимся. Ограничение на выбор вектора r_tilda:
* главное чтобы скалярное произведение Scal(r,r_tilda,n) != 0.0.
*/
Real *BiSoprGrad(IMatrix *xO, SIMPLESPARSE &M,  Real *dV, Real *x, int n){
	printf("\nBiConjugate Gradients Method...:\n");

	// Разреженная матрица СЛАУ
	// в CRS формате.
    Real *val;
    int* col_ind, *row_ptr;

	// преобразование из SIMPLESPARSE формата в CRS формат хранения.
	simplesparsetoCRS(M, val, col_ind, row_ptr, n);

	// ILU предобуславливатель:
    Real *U_val, *L_val;
	int  *U_ind, *U_ptr, *L_ind, *L_ptr;

	printf("Incoplete LU Decomposition begin...\n");
    convertIMatrixtoCSIR_ILU_ITL(xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
	ILU0_Decomp_ITL(U_val, U_ind, U_ptr, L_val, L_ind, L_ptr, n);
	printf("Incoplete LU Decomposition finish...\n");


	Real *r=new Real[n], *r_tilda=new Real[n];
	Real *p=new Real[n], *f=new Real[n], *p_tilda=new Real[n];
	Real nz; // невязка
	Real *ap=new Real[n], *vcopy=new Real[n];
	Real a,b,dold, dnew;

	int i; // счётчик цикла for
	int k=0; // номер итерации.

	// Начальное приближение:
    //X0==
	if (x==NULL) {
        x=new Real[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	Real e = dterminatedTResudual;

	MatrixCRSByVector(val,col_ind,row_ptr,x,ap,n);
	for (i=0; i<n; i++) {
		r[i]=dV[i]-ap[i];
		r_tilda[i]=r[i];
	}

	 // p==M^(-1)*r;
    for (i=0; i<n; i++) vcopy[i]=r[i];
    inverseL_ITL(vcopy, L_val, L_ind, L_ptr, p, n);
    for (i=0; i<n; i++) vcopy[i]=p[i];  
	inverseU_ITL(vcopy, U_val, U_ind, U_ptr, p, n);

    // p_tilda==M^(-T)*r_tilda;
	for (i=0; i<n; i++) vcopy[i]=r_tilda[i];
    inverseL_ITL(vcopy, U_val, U_ind, U_ptr, p_tilda, n);
    for (i=0; i<n; i++) vcopy[i]=p_tilda[i];  
	inverseU_ITL(vcopy, L_val, L_ind, L_ptr, p_tilda, n);
	   


	nz=NormaV(r,n); // начальное значение невязки

	for (i=0; i<n; i++) vcopy[i]=r[i];
    inverseL_ITL(vcopy, L_val, L_ind, L_ptr, f, n);
    for (i=0; i<n; i++) vcopy[i]=f[i];  
	inverseU_ITL(vcopy, U_val, U_ind, U_ptr, f, n);
	// f==M^(-1)*r;
	dold=Scal(f,r_tilda,n); 

    while ((nz>e) && (k<2000)) { 
		MatrixCRSByVector(val,col_ind,row_ptr,p,ap,n);

		a=dold/Scal(ap,p_tilda,n);
		for (i=0; i<n; i++) {
           x[i]+=a*p[i];
		   r[i]-=a*ap[i];
		}
		delete[] ap;
		ap=MatrixTransposeCRSByVector(val,col_ind,row_ptr,p_tilda,n);
        for (i=0; i<n; i++) {
			r_tilda[i]-=a*ap[i];
		}

        for (i=0; i<n; i++) vcopy[i]=r[i];
        inverseL_ITL(vcopy, L_val, L_ind, L_ptr, f, n);
        for (i=0; i<n; i++) vcopy[i]=f[i];  
	    inverseU_ITL(vcopy, U_val, U_ind, U_ptr, f, n);
	    // f==M^(-1)*r;
		dnew=Scal(f,r_tilda,n);
		b=dnew/dold;
		dold=dnew;
		// вычисление невязки.
        nz=NormaV(r,n);
		if (k%10==0) printf("iter residual\n");
		printf(" %d %e\n", k, nz);

		if ((fabs(b) < 1e-60) || (fabs(nz)>1e7)) {
			// метод Бисопряжённых градиентов иногда расходится.
			printf("\nBiCG divergence detected...\n");
            system("pause");
			exit(0); // выход из приложения.
			break; // выход из цикла while
		}

        for (i=0; i<n; i++) {
			p[i]=f[i]+b*p[i];
		}

		
		for (i=0; i<n; i++) vcopy[i]=r_tilda[i];
        inverseL_ITL(vcopy, U_val, U_ind, U_ptr, f, n);
        for (i=0; i<n; i++) vcopy[i]=f[i];  
	    inverseU_ITL(vcopy, L_val, L_ind, L_ptr, f, n);
	    // f==M^(-T)*r_tilda;
        for (i=0; i<n; i++) {
		    p_tilda[i]=f[i]+b*p_tilda[i];
		}

		
		k++; // переход к следующей итерации.
	}

	// Освобождение памяти
	delete[] r; delete[] r_tilda; 
	delete[] p; delete[] p_tilda;
	delete[] ap; delete[] f;

	return x;

} // BiSoprGrad

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
// Добавлен ILU предобуславливатель.
void SoloveichikAlg( IMatrix *xO, SIMPLESPARSE &M,// Разреженная матрица СЛАУ
                         Real *dV,  // вектор правой части
                         Real* &dX0, // вектор начального приближения
                         bool bconsole_message, // выводить ли значения невязки на консоль ?
						 int imaxiter) // максимально допустимое кол-во итераций
{
    

	int isize = xO->n;// размер квадратной матрицы
	 // Разреженная матрица СЛАУ
	 // в CRS формате.
     Real *val;
     int* col_ind, *row_ptr;

	 // преобразование из SIMPLESPARSE формата в CRS формат хранения.
	 simplesparsetoCRS(M, val, col_ind, row_ptr, isize);

	 // ILU предобуславливатель:
     Real *U_val, *L_val;
	 int  *U_ind, *U_ptr, *L_ind, *L_ptr;

	 if (DEBUG) printf("Incoplete LU Decomposition begin...\n");
     convertIMatrixtoCSIR_ILU_ITL(xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
	 ILU0_Decomp_ITL(U_val, U_ind, U_ptr, L_val, L_ind, L_ptr, isize);
	 if (DEBUG) printf("Incoplete LU Decomposition finish...\n");

	 /* // debug проверка ILU decomposition
	 for (i=0; i<U_ptr[isize]; i++) printf("%e ",U_val[i]);
	 printf("\n");
	 for (i=0; i<U_ptr[isize]; i++) printf("%d ",U_ind[i]);
	 printf("\n");
	 for (i=0; i<=isize; i++) printf("%d ",U_ptr[i]);
	 printf("\n");

     system("pause");

	 for (i=0; i<L_ptr[isize]; i++) printf("%e ",L_val[i]);
	 printf("\n");
	 for (i=0; i<L_ptr[isize]; i++) printf("%d ",L_ind[i]);
	 printf("\n");
	 for (i=0; i<=isize; i++) printf("%d ",L_ptr[i]);
	 printf("\n");

	 system("pause");
	 */


     
     Real *dx, *dax, *dr, *dz, *dp, *dar1, *dres, *f, *vcopy;
     Real dar, dbr, dnz, dscalp;
	 Real kend=imaxiter; // ограничение на максимальное число итераций
	 Real epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new Real[isize]; dax=new Real[isize]; dr= new Real[isize];
    dar1=new Real[isize]; vcopy=new Real[isize];dp= new Real[isize];
	dres=new Real[isize]; f=new Real[isize]; dz=new Real[isize];// вектор результата
   

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==NULL) {
	   dX0=new Real[isize];
#pragma omp parallel for
	   for (int i=0; i<isize; ++i) {
		   dx[i]=0.0;
		   dX0[i]=0.0;
	   }
   }
   else {
#pragma omp parallel for
	   for (int i=0; i<isize; ++i) dx[i]=dX0[i];
   }

   
   MatrixCRSByVector(val,col_ind,row_ptr,dx, dax, isize); // результат занесён в  dax
#pragma omp parallel for
   for (int i=0; i<isize; ++i) dr[i]= dV[i] - dax[i];  // начальная невязка
   // dr=L^(-1)*(dV-A*dx);
#pragma omp parallel for
   for (int i=0; i<isize; ++i) vcopy[i]=dr[i]; 
   inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dr, isize);
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   // dz=U^(-1)*dr;
#pragma omp parallel for
   for (int i=0; i<isize; ++i) vcopy[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   inverseU_ITL(vcopy, U_val, U_ind, U_ptr, dz, isize);
   // dp=L^(-1)*A*dz;
   MatrixCRSByVector(val,col_ind,row_ptr,dz,dp, isize);// результат занесён в dp
#pragma omp parallel for
   for (int i=0; i<isize; ++i) vcopy[i]=dp[i]; 
   inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dp, isize);

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      int k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (dnz > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
#pragma omp parallel for
         for (int i=0; i<isize; ++i)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         dnz=dnz-dar*dar/dscalp; // норма невязки
         
         if (DEBUG) if (bconsole_message)
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
            printf("%d %e \n",k,dnz);
		 } 
		 
         // f=U^(-1)*dr;
#pragma omp parallel for
         for (int i=0; i<isize; ++i) vcopy[i]=dr[i];  
         inverseU_ITL(vcopy, U_val, U_ind, U_ptr, f, isize);
#pragma omp parallel for
         for (int i=0; i<isize; ++i) vcopy[i]=f[i]; 
		 MatrixCRSByVector(val,col_ind,row_ptr,vcopy, dar1, isize);// результат занесён в dar1=A*U^(-1)*dr
#pragma omp parallel for
		 for (int i=0; i<isize; ++i) vcopy[i]=dar1[i]; 
		 // dar1=L^(-1)*A*U^(-1)*dr;
         inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dar1, isize);

         dbr=-Scal(dp,dar1,isize)*dscalp;
#pragma omp parallel for
         for (int i=0; i<isize; ++i)
		 {
            dz[i]=f[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }


         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e7) 
		 {
            // восстановление начального приближения
#pragma omp parallel for
            for (int i=0; i<isize; ++i) if (dX0==NULL) dx[i]=0.0; else dx[i]=dX0[i];
            printf("\n divergence Soloveichik solver \n");
            bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
      // возвращение результата
#pragma omp parallel for
      for (int i=0; i<isize; ++i) dres[i]=dx[i];
   }
   else
   {
      // возвращает начальное приближение
#pragma omp parallel for
	  for (int i=0; i<isize; ++i) dres[i]=dX0[i];
   }

   // освобождение памяти выделенной под динамические массивы
   delete[] dx; delete[] dax; delete[] dr; delete[] vcopy;
   delete[] dz; delete[] dp; delete[] dar1; delete[] f;
   delete[] U_val; delete[] U_ind; delete[] U_ptr;
   delete[] L_val; delete[] L_ind; delete[] L_ptr;
   delete[] val; delete[] col_ind; delete[] row_ptr;

#pragma omp parallel for
   for (int i=0; i<isize; ++i) {
		 dX0[i]=dres[i];
	   }
   delete[] dres; 

} // SoloveichikAlg

// Метод Ван Дер Ворста Bi-CGStab
// работает для возможно несимметричных вещественных матриц.
// Несимметричная матрица СЛАУ передаётся в CRS формате
// A (val, col_ind, row_ptr).
// Метод является комбинацией методов BiCG и GMRES(1). 
Real  *Bi_CGStab(int n, Real *val, int* col_ind, int* row_ptr, Real *dV, Real *dX0, int maxit)
{

	int iflag=1, icount=0;
	Real delta0, deltai;
	Real bet, roi;
	Real roim1=1.0, al=1.0, wi=1.0;
	Real *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	Real epsilon=dterminatedTResudual;  // точность вычисления
	int i;

	ri=new Real[n]; roc=new Real[n]; s=new Real[n]; t=new Real[n];
	vi=new Real[n]; pi=new Real[n]; dx=new Real[n]; dax=new Real[n];

	for (i=0; i<n; i++) {
		s[i]=0.0;
		t[i]=0.0;
		vi[i]=0.0;
		pi[i]=0.0;
	}

    // начальное приближение
    // X0 ==
    // под X0 понимается вектор поля температур к примеру.
    if (dX0==NULL) {
	   for (i=0; i<n; i++) dx[i]=0.0;
    }
    else {
	   for (i=0; i<n; i++) dx[i]=dX0[i];
    }

    MatrixCRSByVector(val,col_ind,row_ptr,dx,dax, n); // результат занесён в  dax
	for (i=0; i<n; i++) {
		ri[i]=dV[i]-dax[i];
		roc[i]=ri[i];
	}
	delta0=NormaV(ri,n);

	while ( iflag != 0 && icount < maxit) {

		icount++;

		roi=Scal(roc,ri,n);
		bet=(roi/roim1)*(al/wi);
		for (i=0; i<n; i++) {
			pi[i]=ri[i]+(pi[i]-vi[i]*wi)*bet;
		}
	
		MatrixCRSByVector(val,col_ind,row_ptr,pi,vi, n);
		al=roi/Scal(roc,vi,n);
        for (i=0; i<n; i++) {
			s[i]=ri[i]-al*vi[i];
		}
		
        MatrixCRSByVector(val,col_ind,row_ptr,s,t, n);
		wi=Scal(t,s,n)/Scal(t,t,n);
		for (i=0; i<n; i++) {
			dx[i]+=al*pi[i]+wi*s[i];
			ri[i]=s[i]-wi*t[i];
		}
		deltai=NormaV(ri,n);
		// печать невязки на консоль
		if (DEBUG) {
			if ((icount % 10) == 0)  printf("iter  residual\n");
			printf("%d %e \n", icount, deltai);
		}

		if (deltai <epsilon) iflag=0; // конец вычисления
		else roim1=roi;
	}

	return dx;

} // Bi_CGStab

#endif