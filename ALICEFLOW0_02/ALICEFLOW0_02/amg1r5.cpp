/* amg1r5new1.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

typedef long int integer;
typedef double doublereal;

bool yes_print_amg=true;

//#include "f2c.h"  // этот заголовочный файл просто ненужен.
#include <time.h>
//#include <math.h>

// недостающий функционал
// эта процедура определена выше по коду в mysolverv0_03.c, поэтому здесь её дефениция излишна.
/*integer min(integer ia, integer ib)
{
	integer ir;
	if (ia<ib) ir=ia;
	else ir=ib;
	return ir;
} // min
*/



/*integer max(integer ia, integer ib)
{
	integer ir;
	if (ia<ib) ir=ib;
	else ir=ia;
	return ir;
} // max
*/

/*doublereal max(doublereal da, doublereal db)
{
	doublereal dr;
	if (da<db) dr=db;
	else dr=da;
	return dr;
} // max
*/

doublereal d_lg10(doublereal * x) {
    // десятичный логарифм вещественого числа.
	return log10(*x);
}

doublereal pow_dd(doublereal * xa, doublereal * xb) {
	// возведение вещественного числа *xa в вещественную степень *xb.
	return pow(*xa,*xb);
}

integer pow_ii(integer *ia, integer *ib) {
	// Возведение целого числа ia в целую степень ib.
	// на основе вещественной функции класса <math.h>
	doublereal ra = (doublereal)(ia[0]);
	doublereal rb = (doublereal)(ib[0]);
	integer ii = (integer)(pow(ra, rb));
	return ii;
}

integer i_sign(integer *ia, integer *ib) {
	// Возвращает abs(ia)*s, где s=+1 если ib больше либо равно нулю.
	// -1 наоборот.
	integer s=-1;
	if (ib[0]>=0) {
		s=1;
	}
	return (abs(ia[0])*s);
}

/* Table of constant values */

static integer c__2 = 2;
static integer c__4 = 4;
static integer c__25 = 25;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__10 = 10;

/*     Last change:  ERB  22 Aug 2000   10:31 am */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5                                        MAIN SUBROUTINE */

/*     RELEASE 1.5, OCTOBER 1990 */

/*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

/* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
/*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
/*     FOLLOWING LINE IN SUBROUTINE RESC: */
/*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
/*     HAS BEEN CHANGED TO: */
/*     IW(IMAXW(KC-1)+1) = IAUX */

/* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
/*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
/*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
/*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
/*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
/*     IBLCK1 = IMINW(K). */

/* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
/*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
/*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
/*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
/*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
/*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
/*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

/* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
/*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
/*     TIVELY. */

/* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
/*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
/*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

/* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

/* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
/*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

/* ----------------------------------------------------------------------- */

/*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

/* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
/*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
/*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
/*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */

/*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
/*     DO 450 J=IA(I)+1,ICG(I)-1 */

/*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
/*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

/*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
/*     DO 550 J=IA(I)+1,ICG(I)-1 */

/*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
/*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

/* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
/*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */

/* ----------------------------------------------------------------------- */

/*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

/* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
/*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
/*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
/*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
/*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
/*     ENTRIES. */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* Subroutine */ int amg1r5_(doublereal *a, integer *ia, integer *ja, 
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig, 
	integer *nnu, integer *matrix, integer *iswtch, integer *iout, 
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc, 
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco, 
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2, 
	integer *nwt, integer *ntr, integer *ierr)
{
    /* Format strings */
   

    /* Builtin functions */
    

    /* Local variables */
    static integer mda, mdf, mdu;
    static doublereal res;
    static integer ium, iup;
    static doublereal res0;
    extern /* Subroutine */ int idec_(integer *, integer *, integer *, 
	    integer *);
    static integer mdia, mdja, mdig, iarr[25];
    //static real time[20];
	static unsigned int time[20];
    static integer imin[25], imax[25];
    static doublereal resi[25];
    static integer kout, ncyc0, irow0, ndicg, icgst, iminw[25], imaxw[25];
    extern /* Subroutine */ int first_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *), solve_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *), setup_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer ndigit, levels, kevelx, nstcol[25], kswtch;
    extern /* Subroutine */ int wrkcnt_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
  



/*         ----------------------------------------------- */
/*         | AMG-MODULE FOR SOLVING LINEAR SYSTEMS L*U=F | */
/*         ----------------------------------------------- */

/*         -------------------------------------------------------------- */

/*     ASSUMPTIONS ON L: */

/*         THE PROGRAM REQUIRES: */ // Требования к входным данным:

/*             - DIAGONAL ENTRIES ARE ALWAYS POSITIVE (ON ALL GRIDS); */
/*             - L IS A SQUARE MATRIX WHICH IS EITHER REGULAR OR SINGULAR */
/*               WITH ROWSUMS=0. */
	// Матрица L это квадратная матрица с всегда положительными диагональными элементами с нулевой суммой коэффициентов в каждой строке (положительная определённость).

/*         FOR THEORETICAL REASONS THE FOLLOWING SHOULD HOLD: */

/*             - L POSITIVE DEFINITE (OR SEMI-DEFINITE WITH ROWSUM=0) */
/*             - L "ESSENTIALLY" POSITIVE TYPE, I.E., */

/*                  -- DIAGONAL ENTRIES MUST BE > 0 ; */
/*                  -- MOST OF THE OFF-DIAGONAL ENTRIES <= 0 ; */
/*                  -- ROWSUMS SHOULD BE >= 0 . */

// Матрица L положительно определённая : диагональные элементы  строго больше нуля, большинство внедиагональных элементов <= 0.
// Сумма коэффициентов в строке больше либо равна нулю - диагональное преобладание.


/*     THE USER HAS TO PROVIDE THE MATRIX L, THE RIGHT HAND SIDE F AND */
/*     CERTAIN POINTER VECTORS IA AND JA. */

// Пользователь задаёт матрицу коэффициентов L, правую часть F, а также информацию о связях между элементами матрицы в специальных столбцах IA и JA.
// IA - позиция первого элемента в строке. JA - номер столбца для каждого элемента, первым записывается диагональный элемент.
// В фортране нумерация начинается с единицы.

/*         -------------------------------------------------------------- */

/*     STORAGE OF L: */ // Требования к хранению матрицы L.

/*         THE NON-ZERO ENTRIES OF THE MATRIX L ARE STORED IN */
/*         "COMPRESSED" SKY-LINE FASHION IN A 1-D VECTOR A, I.E., ROW */
/*         AFTER ROW, EACH ROW STARTING WITH ITS DIAGONAL ELEMENT. THE */
/*         OTHER NON-ZERO ROW ENTRIES FOLLOW THEIR DIAGONAL ENTRY IN ANY */
/*         ORDER. */


/*         IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO */
/*         PROVIDE TWO POINTER ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL */
/*         NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L */
/*         (1.LE.I.LE.NNU) ARE STORED IN A(J) WHERE THE RANGE OF J */
/*         IS GIVEN BY */

/*                     IA(I) .LE. J .LE. IA(I+1)-1. */

/*         THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF */
/*         ROW I WITHIN THE VECTOR A. IN PARTICULAR, */

/*                     IA(1) = 1 ,  IA(NNU+1) = 1 + NNA */

/*         WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED. */
/*         THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT */
/*         ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E., */
/*         JA(J) POINTS TO THE COLUMN INDEX OF A(J). */
/*         IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I */
/*         AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I. */

/*         IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS */
/*         (FOR ANY I WITH  1.LE.I.LE.NNU): */

/*                  F(I) =        SUM      A(J) * U(JA(J)) */
/*                           J1.LE.J.LE.J2 */

/*         WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND */
/*         SIDE AND */

/*                     J1 = IA(I) ,  J2 = IA(I+1)-1. */

/*         NOTES: THE ENTRY IA(NNU+1) HAS TO POINT TO THE FIRST FREE */
/*                ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, */
/*                AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW. */

/*                THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5. */
/*                SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT */
/*                BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED */
/*                THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE */
/*                SETUP PHASE WILL FAIL. */
/*                  ON THE OTHER HAND, RUNNING AMG A SECOND TIME ON THE */
/*                SAME INPUT DATA WITH ISWTCH=4 HAS NO SENSE, BECAUSE */
/*                THE RESULTS OF THE FIRST SETUP PHASE ARE STILL STORED */
/*                AND THUS THIS PHASE CAN BE SKIPPED IN A SECOND CALL. */
/*                IN ORDER TO DO THIS, SET ISWTCH TO 1, 2 OR 3. */

/* ----------------------------------------------------------------------- */

/*         THE FORM OF THE CALLING PROGRAM HAS TO BE AS FOLLOWS: */

/*               PROGRAM DRIVER */
/*         C */
/*               DOUBLE PRECISION A(#NDA),U(#NDU),F(#NDF) */
/*               INTEGER IA(#NDIA),JA(#NDJA),IG(#NDIG) */
/*         C */
/*               NDA  = #NDA */
/*               NDU  = #NDU */
/*               NDF  = #NDF */
/*               NDIA = #NDIA */
/*               NDJA = #NDJA */
/*               NDIG = #NDIG */
/*         C */
/*         C     SET UP A, F, IA, JA AND SPECIFY NECESSARY PARAMETERS */
/*         C */
/*               .... */
/*               .... */
/*         C */
/*               CALL AMG1R5(A,IA,JA,U,F,IG, */
/*        +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX, */
/*        +                  ISWTCH,IOUT,IPRINT, */
/*        +                LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU, */
/*        +                  ECG1,ECG2,EWT2,NWT,NTR, */
/*        +                  IERR) */
/*         C */
/*               .... */
/*               .... */
/*         C */
/*               STOP */
/*               END */

/* ----------------------------------------------------------------------- */

/*     INPUT VIA ARRAYS (SEE ABOVE): */

/*     A        -   MATRIX L */

/*     IA       -   POINTER VECTOR */

/*     JA       -   POINTER VECTOR */

/*     U        -   FIRST APPROXIMATION TO SOLUTION */

/*     F        -   RIGHT HAND SIDE */


/* ----------------------------------------------------------------------- */


/*     SCALAR INPUT PARAMETERS OF AMG1R5: */

/*     THE INPUT PARAMETERS OF AMG1R5 IN THE LIST BELOW ARE ARRANGED */
/*     ACCORDING TO THEIR IMPORTANCE TO THE GENERAL USER. THE PARAMETERS */
/*     PRECEEDED BY A * MUST BE SPECIFIED EXPLICITELY. ALL THE OTHER */
/*     PARAMETERS ARE SET TO STANDARD VALUES IF ZERO ON INPUT. */

/*     THERE ARE FOUR CLASSES OF INPUT PARAMETERS WITH DECREASING PRI- */
/*     ORITY: */

/*     1. PARAMETERS DESCRIBING THE USER-DEFINED PROBLEM AND DIMENSIONING */
/*        OF VECTORS IN THE CALLING PROGRAM */

/*     2. PARAMETERS SPECIFYING SOME GENERAL ALGORITHMIC ALTERNATIVES AND */
/*        THE AMOUNT OF OUTPUT DURING SOLUTION */

/*     3. PARAMETERS CONTROLLING THE MULTIGRID CYCLING DURING THE SOLU- */
/*        TION PHASE */

/*     4. PARAMETERS CONTROLLING THE CREATION OF COARSER GRIDS AND INTER- */
/*        POLATION FORMULAS. */

/*     ONLY THE CLASS 1 - PARAMETERS MUST BE SPECIFIED EXPLICITELY BY */
/*     THE USER. CLASS 2 - PARAMETERS CONTROL THE GENERAL PERFORMANCE OF */
/*     AMG1R5. CHANGING THEM DOESN'T REQUIRE UNDERSTANDING THE AMG - */
/*     ALGORITHM. SPECIFYING NON-STANDARD-VALUES FOR CLASS 3 - PARAMETERS */
/*     PRESUPPOSES A GENERAL KNOWLEDGE OF MULTIGRID METHODS, WHEREAS THE */
/*     FUNCTION OF CLASS 4 - PARAMETERS IS ONLY UNDERSTANDABLE AFTER */
/*     STUDYING THE AMG-ALGORITHM IN DETAIL. FORTUNATELY IN MOST CASES */
/*     THE CHOICE OF CLASS 3 AND 4 - PARAMETERS ISN'T CRITICAL AND USING */
/*     THE AMG1R5 - SUPPLIED STANDARD VALUES SHOULD GIVE SATISFACTORY */
/*     RESULTS. */

/*         -------------------------------------------------------------- */

/*     CLASS 1 - PARAMETERS: */

/*  *  NDA      -   DIMENSIONING OF VECTOR A IN CALLING PROGRAM */

/*  *  NDIA     -   DIMENSIONING OF VECTOR IA IN CALLING PROGRAM */

/*  *  NDJA     -   DIMENSIONING OF VECTOR JA IN CALLING PROGRAM */

/*  *  NDU      -   DIMENSIONING OF VECTOR U IN CALLING PROGRAM */

/*  *  NDF      -   DIMENSIONING OF VECTOR F IN CALLING PROGRAM */

/*  *  NDIG     -   DIMENSIONING OF VECTOR IG IN CALLING PROGRAM */

/*  *  NNU      -   NUMBER OF UNKNOWNS */

/*  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L. */

/*                  1ST DIGIT OF MATRIX  --  ISYM: */
/*                    =1: L IS SYMMETRIC; */
/*                    =2: L IS NOT SYMMETRIC. */

/*                  2ND DIGIT OF MATRIX  --  IROW0: */
/*                    =1: L HAS ROWSUM ZERO; */
/*                    =2: L DOES NOT HAVE ROWSUM ZERO. */

/*         -------------------------------------------------------------- */

/*     CLASS 2 - PARAMETERS: */

/*     ISWTCH   -   PARAMETER CONTROLLING WHICH MODULES OF AMG1R5 ARE TO */
/*                  BE USED. */
/*                    =1:   CALL FOR -----, -----, -----, WRKCNT. */
/*                    =2:   CALL FOR -----, -----, SOLVE, WRKCNT. */
/*                    =3:   CALL FOR -----, FIRST, SOLVE, WRKCNT. */
/*                    =4:   CALL FOR SETUP, FIRST, SOLVE, WRKCNT. */
/*                  SETUP DEFINES THE OPERATORS NEEDED IN THE SOLUTION */
/*                         PHASE. */
/*                  FIRST INITIALIZES THE SOLUTION VECTOR (SEE PARAMETER */
/*                         IFIRST). */
/*                  SOLVE COMPUTES THE SOLUTION BY AMG CYCLING (SEE */
/*                         PARAMETER NCYC). */
/*                  WRKCNT PROVIDES THE USER WITH INFORMATION ABOUT */
/*                         RESIDUALS, STORAGE REQUIREMENTS AND CP-TIMES */
/*                         (SEE PARAMETER IOUT). */
/*                  IF AMG1R5 IS CALLED THE FIRST TIME, ISWTCH HAS TO */
/*                  BE =4. INDEPENDENT OF ISWTCH, SINGLE MODULES CAN BE */
/*                  BYPASSED BY A PROPER CHOICE OF THE CORRESPONDING */
/*                  PARAMETER. */

/*     IOUT     -   PARAMETER CONTROLLING THE AMOUNT OF OUTPUT DURING */
/*                  SOLUTION PHASE: */

/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT: */
/*                    =0: NO OUTPUT (EXCEPT FOR MESSAGES) */
/*                    =1: RESIDUAL BEFORE AND AFTER SOLUTION PROCESS */
/*                    =2: ADD.: STATISTICS ON CP-TIMES AND STORAGE REQUI- */
/*                        REMENTS */
/*                    =3: ADD.: RESIDUAL AFTER EACH AMG-CYCLE */

/*     IPRINT   -   PARAMETER SPECIFYING THE FORTRAN UNIT NUMBERS FOR */
/*                  OUTPUT: */

/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO */

/*                  2ND AND 3RD DIGIT  --  IUP: UNIT NUMBER FOR RESULTS */

/*                  4TH AND 5TH DIGIT  --  IUM: UNIT NUMBER FOR MESSAGES */

/*         -------------------------------------------------------------- */

/*     CLASS 3 - PARAMETERS: */

/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
/*                        BY THE FOLLWING DIGITS. */

/*                  REST OF IFIRST  --  RNDU: */
/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
/*                    IFIRST=1372815) */

/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

/*                  1ST DIGIT OF NCYC  --  IGAM: */
/*                    =1: V -CYCLE, */
/*                    =2: V*-CYCLE, */
/*                    =3: F -CYCLE, */
/*                    =4: W -CYCLE. */
/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

/*                  2ND DIGIT OF NCYC  --  ICGR: */
/*                    =0: NO CONJUGATE GRADIENT, */
/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
/*                    =2: CONJUGATE GRADIENT (FULL CG). */

/*                  3RD DIGIT OF NCYC  --  ICONV: */
/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
/*                    (FINEST GRID): */
/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
/*                        NCYCLE (SEE BELOW) */
/*                    =2: STOP, IF  ||RES|| < EPS */
/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
/*                           EPS     (SEE INPUT PARAMETER EPS) */
/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
/*                           |U|   = SUPREMUM NORM OF SOLUTION */
/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
/*                    AFTER AT MOST NCYCLE CYCLES. */

/*                  REST OF NCYC  --  NCYCLE: */
/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
/*                    NCYCLE=0: NO CYCLING. */

/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
/*                  ARE PERFORMED, REGARDLESS OF EPS. */

/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
/*                  GRID IN CYCLING: */

/*                  1ST DIGIT OF MADAPT  --  MSEL: */
/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
/*                        PHASE ARE USED WITHOUT CHECK. */
/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
/*                        (SEE BELOW). */

/*                  REST OF MADAPT  --  FAC */
/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
/*                        BY DEFAULT. */


/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF NRD  --  NRDX: */
/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
/*                    =2: FULL GS SWEEP */
/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

/*                  1ST DIGIT  --  NSC: */
/*                    =1: GAUSS-SEIDEL METHOD */
/*                    =2: DIRECT SOLVER (YALE SMP) */

/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

/*         -------------------------------------------------------------- */

/*     CLASS 4 - PARAMETERS: */

/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
/*                        OPERATORS */
/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */


/* ----------------------------------------------------------------------- */

/*     OUTPUT: */

/*     U        -   CONTAINS THE COMPUTED SOLUTION */


/*     IERR     -   ERROR PARAMETER: */

/*                    >0: FATAL ERROR (ABNORMAL TERMINATION OF AMG1R5) */
/*                    <0: NON-FATAL ERROR (EXECUTION OF AMG1R5 CONTINUES) */

/*                  ERROR CODES IN DETAIL: */

/*                  1. DIMENSIONING TOO SMALL FOR VECTOR */
/*                        A      (IERR = 1) */
/*                        IA     (IERR = 2) */
/*                        JA     (IERR = 3) */
/*                        U      (IERR = 4) */
/*                        F      (IERR = 5) */
/*                        IG     (IERR = 6) */

/*                     NO YALE-SMP BECAUSE OF STORAGE (NDA TOO SMALL): */
/*                               (IERR = -1) */
/*                     NO YALE-SMP BECAUSE OF STORAGE (NDJA TOO SMALL): */
/*                               (IERR = -3) */
/*                     NO CG BECAUSE OF STORAGE (NDU TOO SMALL): */
/*                               (IERR = -4) */
/*                     NO SPACE FOR TRANSPOSE OF INTERPOLATION (NDA OR */
/*                                                     NDJA TOO SMALL): */
/*                               (IERR = -1) */

/*                  2. INPUT DATA ERRONEOUS: */

/*                     A-ENTRY MISSING, ISYM = 1:           (IERR = -11) */
/*                     PARAMETER MATRIX MAY BE ERRONEOUS:   (IERR = -12) */
/*                     DIAGONAL ELEMENT NOT STORED FIRST:   (IERR =  13) */
/*                     DIAGONAL ELEMENT NOT POSITIV:        (IERR =  14) */
/*                     POINTER IA ERRONEOUS:                (IERR =  15) */
/*                     POINTER JA ERRONEOUS:                (IERR =  16) */
/*                     PARAMETER ISWTCH ERRONEOUS:          (IERR =  17) */
/*                     PARAMETER LEVELX ERRONEOUS:          (IERR =  18) */

/*                  3. ERRORS OF THE AMG1R5-SYSTEM (SHOULD NOT OCCUR): */

/*                     TRANSPOSE A-ENTRY MISSING:           (IERR =  21) */
/*                     INTERPOLATION ENTRY MISSING:         (IERR =  22) */

/*                  4. ALGORITHMIC ERRORS: */

/*                     CG-CORRECTION NOT DEFINED:           (IERR =  31) */
/*                     NO YALE-SMP BECAUSE OF ERROR IN */
/*                     FACTORIZATION:                       (IERR = -32) */

/* ----------------------------------------------------------------------- */

/*     WORK SPACE: */

/*     THE INTEGER VECTOR IG HAS TO BE PASSED TO AMG1R5 AS WORK SPACE. */

/* ----------------------------------------------------------------------- */

/*     DIMENSIONING OF INPUT VECTORS AND WORK SPACE: */

/*     IT'S IMPOSSIBLE TO TELL IN ADVANCE THE EXACT STORAGE REQUIREMENTS */
/*     OF AMG. THUS, THE FOLLOWING FORMULAS GIVE ONLY REASONABLE GUESSES */
/*     FOR THE VECTOR LENGTHS WHICH HAVE TO BE DECLARED IN THE CALLING */
/*     PROGRAM. IN THESE FORMULAS NNA DENOTES THE NUMBER OF NON-ZERO */
/*     ENTRIES IN THE INPUT-MATRIX L AND NNU IS THE NUMBER OF UNKNOWNS. */

/*     VECTOR         NEEDED LENGTH (GUESS) */
/*       A               3*NNA + 5*NNU */
/*       JA              3*NNA + 5*NNU */
/*       IA              2.2*NNU */
/*       U               2.2*NNU */
/*       F               2.2*NNU */
/*       IG              5.4*NNU */

/* ----------------------------------------------------------------------- */


/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

/*          ISWTCH = 4 */
/*          IOUT   = 12 */
/*          IPRINT = 10606 */

/*          LEVELX = 25 */
/*          IFIRST = 13 */
/*          NCYC   = 10110 */
/*          EPS    = 1.D-12 */
/*          MADAPT = 27 */
/*          NRD    = 1131 */
/*          NSOLCO = 110 */
/*          NRU    = 1131 */

/*          ECG1   = 0. */
/*          ECG2   = 0.25 */
/*          EWT2   = 0.35 */
/*          NWT    = 2 */
/*          NTR    = 0 */

/*     IF ANY ONE OF THESE PARAMETERS IS 0 ON INPUT, ITS CORRESPONDING */
/*     STANDARD VALUE IS USED BY AMG1R5. */

/* ----------------------------------------------------------------------- */

/*     PORTABILITY RESTRICTIONS: */

/*     1. ROUTINE CTIME IS MACHINE DEPENDENT AND HAS TO BE ADAPTED TO */
/*        YOUR COMPUTER INSTALLATION OR REPLACED BY A DUMMY ROUTINE. */

/*     2. MOST INPUT PARAMETERS ARE COMPOSED OF SEVERAL DIGITS, THEIR */
/*        SIGNIFICANCE HAVING BEEN DESCRIBED ABOVE. BE SURE NOT TO ENTER */
/*        MORE DIGITS THAN YOUR COMPUTER CAN STORE ON AN INTEGER VARI- */
/*        ABLE. */

/*     3. APART FROM FORTRAN INTRINSIC FUNCTIONS AND SERVICE ROUTINES, */
/*        THERE IS ONLY ONE EXTERNAL REFERENCE TO A PROGRAM NOT CONTAINED */
/*        IN THE AMG1R5 - SYSTEM, I.E. THE LINEAR SYSTEM SOLVER NDRV OF */
/*        THE YALE SPARSE MATRIX PACKAGE. IF YOU HAVN'T ACCESS TO THIS */
/*        PACKAGE, ENTER A DUMMY ROUTINE NDRV AND AVOID CHOOSING NSC=2 */
/*        (SUBPARAMETER OF NSOLCO). THEN NDRV ISN'T CALLED BY AMG1R5. */
/*        IN THIS CASE, HOWEVER, INDEFINITE PROBLEMS WILL NOT BE SOLV- */
/*        ABLE. */
/*          THE YALE SPARSE MATRIX PACKAGE IS FREELY AVAILABLE FOR NON- */
/*        PROFIT PURPOSES. CONTACT THE DEPARTMENT OF COMPUTER SCIENCE, */
/*        YALE UNITVERSITY. */

/*     4. IN AMG1R5 THERE IS THE PARAMETER LRATIO, DENOTING THE RATIO */
/*        OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
/*        THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
/*        TO 2. CHANGE THIS VALUE IF NECESSARY. (THE SAME HAS TO BE */
/*        DONE WITH THE YALE SMP-ROUTINE NDRV.) */


/* ----------------------------------------------------------------------- */

/*     AUTHORS: */

/*          JOHN RUGE, FORT COLLINS (USA), */
/*              INSTITUTE FOR COMPUTATIONAL STUDIES AT CSU; */

/*          KLAUS STUEBEN, D-5205 ST. AUGUSTIN (W.-GERMANY), */
/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

/*          ROLF HEMPEL, D-5205 ST. AUGUSTIN (W.-GERMANY), */
/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

/* ----------------------------------------------------------------------- */


/* ===> LRATIO HAS TO BE SET TO THE NUMBER OF INTEGERS OCCUPYING THE SAME */
/*     AMOUNT OF STORAGE AS ONE DOUBLE PRECISION REAL. */


/* ===> MAXGR IS THE MAXIMAL NUMBER OF GRIDS. CHANGING THIS UPPER LIMIT */
/*     JUST REQUIRES CHANGING THE PARAMETER STATEMENT. */


    /* Parameter adjustments */
    --ig;
    --f;
    --u;
    --ja;
    --ia;
    --a;

    /* Function Body */
    *ierr = 0;

/* ===> SET PARAMETERS TO STANDARD VALUES, IF NECCESSARY */


    if (*iout != 0) {
	idec_(iout, &c__2, &ndigit, iarr);
	kout = iarr[1];
    } else {
	kout = 2;
    }

    if (*iswtch != 0) {
	kswtch = *iswtch;
    } else {
	kswtch = 4;
    }

    if (*levelx > 0) {
	kevelx = min(*levelx,25);
    } else if (*levelx < 0) {
	goto L70;
    } else {
	kevelx = 25;
    }

    if (*iprint != 0) {
	idec_(iprint, &c__4, &ndigit, iarr);
	iup = iarr[1] * 10 + iarr[2];
	ium = iarr[3];
    } else {
	iup = 6;
	ium = 6;
    }
    icgst = *nnu + 3;
    ndicg = (*ndig - icgst + 1) / 2;
    if (ndicg <= 0) {
	goto L60;
    }

    switch (kswtch) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }
    //io___10.ciunit = ium;
    //s_wsfe(&io___10);
    //e_wsfe();

	
	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER I.  SWTCH ***\n");
	
    *ierr = 17;
    return 0;

L40:
    setup_(nnu, matrix, &kevelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[
	    1], &ia[1], &ja[1], &ig[1], imin, imax, iminw, imaxw, &ig[icgst], 
	    &ig[icgst + ndicg], nstcol, iarr, time, &levels, &irow0, nda, 
	    ndia, ndja, ndu, ndf, &ndicg, &ium, &mda, &mdia, &mdja, &mdu, &
	    mdf, &mdig, &c__25, &c__2);
    if (*ierr > 0) {
	return 0;
    }
L30:
    first_(ifirst, &u[1], imin, imax, iarr, &irow0);
L20:
    solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
	    ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst], 
	    &ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels, 
	    nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
	    res0, &res);
    if (*ierr > 0) {
	return 0;
    }
L10:
    wrkcnt_(&kout, &ia[1], &ig[1], imin, imax, iminw, &levels, time, &ncyc0, &
	    iup, &mda, &mdia, &mdja, &mdu, &mdf, &mdig, &res0, &res);
    return 0;
L60:
    //io___29.ciunit = ium;
    //s_wsfe(&io___29);
    //e_wsfe();
	printf("*** ERROR IN AMG1R5: NDIG TOO SMALL ***\n");

    *ierr = 6;
    return 0;
L70:
    //io___30.ciunit = ium;
    //s_wsfe(&io___30);
    //e_wsfe();

	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER LEVELX ***\n");

    *ierr = 18;
    return 0;
} /* amg1r5_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     GENERAL AMG1R5 SETUP-SUBROUTINES */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     SETUP                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int setup_(integer *nnu, integer *matrix, integer *levelx, 
	doublereal *ecg1, doublereal *ecg2, doublereal *ewt2, integer *nwt, 
	integer *ntr, integer *ierr, doublereal *a, doublereal *u, integer *
	ia, integer *ja, integer *iw, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *icg, integer *ifg, integer *nstcol, 
	integer *iarr, /*real*/ unsigned int *time, integer *levels, integer *irow0, integer *
	nda, integer *ndia, integer *ndja, integer *ndu, integer *ndf, 
	integer *ndicg, integer *ium, integer *mda, integer *mdia, integer *
	mdja, integer *mdu, integer *mdf, integer *mdig, integer *maxgr, 
	integer *lratio)
{
    static integer i__;
    extern /* Subroutine */ int idec_(integer *, integer *, integer *, 
	    integer *);
    static integer isym;
    extern /* Subroutine */ int check_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/ unsigned int *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     crsng_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    /*real*/ unsigned int *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer ndigit;


/*     PREPARATION PHASE OF AMG1R5 (GENERAL PART) */


/* ===> DECOMPOSE "MATRIX" */

    /* Parameter adjustments */
    --time;
    --iarr;
    --nstcol;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --u;
    --a;

    /* Function Body */
    idec_(matrix, &c__2, &ndigit, &iarr[1]);
    isym = iarr[1];
    *irow0 = iarr[2];

/* ===> PREPARATION (IGNORED IN TIMING) */

    imin[1] = 1;
    imax[1] = *nnu;
    check_(ierr, &a[1], &ia[1], &ja[1], &imin[1], &imax[1], &icg[1], &ifg[1], 
	    &time[1], &isym, irow0, nda, ndicg, ndja, ium);
    if (*ierr > 0) {
	return 0;
    }

/* ===> RESET TIME COUNTERS */

    for (i__ = 1; i__ <= 20; ++i__) {
	//time[i__] = 0.f;
		time[i__] = 0;
/* L30: */
    }

/* ===> DEFINE COARSER GRIDS + OPERATORS. RESET LEVELS IF NECESSARY. */

    crsng_(levelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[1], &ia[1], &
	    ja[1], &iw[1], &imin[1], &imax[1], &iminw[1], &imaxw[1], &icg[1], 
	    &ifg[1], &nstcol[1], levels, irow0, nda, ndja, ndia, ndu, ndf, 
	    ndicg, &time[1], ium, mda, mdia, mdja, mdu, mdf, mdig, maxgr, 
	    lratio);
    return 0;
} /* setup_ */


/* ....................................................................... */

/*     CHECK                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int check_(integer *ierr, doublereal *a, integer *ia, 
	integer *ja, integer *imin, integer *imax, integer *icg, integer *ifg,
	 /*real*/ unsigned int *time, integer *isym, integer *irow0, integer *nda, integer *
	ndicg, integer *ndja, integer *ium)
{
    /* Format strings */
    //static char fmt_9000[] = "(\002 CHECK: A PROBABLY SYMMETRIC\002)";
    //static char fmt_9005[] = "(\002 CHECK: A PROBABLY NOT SYMMETRIC. MEASU"
	//    "RE:\002,d11.3)";
    //static char fmt_9010[] = "(\002 CHECK: A PROBABLY NOT POS. TYPE:\002,i6"
	  //  ",\002 OFF-DIAGONAL\002,\002 ELEMENTS POSITIVE\002)";
    //static char fmt_9020[] = "(\002 CHECK: A PROBABLY NOT POS. TYPE:\002,i6"
	  //  ",\002 ROWSUMS\002,\002 NEGATIVE\002)";
    //static char fmt_9025[] = "(\002 CHECK: A PROBABLY SINGULAR - ROWSUMS ARE"
	//    " ZERO\002)";
   // static char fmt_9030[] = "(\002 CHECK: A PROBABLY POSITIVE TYPE\002)";
    //static char fmt_9100[] = "(\002 --- WARNG IN CHECK: PARAM MATRIX MAY BE "
	//    "BAD ---\002)";
    //static char fmt_9600[] = "(\002 CHECK: MATRIX A WAS SYMMETRICALLY STORE"
	//    "D\002)";
    //static char fmt_9500[] = "(\002 CHECK: STORAGE OF A HAS BEEN SYMMETRIZED"
	//    " BY\002,\002 INTRODUCING\002,i6,\002 ZEROES\002)";
   // static char fmt_9220[] = "(\002 *** WARNG IN CHECK:\002,i5,\002 A-ENTRIE"
	//    "S MISSING ***\002/\002     MISSING TRANSPOSE CONNECTIONS WILL BE"
	 //   " FILLED IN\002)";
    //static char fmt_9200[] = "(\002 *** ERROR IN CHECK: NDIG TOO SMALL **"
	//    "*\002)";
    //static char fmt_9210[] = "(\002 *** ERROR IN CHECK: NDA TOO SMALL ***"
	  //  "\002)";
    //static char fmt_9230[] = "(\002 *** ERROR IN CHECK: NDJA TOO SMALL **"
	//    "*\002)";
    //static char fmt_9300[] = "(\002 *** ERROR IN CHECK: POINTER IA ERRONEOUS"
	//    " ***\002)";
    //static char fmt_9310[] = "(\002 *** ERROR IN CHECK: POINTER JA ERRONEOUS"
	//    " ***\002)";
    //static char fmt_9320[] = "(\002 *** ERROR IN CHECK: DIAGONAL IS NOT STOR"
	//    "ED FIRST ***\002)";
    //static char fmt_9330[] = "(\002 *** ERROR IN CHECK: DIAGONAL IS NON-POSI"
	//    "TIVE ***\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    //integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, i1, j1, j2, nna, new__, nnu;
    static doublereal deps;
    static integer jnew;
    static doublereal asym;
    static integer naneg, naoff, nazer, napos;
    extern /* Subroutine */ int trunc_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, 
	    integer *);
    static integer ishift;
    static doublereal anormm, anormp, rowsum;

    /* Fortran I/O blocks */
    //static cilist io___52 = { 0, 0, 0, fmt_9000, 0 };
    //static cilist io___53 = { 0, 0, 0, fmt_9005, 0 };
    //static cilist io___54 = { 0, 0, 0, fmt_9010, 0 };
    //static cilist io___55 = { 0, 0, 0, fmt_9020, 0 };
    //static cilist io___56 = { 0, 0, 0, fmt_9025, 0 };
    //static cilist io___57 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___58 = { 0, 0, 0, fmt_9100, 0 };
    //static cilist io___59 = { 0, 0, 0, fmt_9600, 0 };
    //static cilist io___62 = { 0, 0, 0, fmt_9500, 0 };
    //static cilist io___63 = { 0, 0, 0, fmt_9220, 0 };
    //static cilist io___64 = { 0, 0, 0, fmt_9200, 0 };
    //static cilist io___65 = { 0, 0, 0, fmt_9210, 0 };
    //static cilist io___66 = { 0, 0, 0, fmt_9230, 0 };
    //static cilist io___67 = { 0, 0, 0, fmt_9300, 0 };
    //static cilist io___68 = { 0, 0, 0, fmt_9310, 0 };
    //static cilist io___69 = { 0, 0, 0, fmt_9320, 0 };
    //static cilist io___70 = { 0, 0, 0, fmt_9330, 0 };



/*     CHECKS FOR SEVERAL PROPERTIES OF A, IA, JA. IN PARTICULAR, */
/*     CHECKS FOR SYMMETRIC STORAGE OF GIVEN MATRIX (I.E. L(I,J) IS */
/*     STORED IN A IFF L(J,I) IS STORED IN A). */

/*     IF STORAGE IS SYMMETRIC, PAIRS OF ZEROES ARE REMOVED (IF THERE */
/*     ARE ANY) AND PROGRAM EXECUTION CONTINUES. */

/*     IF, HOWEVER, STORAGE OF A IS NOT SYMMETRIC, IT IS SYMMETRIZED. */
/*     IF ISYM.EQ.1, THE MISSING TRANSPOSE CONNECTIONS ARE COPIED. */
/*     IF ISYM.NE.1, THE STORAGE OF THE MATRIX A IS SYMMETRIZED (BY */
/*     ADDING CERTAIN ZERO ELEMENTS). THEN PAIRS OF ZEROES ARE REMOVED */
/*     (IF THERE ARE ANY) AND PROGRAM EXECUTION CONTINUES. */

/*     ARRAYS USED FOR TEMPORARY STORAGE: ICG, IFG */


/* ===> CHECK IF POINTERS IA, JA ARE REASONABLE */

    /* Parameter adjustments */
    --time;
    --ifg;
    --icg;
    --imax;
    --imin;
    --ja;
    --ia;
    --a;

    /* Function Body */
    if (ia[1] != 1) {
	goto L3300;
    }
    nnu = imax[1];
    if (nnu >= *ndicg) {
	goto L3200;
    }
    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icg[i__] = 0;
/* L2: */
    }
    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ia[i__];
	j2 = ia[i__ + 1] - 1;
	if (j2 < j1 || j2 > *nda) {
	    goto L3300;
	}
	if (ja[j1] != i__) {
	    goto L3320;
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i1 = ja[j];
	    if (i1 < 1 || i1 > nnu || icg[i1] == 1) {
		goto L3310;
	    }
	    icg[i1] = 1;
/* L5: */
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    icg[ja[j]] = 0;
/* L7: */
	}
/* L10: */
    }
    nna = ia[nnu + 1] - 1;

/* ===> CHECK FOR PROPERTIES OF A. IN PARTICULAR, COUNT */
/* ===> MISSING STORAGE PLACES ("NEW"). RETURN IF NEW=0 */

    anormm = 0.;
    anormp = 0.;
    naoff = 0;
    napos = 0;
    naneg = 0;
    nazer = 0;
    new__ = 0;
    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icg[i__] = ia[i__ + 1] - ia[i__];
/* L100: */
    }

    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__ = a[ia[i__]];
	if (d__ <= 0.) {
	    goto L3330;
	}
	deps = d__ * 1e-12;
	rowsum = d__;
/* Computing 2nd power */
	d__1 = d__;
	anormp += d__1 * d__1 * 2.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    rowsum += a[j];
	    if (a[j] >= deps) {
		++naoff;
	    }
	    i1 = ja[j];
	    if (i1 > 0) {
		goto L140;
	    }
	    ja[j] = -i1;
	    goto L150;
L140:
	    if (i1 < i__) {
		goto L135;
	    }
	    i__3 = ia[i1 + 1] - 1;
	    for (j1 = ia[i1] + 1; j1 <= i__3; ++j1) {
		if (ja[j1] != i__) {
		    goto L130;
		}
		ja[j1] = -ja[j1];
/* Computing 2nd power */
		d__1 = a[j] - a[j1];
		anormm += d__1 * d__1;
/* Computing 2nd power */
		d__1 = a[j] + a[j1];
		anormp += d__1 * d__1;
		goto L150;
L130:
		;
	    }
L135:
	    ja[j] = -ja[j];
/* Computing 2nd power */
	    d__1 = a[j];
	    anormm += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = a[j];
	    anormp += d__1 * d__1;
	    ++new__;
	    ++icg[i1];
L150:
	    ;
	}
	if (rowsum > deps) {
	    ++napos;
	} else if (rowsum < -deps) {
	    ++naneg;
	} else {
	    ++nazer;
	}
/* L200: */
    }
    anormm = sqrt(anormm);
    anormp = sqrt(anormp);
    asym = anormm / anormp;
    if (asym <= 1e-12) {
	asym = 0.;
    }

/* ===> MESSAGES ON A */

    if (asym == 0.) {
	   if (yes_print_amg) {
	      printf("CHECK: A PROBABLY SYMMETRIC\n");
	   }
    }
    if (asym != 0.) {
		 if (yes_print_amg) {
	         printf("CHECK: A PROBABLY NOT SYMMETRIC. MEASURE %1.5f\n",asym);
		 }
    }
    if (naoff > 0) {
	 if (yes_print_amg) {
	     printf("CHECK: A PROBABLY NOT POS. TYPE: %d OFF-DIAGONAL ELEMENTS POSITIVE\n",naoff);
	 }
    }
    if (naneg > 0) {
	 if (yes_print_amg) {
	      printf("CHECK: A PROBABLY NOT POS. TYPE: %d ROWSUMS NEGATIVE\n",naneg);
	 }
    }
    if (nazer == nnu) {
	 if (yes_print_amg) {
	      printf("CHECK: A PROBABLY SINGULAR - ROWSUMS ARE ZERO\n");
	 }
    }
    if (naoff == 0 && naneg == 0) {
	 if (yes_print_amg) {
	     printf("CHECK: A PROBABLY POSITIVE TYPE\n");
	 }
    }

/* ===> WARNINGS */

    if (*isym == 1 && asym != 0. || *isym > 1 && asym == 0. || *irow0 == 1 && 
	    nazer != nnu || *irow0 > 1 && nazer == nnu) {
	 if (yes_print_amg) {
         printf("--- WARNG IN CHECK: PARAM MATRIX MAY BE BAD ---");
	 }
	*ierr = -12;
    }

    if (new__ > 0) {
	goto L220;
    }
    if (yes_print_amg) {
	     printf("CHECK: MATRIX A WAS SYMMETRICALLY STORED\n");
	}
    goto L600;

/* ===> REPLACE A BY SYMMETRIZED VERSION */

L220:
    if (nna + new__ >= *nda) {
	goto L3210;
    }
    if (nna + new__ >= *ndja) {
	goto L3220;
    }

/* ===> EXTEND MATRIX A IN SITU */

    ifg[1] = 1;
    i__1 = nnu + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ifg[i__] = ifg[i__ - 1] + icg[i__ - 1];
/* L230: */
    }
    for (i__ = nnu; i__ >= 1; --i__) {
	ishift = ifg[i__] - ia[i__];
	i__1 = ia[i__];
	for (j = ia[i__ + 1] - 1; j >= i__1; --j) {
	    a[j + ishift] = a[j];
	    ja[j + ishift] = ja[j];
/* L240: */
	}
	icg[i__] = ia[i__ + 1] + ishift;
	ia[i__ + 1] = ifg[i__ + 1];
/* L250: */
    }

/* ===> SYMMETRIZE MATRIX A:  ISYM=1: COPY MISSING TRANSPOSE ENTRIES */
/*                           ISYM=2: FILL IN ZEROES */

    if (*isym != 1) {
	i__1 = nnu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j2 = icg[i__] - 1;
	    i__2 = j2;
	    for (j = ia[i__] + 1; j <= i__2; ++j) {
		i1 = ja[j];
		if (i1 > 0) {
		    goto L450;
		}
		ja[j] = -i1;
		jnew = icg[ja[j]];
		a[jnew] = 0.;
		ja[jnew] = i__;
		icg[ja[j]] = jnew + 1;
L450:
		;
	    }
/* L400: */
	}
	 if (yes_print_amg) {
	     printf("CHECK: STORAGE OF A HAS BEEN \n");
	     printf("SYMMETRIZED BY INTRODUCING %d ZEROES\n", new__);
	 }
    } else {
	 if (yes_print_amg) {
	    printf("*** WARNG IN CHECK: %d A-ENTRIES MISSING ***\n",new__);
	    printf("MISSING TRANSPOSE CONNECTIONS WILL BE FILLED IN\n");
	 }
	*ierr = -11;
	i__1 = nnu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j2 = icg[i__] - 1;
	    i__2 = j2;
	    for (j = ia[i__] + 1; j <= i__2; ++j) {
		i1 = ja[j];
		if (i1 > 0) {
		    goto L550;
		}
		ja[j] = -i1;
		jnew = icg[ja[j]];
		a[jnew] = 0.;
		ja[jnew] = i__;
		icg[ja[j]] = jnew + 1;
L550:
		;
	    }
/* L500: */
	}
    }

/* ===> REMOVE PAIRS OF ZEROES */

L600:
    trunc_(&c__1, &c__0, &a[1], &ia[1], &ja[1], &imin[1], &imax[1], &time[1], 
	    ierr, ium);
    return 0;

/* ===> ERROR MESSAGES */

L3200:
    //io___64.ciunit = *ium;
    //s_wsfe(&io___64);
    //e_wsfe();
	printf("*** ERROR IN CHECK: NDIG TOO SMALL ***\n");
    *ierr = 6;
    return 0;

L3210:
    //io___65.ciunit = *ium;
    //s_wsfe(&io___65);
    //e_wsfe();
	printf("*** ERROR IN CHECK: NDA TOO SMALL ***\n");
    *ierr = 1;
    return 0;

L3220:
    //io___66.ciunit = *ium;
    //s_wsfe(&io___66);
    //e_wsfe();
	printf("*** ERROR IN CHECK: NDJA TOO SMALL ***\n");
    *ierr = 3;
    return 0;

L3300:
    //io___67.ciunit = *ium;
    //s_wsfe(&io___67);
    //e_wsfe();
	printf("*** ERROR IN CHECK: POINTER IA ERRONEOUS***\n");
    *ierr = 15;
    return 0;

L3310:
    //io___68.ciunit = *ium;
    //s_wsfe(&io___68);
    //e_wsfe();
	printf("*** ERROR IN CHECK: POINTER JA ERRONEOUS***\n");
    *ierr = 16;
    return 0;

L3320:
    //io___69.ciunit = *ium;
    //s_wsfe(&io___69);
    //e_wsfe();
	printf("*** ERROR IN CHECK: DIAGONAL IS NOT STORED FIRST ***\n");
    *ierr = 13;
    return 0;

L3330:
    //io___70.ciunit = *ium;
    //s_wsfe(&io___70);
    //e_wsfe();
	printf("*** ERROR IN CHECK: DIAGONAL IS NON-POSITIVE ***\n");
    *ierr = 14;
    return 0;


} /* check_ */


/* ....................................................................... */

/*     TRUNC                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int trunc_(integer *k, integer *ntr, doublereal *a, integer *
	ia, integer *ja, integer *imin, integer *imax, /*real*/unsigned int *time, integer *
	ierr, integer *ium)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 *** ERROR IN TRUNC: TRANSPOSE A-ENTRY MI"
	//    "SSING ON GRID\002,i3,\002 ***\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
   // integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, i1, j1, j2, j3;
    static doublereal at;
    static integer jt, nna, imn, imx;
	static unsigned int told;
    static integer jpos;
	static unsigned int tnew;

    /* Fortran I/O blocks */
    //static cilist io___85 = { 0, 0, 0, fmt_9000, 0 };



/*     TRUNCATES OPERATOR ON GRID K CORRESPONDING TO THE VALUE OF NTR. */
/*     NTR HAS TO BE 0 OR 1: */

/*       =0:    PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID OPERATORS; */
/*       =1:    COARSE GRID OPERATORS REMAIN UNCHANGED. */



    /* Parameter adjustments */
    --time;
    --imax;
    --imin;
    --ja;
    --ia;
    --a;

    /* Function Body */
    if (*ntr == 1) {
	return 0;
    }

	told=clock();
    imn = imin[*k];
    imx = imax[*k];
    nna = ia[imx + 1] - ia[imn];
    jpos = ia[imn];

    i__1 = imx;
    for (i__ = imn; i__ <= i__1; ++i__) {
	j1 = ia[i__] + 1;
	j2 = ia[i__ + 1] - 1;
	a[jpos] = a[ia[i__]];
	ja[jpos] = i__;
	ia[i__] = jpos;
	++jpos;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i1 = ja[j];
	    if (i1 < 0) {
		goto L250;
	    }
	    if (i1 < i__) {
		goto L230;
	    }
	    if (a[j] != 0.) {
		goto L230;
	    }
	    i__3 = ia[i1 + 1] - 1;
	    for (j3 = ia[i1]; j3 <= i__3; ++j3) {
		if (ja[j3] != i__) {
		    goto L210;
		}
		jt = j3;
		at = a[j3];
		goto L215;
L210:
		;
	    }
	    goto L1000;
L215:
	    if (at != 0.) {
		goto L230;
	    }
	    ja[jt] = -ja[jt];
	    goto L250;
L230:
	    a[jpos] = a[j];
	    ja[jpos] = ja[j];
	    ++jpos;
L250:
	    ;
	}
/* L205: */
    }
    ia[imx + 1] = jpos;

/* ===> EXIT */

	tnew=clock();
    time[7] = time[7] + tnew - told;
    return 0;

/* ===> ERROR MESSAGE */

L1000:
    //io___85.ciunit = *ium;
    //s_wsfe(&io___85);
    //do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
    //e_wsfe();
	printf("*** ERROR IN TRUNC: TRANSPOSE A-ENTRY MISSING ON GRID %d ***\n",k[0]);
    *ierr = 21;
    return 0;

} /* trunc_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5 SETUP ROUTINES (VERSION 3) */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     CRSNG                                           SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int crsng_(integer *levelx, doublereal *eecg1, doublereal *
	eecg2, doublereal *eewt2, integer *nnwt, integer *ntr, integer *ierr, 
	doublereal *a, doublereal *u, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, integer *imaxw, integer 
	*icg, integer *ifg, integer *nstcol, integer *levels, integer *irow0, 
	integer *nda, integer *ndja, integer *ndia, integer *ndu, integer *
	ndf, integer *ndicg, /*real*/ unsigned int *time, integer *ium, integer *mda, integer *
	mdia, integer *mdja, integer *mdu, integer *mdf, integer *mdig, 
	integer *maxgr, integer *lratio)
{
    /* Format strings */
    //static char fmt_9000[] = "(/\002 **************** SPACE REQUIREMENTS ***"
	//    "*************\002//\002 VECTOR          NEEDED                  "
	//    "            \002/\002 ------------------------------------------"
	//    "----------\002/\002    A \002,i16,\002   ADJUST THE DIMENSIONING"
	//    " OF \002/\002    JA\002,i16,\002   VECTORS A - IG IN THE     "
	//    " \002/\002    IA\002,i16,\002   CALLING PROGRAM ACCORDING  \002"
	//    "/\002    U \002,i16,\002   TO THE CALCULATED SPACE RE-\002/\002 "
	//    "   F \002,i16,\002   QUIREMENTS AND RERUN THE   \002/\002    I"
	//    "G\002,i16,\002   PROGRAM.                   \002/\002 ----------"
	//    "------------------------------------------\002)";
   // static char fmt_9010[] = "(/\002 NOTE: IF YOU WANT TO USE CG-CORRECTIONS"
   //	    " IN THE SOLU-\002/\002       TION PROCESS (NCYC-SUBPARAMETER ICG"
   //	    "R=1 OR =2),\002/\002       PROVIDE FOR ADDITIONAL\002,i6,\002 ST"
   //	    "ORAGE LOCATIONS\002/\002       IN VECTORS U AND F.              "
   //	    "             \002/\002         SIMILARLY, USAGE OF THE YALE-SMP "
   //	    "SOLVER ON  \002/\002       THE COARSEST GRID (NSOLCO=2) WILL REQ"
   //	    "UIRE     \002/\002       ADDITIONAL SPACE IN VECTOR A DURING THE"
   //	    " SOLU- \002/\002       TION PHASE. IN THIS CASE, HOWEVER, ITS EX"
   //	    "ACT  \002/\002       AMOUNT ISN'T PREDICTABLE.                  "
   //	    "  \002/)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
   // integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, k, nwt;
    static doublereal ecg1, ecg2, ewt2;
    static integer ichk, iias, isia, mdir;
    extern /* Subroutine */ int pcol_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal /*integer*/ *, doublereal /*integer*/ *, /*real*/ unsigned int *, integer *, integer *,
	     integer *, integer *, integer *, integer *);
    static integer mdiw, mmax, kerr, ndir, iirs, irst, mdicg, iajas, isaja;
    extern /* Subroutine */ int opdfn_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal /*integer*/ *, 
	    integer *, integer *, /*real*/unsigned int *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer mdjtr, ndjtr, ncolx;
    extern /* Subroutine */ int trunc_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, /*real*/unsigned int *, integer *, 
	    integer *), pwint_(integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal /*integer*/ *, /*real*/ unsigned int 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *);
    static integer jtrst;
    extern /* Subroutine */ int rwsrt_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal /*integer*/ *, /*real*/ unsigned int *, integer 
	    *, integer *, integer *), setifg_(integer *, integer *, integer *,
	     integer *, integer *, integer *, /*real*/ unsigned int *);
    static integer kfirst;

    /* Fortran I/O blocks */
   // static cilist io___110 = { 0, 0, 0, fmt_9000, 0 };
   // static cilist io___111 = { 0, 0, 0, fmt_9010, 0 };



/*     PERFORMS COARSENING, DEFINES INTERPOLATIONS AND CGC-OPERATORS */

/*     ============== STANDARD VALUES OF PARAMETERS ===================== */

/*             ECG1=0.,    ECG2=0.25,    EWT2=0.35,   NWT=2 */

/*     ================ DESCRIPTION OF PARAMETERS ======================= */

/*     ECG1 --    DEFINES CRITERION FOR DETERMINING DIAGONAL */
/*                DOMINANCE IN RWSRT. I.E., IF THE ABSOLUTE */
/*                VALUE OF THE SUM OF THE OFF-DIAGONALS OF ROW */
/*                I IS SMALLER THAN ECG1 TIMES THE ABSOLUTE */
/*                VALUE OF THE DIAGONAL ENTRY, THEN POINT I IS */
/*                IMMEDIATELY FORCED TO BE AN F-POINT IN THE */
/*                PRE-COLORING ALGORITHM PCOL. */
/*                IN THE SECOND PART (WINT), NO INTERPOLATION */
/*                IS DEFINED FOR POINT I, AND NO POINTS USE */
/*                I FOR INTERPOLATION. IN ADDITION, THE WEIGHT */
/*                FOR POINT I IS NOT DISTRIBUTED TO OTHER POINTS */
/*                WHEN DEFINING THE INTERPOLATION WEIGHTS FOR */
/*                POINTS WHICH DEPEND ON POINT I. (THIS IS */
/*                EQUIVALENT TO COMPLETELY IGNORING SUCH POINTS */
/*                IN DETERMINING THE COARSE GRID AND INTERPOLATION */
/*                WEIGHTS. */

/*     ECG2 --    DEFINES STRONG CONNECTIONS (ALPHA IN THE PAPER) */

/*     EWT2 --    DEFINES STRONG DEPENDENCE ON A SET (BETA IN THE PAPER) */

/*     NWT  --    PARAMETER CONTROLLING THE DEFINITION OF INTERPOLATION */
/*                FORMULAS: */
/*                  =1 - CHECKING OF INT-FORMULA: OFF */
/*                  =2 - CHECKING OF INT-FORMULA: ON */


/* ===> ASSIGN DEFAULT VALUES IF ZERO */

    /* Parameter adjustments */
    --time;
    --nstcol;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --u;
    --a;

    /* Function Body */
    nwt = *nnwt;
    ecg1 = *eecg1;
    ecg2 = *eecg2;
    ewt2 = *eewt2;
    if (nwt == 0) {
	nwt = 2;
    }
    if (ecg2 == 0.) {
	ecg2 = .25;
    }
    if (ewt2 == 0.) {
	ewt2 = .35;
    }

/* ===> DECODE PARAMETER NWT */

/* L5: */
    ichk = nwt;

/* ===> COARSENING */

/* L10: */
    *levels = min(*levelx,*maxgr);
    mmax = *maxgr;

/* ===> INITIALIZE PARAMETERS MDA - MDIW, LATER SET TO ACTUAL STORAGE */
/*     REQUIREMENTS OF CORRESPONDING VECTORS */

    *mda = 0;
    *mdja = 0;
    mdicg = 0;
    mdir = 0;
    mdiw = imax[1];

/* ===> KFIRST IS THE NUMBER OF THE FIRST STORED GRID, IAJAS, IIAS AND */
/*     IIRS ARE THE SHIFTS IN VECTORS A, JA, IA AND IR, RESPECTIVELY, AS */
/*     COMPARED TO FULL STORAGE OF ALL GRIDS */

    kerr = 0;
    kfirst = 1;
    iajas = 0;
    iias = 0;
    iirs = 0;

    i__1 = *levels;
    for (k = 2; k <= i__1; ++k) {

/* ===> JTRST: INITIAL POINTER FOR WORK SPACE IN VECTOR A, TO CONTAIN */
/*            THE STRONG TRANSPOSE CONNECTIONS */
/*     NDJTR: AVAILABLE WORK SPACE */

L20:
	jtrst = ia[imax[k - 1] + 1];
	ndjtr = *lratio * (*nda - jtrst + 1);
	i__2 = k - 1;

	

	rwsrt_(&i__2, &ecg1, &ecg2, ierr, &a[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &imaxw[1], &ifg[1],  &a[jtrst], &time[1], &
		ndjtr, ium, &mdjtr);

	

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}
/* Computing MAX */
	i__2 = *mda, i__3 = jtrst + (mdjtr - 1) / *lratio + iajas;
	*mda = max(i__2,i__3);

/* ===>   IRST: INITIAL POINTER FOR WORK SPACE IN VECTOR U, TO CONTAIN */
/*             RESET STACK */
/*       NDIR: AVAILABLE WORK SPACE */

	irst = mdiw + 1 - iirs;
	ndir = *lratio * (*ndu - irst + 1);
	i__2 = k - 1;

	

	

	pcol_(&i__2, ierr, &ia[1], &ja[1], &iw[1], &imin[1], &imax[1], &imaxw[
		1], &icg[1], &ifg[1],  &u[irst] , &a[jtrst], &time[1], ndicg, &
		ndir, ium, &mdicg, &mdir, &iias);

   

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}
	i__2 = k - 1;

	

	pwint_(&i__2, &ewt2, &ichk, &mmax, &a[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &iminw[1], &imaxw[1], &ifg[1], &icg[1], &u[irst],
		&time[1], ierr, irow0, &ncolx, nda, ndja, ndicg, ndia, 
		ium, mda, mdja, &iajas);

		

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}

	

	opdfn_(&k, ierr, &mmax, &a[1], &ia[1], &ja[1], &iw[1], &imin[1], &
		imax[1], &iminw[1], &imaxw[1], &icg[1], &ifg[1],  &u[irst], &
		nstcol[1], &ncolx, &time[1], nda, ndja, ium, mda, mdja, &
		iajas);

	

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}
	if (k > mmax) {
	    goto L100;
	}
	trunc_(&k, ntr, &a[1], &ia[1], &ja[1], &imin[1], &imax[1], &time[1], 
		ierr, ium);
	if (*ierr > 0) {
	    return 0;
	}
	iw[iminw[k]] = ia[imax[k] + 1];
	if (k >= mmax) {
	    goto L100;
	}
	goto L50;
L30:
	if (k <= kfirst + 1 && iirs != 0) {
	    return 0;
	}
	kerr = *ierr;
	*ierr = 0;
	kfirst = k - 1;
	iirs = mdiw;
	isia = imin[k - 1] - 1;
	iias += isia;
	isaja = ia[imin[k - 1]] - 1;
	iajas += isaja;
	i__2 = ia[imax[k - 1] + 1] - 1;
	for (i__ = ia[imin[k - 1]]; i__ <= i__2; ++i__) {
	    ja[i__ - isaja] = ja[i__] - isia;
	    a[i__ - isaja] = a[i__];
/* L35: */
	}
	i__2 = imax[k - 1] + 1;
	for (i__ = imin[k - 1]; i__ <= i__2; ++i__) {
	    ia[i__ - isia] = ia[i__] - isaja;
/* L40: */
	}
	imin[k - 1] = 1;
	imax[k - 1] -= isia;
	iw[iminw[k - 1]] = ia[imax[k - 1] + 1];
	goto L20;
L50:
	;
    }
    goto L200;
L100:
    *levels = mmax;
L200:
    *mdia = imax[*levels] + 1 + iias;
/* Computing MAX */
    i__1 = imax[*levels] + iias, i__2 = mdiw + 1 + (mdir - 1) / *lratio, i__1 
	    = max(i__1,i__2), i__2 = mdiw + 1 + (mdiw - 1) / *lratio;
    *mdu = max(i__1,i__2);
    *mdf = imax[*levels] + iias;
    *mdig = mdiw + 2 + (mdicg << 1);
    if (kerr != 0 || *mdu > *ndu || *mdf > *ndf) {
	//io___110.ciunit = *ium;
	//s_wsfe(&io___110);
	//do_fio(&c__1, (char *)&(*mda), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdja), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdia), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdu), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdf), (ftnlen)sizeof(integer));
	//do_fio(&c__1,(char *)&(*mdig) , (ftnlen)sizeof(integer));
	//e_wsfe();

		if (yes_print_amg) 
		{
	         printf("( **************** SPACE REQUIREMENTS ***\n");
	         printf("************* VECTOR          NEEDED                  \n");
	         printf("             -----------------------------------------\n");
	         printf("----------    A ,%d,   ADJUST THE DIMENSIONING \n",mda[0]);
	         printf(" OF     JA ,%d,   VECTORS A - IG IN THE     \n",mdja[0]);
	         printf("     IA,%d,   CALLING PROGRAM ACCORDING  \n",mdia[0]);
	         printf("    U ,%d,   TO THE CALCULATED SPACE  \n",mdu[0]);
	         printf("  REF ,%d,   QUIREMENTS AND RERUN THE      \n",mdf[0]);
	         printf("IG %d ,   PROGRAM.                    ----------\n",mdig[0]);
	         printf("------------------------------------------)\n");
		}
	
	//io___111.ciunit = *ium;
	//s_wsfe(&io___111);
	//do_fio(&c__1, (char *)&mdiw, (ftnlen)sizeof(integer));
	//e_wsfe();
		if (yes_print_amg) {
	        printf("NOTE: IF YOU WANT TO USE CG-CORRECTIONS\n");
        	printf(" IN THE SOLUTION PROCESS (NCYC-SUBPARAMETER ICG\n");
	        printf("R=1 OR =2),  PROVIDE FOR ADDITIONAL %d ST\n",mdiw);
	        printf("ORAGE LOCATIONS IN VECTORS U AND F. \n");
	        printf("SIMILARLY, USAGE OF THE YALE-SMP \n");
	        printf("SOLVER ON  THE COARSEST GRID (NSOLCO=2) WILL REQUIRE\n");
	        printf(" ADDITIONAL SPACE IN VECTOR A DURING THE\n");
        	printf(" SOLUTION PHASE. IN THIS CASE, HOWEVER, ITS \n");
	        printf("EXACT  AMOUNT ISN'T PREDICTABLE.   \n");
		}    

	*ierr = kerr;
	return 0;
    }
    setifg_(&imin[1], &imax[1], &icg[1], &ifg[1], &nstcol[1], levels, &time[1]
	    );
    return 0;
} /* crsng_ */


/* ....................................................................... */

/*     RWSRT                                             SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int rwsrt_(integer *k, doublereal *ecg1, doublereal *ecg2, 
	integer *ierr, doublereal *a, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *imaxw, integer *ifg, /*integer*/ doublereal *
	jtr, /*real*/ unsigned int *time, integer *ndjtr, integer *ium, integer *mdjtr)
{
    /* Format strings */
    //static char fmt_9910[] = "(\002 *** ERROR IN RWSRT: NDA TOO SMALL ***"
	//    "\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
   // integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, ii;
    static doublereal rs;
    static integer ihi, jhi;
    static doublereal amn;
    static integer ilo, jlo;
    static doublereal amx, ast;
    static integer imx, jmx, iws;
	static unsigned int told;
    static doublereal atmp;
    static integer itmp;
	static unsigned int tnew;

    /* Fortran I/O blocks */
   // static cilist io___130 = { 0, 0, 0, fmt_9910, 0 };



/*     ROW-SORT ALGORITHM FOR ROWS OF A(K). IN DETAIL: */

/*     - SORTS THE ELEMENTS OF EACH ROW OF A(K) SUCH THAT THE STRONG */
/*       CONNECTIONS JUST FOLLOW THE DIAGONAL ELEMENT. ONE OF THE STRONG- */
/*       EST IS ALWAYS FIRST (EVEN IF PARTIAL SORTING IS PERFORMED!). */
/*       NON-NEGATIVE CONNECTIONS ARE ALWAYS DEFINED TO BE WEAK. */
/*       JA(IA(I)) IS RE-DEFINED TO POINT TO THE LAST STRONG CONNECTION */
/*       OF POINT I (OR TO IA(I) IF THERE IS NO STRONG CONNECTION). */

/*     - THE STRONG TRANSPOSE CONNECTIONS ARE LOADED LOGICALLY */
/*       INTO JTR, I.E., ITS POINTERS ARE STORED IN JTR. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - JA(IA(I)) ----- I=IMIN(K),...,IMAX(K) */

/*     IS DEFINED TO POINT TO THE LAST STRONG CONNECTION OF POINT I. */
/*     NOTE THAT THE ORIGINAL CONTENTS OF JA(IA(I)) IS OVERWRITTEN. */
/*     (IT IS PUT BACK IN SUBROUTINE WINT.) */

/*     - JTR(J)--------- J=1,IW(IMAX(K)+IWS+1)-1 */
/*     - IW(I) --------- I=IMIN(K)+IWS,...,IMAX(K)+IWS */
/*     - IFG(I)--------- I=IMIN(K),...,IMAX(K) */

/*     (IWS=0, IF K=1; IWS=IMAXW(K-1)+2-IMIN(K) OTHERWISE) */

/*     JTR IS INITIALIZED TO HAVE SAME FORM AS JA. JTR(J) */
/*     CONTAINS INFORMATION ON STRONG TRANSPOSE CONNECTIONS: */
/*     JTR(J) WITH IW(I+IWS)<=J<=IW(I+IWS+1)-1 POINTS TO THE STRONG */
/*     TRANSPOSE CONNECTIONS OF I. */


/* ===> INITIALIZATION OF WORK SPACE */

    /* Parameter adjustments */
    --time;
    --jtr;
    --ifg;
    --imaxw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --a;

    /* Function Body */
	told=clock();
    ilo = imin[*k];
    ihi = imax[*k];
    if (*k != 1) {
	iws = imaxw[*k - 1] + 2 - ilo;
    } else {
	iws = 0;
    }
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	ifg[i__] = 0;
	ja[ia[i__]] = ia[i__];
/* L5: */
    }
	tnew=clock();
    time[8] = time[8] + tnew - told;
    told = tnew;

/*     ************************* */
/*     * PARTIAL STANDARD SORT * */
/*     ************************* */

/*     NOTE: NON-NEGATIVE CONNECTIONS ARE ALWAYS DEFINED TO BE WEAK! THE */
/*           STRONGEST CONNECTION IS ALWAYS FOLLOWING THE DIAGONAL. */

    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	jlo = ia[i__] + 1;
	jhi = ia[i__ + 1] - 1;
	if (jhi < jlo) {
	    goto L590;
	}

/* ===>   FIND STRONGEST CONNECTION AND SUM OF |OFF-DIAGONALS| OF ROW I */

	amx = a[jlo];
	amn = a[jlo];
	jmx = jlo;
	if (*ecg1 != 0.) {
	    rs = 0.;
	    i__2 = jhi;
	    for (j = jlo + 1; j <= i__2; ++j) {
		rs += (d__1 = a[j], fabs(d__1));
		if (a[j] < amx) {
		    amx = a[j];
		    jmx = j;
		} else if (a[j] > amn) {
		    amn = a[j];
		}
/* L550: */
	    }

/* ===>     TEST FOR POSITIVE OFF-DIAGONALS / DIAGONAL DOMINANCE */

	    if (amx >= 0. || rs <= *ecg1 * a[ia[i__]]) {
		goto L590;
	    }
	} else {
	    i__2 = jhi;
	    for (j = jlo + 1; j <= i__2; ++j) {
		if (a[j] < amx) {
		    amx = a[j];
		    jmx = j;
		} else if (a[j] > amn) {
		    amn = a[j];
		}
/* L555: */
	    }

/* ===>   TEST FOR POSITIVE OFF-DIAGONALS */

	    if (amx >= 0.) {
		goto L590;
	    }
	}

/* ===>   PUT STRONGEST CONNECTION IN FIRST POSITION */

	ast = *ecg2 * amx;
	imx = ja[jmx];
	a[jmx] = a[jlo];
	ja[jmx] = ja[jlo];
	a[jlo] = amx;
	ja[jlo] = imx;
	if (amn <= ast) {
	    goto L580;
	}
	++jhi;

/* ===>   DECREASE JHI UNTIL A STRONG CONNECTION IS FOUND */
/*       (IF JLO >= JHI STOP: ALL CONNECTIONS ARE SORTED) */

L560:
	--jhi;
	if (jlo >= jhi) {
	    goto L580;
	}
	if (a[jhi] > ast) {
	    goto L560;
	}

/* ===>   INCREASE JLO UNTIL A WEAK CONNECTION IS FOUND */
/*       (IF JLO >= JHI STOP: ALL CONNECTIONS ARE SORTED) */

L570:
	++jlo;
	if (jlo >= jhi) {
	    goto L580;
	}
	if (a[jlo] <= ast) {
	    goto L570;
	}

/* ===>   INTERCHANGE A(JHI) AND A(JLO) */

	atmp = a[jhi];
	itmp = ja[jhi];
	a[jhi] = a[jlo];
	ja[jhi] = ja[jlo];
	a[jlo] = atmp;
	ja[jlo] = itmp;
	goto L560;

/* ===>   ROW SORTED --  SET JA(IA(I)) TO LAST STRONG CONNECTION */

L580:
	ja[ia[i__]] = jhi;

/* ===>   COUNT STRONG TRANSPOSE CONNECTIONS IN ROW I */

	i__2 = jhi;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    ++ifg[ja[j]];
/* L585: */
	}
L590:
	;
    }

/* ===> INITIALIZATION OF WORK SPACE FOR STRONG TRANSPOSE CONNECTIONS */

    iw[ilo + iws] = 1;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	iw[i__ + iws + 1] = iw[i__ + iws] + ifg[i__];
	ifg[i__] = iw[i__ + iws];
/* L5010: */
    }
    *mdjtr = iw[ihi + iws + 1] - 1;
    if (*mdjtr > *ndjtr) {
	goto L9901;
    }

/* ===> LOAD POINTERS TO STRONG TRANSPOSE CONNECTIONS INTO JTR */

    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	i__2 = ja[ia[i__]];
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    ii = ja[j];
	   // jtr[ifg[ii]] = i__;
		jtr[ifg[ii]]=(doublereal)(i__);
	    ++ifg[ii];
/* L5020: */
	}
/* L5030: */
    }
	tnew=clock();
    time[1] = time[1] + tnew - told;
    return 0;

/* ===> ERROR MESSAGES */

L9901:
    //io___130.ciunit = *ium;
    //s_wsfe(&io___130);
    //e_wsfe();
	printf("*** ERROR IN RWSRT: NDA TOO SMALL ***\n");
    *ierr = 1;
    return 0;

} /* rwsrt_ */


/* ....................................................................... */

/*     PCOL                                            SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int pcol_(integer *k, integer *ierr, integer *ia, integer *
	ja, integer *iw, integer *imin, integer *imax, integer *imaxw, 
	integer *icg, integer *ifg, /*integer*/ doublereal *ir, /*integer*/ doublereal *jtr, /*real*/ unsigned int *time, 
	integer *ndicg, integer *ndir, integer *ium, integer *mdicg, integer *
	mdir, integer *iias)
{
    /* Format strings */
   // static char fmt_9910[] = "(\002 *** ERROR IN PCOL: NDIG TOO SMALL ***"
	//    "\002)";
   // static char fmt_9920[] = "(\002 *** ERROR IN PCOL: NDU TOO SMALL ***\002)"
	    ;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    //integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, ic, ii, jj, in, jv, iic, ihi, iii, iip, ilo, iws, 
	    ilo1, iiii;
	static unsigned int told;
    static integer nscn, itop;
	static unsigned int tnew;
    static integer npts, jval0, npts1;
    static integer jvalx, jcnbhi, ntrlim, nscnmx;

    /* Fortran I/O blocks */
    //static cilist io___157 = { 0, 0, 0, fmt_9910, 0 };
    //static cilist io___158 = { 0, 0, 0, fmt_9920, 0 };



/*     PRE-COLORING ALGORITHM FOR GRID K. THIS IS THE VERSION AS */
/*     DESCRIBED IN RUGE/STUEBEN (BRISTOL). THE GOAL IS TO OBTAIN QUICKLY */
/*     A TENTATIVE SET OF C-POINTS WITH THE FOLLOWING PROPERTIES: */

/*       - THE C-POINTS ARE ONLY WEAKLY CONNECTED AMONG EACH OTHER; */
/*       - EACH F-POINT HAS (AT LEAST) ONE STRONG CONNECTION TO A C-POINT */
/*         (EXCEPT FOR SOME EXCEPTIONAL F-POINTS, E.G., THOSE WHICH DO */
/*         NOT HAVE A STRONG CONNECTION AT ALL: THE "FORCED F-POINTS"). */

/*     ON EXIT, ICG(I) (IMIN(K)<=I<=IMAX(K)) CONTAINS ALL THE INFORMATION */
/*     ON THE COLORING. IN DETAIL: */

/*        ICG(I) = 1: I IS C-POINT; */
/*        ICG(I) = 0: I IS FORCED F-POINT, I.E., IT IS AN F-POINT WHICH */
/*                    HAS NO STRONG CONNECTION AT ALL; */
/*        ICG(I) =-1: I IS "REGULAR" F-POINT, I.E., I HAS AT LEAST ONE */
/*                    STRONG CONNECTION TO A C-POINT; */
/*        ICG(I) =-2: I IS "EXCEPTIONAL" F-POINT, I.E., I HAS NO STRONG */
/*                    CONNECTION TO A C-POINT. FURTHERMORE, NO REGULAR */
/*                    F-POINT IS STRONGLY CONNECTED TO I AND ALL POINTS */
/*                    WITH ICG=-2 HAVE NO STRONG CONNECTION AMONG EACH */
/*                    OTHER. (NOTE: THESE POINTS DON'T CONTRIBUTE FROM */
/*                    ANY C-POINT. ALSO, IT DOES NOT MAKE SENSE TO MAKE */
/*                    THESE POINTS C-POINTS AS NO F-POINT CONTRIBUTES */
/*                    FROM THEM.) */

/*     ============ COMMENTS ON INPUT =================================== */

/*     - JA(IA(I)) ----- I=IMIN(K),...,IMAX(K) */

/*     USED AS DEFINED IN RWSRT. NOT CHANGED. */

/*     - JTR(J)--------- J=1,IW(IMAX(K)+IWS+1)-1 */
/*     - IW(I) --------- I=IMIN(K)+IWS,...,IMAX(K)+IWS */

/*     (IWS=0, IF K=1; IWS=IMAXW(K-1)+2-IMIN(K) OTHERWISE) */

/*     USED AS DEFINED IN RWSRT. NOT CHANGED. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - ICG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     USED TO DISTINGUISH F-, C- AND U- (UNDECIDED) POINTS. */
/*     ON ENTRY, ICG IS DEFINED TO BE 0 FOR FORCED F-POINTS AND -2 */
/*     FOR U-PNTS. ON EXIT: SEE ABOVE. */

/*     - IFG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     DEFINED TO BE A MEASURE FOR IMPORTANCE OF MAKING THE U-PNT I */
/*     A C-POINT. THIS MEASURE IS A VALUE BETWEEN JVAL0 AND JVALX WITH */

/*       JVAL0=IMAX(K)+NPTS+2,   NPTS=IMAX(K)-IMIN(K)+1; */
/*       JVALX=JVAL0+4*NTRAV+1,  NTRAV=AVVERAGE NUMBER OF STRONG TRANS- */
/*                                     POSE CONNECTIONS */

/*     THE HIGHER THIS VALUE, THE MORE IMPORTANT IS IT TO MAKE U-POINT I */
/*     A C-POINT. ALSO, THE VALUE JV=IFG(I) IS JUST THE "ORIGIN" OF */
/*     A LIST WHICH CONNECTS ALL POINTS HAVING THE SAME MEASURE JV OF */
/*     IMPORTANCE. */

/*     - ICG(II) -------- II=IMAX(K)+2,...,JVAL0-1 */
/*     - IFG(II) -------- II=IMAX(K)+2,...,JVAL0-1 */

/*     LEFT AND RIGHT STACK POINTERS IN A DOUBLY LINKED LIST. THERE */
/*     ARE SEVERAL SUCH LISTS, EACH OF THEM LINKING U-POINTS WHICH */
/*     HAVE THE SAME MEASURE OF IMPORTANCE JV TO BE MADE C-POINTS. */
/*     EACH OF THE VALUES JV MAY BE REGARDED AS THE "ORIGIN" OF SUCH */
/*     A LIST DESCRIBED IN THE FOLLOWING. */

/*     - ICG(JV) -------- JV=JVAL0,...,JVALX */
/*     - IFG(JV) -------- JV=JVAL0,...,JVALX */

/*     USED AS POINTERS INTO THE LISTS TO INDICATE ITS BEGINNING AND ITS */
/*     END: IFG(JV) POINTS TO THE FIRST POINT IN THE LIST AND ICG(JV) */
/*     POINTS TO THE LAST ONE. THE ROUGH PICTURE OF THE LIST WHICH */
/*     CORRESPONDS TO A VALUE JVAL0 <= JV <= JVALX LOOKS AS  FOLLOWS: */


/*                                  IFG */
/*       ---------------------------------------------------------- */
/*       |                                                        | */
/*       |            IFG           IFG           IFG             | */
/*       ---> ------ ----> ------- ----> ------- ----> ------- ---- */
/*            |    |       |     |       |     |       |     | */
/*            | JV |       | II1 |       | II2 |       | II3 | */
/*            |    |  ICG  |     |  ICG  |     |  ICG  |     | */
/*       ---- ------ <---- ------- <---- ------- <---- ------- <--- */
/*       |                                                        | */
/*       |                          ICG                           | */
/*       ---------------------------------------------------------- */


/*     HERE, II1, II2,... LOGICALLY REPRESENT THE "PHYSICAL" POINTS */
/*     I1, I2,... (WHICH ARE VALUES BETWEEN IMIN(K) AND IMAX(K)). */
/*     THE INTERCONNECTION BETWEEN THESE TWO REPRESENTATIONS OF THE */
/*     SAME POINTS IS GIVEN BY THE RELATION */

/*                          II := I+NPTS+1. */

/*     AN EMPTY LIST IS CHARACTERIZED BY ICG(JV)=JV, IFG(JV)=JV. */
/*     OBVIOUSLY, IT IS QUITE EASY TO REMOVE OR ADD POINTS TO THE LIST. */
/*     SUCH RE-ARRANGEMENTS ARE NECESSARY AS (IN THE COLORING PART OF */
/*     THE COARSENING ALGORITHM) UNDECIDED POINTS BECOME C- OR F-POINTS */
/*     (THEY HAVE TO BE REMOVED FROM THEIR LIST) AND THE MEASURE VALUES */
/*     JV OF SOME POINTS CHANGE DURING THE ALGORITHM. SUCH POINTS HAVE */
/*     TO BE MOVED FROM ONE LIST TO ANOTHER ONE. IN ORDER TO KEEP TRACK */
/*     OF THOSE POINTS WHICH HAVE TO BE RE-ARRANGED, A RESET-STACK IS */
/*     USED WHICH CONTAINS ALL THESE POINTS I (I.E. THEIR "LOGICAL" */
/*     NUMBERS II). THE POINTER IR IS USED TO POINT FROM ONE POINT IN */
/*     THE STACK TO THE PREVIOUS ONE. */

/*     - IR(J) --------- J=1,...,NPTS */

/*     USED BELOW AS STACK-POINTER FOR POINTS TO BE RESET IN PERFORMING */
/*     THE COLORING ALGORITHM. ITOP IS THE TOP-OF-STACK POINTER. */
/*     THE STACK IS EMPTY IF ITOP=-1. */

/*     THE GLOBAL PICTURE OF THE POINTERS IS SKETCHED IN THE FOLLOWING */


/*                                          LIST END */
/*                                       <-------------| */
/*                                         LIST START  | */
/*                                       <------------|| */
/*                                                    || */
/*                  |-- II=I+NPTS+1 ->|            IFG||ICG */
/*                  |                 |               || */
/*                  |                 |               || */
/*       IMIN(K)          IMAX(K)           JVAL0             JVALX */
/*          |       I        |        II      |       JV        | */
/*          |=======*========|========*=======|========*========| */
/*          |                |                |                 | */
/*               GRID K          LOGICAL NUM      LIST ORIGINS */
/*                                              (MEASURE VALUES) */
/*                 ||                                  | */
/*                 ||         IFG                      | */
/*                 ||--------------------------------->| */
/*                 | */
/*                 |     |---> 1 (C) */
/*                 |     | */
/*                 | ICG |---> 0 (FF) */
/*                 |-----| */
/*                       |--->-1 (F) */
/*                       | */
/*                       |--->-2 (U) */


/*          1                NPTS */
/*          |        J        | */
/*          |========*========| */
/*          |                 | */
/*           SHIFTED GRID K */
/*                   | */
/*                   | IR */
/*                   |-----> RESET STACK (ITOP = TOP-OF-STACK) */



/* ===> PREPARATION */

    /* Parameter adjustments */
    --time;
    --jtr;
    --ir;
    --ifg;
    --icg;
    --imaxw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;

    /* Function Body */
	told=clock();
    ilo = imin[*k];
    ihi = imax[*k];
    if (*k != 1) {
	iws = imaxw[*k - 1] + 2 - ilo;
    } else {
	iws = 0;
    }
    ilo1 = ilo - 1;
    npts = imax[*k] - imin[*k] + 1;
    npts1 = npts + 1;
    ntrlim = ((iw[ihi + iws + 1] - iw[ilo + iws]) << 1) / npts; // скобки поставлены в соответствии с приоритетом.
    nscnmx = 0;
    jval0 = ihi + npts1 + 1;
    jvalx = jval0 + (ntrlim << 1) + 1;
/* Computing MAX */
    i__1 = *mdicg, i__2 = jvalx + *iias;
    *mdicg = max(i__1,i__2);
    *mdir = max(*mdir,npts);
    if (jvalx > *ndicg) {
	goto L9901;
    }
    if (npts > *ndir) {
	goto L9902;
    }

/* ===> PUT INITIAL "MEASURE" FOR EACH POINT I INTO IFG(I). */

    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	nscn = iw[i__ + iws + 1] - iw[i__ + iws];
	if (nscn <= ntrlim) {
	    ifg[i__] = jval0 + nscn;
	} else {
	    ifg[i__] = jvalx;
	}
/* L2: */
    }

/* ===> SET CIRCULARLY LINKED LISTS AND RESET-STACK TO EMPTY */

    itop = -1;
    i__1 = jvalx;
    for (j = jval0; j <= i__1; ++j) {
	icg[j] = j;
	ifg[j] = j;
/* L10: */
    }

/* ===> PUT ALL U-POINTS OF GRID K INTO LISTS (I.E., NO FORCED F-POINTS). */
/*     ADD POINTS ALWAYS TO THE END OF THEIR CORRESPONDING LIST. */
/*     IN THE FOLLOWING, JCNBHI DENOTES THE ACTUAL HIGHEST MEASURE VALUE */
/*     (AMONG U-POINTS). */

    jcnbhi = 0;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	if (ja[ia[i__]] > ia[i__]) {
	    goto L15;
	}
	icg[i__] = 0;
	goto L20;
L15:
	icg[i__] = -2;
	ii = i__ + npts1;
	//ir[i__ - ilo1] = 0;
	ir[i__ - ilo1] = 0.0;
	jv = ifg[i__];
	if (jv > jcnbhi) {
	    jcnbhi = jv;
	}
	icg[ii] = icg[jv];
	ifg[ii] = jv;
	icg[jv] = ii;
	ifg[icg[ii]] = ii;
L20:
	;
    }

/*     **************** */
/*     * PRE-COLORING * */
/*     **************** */

/*     PICK A U-POINT IC WITH MAXIMAL MEASURE AS GIVEN BY JCNBHI. */
/*     (TAKE THE FIRST FROM CORRESPONDING LIST: FIRST-IN/FIRST-OUT) */
/*     MAKE THAT POINT A C-POINT AND REMOVE IT FROM LISTS. THEN MAKE */
/*     ALL STRONG TRANSPOSE U-CONNECTIONS F-POINTS AND REMOVE THEM FROM */
/*     THEIR LISTS. UPDATE THE MEASURE OF IMPORTANCE FOR U-POINTS TO */
/*     BECOME C-POINTS AND ADD THESE POINTS TO THE RESET-STACK. FINALLY, */
/*     RE-ARRANGE THE LISTS BY SWEEPING THROUGH THE RESET STACK. THEN */
/*     PICK ANOTHER U-POINT IC. */

/*     IF LIST CORRESPONDING TO JCNBHI IS EMPTY, GO TO NEXT LOWER VALUE. */
/*     PRE-COLOURING IS FINISHED IF JCNBHI<=JVAL0. ALL U-POINTS LEFT AT */
/*     THAT TIME, WILL BE REGARDED AS F-POINTS LATER. */

L30:
    if (jcnbhi <= jval0) {
	goto L100;
    }
    iic = ifg[jcnbhi];
    if (iic != jcnbhi) {
	goto L40;
    }
    --jcnbhi;
    goto L30;

/* ===> CREATE C-POINT */

L40:
    ic = iic - npts1;
    icg[ic] = 1;
    icg[ifg[iic]] = icg[iic];
    ifg[icg[iic]] = ifg[iic];

/* ===> FOR POINT IC WITH ECCESSIVE NUMBER OF STRONG TRANSPOSE */
/*     CONNECTIONS: MAKE IT A C-POINT BUT LET ITS CONNECTED POINTS */
/*     REMAIN UNDECIDED. */

    if (jcnbhi == jvalx) {
	goto L78;
    }

/* ===> CREATE F-POINTS AROUND ABOVE C-POINT */

    i__1 = iw[ic + iws + 1] - 1;
    for (j = iw[ic + iws]; j <= i__1; ++j) {
	//i__ = jtr[j];
    i__ = (integer)(jtr[j]);
	if (icg[i__] != -2) {
	    goto L77;
	}
	icg[i__] = -1;
	ii = i__ + npts1;
	icg[ifg[ii]] = icg[ii];
	ifg[icg[ii]] = ifg[ii];

/* ===>   INCREMENT MEASURE FOR ALL STRONG U-CONNECTIONS III OF I */
/*       (IF NOT YET MEASURE=JVALX) AND PUT THEM ON RESET STACK (IF NOT */
/*       YET THERE) */

	i__2 = ja[ia[i__]];
	for (jj = ia[i__] + 1; jj <= i__2; ++jj) {
	    iii = ja[jj];
	    if (icg[iii] != -2 || ifg[iii] >= jvalx) {
		goto L76;
	    }
	    ++ifg[iii];
	    iiii = iii - ilo1;
	   // if (ir[iiii] != 0) {
		//goto L76;
	    //}
		if (((integer)(ir[iiii])) != 0) {
		goto L76;
	    }
	   // ir[iiii] = itop;
		 ir[iiii] = (doublereal)(itop);
	    itop = iiii;
L76:
	    ;
	}
L77:
	;
    }

/* ===> DECREMENT MEASURE FOR ALL STRONG U-CONNECTIONS I OF IC */
/*     AND PUT THEM ON RESET-STACK (IF NOT YET THERE) */

L78:
    i__1 = ja[ia[ic]];
    for (j = ia[ic] + 1; j <= i__1; ++j) {
	i__ = ja[j];
	if (icg[i__] != -2) {
	    goto L87;
	}
	--ifg[i__];
	ii = i__ - ilo1;
	//if (ir[ii] != 0) {
	  //  goto L87;
	//}
	if (((integer)(ir[ii])) != 0) {
	    goto L87;
	}
	//ir[ii] = itop;
	ir[ii] = (doublereal)(itop);
	itop = ii;
L87:
	;
    }

/* ===> REARRANGE THE LISTS BY SWEEPING THROUGH RESET-STACK. */
/*     THEN GO BACK TO PICK ANOTHER U-POINT IC. */

    in = itop;
    itop = -1;

L90:
    if (in <= 0) {
	goto L30;
    }
    i__ = in + ilo1;
    ii = i__ + npts1;
    if (icg[i__] != -2) {
	goto L95;
    }
    ifg[icg[ii]] = ifg[ii];
    icg[ifg[ii]] = icg[ii];
    jv = ifg[i__];
    if (jv > jcnbhi) {
	jcnbhi = jv;
    }
    icg[ii] = icg[jv];
    ifg[ii] = jv;
    icg[jv] = ii;
    ifg[icg[ii]] = ii;
L95:
    iip = in;
   // in = ir[in];
	in =(integer)(ir[in]);
    //ir[iip] = 0;
	ir[iip] = 0.0;
    goto L90;

L100:
    //ctime_(&tnew);
	tnew=clock();
    time[2] = time[2] + tnew - told;
    return 0;

/* ===> ERROR MESSAGES */

L9901:
    //io___157.ciunit = *ium;
    //s_wsfe(&io___157);
    //e_wsfe();
	printf("*** ERROR IN PCOL: NDIG TOO SMALL ***\n");
    *ierr = 6;
    return 0;

L9902:
    //io___158.ciunit = *ium;
    //s_wsfe(&io___158);
    //e_wsfe();
	printf("*** ERROR IN PCOL: NDU TOO SMALL ***\n");
    *ierr = 4;
    return 0;

} /* pcol_ */


/* ....................................................................... */

/*     PWINT                                             SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int pwint_(integer *k, doublereal *ewt2, integer *ichk, 
	integer *mmax, doublereal *a, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, integer *imaxw, integer 
	*ifg, integer *icg, /*integer*/ doublereal *ncolor, /*real*/ unsigned int *time, integer *ierr, 
	integer *irow0, integer *ncolx, integer *nda, integer *ndja, integer *
	ndicg, integer *ndia, integer *ium, integer *mda, integer *mdja, 
	integer *iajas)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 INTERPOLATION OPERATOR NO.\002,i3,\002 C"
   //	    "OMPLETED. C-POINTS\002,\002 ADDED IN PWINT:\002,i4)";
   // static char fmt_9030[] = "(\002 *** ERROR IN PWINT: NDIA TOO SMALL **"
	//    "*\002)";
    //static char fmt_9040[] = "(\002 *** ERROR IN PWINT: NDIG TOO SMALL **"
	//    "*\002)";
    //static char fmt_9020[] = "(\002 *** ERROR IN PWINT: NDA TOO SMALL ***"
	 //   "\002)";
   // static char fmt_9050[] = "(\002 *** ERROR IN PWINT: NDJA TOO SMALL **"
	//    "*\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    //integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static integer j1, j2, ic, nc, ii, jj, ip;
    static doublereal si;
    static integer is;
    static doublereal ww;
    static integer jw0, ihi, jhi, ilo, jlo, jwx, icgp, jjhi, jjlo;
	static unsigned int told, tnew;
    static integer npts;
    static doublereal ewt2i;
    static integer ndaja, iblck;
    static doublereal scale;
    static integer nptsc, jwpos, iblck1, ncondc, ncount;

    /* Fortran I/O blocks */
   // static cilist io___192 = { 0, 0, 0, fmt_9000, 0 };
   // static cilist io___194 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___195 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___196 = { 0, 0, 0, fmt_9020, 0 };
    //static cilist io___197 = { 0, 0, 0, fmt_9050, 0 };



/*     SET UP FINAL COARSER GRID K+1 AND INTERPOLATION FORMULA FROM */
/*     GRID K+1 TO GRID K. THIS IS THE VERSION AS DESCRIBED IN RUGE/ */
/*     STUEBEN (BRISTOL). PWINT ASSUMES THE GRID TO BE PRE-COLORED */
/*     BY SUBROUTINE PCOL. */

/*     ON EXIT, A, JA AND IW ARE SET TO CONTAIN THE INTERPOLATION */
/*     WEIGHTS AND CORRESPONDING POINTERS AS REQUIRED IN THE SOLUTION */
/*     PHASE OF AMG1R5. ALSO, ICG(I) (IMIN(K)<=I<=IMAX(K)) ARE SET TO */
/*     THEIR FINAL VALUES, EXCEPT FOR THOSE I WITH ICG(I)<0. */

/*     ============ COMMENTS ON INPUT =================================== */

/*     - JA(IA(I)) ----- I=IMIN(K),...,IMAX(K) */

/*     ASSUMED TO POINT TO THE LAST STRONG CONNECTION OF POINT I (OR TO */
/*     IA(I) IF THERE IS NO SUCH CONNECTION). ON EXIT, RESET TO ORIGINAL */
/*     VALUES (I.E. JA(IA(I))=I). */

/*     - ICG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     ASSUMED TO CONTAIN INFORMATION ON PRE-COLORING: */

/*       ICG(I)>0: I IS C-POINT */
/*       ICG(I)=0: I IS FORCED F-POINT (I.E. I HAS NO STRONG CONNNECTION) */
/*       ICG(I)<0: I IS F-POINT WITH (AT LEAST) ONE STRONG CONNECTION. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - IFG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     IS USED FOR SEVERAL PURPOSES. IN PARTICULAR, TO DISTINGUISH */
/*     INTERPOLATORY AND NON-INTERPOLATORY POINTS. */

/*     - NCOLOR(I) ----- I=1,...,#POINTS ON GRID K */

/*     IS SET TO F-POINT-COLORS TO BE USED LATER IN SUBROUTINE OPDFN. */


    /* Parameter adjustments */
    --time;
    --ncolor;
    --icg;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --a;

    /* Function Body */
	told=clock();
    ndaja = min(*nda,*ndja);
    ncount = 0;
    if (*k == 1) {
	iminw[1] = 1;
	iw[1] = ia[imax[1] + 1];
    }

    ilo = imin[*k];
    ihi = imax[*k];
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	ifg[i__] = 0;
/* L12: */
    }

/* ===> SWEEP OVER F-POINTS I WHICH HAVE AT LEAST ONE STRONG CONNECTION */

    iblck = iminw[*k];
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	if (icg[i__] >= 0) {
	    goto L400;
	}
	jlo = ia[i__] + 1;
	jhi = ja[ia[i__]];
	ewt2i = *ewt2 / a[jlo];
	ncondc = 0;

/* ===>   INITIALIZE "BLOCK" OF INTERPOLATION WEIGHTS OF POINT I */

L30:
	jw0 = iw[iblck];
	jwx = jw0;
	if (jwx > ndaja) {
	    goto L2000;
	}
	i__2 = jhi;
	for (j = jlo; j <= i__2; ++j) {
	    ii = ja[j];
	    if (icg[ii] <= 0) {
		goto L20;
	    }
	    a[jwx] = a[j];
	    ja[jwx] = ii;
	    ifg[ii] = jwx;
	    ++jwx;
	    if (jwx > ndaja) {
		goto L2000;
	    }
L20:
	    ;
	}
	a[jwx] = a[ia[i__]];
	ja[jwx] = i__;
/* Computing MAX */
	i__2 = *mda, i__3 = jwx + *iajas;
	*mda = max(i__2,i__3);
/* Computing MAX */
	i__2 = *mdja, i__3 = jwx + *iajas;
	*mdja = max(i__2,i__3);
	ifg[i__] = jwx;

/* ===>   SWEEP OVER STRONGLY CONNECTED F-POINTS. THESE MUST BE "COVERED" */
/*       BY A TOTAL WEIGHT DEFINED BY EWT2. IF AN F-POINT HAS NO STRONG */
/*       CONNECTIONS, REGARD IT TO BE COVERED, BUT DO NOT DISTRIBUTE */
/*       THE CORRESPONDING WEIGHT (ERROR AT SUCH A POINT CAN BE ASSUMED */
/*       TO BE VERY SMALL!). */

	i__2 = jhi;
	for (j = jlo; j <= i__2; ++j) {
	    ii = ja[j];
	    if (icg[ii] >= 0) {
		goto L150;
	    }

/* ===>     COMPUTE DEPENDENCE ON SET OF INTERPOLATION POINTS */

	    s = 0.;
	    si = 0.;
	    jjlo = ia[ii] + 1;
	    jjhi = ia[ii + 1] - 1;
	    i__3 = jjhi;
	    for (jj = jjlo; jj <= i__3; ++jj) {
		if (ifg[ja[jj]] < jw0) {
		    goto L110;
		}
		if (ja[jj] == i__) {
		    si = a[jj];
		}
		s += a[jj];
L110:
		;
	    }
	    if (*ichk == 2) {
		goto L111;
	    }
	    if (s == 0.) {
		a[jwx] += a[j];
		goto L150;
	    } else {
		goto L135;
	    }

/* ===>     CHECK DEPENDENCE ON SET OF INTERPOLATION POINTS */

L111:
	    if (s - si <= ewt2i * a[j] * a[jjlo]) {
		goto L135;
	    }

/* ===>     DEPENDENCE TOO SMALL: IF THERE IS NOT YET A CONDITIONAL */
/*         C-POINT, MAKE II SUCH A POINT AND RESTART THE */
/*         PROCESS FOR DEFINING INTERPOLATION WEIGHTS FOR POINT I. */
/*         OTHERWISE MAKE I ITSELF A C-POINT AND LEAVE II AN F-POINT. */

	    if (ncondc == 0) {
		++ncount;
		ncondc = 1;
		ip = ii;
		icgp = icg[ii];
		icg[ii] = 1;
		goto L30;
	    } else {
		icg[i__] = 1;
		icg[ip] = icgp;
		i__3 = jwx;
		for (jj = jw0; jj <= i__3; ++jj) {
		    ifg[ja[jj]] = 0;
/* L120: */
		}
		goto L400;
	    }

/* ===>     DISTRIBUTE THE WEIGHT OF POINT II */

L135:
	    ww = a[j] / s;
	    i__3 = jjhi;
	    for (jj = jjlo; jj <= i__3; ++jj) {
		if (ifg[ja[jj]] >= jw0) {
		    a[ifg[ja[jj]]] += a[jj] * ww;
		}
/* L140: */
	    }
L150:
	    ;
	}

/* ===>   ALL NECESSARY POINTS ARE COVERED. NOW DISTRIBUTE WEIGHTS FROM */
/*       WEAK CONNECTIONS OF POINT I (ANALOGOUS AS ABOVE) */

	i__2 = ia[i__ + 1] - 1;
	for (j = jhi + 1; j <= i__2; ++j) {
	    ii = ja[j];
	    if (icg[ii] == 0) {
		goto L190;
	    }
	    s = 0.;
	    jjlo = ia[ii] + 1;
	    jjhi = ia[ii + 1] - 1;
	    i__3 = jjhi;
	    for (jj = jjlo; jj <= i__3; ++jj) {
		if (ifg[ja[jj]] >= jw0) {
		    s += a[jj];
		}
/* L160: */
	    }
	    if (s == 0.) {
		a[jwx] += a[j];
	    } else {
		ww = a[j] / s;
		i__3 = jjhi;
		for (jj = jjlo; jj <= i__3; ++jj) {
		    if (ifg[ja[jj]] >= jw0) {
			a[ifg[ja[jj]]] += a[jj] * ww;
		    }
/* L170: */
		}
	    }
L190:
	    ;
	}

	icg[i__] = -iblck;
	++iblck;
	iw[iblck] = jwx + 1;
L400:
	;
    }

/* ===> SET ICG; RESET JA(IA(I)); CHECK SIZE OF COARSEST GRID */

    ic = ihi;
    imin[*k + 1] = ic + 1;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	ja[ia[i__]] = i__;
	if (icg[i__] <= 0) {
	    goto L900;
	}
	++ic;
	icg[i__] = ic;
	if (ic < *ndicg) {
	    icg[ic] = 0;
	}
L900:
	;
    }
    imax[*k + 1] = ic;

    npts = ihi - ilo + 1;
    nptsc = imax[*k + 1] - imin[*k + 1] + 1;
    if (nptsc == 1) {
	*mmax = *k + 1;
    }
    if (nptsc == 1 && *irow0 < 2) {
	*mmax = *k;
    }
    if (nptsc == npts || nptsc == 0) {
	*mmax = *k;
    }
    if (*k >= *mmax) {
	goto L1000;
    }
    if (ic >= *ndia) {
	goto L1700;
    }
    if (ic >= *ndicg) {
	goto L1800;
    }

/* ===> RE-ARRANGE A */

L1000:
    iblck1 = iminw[*k];
    jwpos = iw[iblck1];
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	if (icg[i__] >= 0) {
	    goto L950;
	}
	iblck = -icg[i__];
	j1 = iw[iblck];
	j2 = iw[iblck + 1] - 1;
	if (j2 <= j1) {
	    icg[i__] = 0;
	} else {
	    icg[i__] = -iblck1;
	    iw[iblck1] = jwpos;
	    scale = -1. / a[j2];
	    i__2 = j2 - 1;
	    for (j = j1; j <= i__2; ++j) {
		a[jwpos] = a[j] * scale;
		ja[jwpos] = icg[ja[j]];
		++jwpos;
/* L920: */
	    }
	    ++iblck1;
	}
L950:
	;
    }
    imaxw[*k] = iblck1 - 1;
    iw[iblck1] = jwpos;

/* ===> STORE TYPE OF POINTS (I.E. C, F OR FF)  ON VECTOR NCOLOR */

    is = 1 - ilo;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	nc = icg[i__];
	if (nc < 0) {
	   // ncolor[i__ + is] = 1;
		ncolor[i__ + is]=1.0;
	} else if (nc > 0) {
	   // ncolor[i__ + is] = 2;
		ncolor[i__ + is] = 2.0;
	} else {
	   // ncolor[i__ + is] = 3;
		ncolor[i__ + is] = 3.0;
	}
/* L960: */
    }
    *ncolx = 1;

/* ===> EXIT */

    //io___192.ciunit = *ium;
    //s_wsfe(&io___192);
    //do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
    //do_fio(&c__1, (char *)&ncount, (ftnlen)sizeof(integer));
    //e_wsfe();
	if (yes_print_amg) {
	 printf("( INTERPOLATION OPERATOR NO.,%d, \n", k[0]);
	 printf("COMPLETED. C-POINTS, ADDED IN PWINT:,%d)\n",ncount);
	}
    //ctime_(&tnew);
	tnew=clock();
    time[4] = time[4] + tnew - told;
    return 0;

/* ===> ERROR MESSAGES */

L1700:
    //io___194.ciunit = *ium;
    //s_wsfe(&io___194);
    //e_wsfe();
	printf("( *** ERROR IN PWINT: NDIA TOO SMALL ***)\n");
    *ierr = 2;
    return 0;

L1800:
    //io___195.ciunit = *ium;
    //s_wsfe(&io___195);
    //e_wsfe();
	printf( "( *** ERROR IN PWINT: NDIG TOO SMALL ***)\n");
    *ierr = 6;
    return 0;

L2000:
    if (*nda <= *ndja) {
	//io___196.ciunit = *ium;
	//s_wsfe(&io___196);
	//e_wsfe();
	printf("( *** ERROR IN PWINT: NDA TOO SMALL ***)\n");
	*ierr = 1;
    } else {
	//io___197.ciunit = *ium;
	//s_wsfe(&io___197);
	//e_wsfe();
	printf("( *** ERROR IN PWINT: NDJA TOO SMALL ***)\n");
	*ierr = 3;
    }
    return 0;

} /* pwint_ */


/* ....................................................................... */

/*     OPDFN                                             SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int opdfn_(integer *k, integer *ierr, integer *mmax, 
	doublereal *a, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *imaxw, integer *icg, integer *
	ifg, /*integer*/ doublereal *ncolor, integer *nstcol, integer *ncolx, /*real*/ unsigned int *time, 
	integer *nda, integer *ndja, integer *ium, integer *mda, integer *
	mdja, integer *iajas)
{
    /* Format strings */
    //static char fmt_9940[] = "(\002 *** ERROR IN OPDFN: INTERPOLATION ENTRY "
	//    "MISSING ON GRID\002,i3)";
   // static char fmt_9930[] = "(\002 --- WARNG: UNABLE TO STORE TRANSPOSE OF "
	//    "INTERPOLATION\002,\002 ON GRID \002,i2,\002 DURING \002/\002    "
	//    "        EXECUTION OF OPDFN, BECAUSE NDA OR NDJA \002,\002TOO SMA"
	 //   "LL.\002/\002            SETUP COMPUTATION IS SLOWING DOWN.\002)";
    //static char fmt_9000[] = "(\002 COARSE  GRID  OPERATOR NO.\002,i3,\002 C"
	//    "OMPLETED\002)";
    //static char fmt_9910[] = "(\002 *** ERROR IN OPDFN: NDA TOO SMALL ***"
	  //  "\002)";
   // static char fmt_9920[] = "(\002 *** ERROR IN OPDFN: NDJA TOO SMALL **"
	//    "*\002)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    //integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k1, ic, jb, if__;
    static doublereal ww;
    static integer ic1, ic2, ic3, if1, jf1, jf2, jc3, if2, jf3, ibl, ihi, ilo,
	     ist, ilo1;
    static doublereal wjf1;
    static integer icol;
	static unsigned int told;
    static integer jpos;
	static unsigned int tnew;
    static integer npts, ndaja;
    static integer iadrs, ialow, istti;

    /* Fortran I/O blocks */
   // static cilist io___216 = { 0, 0, 0, fmt_9940, 0 };
   // static cilist io___223 = { 0, 0, 0, fmt_9930, 0 };
    //static cilist io___224 = { 0, 0, 0, fmt_9000, 0 };
    //static cilist io___232 = { 0, 0, 0, fmt_9910, 0 };
    //static cilist io___233 = { 0, 0, 0, fmt_9920, 0 };



/*     THIS SUBROUTINE CONSTRUCTS THE CG-OPERATOR A(K) (K>1). */

/*     ONE ROW IS CONSTRUCTED AT A TIME, AND THE ROW STRUCTURE */
/*     (I.E., WHICH CONNECTION GOES WHERE) MUST BE DETERMINED */
/*     AT THE SAME TIME.  IN ORDER TO AVOID SEARCHES THROUGH */
/*     THE CURRENT ROW TO DETERMINE IF A POSITION FOR A CONNECTION */
/*     HAS ALREADY BEEN DEFINED, ICG (FOR LEVEL K) AND IFG (FOR LEVEL */
/*     K-1) ARE USED AUXILIARILY. FURTHERMORE, TO SPEED UP COMPUTATION, */
/*     THE TRANSPOSE OF INTERPOLATION IS TEMPORARILY STORED ON A (AT */
/*     THE END OF THE AVAILABLE SPACE). IFG (FOR LEVEL K) IS USED */
/*     AS POINTER TO THE CORRESPONDING ROWS. MORE DETAILS: SEE BELOW. */

/*     ============ COMMENTS ON INPUT =================================== */

/*     - IMIN(K), IMAX(K) */

/*     FIRST/LAST NUMBER OF GRID POINTS ON GRID K. */

/*     - IMINW(K-1), IMAXW(K-1) */

/*     FIRST/LAST "BLOCK" NUMBER OF INTERPOLATION WEIGHTS CONTRIBUTING */
/*     IN INTERPOLATION TO GRID K-1. */

/*     - IW(IBL) ----- IBL=IMINW(K-1),IMAXW(K-1)+1 */
/*     - A(J) -------- J=JLO,JHI      (JLO=IW(IMINW(K-1))) */
/*     - JA(J) ------- J=JLO,JHI      (JHI=IW(IMAXW(K-1)+1)-1) */

/*     ARE ASSUMED TO CONTAIN THE WEIGHTS OF INTERPOLATION ALONG WITH */
/*     THE NECESSARY POINTERS. IW(IBL) POINTS TO THE FIRST ENTRY OF */
/*     "BLOCK" NUMBER IBL IN A (CF. BELOW). */

/*     - ICG(IF) ----- IF=IMIN(K-1),IMAX(K-1) */

/*     ARE ASSUMED TO BE SET TO THEIR FINAL VALUES: */

/*       ICG(IF)>0: IF IS C-POINT OF GRID K-1 AND ICG(IF) POINTS JUST */
/*                  TO THE CORRESPONDING POINT ON GRID K; */
/*       ICG(IF)=0: IF IS F-POINT OF GRID K-1 WITHOUT ANY CONTRIBUTION */
/*                  IN INTERPOLATION FROM GRID K; */
/*       ICG(IF)<0: IF IS F-POINT OF GRID K-1; IBL=-ICG(IF) POINTS JUST */
/*                  TO THE "BLOCK" OF INTERPOLATION WEIGHTS. THAT MEANS, */
/*                  JW(J) (IW(IBL)<=J<=IW(IBL+1)-1) POINTS TO THE POINTS */
/*                  (ON GRID K) WHICH CONTRIBUTE IN INTERPOLATION TO IF. */
/*                  THE CORRESPONDING WEIGHTS ARE STORED IN A(J). */

/*     - IFG(IF) ----- IF=IMIN(K-1),IMAX(K-1) */

/*     WORK SPACE (SEE BELOW). IFG(IF) IS ASSUMED TO BE > -IMIN(K). */

/*     - IFG(IC) ----- IC=IMIN(K),IMAX(K)+1 */

/*     WORK SPACE (SEE BELOW). THE CONTENTS OF IFG(IC) IS ARBITRARY. */

/*     - ICG(IC) ----- IC=IMIN(K),IMAX(K) */

/*     WORK SPACE (SEE BELOW). ICG(IC) IS ASSUMED TO BE ZERO. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - A(J) -------- J=JJLO,JJHI */
/*     - JA(J) ------- J=JJLO,JJHI */

/*     (JJLO=MIN(NDA,NDJA)-IW(IMAXW(K-1)+1)+IW(IMINW(K-1)+1), */
/*      JJHI=MIN(NDA,NDJA), I.E. THE SPACE OF LENGTH OF INTERPOLATION */
/*      (K-1) AT THE END OF A AND JA, RESPECTIVELY.) */

/*     IS USED TO STORE THE TRANSPOSE OF INTERPOLATION W(K-1). THIS IS */
/*     DONE IN ORDER TO SPEED UP THE COMPUTATION OF THE OPERATOR A(K). */
/*     IF DURING ASSEMBLAGE OF A(K) THIS WORK SPACE IS REQUIRED BY A(K), */
/*     PROCESSING CONTINUES USING THE NOT YET REDEFINED PART OF THE WORK */
/*     SPACE AND RECALCULATING THE LOST ENTRIES EACH TIME THEY ARE */
/*     NEEDED. THUS THE CALCULATION SLOWS DOWN. */

/*     - IFG(IC) ----- IC=IMIN(K),IMAX(K)+1 */

/*     IS USED AS POINTER FOR THE TRANSPOSE OF INTERPOLATION: THE COARSE- */
/*     GRID POINT IC CONTRIBUTES IN INTERPOLATION TO THE FINE-GRID POINTS */
/*     IF=JA(J) (IFG(IC)<=J<=IFG(IC+1)-1) WITH WEIGHT A(J). THE */
/*     CONTRIBUTION TO ITSELF (BY THE WEIGHT 1.0) IS NOT CONTAINED. */

/*     - ICG(IC) ----- IC=IMIN(K),IMAX(K) */

/*     IS USED FOR SEVERAL PURPOSES. IN ASSEMBLING THE CG-OPERATOR A(K), */
/*     IT SERVES AS POINTER TO POSITIONS IN A(K) WHICH HAVE ALREADY BEEN */
/*     DEFINED: IF THE CURRENT ROW CORRESPONDS TO POINT IC1, AND A */
/*     CONNECTION TO ANOTHER CG-POINT IC2 HAS JUST BEEN FOUND, THEN IF */
/*     ICG(IC2)<IA(IC1), THE CORRESPONDING ENTRY IN ROW IC1 OF A(K) */
/*     HAS NOT YET BEEN DEFINED. OTHERWISE, ICG(IC2) POINTS TO THE */
/*     LOCATION FOR THAT ENTRY. (ALSO SEE IFG BELOW.) */

/*     - IFG(IF) ----- IF=IMIN(K-1),IMAX(K-1) */

/*     IN ASSEMBLING THE CG-OPERATOR A(K), THIS VECTOR CONTAINS INFORMAT- */
/*     ION ON WHETHER THE EXISTENCE OF ENTRIES IN A(K) HAS TO BE CHECKED */
/*     OR NOT. */


/* ===> EXTEND A, JA TO STORE TRANSPOSE OF INTERPOLATION */

    /* Parameter adjustments */
    --time;
    --nstcol;
    --ncolor;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --a;

    /* Function Body */
	told=clock();
    ndaja = min(*nda,*ndja);
    iminw[*k] = imaxw[*k - 1] + 1;
    if (*k > *mmax) {
	goto L8000;
    }
    jpos = iw[iminw[*k]];
    i__1 = iw[imaxw[*k - 1] + 1] - 1;
    for (j = iw[iminw[*k - 1]]; j <= i__1; ++j) {
	++icg[ja[j]];
/* L200: */
    }

/* Computing MAX */
    i__1 = *mda, i__2 = jpos + jpos - iw[iminw[*k - 1]] + *iajas;
    *mda = max(i__1,i__2);
/* Computing MAX */
    i__1 = *mdja, i__2 = jpos + jpos - iw[iminw[*k - 1]] + *iajas;
    *mdja = max(i__1,i__2);
    ifg[imin[*k]] = ndaja - jpos + iw[iminw[*k - 1]] + 1;
    if (ifg[imin[*k]] <= jpos) {
	goto L9900;
    }
    i__1 = imax[*k];
    for (ic = imin[*k]; ic <= i__1; ++ic) {
	ifg[ic + 1] = ifg[ic] + icg[ic];
	icg[ic] = ifg[ic];
/* L220: */
    }

    i__1 = imax[*k - 1];
    for (if__ = imin[*k - 1]; if__ <= i__1; ++if__) {
	if (icg[if__] >= 0) {
	    goto L250;
	}
	ibl = -icg[if__];
	i__2 = iw[ibl + 1] - 1;
	for (j = iw[ibl]; j <= i__2; ++j) {
	    ic = ja[j];
	    a[icg[ic]] = a[j];
	    ja[icg[ic]] = if__;
	    ++icg[ic];
/* L260: */
	}
L250:
	;
    }

    i__1 = imax[*k];
    for (ic = imin[*k]; ic <= i__1; ++ic) {
	icg[ic] = 0;
/* L300: */
    }

/* ===> SWEEP OVER ALL CG-POINTS IC TO ASSEMBLE ROWS OF CG MATRIX */

    istti = 1;
    i__1 = imax[*k - 1];
    for (if__ = imin[*k - 1]; if__ <= i__1; ++if__) {
	if (icg[if__] <= 0) {
	    goto L100;
	}
	ic = icg[if__];
	if (jpos > ndaja) {
	    goto L9900;
	}
	ialow = jpos;
	icg[ic] = jpos;
	a[jpos] = a[ia[if__]];
	ja[jpos] = ic;
	++jpos;

/*       ---------------------------------------------- */
/*       | SEARCH FOR C-C-C-C AND C-C-F-C CONNECTIONS | */
/*       ---------------------------------------------- */

	i__2 = ia[if__ + 1] - 1;
	for (jf1 = ia[if__] + 1; jf1 <= i__2; ++jf1) {
	    if1 = ja[jf1];
	    ic1 = icg[if1];
	    if (ic1 < 0) {
		goto L11;
	    } else if (ic1 == 0) {
		goto L25;
	    } else {
		goto L20;
	    }

/* ===>     IF1 IS F-POINT: SWEEP OVER C-C-F-C CONNECTIONS */

L11:
	    ifg[if1] = -ic;
	    i__3 = iw[-ic1 + 1] - 1;
	    for (jf2 = iw[-ic1]; jf2 <= i__3; ++jf2) {
		ic2 = ja[jf2];
		if (icg[ic2] >= ialow) {
		    goto L10;
		}
		if (jpos > ndaja) {
		    goto L9900;
		}
		icg[ic2] = jpos;
		a[jpos] = a[jf1] * a[jf2];
		ja[jpos] = ic2;
		++jpos;
		goto L15;
L10:
		a[icg[ic2]] += a[jf1] * a[jf2];
L15:
		;
	    }
	    goto L25;

/* ===>     IF1 IS C-POINT: C-C-C-C CONNECTION */

L20:
	    if (icg[ic1] >= ialow) {
		goto L23;
	    }
	    if (jpos > ndaja) {
		goto L9900;
	    }
	    icg[ic1] = jpos;
	    a[jpos] = a[jf1];
	    ja[jpos] = ic1;
	    ++jpos;
	    goto L25;
L23:
	    a[icg[ic1]] += a[jf1];
L25:
	    ;
	}
	ist = imin[*k - 1] - 1;

/*        ---------------------------------------------- */
/*        | SEARCH FOR C-F-C-C AND C-F-F-C CONNECTIONS | */
/*        ---------------------------------------------- */

	i__2 = ifg[ic + 1] - 1;
	for (jf1 = ifg[ic]; jf1 <= i__2; ++jf1) {
	    if (jf1 >= jpos) {
		if1 = ja[jf1];
		wjf1 = a[jf1];
		ist = if1;
	    } else {
		istti = 0;
		i__3 = imax[*k - 1];
		for (jf2 = ist + 1; jf2 <= i__3; ++jf2) {
		    if (icg[jf2] >= 0) {
			goto L120;
		    }
		    ibl = -icg[jf2];
		    i__4 = iw[ibl + 1] - 1;
		    for (jb = iw[ibl]; jb <= i__4; ++jb) {
			jc3 = ja[jb];
			if (jc3 != ic) {
			    goto L110;
			}
			if1 = jf2;
			wjf1 = a[jb];
			goto L130;
L110:
			;
		    }
L120:
		    ;
		}

/* ===>       ERROR EXIT */

		//io___216.ciunit = *ium;
		//s_wsfe(&io___216);
		i__3 = *k - 1;
		//do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		//e_wsfe();

		printf("( *** ERROR IN OPDFN: INTERPOLATION ENTRY MISSING ON GRID,%d)\n",i__3);
		*ierr = 22;
		return 0;
L130:
		ist = if1;
	    }
	    i__3 = ia[if1 + 1] - 1;
	    for (jf2 = ia[if1]; jf2 <= i__3; ++jf2) {
		if2 = ja[jf2];
		ic2 = icg[if2];
		ww = wjf1 * a[jf2];
		if (ic2 < 0) {
		    goto L35;
		} else if (ic2 == 0) {
		    goto L90;
		} else {
		    goto L50;
		}

/* ===>       IF2 IS F-POINT: SWEEP OVER C-F-F-C CONNECTIONS */

L35:
		if (ifg[if2] == -ic) {
		    goto L70;
		}
		ifg[if2] = -ic;
		i__4 = iw[-ic2 + 1] - 1;
		for (jf3 = iw[-ic2]; jf3 <= i__4; ++jf3) {
		    ic3 = ja[jf3];
		    if (icg[ic3] >= ialow) {
			goto L30;
		    }
		    if (jpos > ndaja) {
			goto L9900;
		    }
		    icg[ic3] = jpos;
		    a[jpos] = ww * a[jf3];
		    ja[jpos] = ic3;
		    ++jpos;
		    goto L40;
L30:
		    a[icg[ic3]] += ww * a[jf3];
L40:
		    ;
		}
		goto L90;

/* ===>       IF2 HAS BEEN ENCOUNTERED BEFORE; DO NOT CHECK POSITIONS! */

L70:
		i__4 = iw[-ic2 + 1] - 1;
		for (jf3 = iw[-ic2]; jf3 <= i__4; ++jf3) {
		    iadrs = icg[ja[jf3]];
		    a[iadrs] += ww * a[jf3];
/* L80: */
		}
		goto L90;

/* ===>       IF2 IS C-POINT: C-F-C-C CONNECTION */

L50:
		if (icg[ic2] >= ialow) {
		    goto L60;
		}
		if (jpos > ndaja) {
		    goto L9900;
		}
		icg[ic2] = jpos;
		a[jpos] = ww;
		ja[jpos] = ic2;
		++jpos;
		goto L90;
L60:
		a[icg[ic2]] += ww;
L90:
		;
	    }
/* L95: */
	}
	ia[ic + 1] = jpos;
L100:
	;
    }
/* Computing MAX */
    i__1 = *mda, i__2 = ia[imax[*k] + 1] - 1 + *iajas;
    *mda = max(i__1,i__2);
/* Computing MAX */
    i__1 = *mdja, i__2 = ia[imax[*k] + 1] - 1 + *iajas;
    *mdja = max(i__1,i__2);

/* ===> WARNING FOR STORAGE SHORTAGE */

    if (istti != 1) {
	k1 = *k - 1;
	//io___223.ciunit = *ium;
	//s_wsfe(&io___223);
	//do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	//e_wsfe();

	    printf("( --- WARNG: UNABLE TO STORE TRANSPOSE OF \n");
	    printf("INTERPOLATION ON GRID ,%d, DURING     \n",k1);
	    printf("  EXECUTION OF OPDFN, BECAUSE NDA OR NDJA TOO SMALL.\n");
	    printf("  SETUP COMPUTATION IS SLOWING DOWN.)\n");

	*ierr = -1;
    }
    //io___224.ciunit = *ium;
    //s_wsfe(&io___224);
    //do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
    //e_wsfe();
	if (yes_print_amg) {
	    printf("( COARSE  GRID  OPERATOR NO.,%d, COMPLETED)\n", k[0]);
	}

/* ===> SET UP LINKED LIST FOR RELAXATION ON GRID K-1 */

L8000:
    ia[imin[*k]] = iw[iminw[*k]];
    ilo = imin[*k - 1];
    ihi = imax[*k - 1];
    ilo1 = ilo - 1;
    npts = ihi - ilo1;
    ist = 100000000;
    i__1 = *ncolx;
    for (icol = 1; icol <= i__1; ++icol) {
	for (i__ = npts; i__ >= 1; --i__) {
	    if (((integer)(ncolor[i__])) != icol) {
		goto L8040;
	    }
	    icg[i__ + ilo1] = -ist;
	    ist = i__ + ilo1;
L8040:
	    ;
	}
/* L8100: */
    }
    nstcol[*k - 1] = ist;

/* ===> EXIT / ERROR MESSAGES */

    //ctime_(&tnew);
	tnew=clock();
    time[6] = time[6] + tnew - told;
    return 0;

L9900:
    if (*nda <= *ndja) {
	//io___232.ciunit = *ium;
	//s_wsfe(&io___232);
	//e_wsfe();
	printf( "( *** ERROR IN OPDFN: NDA TOO SMALL ***)\n");
	*ierr = 1;
    } else {
	//io___233.ciunit = *ium;
	//s_wsfe(&io___233);
	//e_wsfe();
	printf("( *** ERROR IN OPDFN: NDJA TOO SMALL ***)\n");
	*ierr = 3;
    }
    return 0;

} /* opdfn_ */


/* ....................................................................... */

/*     SETIFG                                              SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int setifg_(integer *imin, integer *imax, integer *icg, 
	integer *ifg, integer *nstcol, integer *levels, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, ib, ist;
	static unsigned int told, tnew;


/*     SET "INVERSE" POINTER IFG */


    /* Parameter adjustments */
    --time;
    --nstcol;
    --ifg;
    --icg;
    --imax;
    --imin;

    /* Function Body */
	told=clock();
    i__1 = imax[*levels - 1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	if (icg[i__] > 0) {
	    ifg[icg[i__]] = i__;
	}
/* L10: */
    }
    ib = 1;
    i__1 = *levels - 1;
    for (k = 1; k <= i__1; ++k) {
	ist = nstcol[k];
L20:
	if (ist >= 100000000) {
	    goto L30;
	}
	ifg[ib] = ist;
	ist = -icg[ist];
	++ib;
	goto L20;
L30:
	;
    }
	tnew=clock();
    time[8] = time[8] + tnew - told;
    return 0;
} /* setifg_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5 SOLUTION-SUBROUTINES */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     SOLVE                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int solve_(integer *madapt, integer *ncyc, integer *nrd, 
	integer *nsolco, integer *nru, integer *iout, integer *ierr, 
	doublereal *a, doublereal *u, doublereal *f, integer *ia, integer *ja,
	 integer *iw, doublereal *eps, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *icg, integer *ifg, integer *nstcol, 
	integer *iarr, /*real*/unsigned int *time, integer *ncyc0, integer *irow0, integer *
	levels, integer *nda, integer *ndja, integer *ndu, integer *ndf, 
	integer *mda, integer *mdja, integer *mdu, integer *mdf, integer *iup,
	 integer *ium, doublereal *resi, doublereal *res0, doublereal *res)
{
    /* Format strings */
   // static char fmt_9050[] = "(\002 *** ERROR IN SOLVE: NDU TOO SMALL ***"
	//    "\002)";
    //static char fmt_9060[] = "(\002 *** ERROR IN SOLVE: NDF TOO SMALL ***"
	//    "\002)";
    //static char fmt_9005[] = "(/\002 ************* CYCLING..... **********"
	//    "***\002/)";
   // static char fmt_9000[] = "(\002 CYCLE  0:\002,3x,\002RES=\002,d9.3)";
    //static char fmt_9040[] = "(/\002 CYCLING BETWEEN GRIDS 1 AND\002,i3"
	  //  ",\002:\002/)";
   // static char fmt_9010[] = "(\002 CYCLE \002,i2,\002:   RESCG=\002,d9.3"
	//    ",\002   RES=\002,d9.3,\002   CFAC=\002,d9.3)";
    //static char fmt_9020[] = "(\002 CYCLE \002,i2,\002:   RES=\002,d9.3,\002"
	//    "   CFAC=\002,d9.3)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
   // integer s_wsfe(cilist *), e_wsfe(void), i_sign(integer *, integer *), 
	 //   do_fio(integer *, char *, ftnlen);
	integer i_sign(integer *, integer *);

    /* Local variables */
    static integer i__, l, m, n;
    extern /* Subroutine */ int cg_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/unsigned int *, 
	    integer *, integer *);
    static doublereal fac, ama;
    extern /* Subroutine */ int cyc_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    static integer nsc;
    static doublereal cfac;
    extern /* Subroutine */ int idec_(integer *, integer *, integer *, 
	    integer *);
    static integer igam, icgr;
    static doublereal fmax, epsi;
    static integer msel, iter, nrcx, nrdx;
    static doublereal umax;
    static integer nrux;
    static doublereal rescg;
    extern /* Subroutine */ int resid_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static doublereal epsil;
    static integer iconv;
    extern /* Subroutine */ int usave_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, /*real*/unsigned int *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer ncycle, ndigit, nrdlen;
    static doublereal resold;
    static integer nrulen, mfirst, nrdtyp[10], nrutyp[10];

    /* Fortran I/O blocks */
   // static cilist io___240 = { 0, 0, 0, fmt_9050, 0 };
    //static cilist io___241 = { 0, 0, 0, fmt_9060, 0 };
   // static cilist io___264 = { 0, 0, 0, fmt_9005, 0 };
   // static cilist io___265 = { 0, 0, 0, fmt_9000, 0 };
  //  static cilist io___269 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___270 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___273 = { 0, 0, 0, fmt_9010, 0 };
    //static cilist io___274 = { 0, 0, 0, fmt_9020, 0 };



/*     SOLUTION PHASE OF AMG1R5 */


/* ===> TEST OF AVAILABLE STORAGE */

    /* Parameter adjustments */
    --resi;
    --time;
    --iarr;
    --nstcol;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*ndu < imax[*levels]) {
	//io___240.ciunit = *ium;
	//s_wsfe(&io___240);
	//e_wsfe();
	printf("( *** ERROR IN SOLVE: NDU TOO SMALL ***)\n");
	*ierr = 4;
	return 0;
    }
    if (*ndf < imax[*levels]) {
	//io___241.ciunit = *ium;
	//s_wsfe(&io___241);
	//e_wsfe();
	printf("( *** ERROR IN SOLVE: NDF TOO SMALL ***)\n");
	*ierr = 5;
	return 0;
    }

    m = *levels;
    *ncyc0 = 0;
    for (n = 11; n <= 20; ++n) {
	//time[n] = 0.f;
		time[n]=0;
/* L5: */
    }
    if (*eps != 0.) {
	epsi = *eps;
    } else {
	epsi = 1e-12;
    }

/* ===> DECOMPOSE MADAPT */

    if (*madapt != 0) {
	idec_(madapt, &c__2, &ndigit, &iarr[1]);
	msel = iarr[1];
	if (msel == 2) {
	    if (iarr[2] != 0) {
		fac = (doublereal) iarr[2];
		for (i__ = 1; i__ <= 100; ++i__) {
		    fac /= 10.;
		    if (fac <= 1.) {
			goto L9;
		    }
/* L8: */
		}
	    } else {
		fac = .7;
	    }
	}
    } else {
	msel = 2;
	fac = .7;
    }

/* ===> DECOMPOSE NCYC */

L9:
    if (*ncyc != 0) {
	i__1 = abs(*ncyc);
	idec_(&i__1, &c__4, &ndigit, &iarr[1]);
	igam = i_sign(&iarr[1], ncyc);
	icgr = iarr[2];
	iconv = iarr[3];
	ncycle = iarr[4];
	if (ncycle == 0) {
	    return 0;
	}
    } else {
	igam = 1;
	icgr = 0;
	iconv = 1;
	ncycle = 10;
    }

/* ===> SET EPSI ACCORDING TO CONVERGENCE CRITERION GIVEN BY ICONV */

    if (iconv != 3) {
	if (iconv == 4) {
	    ama = 0.;
	    i__1 = imax[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__1 = ama, d__2 = a[ia[i__]];
		ama = max(d__1,d__2);
/* L6: */
	    }
	    epsi *= ama;
	}
    } else {
	fmax = 0.;
	i__1 = imax[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = fmax, d__3 = (d__1 = f[i__], fabs(d__1));
	    fmax = max(d__2,d__3);
/* L7: */
	}
	epsi *= fmax;
    }

/* ===> DECOMPOSE NRD */

    if (*nrd != 0) {
	idec_(nrd, &c__9, &ndigit, nrdtyp);
	nrdx = nrdtyp[1];
	nrdlen = ndigit - 2;
	i__1 = nrdlen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nrdtyp[i__ - 1] = nrdtyp[i__ + 1];
/* L10: */
	}
    } else {
	nrdx = 1;
	nrdlen = 2;
	nrdtyp[0] = 3;
	nrdtyp[1] = 1;
    }

/* ===> DECOMPOSE NRU */

    if (*nru != 0) {
	idec_(nru, &c__9, &ndigit, nrutyp);
	nrux = nrutyp[1];
	nrulen = ndigit - 2;
	i__1 = nrulen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nrutyp[i__ - 1] = nrutyp[i__ + 1];
/* L40: */
	}
    } else {
	nrux = 1;
	nrulen = 2;
	nrutyp[0] = 3;
	nrutyp[1] = 1;
    }

/* ===> DECOMPOSE NSOLCO */

    if (*nsolco != 0) {
	idec_(nsolco, &c__2, &ndigit, &iarr[1]);
	nsc = iarr[1];
	nrcx = iarr[2];

/* ===> IN CASE OF YALE-SMP COARSE GRID SOLUTION, DON'T USE COARSEST */
/*     GRID WITH LESS THAN 10 POINTS */

	if (nsc == 2) {
	    for (i__ = m; i__ >= 1; --i__) {
		l = i__;
		if (imax[i__] - imin[i__] >= 9) {
		    goto L60;
		}
/* L50: */
	    }
L60:
	    m = i__;
	    *levels = i__;
	}
    } else {
	nsc = 1;
	nrcx = 0;
    }

/* ===> CYCLING */

/* L100: */
    if (*iout != 0) {
	resid_(&c__1, res0, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &iminw[1]);
	if (*iout == 3) {
	    
		if (yes_print_amg) {
		     printf("( ************* CYCLING..... *************)\n");
		     printf("( CYCLE  0:,3x,RES=,%1.4f)\n",*res0);
		}
	}
	resold = *res0;
    }

    i__1 = ncycle;
    for (iter = 1; iter <= i__1; ++iter) {
	usave_(&c__1, &icgr, &u[1], &imin[1], &imax[1], ndu, &m, &time[1], 
		ierr, ium, mdu, ndf, mdf);
	cyc_(&c__1, &nrdx, nrdtyp, &nrdlen, &nrcx, &nrux, nrutyp, &nrulen, &
		igam, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &
		imax[1], &iminw[1], &imaxw[1], &ifg[1], &icg[1], &nstcol[1], &
		iarr[1], &time[1], irow0, &m, ium, ierr, &iter, &nsc, nda, 
		ndja, mda, mdja, &msel, &fac, &resi[1], levels);
	cg_(&c__1, &icgr, &iter, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], 
		&imin[1], &imax[1], &iminw[1], &m, &time[1], ierr, ium);
	if (*ierr > 0) {
	    return 0;
	}
	if (iter == 1) {
	    mfirst = m;
	    if (*iout == 3) {
		
			if (yes_print_amg) {
		         printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%d,)\n", m);
			}
	    }
	} else if (*iout == 3 && m != mfirst) {
	    mfirst = m;
	  
		if (yes_print_amg) {
                printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%d,)\n", m);
		}
		
		
	}
	if (*iout == 3 || iconv != 1) {
	    resid_(&c__1, res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		    imin[1], &imax[1], &iminw[1]);
	}
	*ncyc0 = iter;
	if (*iout != 3) {
	    goto L110;
	}
	cfac = *res / (resold + 1e-40);
	resold = *res;
	if (1 == m) {
	    goto L150;
	}
	resid_(&c__2, &rescg, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &iminw[1]);
	
	if (yes_print_amg) {
	     printf("( CYCLE ,%d,:   RESCG=,%1.4f\n", iter, rescg);
	     printf(",   RES=,%1.4f,   CFAC=,%1.4f)\n", *res, cfac);
	}
	goto L110;
L150:
	if (yes_print_amg) {
	     printf("( CYCLE ,%d,:   RES=,%1.4f,\n",iter,*res);
	     printf("   CFAC=,%1.4f)\n",cfac);
	}

L110:
	if (iconv == 1) {
	    goto L120;
	}
	epsil = epsi;
	if (iconv != 4) {
	    goto L115;
	}
	umax = 0.;
	i__2 = imax[1];
	for (i__ = imin[1]; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = umax, d__3 = (d__1 = u[i__], fabs(d__1));
	    umax = max(d__2,d__3);
/* L160: */
	}
	epsil = epsi * umax;
L115:
	if (*res < epsil) {
	    goto L170;
	}
L120:
	;
    }
L170:
    if (*iout != 3 && *iout != 0) {
	resid_(&c__1, res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[
		1], &imax[1], &iminw[1]);
    }
    return 0;

/* L9030: */
} /* solve_ */


/* ....................................................................... */

/*     CYC                                                    SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int cyc_(integer *l, integer *nrdx, integer *nrdtyp, integer 
	*nrdlen, integer *nrcx, integer *nrux, integer *nrutyp, integer *
	nrulen, integer *igam, doublereal *a, doublereal *u, doublereal *f, 
	integer *ia, integer *ja, integer *iw, integer *imin, integer *imax, 
	integer *iminw, integer *imaxw, integer *ifg, integer *icg, integer *
	nstcol, integer *ng, /*real*/ unsigned int *time, integer *irow0, integer *m, integer *
	ium, integer *ierr, integer *iter, integer *nsc, integer *nda, 
	integer *ndja, integer *mda, integer *mdja, integer *msel, doublereal 
	*fac, doublereal *resi, integer *levels)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, n, nl, ifi;
    static doublereal res;
    static integer ifac;
    extern /* Subroutine */ int inta_(integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, /*real*/unsigned int *), resc_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/unsigned int *);
	static unsigned int told;
    static integer mink;
    extern /* Subroutine */ int relx_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/ unsigned int *);
	static unsigned int tnew;
    extern /* Subroutine */ int nrmu_(integer *, doublereal *, integer *, 
	    integer *, integer *), putz_(integer *, doublereal *, integer *, 
	    integer *, /*real*/unsigned int *),/* ctime_(real *),*/ resid_(integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer nptsf;
    extern /* Subroutine */ int coarse_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *, 
	    integer *, integer *, integer *), vscale_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, /*real*/unsigned int *);
    static integer ivstar;


/*     PERFORMS ONE AMG CYCLE WITH GRID L AS FINEST GRID */



/* ===> DURING FIRST CYCLE: INITIALIZE PARAMETERS CONTROLLING YALE-SMP */
/*     FACTORIZATION AND ADAPTIVE DETERMINATION OF COARSEST GRID: */
/*     IFAC=1: ON NEXT CALL OF YALE-SMP FACTORIZE MATRIX */
/*     IFI =1: ON FIRST RETURN TO NEXT TO COARSEST GRID AFTER COARSE */
/*             GRID SOLUTION COMPARE RESIDUAL WITH THE RESIDUAL ON THE */
/*             SAME GRID BEFORE COARSE GRID SOLUTION. IF REDUCTION OF */
/*             RESIDUAL NOT SATISFYING, REDUCE NUMBER OF GRIDS USED IN */
/*             CYCLING BY ONE AND REPETE THE PROCESS WITH THE NOW */
/*             COARSEST GRID. */
/*     NPTSF:  NUMBER OF GRID POINTS ON FINEST GRID USED IN CYCLE, */
/*             DIVIDED BY 10. ONLY GRIDS WITH LESS THEN NPTSF POINTS */
/*             ARE ALLOWED TO BECOME COARSEST GRID DURING ADAPTIVE */
/*             COARSE GRID DETERMINATION. */

    /* Parameter adjustments */
    --resi;
    --time;
    --ng;
    --nstcol;
    --icg;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;
    --nrutyp;
    --nrdtyp;

    /* Function Body */
    if (*iter == 1) {
	ifac = 1;
    }
    if (*msel != 2) {
	ifi = 0;
	nptsf = 0;
    } else {
	mink = 1000000;
	ifi = 1;
	nptsf = (imax[*l] - imin[*l] + 1) / 10;
    }
    if (*l < *m) {
	goto L100;
    }

/* ===> ONE GRID ONLY */

    coarse_(m, &ifac, nsc, nrcx, ium, ierr, &a[1], &u[1], &f[1], &ia[1], &ja[
	    1], &iw[1], &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
	    time[1], nda, ndja, mda, mdja, irow0);
    goto L1000;

/* ===> MORE THAN ONE GRID */

L100:
    i__1 = *m;
    for (k = *l; k <= i__1; ++k) {
	ng[k] = 0;
/* L110: */
    }
    ivstar = 3 - *igam;
    k = *l;

/* ===> RELAX (DOWNWARDS) */

L150:
    i__1 = *nrdx;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *nrdlen;
	for (nl = 1; nl <= i__2; ++nl) {
	    relx_(&k, &nrdtyp[nl], &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1]
		    , &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
		    time[1]);
/* L160: */
	}
/* L170: */
    }
    if (ifi != 1 || imax[k] - imin[k] >= nptsf) {
	goto L190;
    }

/* ===> MINK: LOWEST GRID NUMBER FOR WHICH RESIDUAL IS STORED DURING */
/*           FIRST DOWNWARDS RELAXATION */

    if (mink == 1000000) {
	mink = k;
    }
	told=clock();
    resid_(&k, &resi[k], &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1]);
	tnew=clock();
    time[15] = time[15] + tnew - told;
L190:
    ++ng[k];
    ++k;
L195:
    putz_(&k, &u[1], &imin[1], &imax[1], &time[1]);
    resc_(&k, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &imax[1],
	     &iminw[1], &imaxw[1], &ifg[1], &time[1]);
    if (k < *m) {
	goto L150;
    }

/* ===> SOLVE ON COARSEST GRID */

    coarse_(m, &ifac, nsc, nrcx, ium, ierr, &a[1], &u[1], &f[1], &ia[1], &ja[
	    1], &iw[1], &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
	    time[1], nda, ndja, mda, mdja, irow0);

/* ===> RELAX (UPWARDS) */

L200:
    vscale_(&k, &ivstar, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1], &time[1]);
    --k;
    inta_(&k, &a[1], &u[1], &ia[1], &ja[1], &iw[1], &imin[1], &imax[1], &
	    iminw[1], &imaxw[1], &ifg[1], &time[1]);
    i__1 = *nrux;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *nrulen;
	for (nl = 1; nl <= i__2; ++nl) {
	    relx_(&k, &nrutyp[nl], &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1]
		    , &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
		    time[1]);
/* L210: */
	}
/* L215: */
    }
    if (ifi != 1 || k < mink) {
	goto L219;
    }

/* ===> ON FIRST RETURN TO NEXT TO COARSEST GRID COMPARE RESIDUAL WITH */
/*     THE PREVIOUS ONE */

	told=clock();
    resid_(&k, &res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &
	    imax[1], &iminw[1]);
	tnew=clock();
    time[15] = time[15] + tnew - told;

/* ===> IF RESIDUAL REDUCTION SATISFYING: COARSE GRID ADAPTION FINISHED */

    if (res < resi[k] * *fac) {
	ifi = 0;
    } else {
	if (*nsc == 2) {
	    *levels = k;
	}
	*m = k;
	ifac = 1;
	goto L195;
    }
L219:
    if (k == *l) {
	goto L1000;
    }

/* ===> GRID SWITCHING CORRESPONDING TO IGAM */

    if (*igam >= 3) {
	goto L220;
    }
    if (*igam >= 0) {
	goto L200;
    }
    if (k == *l + 1 && ng[k] < abs(*igam)) {
	goto L150;
    }
    goto L200;
L220:
    if (ng[k] < 2) {
	goto L150;
    }
    if (*igam == 4) {
	ng[k] = 0;
    }
    goto L200;

/* ===> RETURN */

L1000:
    nrmu_(l, &u[1], &imin[1], &imax[1], irow0);
    return 0;
} /* cyc_ */


/* ....................................................................... */

/*     COARSE                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int coarse_(integer *m, integer *ifac, integer *nsc, integer 
	*nrcx, integer *ium, integer *ierr, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *icg, integer *nstcol, /*real*/ unsigned int *
	time, integer *nda, integer *ndja, integer *mda, integer *mdja, 
	integer *irow0)
{
    /* Format strings */
    //static char fmt_9020[] = "(\002 --- WARNG IN COARSE: NO YALE-SMP BECAUSE"
	//    " NDJA TOO \002,\002SMALL\002)";
    //static char fmt_9030[] = "(\002 --- WARNG IN COARSE: NO YALE-SMP BECAUSE"
	//    " NDA TOO SMALL\002)";
    //static char fmt_9040[] = "(\002 --- WARNG IN COARSE: NO YALE-SMP BECAUSE"
	//    " OF ERROR IN \002,\002FACTORIZATION,\002/\002     CODE=\002,i8)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;
	


    /* Builtin functions */
    //integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, ii, jj, is, js, np, ihi, jhi, ilo, jlo, esp, nsp, 
	    flag__;
    static doublereal fmax;
    static integer path;
    static doublereal aaux;
	static unsigned int told;
    static integer iter, iaux;
    extern /* Subroutine */ int ndrv_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal /*integer*/ *, doublereal *, integer *, 
	    integer *, integer *);
	 extern /* Subroutine */ int relx_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/ unsigned int *);
    static integer jpos;
	static unsigned int tnew;
    extern /* Subroutine */ int  resid_(integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static doublereal resold, resnew;
    static integer npoint;

    /* Fortran I/O blocks */
    //static cilist io___296 = { 0, 0, 0, fmt_9020, 0 };
    //static cilist io___309 = { 0, 0, 0, fmt_9030, 0 };
   // static cilist io___310 = { 0, 0, 0, fmt_9040, 0 };



/*     SOLVES ON COARSEST GRID, EITHER WITH GAUSS-SEIDEL RELAXATION */
/*     (NSC=1) OR WITH THE YALE-SMP DIRECT SOLVER NDRV (NSC=2) */


/* ===> CONV: IF COARSE GRID SOLUTION IS DONE WITH GS-RELAXATION AND */
/*     NRCX=0, AS MANY GS-SWEEPS ARE PERFORMED AS ARE NECESSARY TO RE- */
/*     DUCE THE RESIDUAL BY THE FACTOR CONV */


    /* Parameter adjustments */
    --time;
    --nstcol;
    --icg;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*nsc != 2) {
	goto L400;
    }

/* ===> SOLUTION WITH YALE-SMP */

	told=clock();
    ilo = imin[*m];
    jlo = ia[ilo];
    if (*ifac == 1) {

/* ===> FIRST CALL ON GRID M, FIRST FACTORIZE MATRIX */

	ihi = imax[*m];
	jhi = iw[iminw[*m]] - 1;
	np = ihi - ilo + 1;
	is = ilo - 1;
	js = jlo - 1;

/* ===>   TEST OF AVAILABLE WORK SPACE */

	if (jhi + np * 3 > *ndja) {
	    *nsc = 1;
	    //io___296.ciunit = *ium;
	    //s_wsfe(&io___296);
	    //e_wsfe();

		printf("( --- WARNG IN COARSE: NO YALE-SMP BECAUSE\n");
	    printf(" NDJA TOO ,SMALL)\n");

	    *ierr = -3;
	    goto L400;
	}

/* ===>   INITIALISATION OF YALE-SMP POINTER VECTORS */

	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ja[jhi + i__] = i__;
	    ja[jhi + np + i__] = i__;
	    ja[jhi + (np << 1) + i__] = i__;
/* L310: */
	}
	if (*irow0 != 1) {

/* ===>   COARSE GRID OPERATOR REGULAR, SHIFT CONTENTS OF POINTER */
/*       VECTORS IA AND JA */

	    i__1 = ihi;
	    for (i__ = ilo; i__ <= i__1; ++i__) {
		ia[i__] -= js;
/* L320: */
	    }
	    iaux = ia[ihi + 1];
	    ia[ihi + 1] = iw[iminw[*m]] - js;
	    i__1 = jhi;
	    for (i__ = jlo; i__ <= i__1; ++i__) {
		ja[i__] -= is;
/* L300: */
	    }
	    npoint = np;
	} else {

/* ===>   COARSE GRID OPERATOR HAS ROWSUM ZERO, ELIMINATE LAST SOLUTION */
/*       COMPONENT BY SETTING U(IHI) TO ZERO AND CANCELLING THE COR- */
/*       RESPONDING ENTRIES IN A AND JA, RESPECTIVELY. THE CANCELLED */
/*       ENTRIES ARE STORED BEFORE THE LAST ROW IN A AND JA. */

	    iaux = ia[ihi];

/* ===>     JPOS: POINTER TO POSITION IN A AND JA TO CONTAIN NEXT */
/*               ELIMINATED ENTRY */

	    jpos = iaux - 1;
	    j = ia[ilo];
	    i__1 = ihi - 1;
	    for (i__ = ilo; i__ <= i__1; ++i__) {
		ia[i__] = j - js;
		ja[j] = i__ - is;
		++j;
L13:
		if (j == ia[i__ + 1]) {
		    goto L20;
		}
		if (ja[j] != ihi) {
		    ja[j] -= is;
		    ++j;
		} else {
		    aaux = a[j];
		    i__2 = jpos - 1;
		    for (jj = j; jj <= i__2; ++jj) {
			a[jj] = a[jj + 1];
			ja[jj] = ja[jj + 1];
/* L11: */
		    }
		    i__2 = ihi;
		    for (ii = i__ + 1; ii <= i__2; ++ii) {
			--ia[ii];
/* L12: */
		    }
		    a[jpos] = aaux;
		    ja[jpos] = i__;
		    --jpos;
		}
		goto L13;
L20:
		;
	    }
	    ia[ihi] -= js;

/* ===>     DECREASE NUMBER OF POINTS BY ONE AND SET LAST SOLUTION */
/*         COMPONENT TO ZERO */

	    npoint = np - 1;
	    u[ihi] = 0.;
	}
	nsp = *nda - jhi;
	path = 1;

	


	ndrv_(&npoint, &ja[jhi + 1], &ja[jhi + np + 1], &ja[jhi + (np << 1) + 
		1], &ia[ilo], &ja[jlo], &a[jlo], &f[ilo], &u[ilo], &nsp, 
		&a[jhi + 1], &a[jhi + 1], &esp, &path, &flag__);

		

/* ===>   RESTORE PREVIOUS VALUES FOR IA, JA AND A, IN PARTICULAR PUT */
/*       BACK ELIMINATED MATRIX ENTRIES, IF ROWSUM ZERO */

	i__1 = ihi;
	for (i__ = ilo; i__ <= i__1; ++i__) {
	    ia[i__] += js;
/* L330: */
	}
	if (*irow0 != 1) {
	    ia[ihi + 1] = iaux;
	    i__1 = jhi;
	    for (i__ = jlo; i__ <= i__1; ++i__) {
		ja[i__] += is;
/* L335: */
	    }
	} else {
	    i__1 = jpos;
	    for (i__ = jlo; i__ <= i__1; ++i__) {
		ja[i__] += is;
/* L25: */
	    }
	    i__1 = iaux - 1;
	    for (j = jpos + 1; j <= i__1; ++j) {
		aaux = a[j];
		i__ = ja[j];
		i__2 = ihi;
		for (ii = i__ + 1; ii <= i__2; ++ii) {
		    ++ia[ii];
/* L27: */
		}
		i__2 = ia[i__ + 1];
		for (jj = j; jj >= i__2; --jj) {
		    a[jj] = a[jj - 1];
		    ja[jj] = ja[jj - 1];
/* L30: */
		}
		a[ia[i__ + 1] - 1] = aaux;
		ja[ia[i__ + 1] - 1] = ihi;
/* L40: */
	    }
	}

/* ===>   IF AN ERROR OCCURED DURING EXECUTION OF NDRV, SOLVE WITH */
/*       GAUSS-SEIDEL RELAXATION */

	if (flag__ != 0) {
	    *nsc = 1;
	    if (esp < 0) {
		

		printf("( --- WARNG IN COARSE: NO YALE-SMP BECAUSE\n");
	    printf(" NDA TOO SMALL)\n");

		*ierr = -1;
	    } else {
		

		printf("( --- WARNG IN COARSE: NO YALE-SMP BECAUSE");
	    printf(" OF ERROR IN ,FACTORIZATION,     CODE=,%d)", flag__);

		*ierr = -32;
	    }
	    goto L400;
	} else {

/* ===>     FACTORIZATION SUCCESSFULL, UPDATE LIMITS OF USED STORAGE */

/* Computing MAX */
	    i__1 = *mda, i__2 = *nda - esp;
	    *mda = max(i__1,i__2);
/* Computing MAX */
	    i__1 = *mdja, i__2 = jhi + np * 3;
	    *mdja = max(i__1,i__2);
	    *ifac = 0;
	}
    } else {

/* ===>   FACTORIZATION ALLREADY DONE */

	if (*irow0 != 1) {
	    npoint = np;
	} else {
	    npoint = np - 1;
	    u[ihi] = 0.;
	}
	path = 3;

	

	ndrv_(&npoint, &ja[jhi + 1], &ja[jhi + np + 1], &ja[jhi + (np << 1) + 
		1], &ia[ilo], &ja[jlo], &a[jlo], &f[ilo], &u[ilo], &nsp, &a[jhi + 1], 
		&a[jhi + 1], &esp, &path, &flag__);
    }

	

/* ===> UPDATE TIME COUNTER */

	tnew=clock();
    time[17] = time[17] + tnew - told;
    goto L190;

/* ===> SOLUTION WITH GAUSS-SEIDEL RELAXATION */

L400:
    if (*nrcx != 0) {
	i__1 = *nrcx;
	for (iter = 1; iter <= i__1; ++iter) {
	    relx_(m, &c__2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		    imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &time[
		    1]);
/* L180: */
	}
    } else {

/* ===> IF NRCX=0: REDUCE RESIDUAL ON COARSEST GRID BY FACTOR CONV */
/*                (IF NOT YET IN THE RANGE OF THE TRUNCATION ERROR) */

	told=clock();

/* ===>   CALCULATE SUPREMUM NORM OF RIGHT HAND SIDE */

	fmax = 0.;
	i__1 = imax[*m];
	for (i__ = imin[*m]; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = fmax, d__3 = (d__1 = f[i__], fabs(d__1));
	    fmax = max(d__2,d__3);
/* L181: */
	}
	resid_(m, &resold, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[
		1], &imax[1], &iminw[1]);
	tnew=clock();
	time[15] = time[15] + tnew - told;
/* Computing MAX */
	d__1 = resold * .01, d__2 = fmax * 1e-12;
	resold = max(d__1,d__2);
	for (i__ = 1; i__ <= 10; ++i__) {
	    for (j = 1; j <= 10; ++j) {
		relx_(m, &c__2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
			imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
			time[1]);
/* L185: */
	    }
		told=clock();
	    resid_(m, &resnew, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		    imin[1], &imax[1], &iminw[1]);
		tnew=clock();
	    time[15] = time[15] + tnew - told;
	    if (resnew <= resold) {
		goto L190;
	    }
/* L187: */
	}
    }

/* ===> COMPUTE RESIDUAL AFTER SOLUTION ON GRID M */

L190:
    return 0;

} /* coarse_ */


/* ....................................................................... */

/*     RELX                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int relx_(integer *k, integer *irel, doublereal *a, 
	doublereal *u, doublereal *f, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, integer *icg, integer *
	nstcol, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal s;
	static unsigned int told;
    static integer iaux;
	static unsigned int tnew;


/*     PERFORMS ONE (PARTIAL) GAUSS-SEIDEL SWEEP ON GRIK K: */

/*     IREL = 1 :   PARTIAL GAUSS-SEIDEL SWEEP (ONLY F-POINTS) */
/*          = 2 :   FULL GAUSS-SEIDEL SWEEP (ALL POINTS) */
/*          = 3 :   PARTIAL GAUSS-SEIDEL SWEEP (ONLY C-POINTS) */
/*          = 4 :   FULL SWEEP: FF -- C -- COLORS (HIGHEST FIRST) */


    /* Parameter adjustments */
    --time;
    --nstcol;
    --icg;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
	told=clock();
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    switch (*irel) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
    }
    goto L200;

/* ===> F-RELAXATION */

L100:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] > 0) {
	    goto L120;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L110: */
	}
	u[i__] = s / a[ia[i__]];
L120:
	;
    }
    goto L1000;

/* ===> FULL GS RELAXATION */

L200:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L210: */
	}
	u[i__] = s / a[ia[i__]];
/* L220: */
    }
    goto L1000;

/* ===> C-RELAXATION */

L300:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] <= 0) {
	    goto L320;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L310: */
	}
	u[i__] = s / a[ia[i__]];
L320:
	;
    }
    goto L1000;

/* ===> FF-RELAXATION */

L400:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] != 0) {
	    goto L420;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L410: */
	}
	u[i__] = s / a[ia[i__]];
L420:
	;
    }

/* ===> C-RELAXATION */

    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] <= 0) {
	    goto L440;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L430: */
	}
	u[i__] = s / a[ia[i__]];
L440:
	;
    }

/* ===> MORE-COLOR RELAXATION */

    i__ = nstcol[*k];
L470:
    if (i__ >= 100000000) {
	goto L1000;
    }
    s = f[i__];
    i__1 = ia[i__ + 1] - 1;
    for (j = ia[i__] + 1; j <= i__1; ++j) {
	s -= a[j] * u[ja[j]];
/* L480: */
    }
    u[i__] = s / a[ia[i__]];
    i__ = -icg[i__];
    goto L470;

L1000:
	tnew=clock();
    time[13] = time[13] + tnew - told;
    ia[imax[*k] + 1] = iaux;
    return 0;
} /* relx_ */


/* ....................................................................... */

/*     VSCALE                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int vscale_(integer *k, integer *ivstar, doublereal *a, 
	doublereal *u, doublereal *f, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, /*real*/unsigned int  *time)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal s1, s2, sa, fac;
	static unsigned int told;
    static integer iaux;
	static unsigned int tnew;
   


/*     SCALES ACTUAL APPROXIMATE SOLUTION ON GRID K (V*-CYCLE); SCALING */
/*     IS DONE SUCH THAT ENERGY NORM BECOMES MINIMAL */

/*     NOTE: THIS SCALING MAKES SENSE ONLY FOR SYMMETRIC PROBLEMS */


/* ===> COMPUTATION OF SCALING FACTOR "FAC" */

    /* Parameter adjustments */
    --time;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*ivstar != 1) {
	return 0;
    }
   
	told=clock();
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    s1 = 0.;
    s2 = 0.;
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	sa = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    sa += a[j] * u[ja[j]];
/* L20: */
	}
	s1 += u[i__] * f[i__];
	s2 += u[i__] * sa;
/* L10: */
    }
    fac = 1.;
    if (s2 != 0.) {
	fac = s1 / s2;
    }

/* ===> SCALING */

    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] *= fac;
/* L30: */
    }
    ia[imax[*k] + 1] = iaux;
	tnew=clock();
    time[14] = time[14] + tnew - told;
    return 0;
} /* vscale_ */


/* ....................................................................... */

/*     INTA                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int inta_(integer *kf, doublereal *a, doublereal *u, integer 
	*ia, integer *ja, integer *iw, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *ifg, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ic, if__;
	static unsigned int told;
    static integer iaux;
	static unsigned int tnew;


/*     INTERPOLATES CORRECTION FROM GRID KF+1 TO GRID KF */


/* ===> C->C CONTRIBUTIONS */

    /* Parameter adjustments */
    --time;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --u;
    --a;

    /* Function Body */
	told=clock();
    i__1 = imax[*kf + 1];
    for (ic = imin[*kf + 1]; ic <= i__1; ++ic) {
	if__ = ifg[ic];
	u[if__] += u[ic];
/* L50: */
    }

/* ===> C->F CONTRIBUTIONS */

    iaux = iw[imaxw[*kf] + 1];
    iw[imaxw[*kf] + 1] = ia[imin[*kf + 1]];
    i__1 = imaxw[*kf];
    for (i__ = iminw[*kf]; i__ <= i__1; ++i__) {
	if__ = ifg[i__];
	i__2 = iw[i__ + 1] - 1;
	for (j = iw[i__]; j <= i__2; ++j) {
	    u[if__] += a[j] * u[ja[j]];
/* L150: */
	}
/* L100: */
    }
    iw[imaxw[*kf] + 1] = iaux;
	tnew=clock();
    time[11] = time[11] + tnew - told;
    return 0;
} /* inta_ */


/* ....................................................................... */

/*     RESC                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int resc_(integer *kc, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *imaxw, integer *ifg, /*real*/ unsigned int *
	time)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, ic, if__;
	static unsigned int told;
    static integer iaux;
	static unsigned int tnew;
    static integer iaux1;


/*     RESTRICTS RESIDUALS FROM GRID KC-1 TO GRID KC */


/* ===> TRANSFER OF C-POINT DEFECTS */

    /* Parameter adjustments */
    --time;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
	told=clock();
    iaux = ia[imax[*kc - 1] + 1];
    ia[imax[*kc - 1] + 1] = iw[iminw[*kc - 1]];
    iaux1 = iw[imaxw[*kc - 1] + 1];
    iw[imaxw[*kc - 1] + 1] = iaux;
    i__1 = imax[*kc];
    for (ic = imin[*kc]; ic <= i__1; ++ic) {
	if__ = ifg[ic];
	d__ = f[if__];
	i__2 = ia[if__ + 1] - 1;
	for (j = ia[if__]; j <= i__2; ++j) {
	    d__ -= a[j] * u[ja[j]];
/* L80: */
	}
	f[ic] = d__;
/* L100: */
    }

/* ===> TRANSFER OF F-POINT DEFECTS */

    i__1 = imaxw[*kc - 1];
    for (i__ = iminw[*kc - 1]; i__ <= i__1; ++i__) {
	if__ = ifg[i__];
	d__ = f[if__];
	i__2 = ia[if__ + 1] - 1;
	for (j = ia[if__]; j <= i__2; ++j) {
	    d__ -= a[j] * u[ja[j]];
/* L20: */
	}
	i__2 = iw[i__ + 1] - 1;
	for (j = iw[i__]; j <= i__2; ++j) {
	    f[ja[j]] += a[j] * d__;
/* L250: */
	}
/* L200: */
    }
    ia[imax[*kc - 1] + 1] = iaux;
    iw[imaxw[*kc - 1] + 1] = iaux1;
	tnew=clock();
    time[12] = time[12] + tnew - told;
    return 0;
} /* resc_ */


/* ....................................................................... */

/*     PUTZ                                                   SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int putz_(integer *k, doublereal *u, integer *imin, integer *
	imax, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
	static unsigned int told, tnew;


/*     PUTS ZERO TO U-VALUES OF GRID K */


    /* Parameter adjustments */
    --time;
    --imax;
    --imin;
    --u;

    /* Function Body */
	told=clock();
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] = 0.;
/* L10: */
    }
	tnew=clock();
    time[15] = time[15] + tnew - told;
    return 0;
} /* putz_ */


/* ....................................................................... */

/*     FIRST                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int first_(integer *ifirst, doublereal *u, integer *imin, 
	integer *imax, integer *iarr, integer *irow0)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal s;
    extern /* Subroutine */ int idec_(integer *, integer *, integer *, 
	    integer *);
    static integer ifrst, ndigit;
    extern doublereal random_(doublereal *);


/*     PUTS A FIRST APPROXIMATION TO FINEST GRID */


    /* Parameter adjustments */
    --iarr;
    --imax;
    --imin;
    --u;

    /* Function Body */
    if (*ifirst != 0) {
	idec_(ifirst, &c__3, &ndigit, &iarr[1]);
	ifrst = iarr[2];
    } else {
	ifrst = 3;
    }

    switch (ifrst) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
    }
    return 0;

L100:
    i__1 = imax[1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	u[i__] = 0.;
/* L110: */
    }
    return 0;

L200:
    i__1 = imax[1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	u[i__] = 1.;
/* L210: */
    }
    if (*irow0 < 2) {
	u[imax[1]] = 0.;
    }
    return 0;

L300:
    if (iarr[3] * *ifirst == 0) {
	goto L350;
    }
    s = (doublereal) iarr[3];
    for (i__ = 1; i__ <= 10; ++i__) {
	s *= .1;
	if (s < 1.) {
	    goto L370;
	}
/* L310: */
    }
L350:
    s = .72815;
L370:
    i__1 = imax[1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	u[i__] = random_(&s);
/* L390: */
    }
    if (*irow0 < 2) {
	u[imax[1]] = 0.;
    }
    return 0;
} /* first_ */


/* ....................................................................... */

/*     INJF                                                   SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int injf_(integer *kc, doublereal *f, integer *imin, integer 
	*imax, integer *ifg)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ic;


/*     INJECTS F-VALUES FROM GRID KC-1 TO GRID KC */

    /* Parameter adjustments */
    --ifg;
    --imax;
    --imin;
    --f;

    /* Function Body */
    i__1 = imax[*kc];
    for (ic = imin[*kc]; ic <= i__1; ++ic) {
	f[ic] = f[ifg[ic]];
/* L10: */
    }
    return 0;
} /* injf_ */


/* ....................................................................... */

/*     INJU                                                   SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int inju_(integer *kc, doublereal *u, integer *imin, integer 
	*imax, integer *ifg)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ic;


/*     INJECTS U-VALUES FROM GRID KC-1 TO GRID KC */

    /* Parameter adjustments */
    --ifg;
    --imax;
    --imin;
    --u;

    /* Function Body */
    i__1 = imax[*kc];
    for (ic = imin[*kc]; ic <= i__1; ++ic) {
	u[ic] = u[ifg[ic]];
/* L10: */
    }
    return 0;
} /* inju_ */


/* ....................................................................... */

/*     NRMU                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int nrmu_(integer *k, doublereal *u, integer *imin, integer *
	imax, integer *irow0)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal fac;


/*     NORMALIZES U ON GRID K IF ROWSUM=0 (LAST COMPONENT =0) */


    /* Parameter adjustments */
    --imax;
    --imin;
    --u;

    /* Function Body */
    if (*irow0 > 1) {
	return 0;
    }
    fac = u[imax[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] -= fac;
/* L10: */
    }
    return 0;
} /* nrmu_ */


/* ....................................................................... */

/*     RESID                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int resid_(integer *k, doublereal *res, doublereal *a, 
	doublereal *u, doublereal *f, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static integer iaux;


/*     COMPUTES L2-NORM OF RESIDUAL ON GRID K */


    /* Parameter adjustments */
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    *res = 0.;
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L10: */
	}
	*res += s * s;
/* L20: */
    }
    ia[imax[*k] + 1] = iaux;
    *res = sqrt(*res);
    return 0;
} /* resid_ */


/* ....................................................................... */

/*     CG                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int cg_(integer *k, integer *icgr, integer *iter, doublereal 
	*a, doublereal *u, doublereal *f, integer *ia, integer *ja, integer *
	iw, integer *imin, integer *imax, integer *iminw, integer *m, /*real*/ unsigned int *
	time, integer *ierr, integer *ium)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal s2, alf, eps;
    static integer nnu;
	static unsigned int told, tnew;
    extern doublereal cgalf_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), cgeps_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *);
    static integer ishift;


/*     PERFORMS ONE STEP OF PRECONDITIONED CONJUGATE GRADIENT */


    /* Parameter adjustments */
    --time;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*icgr == 0) {
	return 0;
    }

/* ===> COMPUTE MOST RECENT MG CORRECTION */

	told=clock();
    nnu = imax[*k] - imin[*k] + 1;
    ishift = imax[*m] + 1 - imin[*k];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] -= u[i__ + ishift];
/* L10: */
    }
    if (*icgr == 2 && *iter > 1) {
	goto L100;
    }

/* ===> FIRST CG STEP */

    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	f[i__ + ishift] = u[i__];
/* L20: */
    }
    goto L200;

/* ===> NEXT CG STEPS (IF ICGR=2 ONLY) */

L100:
    alf = cgalf_(k, &s2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1], m);
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	f[i__ + ishift] = u[i__] + alf * f[i__ + ishift];
/* L110: */
    }
L200:
    eps = cgeps_(k, &s2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1], m, ierr, ium);
    if (*ierr > 0) {
	return 0;
    }
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] = u[i__ + ishift] + eps * f[i__ + ishift];
/* L210: */
    }
	tnew=clock();
    time[16] = time[16] + tnew - told;
    return 0;
} /* cg_ */


/* ....................................................................... */

/*     USAVE                                            SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int usave_(integer *k, integer *icgr, doublereal *u, integer 
	*imin, integer *imax, integer *ndu, integer *m, /*real*/ unsigned int *time, integer *
	ierr, integer *ium, integer *mdu, integer *ndf, integer *mdf)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 --- WARNG IN USAVE: NO CG BECAUSE OF STO"
	//    "RAGE ---\002/\002     REQUIRED: NDU =\002,i9/\002               "
	 //   "NDF =\002,i9)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    //integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;
	static unsigned int told, tnew;
    static integer ishift;

    /* Fortran I/O blocks */
    //static cilist io___371 = { 0, 0, 0, fmt_9000, 0 };



/*     MAKES A BACK-UP OF THE CURRENT APPROXIMATION ON LEVEL K */
/*     IF ICGR.NE.0. */


    /* Parameter adjustments */
    --time;
    --imax;
    --imin;
    --u;

    /* Function Body */
    if (*icgr == 0) {
	return 0;
    }
    ishift = imax[*m] + 1 - imin[*k];
/* Computing MAX */
    i__1 = *mdu, i__2 = imax[*k] + ishift;
    *mdu = max(i__1,i__2);
/* Computing MAX */
    i__1 = *mdf, i__2 = imax[*k] + ishift;
    *mdf = max(i__1,i__2);
    if (*mdu <= *ndu && *mdf <= *ndf) {
	goto L10;
    }
   

	printf("( --- WARNG IN USAVE: NO CG BECAUSE OF STORAGE\n");
	printf(" ---    REQUIRED: NDU =,%d,      \n",mdu[0]);
	printf( "NDF =,%d)\n",mdf[0]);
	    
	   

    *ierr = -4;
    *icgr = 0;
    return 0;

L10:
	told=clock();
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__ + ishift] = u[i__];
/* L20: */
    }
	tnew=clock();
    time[16] = time[16] + tnew - told;
    return 0;

} /* usave_ */


/* ....................................................................... */

/*     CGEPS                                            FUNCTION */

/* ....................................................................... */

doublereal cgeps_(integer *k, doublereal *s2, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *m, integer *ierr, integer *
	ium)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 *** ERROR IN CGEPS: CG CORRECTION NOT DE"
	//    "FINED ***\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val=0.0; // TODO инициализация была добавлена позже мной её не было изначально.

    /* Builtin functions */
   // integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal s1, sp, sr;
    static integer iaux, ishift;

    /* Fortran I/O blocks */
    //static cilist io___382 = { 0, 0, 0, fmt_9000, 0 };




    /* Parameter adjustments */
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    ishift = imax[*m] + 1 - imin[*k];
    s1 = 0.;
    *s2 = 0.;
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	sr = f[i__];
	sp = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    sr -= a[j] * u[ja[j] + ishift];
	    sp += a[j] * f[ja[j] + ishift];
/* L40: */
	}
	s1 += sr * f[i__ + ishift];
	*s2 += sp * f[i__ + ishift];
/* L50: */
    }
    ia[imax[*k] + 1] = iaux;
    if (*s2 == 0.) {
	goto L100;
    }
    ret_val = s1 / *s2;
    return ret_val;

/* ===> ERROR EXIT */

L100:
    //io___382.ciunit = *ium;
    //s_wsfe(&io___382);
    //e_wsfe();
	printf("( *** ERROR IN CGEPS: CG CORRECTION NOT DEFINED ***)\n");
    *ierr = 31;
    return ret_val;
} /* cgeps_ */


/* ....................................................................... */

/*     CGALF                                            FUNCTION */

/* ....................................................................... */

doublereal cgalf_(integer *k, doublereal *s2, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *m)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j;
    static doublereal s1, sr;
    static integer iaux, ishift;



    /* Parameter adjustments */
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    ishift = imax[*m] + 1 - imin[*k];
    s1 = 0.;
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	sr = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    sr += a[j] * u[ja[j]];
/* L40: */
	}
	s1 += sr * f[i__ + ishift];
/* L50: */
    }
    ret_val = -s1 / *s2;
    ia[imax[*k] + 1] = iaux;
    return ret_val;
} /* cgalf_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     OUTPUT OF STATISTICAL INFORMATION ON CP-TIMES AND DIMENSIONING */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     WRKCNT                                           SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int wrkcnt_(integer *iout, integer *ia, integer *iw, integer 
	*imin, integer *imax, integer *iminw, integer *levels, /*real*/ unsigned int *time, 
	integer *ncyc0, integer *iup, integer *mda, integer *mdia, integer *
	mdja, integer *mdu, integer *mdf, integer *mdig, doublereal *res0, 
	doublereal *res)
{
    /* Format strings */
    //static char fmt_9110[] = "(//\002 **************** CONVERGENCE *********"
	//    "********\002/\002 L2-NORM OF RESIDUAL BEFORE CYCLING =\002,d10.3/"
	//    "\002 L2-NORM OF RESIDUAL AFTER  CYCLING =\002,d10.3/\002 CONVERG"
	//    "ENCE FACTOR                 =\002,d10.3)";
   // static char fmt_9120[] = "(\002 CONVERGENCE FACTOR PER CYCLE       =\002"
	//    ",d10.3)";
    //static char fmt_9000[] = "(//\002 ************** WORK COUNT ************"
	//    "***\002/)";
    //static char fmt_9020[] = "(\002   PREP       SEC       SOL      SEC/CY"
	  //  "CLE\002)";
    //static char fmt_9030[] = "(\002 --------------------------------------"
	//    "---\002)";
   // static char fmt_9040[] = "(\002 1 RWSRT   \002,f7.2,\002   11 INTADD  "
   //	    " \002,f7.2/\002 2 PRE-COL \002,f7.2,\002   12 RESCAL   \002,f7.2/"
   //	    "\002 3 CHK-COL \002,f7.2,\002   13 RELAX    \002,f7.2/\002 4 INT"
   //	    "ERPOL\002,f7.2,\002   14 V-*      \002,f7.2/\002 5 RESTRICT\002,"
//	    "f7.2,\002   15 OTHERS   \002,f7.2/\002 6 OPDFN   \002,f7.2,\002 "
	//    "  16 CONJ-GRAD\002,f7.2/\002 7 TRUNC   \002,f7.2,\002   17 YALE-"
	 //   "SMP \002,f7.2/\002 8 OTHERS  \002,f7.2,\002   18 ------   \002,f"
	 //   "7.2)";
    //static char fmt_9050[] = "(\002   SUM     \002,f7.2,\002      SUM     "
	//    " \002,f7.2)";
    //static char fmt_9100[] = "(//\002 ********* SPACE REQUIREMENTS ********"
	//    "*\002//\002 VECTOR      NEEDED      THEOR. MINIMUM\002/\002 ----"
	//    "----------------------------------\002/\002   A \002,i14,2x,i16"
	//    "/\002   JA\002,i14,2x,i16/\002   IA\002,i14,2x,i16/\002   U \002"
	//    ",i14,2x,i16/\002   F \002,i14,2x,i16/\002   IG\002,i14,2x,i16"
	//    "/\002 --------------------------------------\002/)";
    //static char fmt_9080[] = "(/\002 ******************* COMPLEXITIES ******"
	//    "**************\002/\002 SPACE OCCUPIED BY ALL OPERATORS / SPACE "
	//    "OF OPERATOR  \002/\002 ON THE FINEST GRID   = \002,f8.2,\002   ("
	//    "A-COMPLEXITY)     \002/\002 TOTAL NUMBER OF GRID POINTS / NUMBER"
	//    " OF POINTS IN    \002/\002 THE  FINEST  GRID    = \002,f8.2,\002"
	//    "   (O-COMPLEXITY)     \002/\002 TOTAL SPACE USED BY AMG1R5 / SPA"
	//    "CE OCCUPIED BY USER- \002/\002 DEFINED  PROBLEM     = \002,f8.2"
	//    ",\002   (S-COMPLEXITY)     \002/\002 SPACE USED DURING SOLUTION "
	//    "PHASE / SPACE OCCUPIED BY \002/\002 USER-DEFINED PROBLEM = \002,"
	//    "f8.2/\002 ****************************************************"
	//    "*\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
   // integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, k;
    //static real t[10];
	static unsigned int t[10];
    static integer nnu;
   // static real sum1, sum2;
	static unsigned int sum1, sum2;
    static doublereal cfac, cfpc;
    static integer mdta, mdtf, mdtu, idima, mdtia, mdtja, mdtig;
    static doublereal acmplx, ocmplx, scmplx, tcmplx;

    /* Fortran I/O blocks */
    //static cilist io___390 = { 0, 0, 0, fmt_9110, 0 };
   // static cilist io___392 = { 0, 0, 0, fmt_9120, 0 };
   // static cilist io___393 = { 0, 0, 0, fmt_9000, 0 };
   // static cilist io___399 = { 0, 0, 0, fmt_9020, 0 };
   // static cilist io___400 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___401 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___402 = { 0, 0, 0, fmt_9030, 0 };
   // static cilist io___403 = { 0, 0, 0, fmt_9050, 0 };
    //static cilist io___404 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___413 = { 0, 0, 0, fmt_9100, 0 };
    //static cilist io___418 = { 0, 0, 0, fmt_9080, 0 };



/*     RESIDUALS / CP-TIMES / COMPLEXITY / DIMENSIONING */


    /* Parameter adjustments */
    --time;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ia;

    /* Function Body */
    if (*iout < 1) {
	return 0;
    }

/* ===> RESIDUALS / CONVERGENCE */

    if (*ncyc0 > 0) {
	cfac = *res / (*res0 + 1e-40);
	
	if (yes_print_amg) {
	    printf("( **************** CONVERGENCE *********\n");
	    printf("******** L2-NORM OF RESIDUAL BEFORE CYCLING =,%1.4f\n",res0[0]);
	    printf(" L2-NORM OF RESIDUAL AFTER  CYCLING =,%1.4f, CONVERG\n",res[0]);
	    printf("ENCE FACTOR                 =,%1.4f)\n",cfac);
	}

	d__1 = 1. / (doublereal) (*ncyc0);
	cfpc = pow_dd(&cfac, &d__1);
	
	if (yes_print_amg) {
	    printf("( CONVERGENCE FACTOR PER CYCLE       =%1.4f", cfpc);
	}
	 
    }
    if (*iout <= 1) {
	return 0;
    }
   
	if (yes_print_amg) {
	    printf("( ************** WORK COUNT ***************)\n");
	}

    nnu = imax[1];

/* ===> COMPUTING TIMES */

    //sum1 = 0.f;
    //sum2 = 0.f;
	sum1=0;
	sum2=0;
    for (i__ = 1; i__ <= 10; ++i__) {
	//t[i__ - 1] = 0.f;
		t[i__ - 1]=0;
	if (*ncyc0 > 0) {
	   // t[i__ - 1] = time[i__ + 10] / (real) (*ncyc0);
		 t[i__ - 1] = (unsigned int)( time[i__ + 10] / (*ncyc0));
	}
	sum1 += time[i__];
	sum2 += t[i__ - 1];
/* L10: */
    }

   
	if (yes_print_amg) {
	    printf("(   PREP       SEC       SOL      SEC/CYCLE)\n");
	}

   
	if (yes_print_amg) {
	    printf("( -----------------------------------------)\n");
	}

   
if (yes_print_amg) {
	printf("( 1 RWSRT   ,%d,   11 INTADD  \n", time[1]);
	    printf(" ,%d 2 PRE-COL ,%d,   12 RESCAL   ,%d\n", t[1 - 1], time[2], t[2 - 1]);
	    printf(" 3 CHK-COL ,%d,   13 RELAX    ,%d 4 INTERPOL\n", time[3], t[3 - 1]);
	    printf(",%d,   14 V-*      ,%d 5 RESTRICT,\n", time[4], t[4 - 1]);
	    printf("%d,   15 OTHERS   ,%d 6 OPDFN   ,%d, \n", time[5], t[5 - 1],time[6]);
	    printf("  16 CONJ-GRAD ,%d 7 TRUNC   ,%d,   17 YALE-\n", t[6 - 1], time[7]);
	    printf("SMP ,%d 8 OTHERS  ,%d,   18 ------   ,%d)\n", t[7 - 1], time[8], t[8 - 1]);

	    


   
	printf("( -----------------------------------------)\n");

   

	printf("(   SUM     ,%d,      SUM      ,%d)\n",sum1,sum2);

   

	printf("( -----------------------------------------)\n");

}

/* ===> SPACE OCCUPIED BY OPERATORS A(1) - A(LEVELS) */

    idima = 0;
    i__1 = *levels;
    for (k = 1; k <= i__1; ++k) {
	idima = idima + iw[iminw[k]] - ia[imin[k]];
/* L20: */
    }

/* ===> THEORETICAL MINIMAL SPACE REQUIREMENTS */

    if (*levels < 2) {
	return 0;
    }
    mdta = iw[iminw[*levels]] - 1;
    mdtja = mdta;
    mdtia = imax[*levels] + 1;
    mdtu = imax[*levels];
    mdtf = mdtu;
    mdtig = (mdtu << 1) + nnu;
   

	if (yes_print_amg) {

	 printf("( ********* SPACE REQUIREMENTS ********");
	    printf("* VECTOR      NEEDED      THEOR. MINIMUM ----\n");
	    printf("----------------------------------   A ,%d,2x,%d\n", mda[0], mdta);
	    printf("/   JA,%d,2x,%d   IA,%d,2x,%d   U \n", mdja[0], mdtja, mdia[0], mdtia);
	    printf(",%d,2x,%d   F ,%d,2x,%d   IG,%d,2x,%d\n", mdu[0], mdtu, mdf[0], mdtf, mdig[0], mdtig);
	    printf("/ --------------------------------------)\n");
	}

/* ===> COMPLEXITIES */

    scmplx = (doublereal) (((*mda + *mdu + *mdf) << 1) + *mdja + *mdia + *mdig) 
	    / (doublereal) (nnu * 5 + 1 + (iw[iminw[1]] - ia[imin[1]]) * 3);
    tcmplx = (doublereal) (((mdta + mdtu + mdtf) << 1) + mdtja + mdtia + mdtig) // скобки поставлены в соответствии с приоритетом.
	    / (doublereal) (nnu * 5 + 1 + (iw[iminw[1]] - ia[imin[1]]) * 3);

    acmplx = (doublereal) idima / (doublereal) (iw[1] - 1);
    ocmplx = (doublereal) mdtu / (doublereal) nnu;
   
	if (yes_print_amg) {

     	printf("( ******************* COMPLEXITIES ******\n");
	    printf("************** SPACE OCCUPIED BY ALL OPERATORS / SPACE \n");
	    printf("OF OPERATOR   ON THE FINEST GRID   = ,%1.4f,   (\n",acmplx);
	    printf("A-COMPLEXITY)      TOTAL NUMBER OF GRID POINTS / NUMBER\n");
	    printf(" OF POINTS IN     THE  FINEST  GRID    = ,%1.4f,\n",ocmplx);
	    printf("   (O-COMPLEXITY)      TOTAL SPACE USED BY AMG1R5 / SPA\n");
	    printf("CE OCCUPIED BY USER-  DEFINED  PROBLEM     = ,%1.4f\n",scmplx);
	    printf(",   (S-COMPLEXITY)      SPACE USED DURING SOLUTION \n");
	    printf("PHASE / SPACE OCCUPIED BY  USER-DEFINED PROBLEM = %1.4f,\n",tcmplx);
	    printf(" *****************************************************)\n");

	}

    return 0;

} /* wrkcnt_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5 AUXILIARY PROGRAMS */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     IDEC                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ int idec_(integer *into, integer *nnum, integer *ndigit, 
	integer *iarr)
{
    /* Initialized data */

    static doublereal eps = .5;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_lg10(doublereal *);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, iq, nrest;


/*     DECOMPOSE NON-NEGATIVE INTEGER INTO INTO NNUM INTEGERS */

/*     INPUT:  INTO   - INTEGER (0.LE. INTO .LE.999999999) */
/*             NNUM   - INTEGER (1.LE. NNUM .LE.9); NUMBER OF INTEGERS */
/*                      TO BE RETURNED ON ARRAY IARR (SEE BELOW) */

/*     OUTPUT: NDIGIT - INTEGER; NUMBER OF DIGITS OF INTO */
/*             IARR   - INTEGER-ARRAY OF LENGTH 10: */
/*                      IARR(1)        = FIRST      DIGIT OF INTO, */
/*                      IARR(2)        = SECOND     DIGIT OF INTO, ..... */
/*                      IARR(NNUM-1)   = (NNUM-1)ST DIGIT OF INTO, */
/*                      IARR(NNUM)     = REST OF INTO */
/*                      IF NNUM > NDIGIT, THE CORRESPONDING COMPONENTS */
/*                      OF IARR ARE PUT TO ZERO. */

/*     WARNING: BE SURE THAT YOUT COMPUTER CAN STORE NNUM DIGITS ON AN */
/*              INTEGER VARIABLE. */

    /* Parameter adjustments */
    --iarr;

    /* Function Body */

    nrest = *into;
    d__1 = eps + (doublereal) (*into);
    *ndigit = (integer) d_lg10(&d__1) + 1;
    if (*nnum >= *ndigit) {
	goto L20;
    }
    i__1 = *ndigit - *nnum + 1;
    iq = pow_ii(&c__10, &i__1);
    iarr[*nnum] = nrest - nrest / iq * iq;
    nrest /= iq;
    for (i__ = *nnum - 1; i__ >= 1; --i__) {
	iarr[i__] = nrest - nrest / 10 * 10;
	nrest /= 10;
/* L10: */
    }
    return 0;

L20:
    for (i__ = *ndigit; i__ >= 1; --i__) {
	iarr[i__] = nrest - nrest / 10 * 10;
	nrest /= 10;
/* L30: */
    }
    i__1 = *nnum;
    for (i__ = *ndigit + 1; i__ <= i__1; ++i__) {
	iarr[i__] = 0;
/* L40: */
    }
    return 0;
} /* idec_ */


/* ....................................................................... */

/*     RANDOM                                                FUNCTION */

/* ....................................................................... */

doublereal random_(doublereal *s)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);


/*     FUNCTION TO CREATE "RANDOM" SEQUENCE OF NUMBERS BETWEEN 0 AND 0.1 */

/*     INPUT:   S      - NUMBER BETWEEN 0 AND 0.1 */

/*     OUTPUT:  RANDOM - NUMBER BETWEEN 0 AND 0.1 */
/*              S      - S=RANDOM */

    ret_val = exp(*s) * 100.;
    ret_val -= (doublereal) ((integer) ret_val);
    *s = ret_val;
    return ret_val;
} /* random_ */

/* *********************************************************************** */
/* $LARGE */
/* $NOFLOATCALLS */
/*                            APPENDIX 1                          1/15/81 */

/*        SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC SYSTEMS */
/*        OF LINEAR EQUATIONS  (UNCOMPRESSED POINTER STORAGE) */

/*        REAL*8 VERSION. NOTE: THE ORIGINAL SUBROUTINES */

/*            NDRV, NSF, NNF, NNS AND NNT */

/*        HAVE BEEN RENAMED TO */

/*            YALE8, NSF8, NNF8, NNS8 AND NNT8, RESPECTIVELY. */

/* *** SUBROUTINE YALE8 (OLD NAME: NDRV) */
/* *** DRIVER FOR SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC SYSTEMS OF */
/*       LINEAR EQUATIONS (UNCOMPRESSED POINTER STORAGE) */

/*       SUBROUTINE  NDRV  (= OLD NAME) */
/*       SUBROUTINE  YALE8 (= NEW NAME) */
/* Subroutine */ int ndrv_(integer *n, integer *r__, integer *c__, integer *
	ic, integer *ia, integer *ja, doublereal *a, doublereal *b, 
	doublereal *z__, integer *nsp, /*integer*/ doublereal *isp, doublereal *rsp, integer 
	*esp, integer *path, integer *flag__)
{
    /* Initialized data */

    static integer lratio = 2;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer d__, j, l, q, u, il, im, jl, iu, ju, max__, tmp, row;
    extern /* Subroutine */ int nnf8_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, integer *, 
	    doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), nsf8_(integer *, integer *
	    , integer *, integer *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, integer *
	    , /*integer*/ doublereal *, /*integer*/ doublereal *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, integer *
	    ), nns8_(integer *, integer *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, 
	    doublereal *, doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), nnt8_(integer *, 
	    integer *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, 
	    doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer lmax, umax, jlmax, jumax, jutmp;


/*    PARAMETERS */
/*    CLASS ABBREVIATIONS ARE -- */
/*       N - INTEGER VARIABLE */
/*       F - REAL VARIABLE */
/*       V - SUPPLIES A VALUE TO THE DRIVER */
/*       R - RETURNS A RESULT FROM THE DRIVER */
/*       I - USED INTERNALLY BY THE DRIVER */
/*       A - ARRAY */

/* CLASS   PARAMETER */
/* ------+---------- */

/*         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED */
/*    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO */
/*    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY */
/*    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES */
/*    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN */
/*    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND */
/*    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF */
/*    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST */
/*    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN */
/*    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION */
/*    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1). */
/*    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY */
/*    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED */
/*    CONSECUTIVELY IN */
/*            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1), */
/*    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN */
/*            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1). */
/*    FOR EXAMPLE, THE 5 BY 5 MATRIX */
/*                ( 1. 0. 2. 0. 0.) */
/*                ( 0. 3. 0. 0. 0.) */
/*            M = ( 0. 4. 5. 6. 0.) */
/*                ( 0. 0. 0. 7. 0.) */
/*                ( 0. 0. 0. 8. 9.) */
/*    WOULD BE STORED AS */
/*                 1  2  3  4  5  6  7  8  9 */
/*            ---+-------------------------- */
/*            IA   1  3  4  7  8 10 */
/*            JA   1  3  2  2  3  4  4  4  5 */
/*             A   1. 2. 3. 4. 5. 6. 7. 8. 9.         . */

/* NV      N     - NUMBER OF VARIABLES/EQUATIONS. */
/* FVA     A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED */
/*                   BY ROWS. */
/*                   SIZE = NUMBER OF NONZERO ENTRIES IN M. */
/* NVA     IA    - POINTERS TO DELIMIT THE ROWS IN A. */
/*                   SIZE = N+1. */
/* NVA     JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A. */
/*                   SIZE = SIZE OF A. */
/* FVA     B     - RIGHT-HAND SIDE B;  B AND Z CAN THE SAME ARRAY. */
/*                   SIZE = N. */
/* FRA     Z     - SOLUTION X;  B AND Z CAN BE THE SAME ARRAY. */
/*                   SIZE = N. */

/*         THE ROWS AND COLUMNS OF THE ORIGINAL MATRIX M CAN BE */
/*    REORDERED (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY) */
/*    BEFORE CALLING THE DRIVER.  IF NO REORDERING IS DONE, THEN SET */
/*    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE SOLUTION Z IS RETURNED */
/*    IN THE ORIGINAL ORDER. */

/* NVA     R     - ORDERING OF THE ROWS OF M. */
/*                   SIZE = N. */
/* NVA     C     - ORDERING OF THE COLUMNS OF M. */
/*                   SIZE = N. */
/* NVA     IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E., */
/*                   IC(C(I)) = I  FOR I=1,...,N. */
/*                   SIZE = N. */

/*         THE SOLUTION OF THE SYSTEM OF LINEAR EQUATIONS IS DIVIDED INTO */
/*    THREE STAGES -- */
/*      NSF -- THE MATRIX M IS PROCESSED SYMBOLICALLY TO DETERMINE WHERE */
/*              FILLIN WILL OCCUR DURING THE NUMERIC FACTORIZATION. */
/*      NNF -- THE MATRIX M IS FACTORED NUMERICALLY INTO THE PRODUCT LDU */
/*              OF A UNIT LOWER TRIANGULAR MATRIX L, A DIAGONAL MATRIX D, */
/*              AND A UNIT UPPER TRIANGULAR MATRIX U, AND THE SYSTEM */
/*              MX = B  IS SOLVED. */
/*      NNS -- THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE LDU */
/*  OR          FACTORIZATION FROM NNF. */
/*      NNT -- THE TRANSPOSED LINEAR SYSTEM  MT X = B  IS SOLVED USING */
/*              THE LDU FACTORIZATION FROM NNF. */
/*    FOR SEVERAL SYSTEMS WHOSE COEFFICIENT MATRICES HAVE THE SAME */
/*    NONZERO STRUCTURE, NSF NEED BE DONE ONLY ONCE (FOR THE FIRST */
/*    SYSTEM);  THEN NNF IS DONE ONCE FOR EACH ADDITIONAL SYSTEM.  FOR */
/*    SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, NSF AND NNF NEED */
/*    BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN NNS OR NNT IS DONE */
/*    ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE. */

/* NV      PATH  - PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -- */
/*                   1  PERFORM NSF AND NNF. */
/*                   2  PERFORM NNF ONLY  (NSF IS ASSUMED TO HAVE BEEN */
/*                       DONE IN A MANNER COMPATIBLE WITH THE STORAGE */
/*                       ALLOCATION USED IN THE DRIVER). */
/*                   3  PERFORM NNS ONLY  (NSF AND NNF ARE ASSUMED TO */
/*                       HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE */
/*                       STORAGE ALLOCATION USED IN THE DRIVER). */
/*                   4  PERFORM NNT ONLY  (NSF AND NNF ARE ASSUMED TO */
/*                       HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE */
/*                       STORAGE ALLOCATION USED IN THE DRIVER). */
/*                   5  PERFORM NSF ONLY. */

/*         VARIOUS ERRORS ARE DETECTED BY THE DRIVER AND THE INDIVIDUAL */
/*    SUBROUTINES. */

/* NR      FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -- */
/*                     0     NO ERRORS DETECTED */
/*                     N+K   NULL ROW IN A  --  ROW = K */
/*                    2N+K   DUPLICATE ENTRY IN A  --  ROW = K */
/*                    3N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K */
/*                    4N+1   INSUFFICIENT STORAGE IN NNF */
/*                    5N+K   NULL PIVOT  --  ROW = K */
/*                    6N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K */
/*                    7N+1   INSUFFICIENT STORAGE IN NNF */
/*                    8N+K   ZERO PIVOT  --  ROW = K */
/*                   10N+1   INSUFFICIENT STORAGE IN NDRV */
/*                   11N+1   ILLEGAL PATH SPECIFICATION */

/*         WORKING STORAGE IS NEEDED FOR THE FACTORED FORM OF THE MATRIX */
/*    M PLUS VARIOUS TEMPORARY VECTORS.  THE ARRAYS ISP AND RSP SHOULD BE */
/*    EQUIVALENCED;  INTEGER STORAGE IS ALLOCATED FROM THE BEGINNING OF */
/*    ISP AND REAL STORAGE FROM THE END OF RSP. */

/* NV      NSP   - DECLARED DIMENSION OF RSP;  NSP GENERALLY MUST */
/*                   BE LARGER THAN  5N+3 + 2K  (WHERE  K = (NUMBER OF */
/*                   NONZERO ENTRIES IN M)). */
/* NVIRA   ISP   - INTEGER WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS */
/*                   NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE */
/*                   EQUIVALENCED. */
/*                   SIZE = LRATIO*NSP */
/* FVIRA   RSP   - REAL WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS */
/*                   NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE */
/*                   EQUIVALENCED. */
/*                   SIZE = NSP. */
/* NR      ESP   - IF SUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE */
/*                   SYMBOLIC FACTORIZATION (NSF), THEN ESP IS SET TO THE */
/*                   AMOUNT OF EXCESS STORAGE PROVIDED (NEGATIVE IF */
/*                   INSUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE */
/*                   NUMERIC FACTORIZATION (NNF)). */


/*  CONVERSION TO DOUBLE PRECISION */

/*    TO CONVERT THESE ROUTINES FOR DOUBLE PRECISION ARRAYS, SIMPLY USE */
/*    THE DOUBLE PRECISION DECLARATIONS IN PLACE OF THE REAL DECLARATIONS */
/*    IN EACH SUBPROGRAM;  IN ADDITION, THE DATA VALUE OF THE INTEGER */
/*    VARIABLE LRATIO MUST BE SET AS INDICATED IN SUBROUTINE NDRV */

/*       REAL  A(1),  B(1),  Z(1),  RSP(1) */

/*  SET LRATIO EQUAL TO THE RATIO BETWEEN THE LENGTH OF FLOATING POINT */
/*  AND INTEGER ARRAY DATA;  E. G., LRATIO = 1 FOR (REAL, INTEGER), */
/*  LRATIO = 2 FOR (DOUBLE PRECISION, INTEGER) */

    /* Parameter adjustments */
    --rsp;
    --isp;
    --z__;
    --b;
    --a;
    --ja;
    --ia;
    --ic;
    --c__;
    --r__;

    /* Function Body */

    if (*path < 1 || 5 < *path) {
	goto L111;
    }
/*  ******  INITIALIZE AND DIVIDE UP TEMPORARY STORAGE  ***************** */
    il = 1;
    iu = il + *n + 1;
    jl = iu + *n + 1;

/*  ******  CALL NSF IF FLAG IS SET  ************************************ */
    if ((*path - 1) * (*path - 5) != 0) {
	goto L2;
    }
    max__ = lratio * *nsp + 1 - jl - (*n + 1) - *n;
    jlmax = max__ / 2;
    q = jl + jlmax;
    im = q + (*n + 1);
    jutmp = im + *n;
    jumax = lratio * *nsp + 1 - jutmp;
    *esp = max__ / lratio;
    if (jlmax <= 0 || jumax <= 0) {
	goto L110;
    }
    nsf8_(n, &r__[1], &ic[1], &ia[1], &ja[1], &isp[il], &isp[jl], &jlmax, &
	    isp[iu], &isp[jutmp], &jumax, &isp[q], &isp[im], flag__);
    if (*flag__ != 0) {
	goto L100;
    }
/*  ******  MOVE JU NEXT TO JL  ***************************************** */
    //jlmax = isp[il + *n] - 1;
	jlmax = ((integer)(isp[il + *n])) - 1;
    ju = jl + jlmax;
    //jumax = isp[iu + *n] - 1;
	jumax = ((integer)(isp[iu + *n])) - 1;
    if (jumax <= 0) {
	goto L2;
    }
    i__1 = jumax;
    for (j = 1; j <= i__1; ++j) {
/* L1: */
	isp[ju + j - 1] = isp[jutmp + j - 1];
    }

/*  ******  CALL REMAINING SUBROUTINES  ********************************* */
L2:
    //jlmax = isp[il + *n] - 1;
	jlmax = ((integer)(isp[il + *n])) - 1;
    ju = jl + jlmax;
    //jumax = isp[iu + *n] - 1;
	jumax = ((integer)(isp[iu + *n])) - 1;
    l = (ju + jumax - 2 + lratio) / lratio + 1;
    lmax = jlmax;
    d__ = l + lmax;
    u = d__ + *n;
    row = *nsp + 1 - *n;
    tmp = row - *n;
    umax = tmp - u;
    *esp = umax - jumax;

    if ((*path - 1) * (*path - 2) != 0) {
	goto L3;
    }
    if (umax <= 0) {
	goto L110;
    }
    nnf8_(n, &r__[1], &c__[1], &ic[1], &ia[1], &ja[1], &a[1], &z__[1], &b[1], 
	    &isp[il], &isp[jl], &rsp[l], &lmax, &rsp[d__], &isp[iu], &isp[ju],
	     &rsp[u], &umax, &rsp[row], &rsp[tmp], flag__);
    if (*flag__ != 0) {
	goto L100;
    }
    return 0;

L3:
    if (*path - 3 != 0) {
	goto L4;
    }
    nns8_(n, &r__[1], &c__[1], &isp[il], &isp[jl], &rsp[l], &rsp[d__], &isp[
	    iu], &isp[ju], &rsp[u], &z__[1], &b[1], &rsp[tmp]);

L4:
    if (*path - 4 != 0) {
	goto L5;
    }
    nnt8_(n, &r__[1], &c__[1], &isp[il], &isp[jl], &rsp[l], &rsp[d__], &isp[
	    iu], &isp[ju], &rsp[u], &z__[1], &b[1], &rsp[tmp]);
L5:
    return 0;

/* ** ERROR:  ERROR DETECTED IN NSF, NNF, NNS, OR NNT */
L100:
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE */
L110:
    *flag__ = *n * 10 + 1;
    return 0;
/* ** ERROR:  ILLEGAL PATH SPECIFICATION */
L111:
    *flag__ = *n * 11 + 1;
    return 0;
} /* ndrv_ */


/*       ---------------------------------------------------------------- */

/*               YALE SPARSE MATRIX PACKAGE - NONSYMMETRIC CODES */
/*                    SOLVING THE SYSTEM OF EQUATIONS MX = B */
/*                        (UNCOMPRESSED POINTER STORAGE) */

/*    I.   CALLING SEQUENCES */
/*         THE COEFFICIENT MATRIX CAN BE PROCESSED BY AN ORDERING ROUTINE */
/*    (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY) BEFORE USING */
/*    THE REMAINING SUBROUTINES.  IF NO REORDERING IS DONE, THEN SET */
/*    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE CALLING SEQUENCE IS -- */
/*        (      (MATRIX ORDERING)) */
/*         NSF   (SYMBOLIC FACTORIZATION TO DETERMINE WHERE FILLIN WILL */
/*                 OCCUR DURING NUMERIC FACTORIZATION) */
/*         NNF   (NUMERIC FACTORIZATION INTO PRODUCT LDU OF UNIT LOWER */
/*                 TRIANGULAR MATRIX L, DIAGONAL MATRIX D, AND UNIT UPPER */
/*                 TRIANGULAR MATRIX U, AND SOLUTION OF LINEAR SYSTEM) */
/*         NNS   (SOLUTION OF LINEAR SYSTEM FOR ADDITIONAL RIGHT-HAND */
/*     OR          SIDE USING LDU FACTORIZATION FROM NNF) */
/*         NNT   (SOLUTION OF TRANSPOSED LINEAR SYSTEM FOR ADDITIONAL */
/*                 RIGHT-HAND SIDE USING LDU FACTORIZATION FROM NNF) */

/*    II.  STORAGE OF SPARSE MATRICES */
/*         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED */
/*    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO */
/*    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY */
/*    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES */
/*    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN */
/*    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND */
/*    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF */
/*    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST */
/*    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN */
/*    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION */
/*    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1). */
/*    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY */
/*    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED */
/*    CONSECUTIVELY IN */
/*            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1), */
/*    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN */
/*            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1). */
/*    FOR EXAMPLE, THE 5 BY 5 MATRIX */
/*                ( 1. 0. 2. 0. 0.) */
/*                ( 0. 3. 0. 0. 0.) */
/*            M = ( 0. 4. 5. 6. 0.) */
/*                ( 0. 0. 0. 7. 0.) */
/*                ( 0. 0. 0. 8. 9.) */
/*    WOULD BE STORED AS */
/*                 1  2  3  4  5  6  7  8  9 */
/*            ---+-------------------------- */
/*            IA   1  3  4  7  8 10 */
/*            JA   1  3  2  2  3  4  4  4  5 */
/*             A   1. 2. 3. 4. 5. 6. 7. 8. 9.         . */

/*         THE STRICT TRIANGULAR PORTIONS OF THE MATRICES L AND U ARE */
/*    STORED IN THE SAME FASHION USING THE ARRAYS  IL, JL, L  AND */
/*    IU, JU, U  RESPECTIVELY.  THE DIAGONAL ENTRIES OF L AND U ARE */
/*    ASSUMED TO BE EQUAL TO ONE AND ARE NOT STORED.  THE ARRAY D */
/*    CONTAINS THE RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D. */

/*    III. ADDITIONAL STORAGE SAVINGS */
/*         IN NSF, R AND IC CAN BE THE SAME ARRAY IN THE CALLING */
/*    SEQUENCE IF NO REORDERING OF THE COEFFICIENT MATRIX HAS BEEN DONE. */
/*         IN NNF, R, C AND IC CAN ALL BE THE SAME ARRAY IF NO REORDERING */
/*    HAS BEEN DONE.  IF ONLY THE ROWS HAVE BEEN REORDERED, THEN C AND IC */
/*    CAN BE THE SAME ARRAY.  IF THE ROW AND COLUMN ORDERINGS ARE THE */
/*    SAME, THEN R AND C CAN BE THE SAME ARRAY.  Z AND ROW CAN BE THE */
/*    SAME ARRAY. */
/*         IN NNS OR NNT, R AND C CAN BE THE SAME ARRAY IF NO REORDERING */
/*    HAS BEEN DONE OR IF THE ROW AND COLUMN ORDERINGS ARE THE SAME.  Z */
/*    AND B CAN BE THE SAME ARRAY;  HOWEVER, THEN B WILL BE DESTROYED. */

/*    IV.  PARAMETERS */
/*         FOLLOWING IS A LIST OF PARAMETERS TO THE PROGRAMS.  NAMES ARE */
/*    UNIFORM AMONG THE VARIOUS SUBROUTINES.  CLASS ABBREVIATIONS ARE -- */
/*       N - INTEGER VARIABLE */
/*       F - REAL VARIABLE */
/*       V - SUPPLIES A VALUE TO A SUBROUTINE */
/*       R - RETURNS A RESULT FROM A SUBROUTINE */
/*       I - USED INTERNALLY BY A SUBROUTINE */
/*       A - ARRAY */

/* CLASS   PARAMETER */
/* ------+---------- */
/* FVA     A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED */
/*                   BY ROWS. */
/*                   SIZE = NUMBER OF NONZERO ENTRIES IN M. */
/* FVA     B     - RIGHT-HAND SIDE B. */
/*                   SIZE = N. */
/* NVA     C     - ORDERING OF THE COLUMNS OF M. */
/*                   SIZE = N. */
/* FVRA    D     - RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D. */
/*                   SIZE = N. */
/* NR      FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -- */
/*                    0     NO ERRORS DETECTED */
/*                    N+K   NULL ROW IN A  --  ROW = K */
/*                   2N+K   DUPLICATE ENTRY IN A  --  ROW = K */
/*                   3N+K   INSUFFICIENT STORAGE FOR JL  --  ROW = K */
/*                   4N+1   INSUFFICIENT STORAGE FOR L */
/*                   5N+K   NULL PIVOT  --  ROW = K */
/*                   6N+K   INSUFFICIENT STORAGE FOR JU  --  ROW = K */
/*                   7N+1   INSUFFICIENT STORAGE FOR U */
/*                   8N+K   ZERO PIVOT  --  ROW = K */
/* NVA     IA    - POINTERS TO DELIMIT THE ROWS IN A. */
/*                   SIZE = N+1. */
/* NVA     IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E., */
/*                   IC(C(I) = I  FOR I=1,...N. */
/*                   SIZE = N. */
/* NVRA    IL    - POINTERS TO DELIMIT THE ROWS IN L. */
/*                   SIZE = N+1. */
/* NVRA    IU    - POINTERS TO DELIMIT THE ROWS IN U. */
/*                   SIZE = N+1. */
/* NVA     JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A. */
/*                   SIZE = SIZE OF A. */
/* NVRA    JL    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF L. */
/*                   SIZE = JLMAX. */
/* NV      JLMAX - DECLARED DIMENSION OF JL;  JLMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER */
/*                   TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF). */
/* NVRA    JU    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF U. */
/*                   SIZE = JUMAX. */
/* NV      JUMAX - DECLARED DIMENSION OF JU;  JUMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER */
/*                   TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF). */
/* FVRA    L     - NONZERO ENTRIES IN THE STRICT LOWER TRIANGULAR PORTION */
/*                   OF THE MATRIX L, STORED BY ROWS. */
/*                   SIZE = LMAX */
/* NV      LMAX  - DECLARED DIMENSION OF L;  LMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER */
/*                   TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF). */
/* NV      N     - NUMBER OF VARIABLES/EQUATIONS. */
/* NVA     R     - ORDERING OF THE ROWS OF M. */
/*                   SIZE = N. */
/* FVRA    U     - NONZERO ENTRIES IN THE STRICT UPPER TRIANGULAR PORTION */
/*                   OF THE MATRIX U, STORED BY ROWS. */
/*                   SIZE = UMAX. */
/* NV      UMAX  - DECLARED DIMENSION OF U;  UMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER */
/*                   TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF). */
/* FRA     Z     - SOLUTION X. */
/*                   SIZE = N. */


/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NSF */
/* *** SYMBOLIC LDU-FACTORIZATION OF A NONSYMMETRIC SPARSE MATRIX */
/*      (UNCOMPRESSED POINTER STORAGE) */

/* Subroutine */ int nsf8_(integer *n, integer *r__, integer *ic, integer *ia,
	 integer *ja, /*integer*/ doublereal *il, /*integer*/ doublereal *jl, integer *jlmax, /*integer*/ doublereal *iu, 
	/*integer*/ doublereal *ju, integer *jumax, /*integer*/ doublereal *q, /*integer*/ doublereal *im, integer *flag__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m, qm, vj, jmin, jmax, jlptr, juptr;


/*       INPUT VARIABLES:   N, R,IC, IA,JA, JLMAX, JUMAX. */
/*       OUTPUT VARIABLES:  IL,JL, IU,JU, FLAG. */

/*       PARAMETERS USED INTERNALLY: */
/* NIA     Q     - SUPPOSE M' IS THE RESULT OF REORDERING M;  IF */
/*                   PROCESSING OF THE KTH ROW OF M' (HENCE THE KTH ROWS */
/*                   OF L AND U) IS BEING DONE, THEN Q(J) IS INITIALLY */
/*                   NONZERO IF M'(K,J) IS NONZERO;  SINCE VALUES NEED */
/*                   NOT BE STORED, EACH ENTRY POINTS TO THE NEXT */
/*                   NONZERO;  FOR EXAMPLE, IF  N=9  AND THE 5TH ROW OF */
/*                   M' IS */
/*                           0 X X 0 X 0 0 X 0, */
/*                   THEN Q WILL INITIALLY BE */
/*                           A 3 5 A 8 A A 10 A 2        (A - ARBITRARY); */
/*                   Q(N+1) POINTS TO THE FIRST NONZERO IN THE ROW AND */
/*                   THE LAST NONZERO POINTS TO  N+1;  AS THE ALGORITHM */
/*                   PROCEEDS, OTHER ELEMENTS OF Q ARE INSERTED IN THE */
/*                   LIST BECAUSE OF FILLIN. */
/*                   SIZE = N+1. */
/* NIA     IM    - AT EACH STEP IN THE FACTORIZATION, IM(I) IS THE LAST */
/*                   ELEMENT IN THE ITH ROW OF U WHICH NEEDS TO BE */
/*                   CONSIDERED IN COMPUTING FILLIN. */
/*                   SIZE = N. */

/*  INTERNAL VARIABLES-- */
/*    JLPTR - POINTS TO THE LAST POSITION USED IN  JL. */
/*    JUPTR - POINTS TO THE LAST POSITION USED IN  JU. */


/*  ******  INITIALIZE POINTERS  **************************************** */
    /* Parameter adjustments */
    --im;
    --q;
    --ju;
    --iu;
    --jl;
    --il;
    --ja;
    --ia;
    --ic;
    --r__;

    /* Function Body */
    jlptr = 0;
   // il[1] = 1;
	 il[1] = 1.0;
    juptr = 0;
   // iu[1] = 1;
	iu[1] = 1.0;

/*  ******  FOR EACH ROW OF L AND U  ************************************ */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*  ******  SET Q TO THE REORDERED ROW OF A  **************************** */
	//q[*n + 1] = *n + 1;
	q[*n + 1] = (doublereal)(*n + 1);
	jmin = ia[r__[k]];
	jmax = ia[r__[k] + 1] - 1;
	if (jmin > jmax) {
	    goto L101;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
	    vj = ic[ja[j]];
	    qm = *n + 1;
L1:
	    m = qm;
	   // qm = q[m];
		 qm = (integer)(q[m]);
	    if (qm < vj) {
		goto L1;
	    }
	    if (qm == vj) {
		goto L102;
	    }
	   // q[m] = vj;
	   // q[vj] = qm;
		 q[m] = (doublereal)(vj);
	    q[vj] = (doublereal)(qm);
/* L2: */
	}

/*  ******  FOR EACH ENTRY IN THE LOWER TRIANGLE  *********************** */
	i__ = *n + 1;
L3:
	//i__ = q[i__];
	i__ = (integer)(q[i__]);
	if (i__ >= k) {
	    goto L7;
	}
/*  ******  L(K,I) WILL BE NONZERO, SO ADD IT TO JL  ******************** */
	++jlptr;
	if (jlptr > *jlmax) {
	    goto L103;
	}
	//jl[jlptr] = i__;
	jl[jlptr] = (doublereal)(i__);
	qm = i__;
/*  ******  INSPECT ITH ROW FOR FILLIN, ADJUST IM IF POSSIBLE  ********** */
	//jmin = iu[i__];
	jmin = (integer)(iu[i__]);
	//jmax = im[i__];
	jmax = (integer)(im[i__]);
	if (jmin > jmax) {
	    goto L6;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
	   // vj = ju[j];
		vj = (integer)(ju[j]);
	    if (vj == k) {
		//im[i__] = j;
			im[i__] = (doublereal)(j);
	    }
L4:
	    m = qm;
	    //qm = q[m];
		qm = (integer)(q[m]);
	    if (qm < vj) {
		goto L4;
	    }
	    if (qm == vj) {
		goto L5;
	    }
	    //q[m] = vj;
	    //q[vj] = qm;
		q[m] = (doublereal)(vj);
	    q[vj] = (doublereal)(qm);
	    qm = vj;
L5:
	    ;
	}
L6:
	goto L3;

/*  ******  CHECK FOR NULL PIVOT  *************************************** */
L7:
	if (i__ != k) {
	    goto L105;
	}
/*  ******  REMAINING ELEMENTS OF Q DEFINE STRUCTURE OF U(K, )  ********* */
L8:
	//i__ = q[i__];
	i__ = (integer)(q[i__]);
	if (i__ > *n) {
	    goto L9;
	}
	++juptr;
	if (juptr > *jumax) {
	    goto L106;
	}
	//ju[juptr] = i__;
	ju[juptr] = (doublereal)(i__);
	goto L8;
/*  ******  GET READY FOR NEXT ROW  ************************************* */
L9:
	//im[k] = juptr;
	im[k] = (doublereal)(juptr);
	//il[k + 1] = jlptr + 1;
	il[k + 1] = (doublereal)(jlptr + 1);
/* L10: */
	//iu[k + 1] = juptr + 1;
	iu[k + 1] = (doublereal)(juptr + 1);
    }

    *flag__ = 0;
    return 0;

/* ** ERROR:  NULL ROW IN A */
L101:
    *flag__ = *n + r__[k];
    return 0;
/* ** ERROR:  DUPLICATE ENTRY IN A */
L102:
    *flag__ = (*n << 1) + r__[k];
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE FOR JL */
L103:
    *flag__ = *n * 3 + k;
    return 0;
/* ** ERROR:  NULL PIVOT */
L105:
    *flag__ = *n * 5 + k;
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE FOR JU */
L106:
    *flag__ = *n * 6 + k;
    return 0;
} /* nsf8_ */


/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NNF */
/* *** NUMERIC LDU-FACTORIZATION OF SPARSE NONSYMMETRIC MATRIX AND */
/*      SOLUTION OF SYSTEM OF LINEAR EQUATIONS (UNCOMPRESSED POINTER */
/*      STORAGE) */

/* Subroutine */ int nnf8_(integer *n, integer *r__, integer *c__, integer *
	ic, integer *ia, integer *ja, doublereal *a, doublereal *z__, 
	doublereal *b, /*integer*/ doublereal *il, /*integer*/ doublereal *jl, doublereal *l, integer *lmax,
	 doublereal *d__, /*integer*/ doublereal *iu, /*integer*/ doublereal *ju, doublereal *u, integer *
	umax, doublereal *row, doublereal *tmp, integer *flag__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal dk, li, sum;
    static integer imin, jmin, imax, jmax;


/*       INPUT VARIABLES:   N, R,C,IC, IA,JA,A, B, IL,JL,LMAX, IU,JU,UMAX */
/*       OUTPUT VARIABLES:  Z, L,D,U, FLAG */

/*       PARAMETERS USED INTERNALLY: */
/* FIA     ROW   - HOLDS INTERMEDIATE VALUES IN CALCULATION OF L, D, U. */
/*                   SIZE = N. */
/* FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE */
/*                   EQUATION  UX = B'. */
/*                   SIZE = N. */

/*       REAL  A(1), Z(1), B(1),  L(1), D(1), U(1), */
/*    *     ROW(1), TMP(1),  LI, SUM, DK */

/*  ******  CHECK STORAGE  ********************************************** */
    /* Parameter adjustments */
    --tmp;
    --row;
    --u;
    --ju;
    --iu;
    --d__;
    --l;
    --jl;
    --il;
    --b;
    --z__;
    --a;
    --ja;
    --ia;
    --ic;
    --c__;
    --r__;

    /* Function Body */
    if (((integer)(il[*n + 1])) - 1 > *lmax) {
	goto L104;
    }
    if (((integer)(iu[*n + 1])) - 1 > *umax) {
	goto L107;
    }

/*  ******  FOR EACH ROW  *********************************************** */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*  ******  SET THE INITIAL STRUCTURE OF ROW  *************************** */
	//jmin = il[k];
	//jmax = il[k + 1] - 1;
	jmin = (integer)(il[k]);
	jmax = (integer)(il[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L2;
	}
/*  ******  IF L(K,M) .NE. 0, ROW(M)=0  ********************************* */
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L1: */
	   // row[jl[j]] = 0.;
		 row[((integer)(jl[j]))] = 0.;
	}
L2:
	row[k] = 0.;
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = (integer)(iu[k]);
	jmax = (integer)(iu[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L4;
	}
/*  ******  IF U(K,M) .NE. 0, ROW(M)=0  ********************************* */
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L3: */
	    //row[ju[j]] = 0.;
		 row[(integer)(ju[j])] = 0.;
	}
L4:
	jmin = ia[r__[k]];
	jmax = ia[r__[k] + 1] - 1;
/*  ******  SET ROW TO KTH ROW OF REORDERED A  ************************** */
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L5: */
	    row[ic[ja[j]]] = a[j];
	}
/*  ******  INITIALIZE SUM  ********************************************* */
	sum = b[r__[k]];

/*  ******  ASSIGN THE KTH ROW OF L AND ADJUST ROW, SUM  **************** */
	//imin = il[k];
	//imax = il[k + 1] - 1;
	imin = (integer)(il[k]);
	imax = (integer)(il[k + 1]) - 1;
	if (imin > imax) {
	    goto L8;
	}
	i__2 = imax;
	for (i__ = imin; i__ <= i__2; ++i__) {
	   // li = -row[jl[i__]];
/*  ******  IF L IS NOT REQUIRED, THEN COMMENT OUT THE FOLLOWING LINE  ** */
	    //l[i__] = -li;
	    //sum += li * tmp[jl[i__]];
	    //jmin = iu[jl[i__]];
	    //jmax = iu[jl[i__] + 1] - 1;

		int i__ind=(integer)(jl[i__]);
		li = -row[i__ind];
/*  ******  IF L IS NOT REQUIRED, THEN COMMENT OUT THE FOLLOWING LINE  ** */
	    l[i__] = -li;
	    sum += li * tmp[i__ind];
	    jmin = (integer)(iu[i__ind]);
	    jmax = (integer)(iu[i__ind + 1]) - 1;

	    if (jmin > jmax) {
		goto L7;
	    }
	    i__3 = jmax;
	    for (j = jmin; j <= i__3; ++j) {
/* L6: */
		//row[ju[j]] += li * u[j];
			row[(integer)(ju[j])] += li * u[j];
	    }
L7:
	    ;
	}

/*  ******  ASSIGN DIAGONAL D AND KTH ROW OF U, SET TMP(K)  ************* */
L8:
	if (row[k] == 0.) {
	    goto L108;
	}
	dk = 1 / row[k];
	d__[k] = dk;
	tmp[k] = sum * dk;
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = (integer)(iu[k]);
	jmax = (integer)(iu[k + 1]) - 1;

	if (jmin > jmax) {
	    goto L10;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L9: */
	   // u[j] = row[ju[j]] * dk;
		u[j] = row[(integer)(ju[j])] * dk;
	}
L10:
	;
    }

/*  ******  SOLVE  UX = TMP  BY BACK SUBSTITUTION  ********************** */
    k = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = tmp[k];
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = (integer)(iu[k]);
	jmax = (integer)(iu[k + 1]) - 1;

	if (jmin > jmax) {
	    goto L12;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L11: */
	   // sum -= u[j] * tmp[ju[j]];
		sum -= u[j] * tmp[(integer)(ju[j])];
	}
L12:
	tmp[k] = sum;
	z__[c__[k]] = sum;
/* L13: */
	--k;
    }

    *flag__ = 0;
    return 0;

/* ** ERROR:  INSUFFICIENT STORAGE FOR L */
L104:
    *flag__ = (*n << 2) + 1;
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE FOR U */
L107:
    *flag__ = *n * 7 + 1;
    return 0;
/* ** ERROR:  ZERO PIVOT */
L108:
    *flag__ = (*n << 3) + k;
    return 0;
} /* nnf8_ */


/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NNS 
/* *** NUMERIC SOLUTION OF A SPARSE NONSYMMETRIC SYSTEM OF LINEAR */
/*      EQUATIONS GIVEN LDU-FACTORIZATION (UNCOMPRESSED POINTER STORAGE) */

/* Subroutine */ int nns8_(integer *n, integer *r__, integer *c__, /*integer*/ doublereal *
	il, /*integer*/ doublereal *jl, doublereal *l, doublereal *d__, /*integer*/ doublereal *iu, /*integer*/ 
	doublereal *ju, doublereal *u, doublereal *z__, doublereal *b, doublereal *tmp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum;
    static integer jmin, jmax;


/*       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B */
/*       OUTPUT VARIABLES:  Z */

/*       PARAMETERS USED INTERNALLY: */
/* FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE */
/*                   EQUATION UX = B'. */
/*                   SIZE = N. */

/*       REAL  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), SUM */

/*  ******  SOLVE LDY = B  BY FORWARD SUBSTITUTION  ********************* */
    /* Parameter adjustments */
    --tmp;
    --b;
    --z__;
    --u;
    --ju;
    --iu;
    --d__;
    --l;
    --jl;
    --il;
    --c__;
    --r__;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	sum = b[r__[k]];
	//jmin = il[k];
	//jmax = il[k + 1] - 1;
	jmin = (integer)(il[k]);
	jmax = (integer)(il[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L2;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L1: */
	    //sum -= l[j] * tmp[jl[j]];
		sum -= l[j] * tmp[(integer)(jl[j])];
	}
L2:
	tmp[k] = sum * d__[k];
    }

/*  ******  SOLVE  UX = Y  BY BACK SUBSTITUTION  ************************ */
    k = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = tmp[k];
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = (integer)(iu[k]);
	jmax = (integer)(iu[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L4;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L3: */
	    //sum -= u[j] * tmp[ju[j]];
		sum -= u[j] * tmp[(integer)(ju[j])];
	}
L4:
	tmp[k] = sum;
	z__[c__[k]] = sum;
/* L5: */
	--k;
    }
    return 0;
} /* nns8_ */

/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NNT */
/* *** NUMERIC SOLUTION OF THE TRANSPOSE OF A SPARSE NONSYMMETRIC SYSTEM */
/*      OF LINEAR EQUATIONS GIVEN LDU-FACTORIZATION (UNCOMPRESSED POINTER */
/*      STORAGE) */

/* Subroutine */ int nnt8_(integer *n, integer *r__, integer *c__, /*integer*/ doublereal *
	il, /*integer*/ doublereal *jl, doublereal *l, doublereal *d__, /*integer*/ doublereal *iu, /*integer*/ doublereal 
	*ju, doublereal *u, doublereal *z__, doublereal *b, doublereal *tmp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, jmin, jmax;
    static doublereal tmpk;


/*       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B */
/*       OUTPUT VARIABLES:  Z */

/*       PARAMETERS USED INTERNALLY: */
/* FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE */
/*                   EQUATION LX = B'. */
/*                   SIZE = N. */

/*       REAL  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), TMPK */

/*  ******  SOLVE  UT Y = B  BY FORWARD SUBSTITUTION  ******************* */
    /* Parameter adjustments */
    --tmp;
    --b;
    --z__;
    --u;
    --ju;
    --iu;
    --d__;
    --l;
    --jl;
    --il;
    --c__;
    --r__;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* L1: */
	tmp[k] = b[c__[k]];
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	tmpk = -tmp[k];
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = (integer)(iu[k]);
	jmax = (integer)(iu[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L3;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L2: */
	   // tmp[ju[j]] += u[j] * tmpk;
		tmp[(integer)(ju[j])] += u[j] * tmpk;
	}
L3:
	;
    }

/*  ******  SOLVE  D LT X = Y  BY BACK SUBSTITUTION  ******************** */
    k = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmpk = -tmp[k] * d__[k];
	//jmin = il[k];
	//jmax = il[k + 1] - 1;
	jmin = (integer)(il[k]);
	jmax = (integer)(il[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L5;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L4: */
	    //tmp[jl[j]] += l[j] * tmpk;
		tmp[(integer)(jl[j])] += l[j] * tmpk;
	}
L5:
	z__[r__[k]] = -tmpk;
/* L6: */
	--k;
    }
    return 0;
} /* nnt8_ */


// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
void swapCSIRamg(integer * &v, doublereal * &dr, integer i, integer j)
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

} // swapCSIRamg

// Вот алгоритм PivotList
int PivotListCSIRamg(integer * &jptr, doublereal * &altr, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	int PivotValue = jptr[first];
	int PivotPoint = first;

	for (int index=(first+1); index<=last; index++) {
		if (jptr[index]<PivotValue) {
			PivotPoint++;
			swapCSIRamg(jptr, altr, PivotPoint, index);
		}
	}

	swapCSIRamg(jptr, altr, first, PivotPoint);

	return PivotPoint;
} // PivotListamg


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
void QuickSortCSIR_amg(integer * &jptr, doublereal * &altr, integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	if (0) {
		// BubbleSort
		int numberOfPairs=last-first+1;
		bool swappedElements=true;
		while (swappedElements) {
			 numberOfPairs--;
			 swappedElements=false;
			 for (int i=first; i<=first+numberOfPairs-1; i++) {
				 if (jptr[i]>jptr[i+1]) {
					 swapCSIRamg(jptr, altr, i, i+1);
					 swappedElements=true;
				 }
			 }
		}
	}
	else
	{
	    int pivot;

	    if (first < last) {
             pivot = PivotListCSIRamg(jptr, altr, first, last);
             QuickSortCSIR_amg(jptr, altr, first, pivot-1);
		     QuickSortCSIR_amg(jptr, altr, pivot+1, last);
	    }
	}
} // QuickSortCSIR_amg

// переменная с глобальной областью видимости необходима для 
// защиты от холостого рестарта, (перезапуск на сошедшмся решении).
// Это должно существенным образом экономить время пользователя.
// Нехолостой рестарт необходим на нелинейных задачах.
// 23 июля 2015.
Real finish_residual=0.0;


// Здесь содержится обвязка вызывающая amg1r5.
// локальное выдление памяти :всё внутри, многократные alloc и free.
/*void amg_loc_memory(equation3D*& sl, equation3D_bon*& slb,
	int maxelm, int maxbound,
	Real *dV, Real* &dX0,
	Real alpharelax, int iVar, bool bLRfree, QuickMemVorst& m) */


/*
void amg_loc_memory(Real** &u_original, Real** rthdsd, MATRNODE** &A,
                    int m, int n, int iVar, Real* &dX0)
{
	// dX0 - внешняя память.

	// iVar : 0 potencial, 1 - n, 2 - nu.

	// Замер времени.
	unsigned int calculation_main_start_time=0; // начало счёта мс.
	unsigned int calculation_main_end_time=0; // окончание счёта мс.
	unsigned int calculation_search_time=0; // время нахождения решения.

	calculation_main_start_time=clock(); // момент начала счёта.

	int maxelm = m*n; // число внутренних узлов.
	int maxbound = 2 * (m + n) + 4; // число граничных узлов.

	 // На случай если память не была выделена.
	 if (dX0==NULL) {
	    dX0=new Real[maxelm+maxbound];	    
	    for (int i=0; i<maxelm+maxbound; i++) {
	        dX0[i]=0.0;
	    }
	 }
	 else {
		 // инициализация.
		 int k0 = 0;
		
		for (int j1 = 0; j1 < n + 2; j1++) {
			 for (int i1 = 0; i1 < m + 2; i1++) {
				 dX0[k0] = u_original[i1][j1];
				 k0++;
			 }
		 }
	 }

	
	const Real nonzeroEPS=1e-37; // для отделения вещественного нуля
	Real res_sum=0.0;
	res_sum=0.0;
	for (int i = 1; i<m + 1; i++) {
		for (int j = 1; j<n + 1; j++) {
			Real buf = 0.0;
			buf = (A[i][j].aw*u_original[i - 1][j] + A[i][j].as*u_original[i][j - 1] - A[i][j].ap*u_original[i][j] + A[i][j].an*u_original[i][j + 1] + A[i][j].ae*u_original[i + 1][j] + rthdsd[i][j]);
			buf *= buf;
			res_sum += buf;
		}
	}
	// Граничные условия :
	// граничные условия.
	for (int i = 1; i<m + 1; i++) {
		Real buf = 0.0;
		// нижняя граница
		buf = (A[i][0].an*u_original[i][1] - A[i][0].ap*u_original[i][0] + rthdsd[i][0]);
		buf *= buf;
		res_sum += buf;
		// верхняя граница
		buf = (A[i][n + 1].as*u_original[i][n] - A[i][n + 1].ap*u_original[i][n + 1] + rthdsd[i][n + 1]);
		buf *= buf;
		res_sum += buf;
	}
	// граничные условия.
	for (int j = 1; j<n + 1; j++) {
		Real buf = 0.0;
		// левая граница
		buf = (A[0][j].ae*u_original[1][j] - A[0][j].ap*u_original[0][j] + rthdsd[0][j]);
		buf *= buf;
		res_sum += buf;
		// правая граница
		buf = (A[m + 1][j].aw*u_original[m][j] - A[m + 1][j].ap*u_original[m + 1][j] + rthdsd[m + 1][j]);
		buf *= buf;
		res_sum += buf;
	}


	{
		Real buf = 0.0;
		// Угловые точки, которые не рассчитываются но визуализируются.
		// среднее арефметическое
		// левый нижний угол.
		buf = (A[0][0].an*u_original[0][1] + A[0][0].ae*u_original[1][0] - A[0][0].ap*u_original[0][0] + rthdsd[0][0]);
		buf *= buf;
		res_sum += buf;
		// правый нижний угол.
		buf = (A[m + 1][0].aw*u_original[m][0] + A[m + 1][0].an*u_original[m + 1][1] - A[m + 1][0].ap*u_original[m + 1][0] + rthdsd[m + 1][0]);
		buf *= buf;
		res_sum += buf;
		// левый верхний угол.
		buf = (A[0][n + 1].as*u_original[0][n] + A[0][n + 1].ae*u_original[1][n + 1] - A[0][n + 1].ap*u_original[0][n + 1] + rthdsd[0][n + 1]);
		buf *= buf;
		res_sum += buf;
		// правый верхний угол.
		buf = (A[m + 1][n + 1].aw*u_original[m][n + 1] + A[m + 1][n + 1].as*u_original[m + 1][n] - A[m + 1][n + 1].ap*u_original[m + 1][n + 1] + rthdsd[m + 1][n + 1]);
		buf *= buf;
		res_sum += buf;
	}

	res_sum=sqrt(res_sum);
	//printf("residual start=%1.4e\n",res_sum);
	//getchar();
	
	// результаты тестирования
	// задача, начальная невязка , значение евклидовой нормы невязки при которой решение является полученным.
	// tgf01 5.4357e-1 1.0209e-11
	// CGHV1J с метализацией 3.3667e-1 5.0712e-12
	// tgf02 7.6872e-11 1.434e-11
	// tgf05 1.0871e+0  2.2895e-11
	// резистор на 1мм поликоре 5.0e-2 4.9174e-14
	//Diamond ZUb 4 4.0016e-1  4.64444e-11
	// DiamondZUB 4.0016e-1 1.1443e-8
	// NXP100 4.3399e+0  7.8347e-11 (для решения хватило 8Гб ОЗУ.)


	//if (res_sum>1.0E-10) 
	if (res_sum>1.05*finish_residual) // защита от повторного холостого запуска экономит время конечного пользователя.
	{

	//yes_print_amg=false;
	yes_print_amg=false;
	 
	

	int id=0;

	integer ierr=0;
	doublereal eps=1.0e-12;

	ierr=0; // изначальное состояние безошибочное.
	// Порог точности решения СЛАУ. Значение 1.0E-12 достаточно что проверено в ANSYS icepak.
	eps=1.0e-12; // рекомендуемое значение которого достаточно. 

// Требования к оперативной памяти.
//     VECTOR         NEEDED LENGTH (GUESS) 
//       A               3*NNA + 5*NNU 
//       JA              3*NNA + 5*NNU 
//       IA              2.2*NNU 
//       U               2.2*NNU 
//       F               2.2*NNU 
//       IG              5.4*NNU 

	
	integer nna=0; // количество ненулевых элементов в матрице СЛАУ.
	

	// подсчёт числа ненулевых элементов в матрице.
	nna=0;
	for (int i = 1; i<m + 1; i++) {
		for (int j = 1; j<n + 1; j++) {
			 if (fabs(A[i][j].aw)>nonzeroEPS) (nna)++;
			 if (fabs(A[i][j].as)>nonzeroEPS) (nna)++;
			 if (fabs(A[i][j].an)>nonzeroEPS) (nna)++;
			 if (fabs(A[i][j].ae)>nonzeroEPS) (nna)++;
			 if (fabs(A[i][j].ap)>nonzeroEPS) (nna)++;
		}
	}

	for (int i = 1; i<m + 1; i++) {
		// нижняя граница
		if (fabs(A[i][0].an)>nonzeroEPS) (nna)++;
		if (fabs(A[i][0].ap)>nonzeroEPS) (nna)++;
		// верхняя граница
		if (fabs(A[i][n + 1].as)>nonzeroEPS) (nna)++;
		if (fabs(A[i][n + 1].ap)>nonzeroEPS) (nna)++;
	}

	// граничные условия.
	for (int j = 1; j<n + 1; j++) {
		// левая граница
		if (fabs(A[0][j].ae)>nonzeroEPS) (nna)++;
		if (fabs(A[0][j].ap)>nonzeroEPS) (nna)++;
		// правая граница
		if (fabs(A[m + 1][j].aw)>nonzeroEPS) (nna)++;
		if (fabs(A[m + 1][j].ap)>nonzeroEPS) (nna)++;
	}
	// Угловые точки, которые не рассчитываются но визуализируются.
	// среднее арефметическое
	// левый нижний угол.
	if (fabs(A[0][0].an)>nonzeroEPS) (nna)++;
	if (fabs(A[0][0].ae)>nonzeroEPS) (nna)++;
	if (fabs(A[0][0].ap)>nonzeroEPS) (nna)++;

	// правый нижний угол.
	if (fabs(A[m + 1][0].aw)>nonzeroEPS) (nna)++;
	if (fabs(A[m + 1][0].an)>nonzeroEPS) (nna)++;
	if (fabs(A[m + 1][0].ap)>nonzeroEPS) (nna)++;

	// левый верхний угол.
	if (fabs(A[0][n + 1].as)>nonzeroEPS) (nna)++;
	if (fabs(A[0][n + 1].ae)>nonzeroEPS) (nna)++;
	if (fabs(A[0][n + 1].ap)>nonzeroEPS) (nna)++;

	// правый верхний угол.
	if (fabs(A[m + 1][n + 1].aw)>nonzeroEPS) (nna)++;
	if (fabs(A[m + 1][n + 1].as)>nonzeroEPS) (nna)++;
	if (fabs(A[m + 1][n + 1].ap)>nonzeroEPS) (nna)++;


	integer nnu=0; // число неизвестных.
	nnu=maxelm+maxbound;

	
	// Рекомендуемые по умолчанию параметры.
	//integer nda=0; // память под вектор значений матрицы слау.
	//nda=3*(nna)+5*(nnu);
	//integer ndia=0;
	//ndia=(integer)(2.2*(nnu));
	//integer ndja=0;
	//ndja=3*(nna)+5*(nnu);
	//integer ndu=0;
	//ndu=(integer)(2.2*(nnu));
	//integer ndf=0;
	//ndf=(integer)(2.2*(nnu));
	//integer ndig=0;
	//ndig=(integer)(5.4*(nnu));
	

	
	// в двое больше памяти чем рекомендовано.
	//integer nda=0; // память под вектор значений матрицы слау.
	//nda=6*(nna)+10*(nnu);
	//integer ndia=0;
	//ndia=(integer)(4.4*(nnu));
	//integer ndja=0;
	//ndja=6*(nna)+10*(nnu);
	//integer ndu=0;
	//ndu=(integer)(4.4*(nnu));
	//integer ndf=0;
	//ndf=(integer)(4.4*(nnu));
	///integer ndig=0;
	//ndig=(integer)(10.8*(nnu));
	

	// данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
	//doublereal rsize=1.51; // 1048416
	// Вынужденные течения достаточно 2.5. 
	doublereal rsize=3.0; // на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.

	integer nda=0; // память под вектор значений матрицы слау.
	nda=(integer)(rsize*(3*(nna)+5*(nnu)));
	integer ndia=0;
	ndia=(integer)(rsize*2.2*(nnu));
	integer ndja=0;
	ndja=(integer)(rsize*(3*(nna)+5*(nnu)));
	integer ndu=0;
	ndu=(integer)(rsize*2.2*(nnu));
	integer ndf=0;
	ndf=(integer)(rsize*2.2*(nnu));
	integer ndig=0;
	ndig=(integer)(rsize*5.4*(nnu));

	//     CLASS 3 - PARAMETERS: 

//     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). 

//     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. 

//                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. 

//                  2ND DIGIT OF IFIRST  --  ITYPU: 
//                    =0: NO SETTING OF FIRST APPROXIMATION, 
//                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, 
//                    =2: FIRST APPROXIMATION CONSTANT TO ONE, 
//                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH 
//                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED 
//                        BY THE FOLLWING DIGITS. 

//                  REST OF IFIRST  --  RNDU: 
//                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN 
//                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO 
//                    IFIRST=1372815) 

//     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE 
//                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. 

//                  1ST DIGIT OF NCYC  --  IGAM: 
//                    =1: V -CYCLE, 
//                    =2: V*-CYCLE, 
//                    =3: F -CYCLE, 
//                    =4: W -CYCLE. 
//                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE 
//                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY 
//                  IGAM V-CYCLES ON THAT PARTICULAR GRID. 

//                  2ND DIGIT OF NCYC  --  ICGR: 
//                    =0: NO CONJUGATE GRADIENT, 
//                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), 
//                    =2: CONJUGATE GRADIENT (FULL CG). 

//                  3RD DIGIT OF NCYC  --  ICONV: 
//                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM 
//                    (FINEST GRID): 
//                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY 
//                        NCYCLE (SEE BELOW) 
//                    =2: STOP, IF  ||RES|| < EPS 
//                    =3: STOP, IF  ||RES|| < EPS * |F| 
//                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| 
//                    WITH ||RES|| = L2-NORM OF RESIDUAL, 
//                           EPS     (SEE INPUT PARAMETER EPS) 
//                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE 
//                           |U|   = SUPREMUM NORM OF SOLUTION 
//                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L 
//                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS 
//                    AFTER AT MOST NCYCLE CYCLES. 

//                  REST OF NCYC  --  NCYCLE: 
//                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR 
//                    NCYCLE=0: NO CYCLING. 

//     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE 
//                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES 
//                  ARE PERFORMED, REGARDLESS OF EPS. 

//     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST 
//                  GRID IN CYCLING: 

//                  1ST DIGIT OF MADAPT  --  MSEL: 
//                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP 
//                        PHASE ARE USED WITHOUT CHECK. 
//                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED 
//                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS 
//                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC 
//                        (SEE BELOW). 

//                  REST OF MADAPT  --  FAC 
//                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART 
//                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. 
//                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT 
//                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 
//                        BY DEFAULT. 


//     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): 

//                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. 

//                  2ND DIGIT OF NRD  --  NRDX: 
//                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED 
//                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS 

//                  FOLLOWING DIGITS  --  ARRAY NRDTYP: 
//                    =1: RELAXATION OVER THE F-POINTS ONLY 
//                    =2: FULL GS SWEEP 
//                    =3: RELAXATION OVER THE C-POINTS ONLY 
//                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST 

//     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: 

//                  1ST DIGIT  --  NSC: 
//                    =1: GAUSS-SEIDEL METHOD 
//                    =2: DIRECT SOLVER (YALE SMP) 

//                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) 
//                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). 
//                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED 
//                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS 
//                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) 

//     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. 

//         -------------------------------------------------------------- 

//     CLASS 4 - PARAMETERS: 

//     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER 
//     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. 
//                  THE CHOICE OF THESE PARAMETERS DEPENDS ON 
//                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) 

//     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER 
//                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. 
//                  THE CHOICE OF THIS PARAMETER DEPENDS ON 
//                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) 

//     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION 
//                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID 
//                        OPERATORS 
//                    =1: NO COARSE-GRID OPERATOR TRUNCATION 



//     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): 

//          ISWTCH = 4 
//          IOUT   = 12 
//          IPRINT = 10606 

//          LEVELX = 25 
//          IFIRST = 13 
//          NCYC   = 10110 
//          EPS    = 1.D-12 
//          MADAPT = 27 
//          NRD    = 1131 
//          NSOLCO = 110 
//          NRU    = 1131 

//          ECG1   = 0. 
//          ECG2   = 0.25 
//          EWT2   = 0.35 
//          NWT    = 2 
//          NTR    = 0 



	// рекомедуемые параметры по дефолту.

	integer iswtch=0;
	iswtch=4;
    integer iout=0;
	iout=12;
	integer iprint=0;
	iprint=10606;
	integer levelx=0;
	levelx=25;
	integer ifirst=0;
	// начальное приближение :
	// 0 - используется из вне.
	// 1 - нулевое.
	// 2 - единицы.
	// 3 - случайная последовательность.
	ifirst=10;//13 по умолчанию.
	//ifirst=11; // нулевое начальное приближение.
	//ifirst=10; // вроде как начальное приближение берётся из dX0.
	// но 10 никоим образом не улучшает сходимость.
	integer ncyc=0;
	ncyc=10110; 
	integer madapt=0;
	madapt=27;
	integer nrd=0;
	nrd=1131;
	integer nsolco=0;
	nsolco=110;
	integer nru=0;
	nru=1131;
	doublereal ecg1=0.0;
	ecg1=0.0;
	doublereal ecg2=0.0;
	ecg2=0.25;
	doublereal ewt2=0.0;
	ewt2=0.35;
	integer nwt=0;
	nwt=2;
	integer ntr=0;
    ntr=0;

	integer matrix=0;
	//matrix=11; // symmetric SPD.
	matrix=22; 

	

	// allocate memory.
	doublereal *a=NULL;
	a=new doublereal[nda+1];
	if (a==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	integer *ia=NULL;
	ia=new integer[ndia+1];
	if (ia==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	integer *ja=NULL;
	ja=new integer[ndja+1];
	if (ja==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	doublereal *u=NULL;
	u=new doublereal[ndu+1];
	if (u==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	doublereal *f=NULL;
	f=new doublereal[ndf+1];
	if (f==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	integer *ig=NULL;
	ig=new integer[ndig+1];
	if (ig==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}

	// Блок инициализации нулём, возможно будет работоспособно и без него.

	for (int k=0; k<=nda; k++) {
		a[k]=0.0;
	}
	for (int k=0; k<=ndia; k++) {
		ia[k]=0;
	}
	for (int k=0; k<=ndja; k++) {
		ja[k]=0;
	}
	for (int k=0; k<=ndu; k++) {
		u[k]=0.0;
	}
	for (int k=0; k<=ndf; k++) {
		f[k]=0.0;
	}
	for (int k=0; k<=ndig; k++) {
		ig[k]=0;
	}


	// обязателная инициализация.
	for (int k=0; k<=nnu+1; k++) ia[k+id]=nna+1; // инициализация.
	if (id==1) ia[nnu+2]=0;
	

	
	   


	// начальное приближение.
	for (int i=0; i<=ndu; i++) {
		u[i]=0.0;
		if (i<maxelm+maxbound) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			u[i+id]=dX0[i]; 
		}
	}

	// правая часть.
    for (int i=0; i<=ndf; i++) {
		f[i]=0.0;
	}

	// правая часть .
	// инициализация.
	{
		int k0 = 0;
		
			for (int j1 = 0; j1 < n + 2; j1++) {
				for (int i1 = 0; i1 < m + 2; i1++) {
				f[k0] = rthdsd[i1][j1];
				k0++;
			}
		}
	}

	// см. equation3DtoCRS.

	    int ik=0; // счётчик ненулевых элементов СЛАУ
		int k = 0;
		
			for (int j1 = 0; j1 < n + 2; j1++) {
				for (int i1 = 0; i1 < m + 2; i1++) {
				if ((i1 == 0) && (j1 == 0)) {
                    // левый нижний угол
					if (fabs(A[0][0].ap) > nonzeroEPS) {
						a[ik + id] = A[0][0].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[0][0].an) > nonzeroEPS) {
						a[ik + id] = -A[0][0].an;
						ja[ik + id] = k + 1 + m + 2;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[0][0].ae) > nonzeroEPS) {
						a[ik + id] = -A[0][0].ae;
						ja[ik + id] = k + 1 + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == m+1) && (j1 == 0)) {
					// правый нижний угол
					if (fabs(A[m + 1][0].ap) > nonzeroEPS) {
						a[ik + id] = A[m+1][0].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][0].aw) > nonzeroEPS) {
						a[ik + id] = -A[m+1][0].aw;
						ja[ik + id] = k + 1 - 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][0].an) > nonzeroEPS) {
						a[ik + id] = -A[m+1][0].an;
						ja[ik + id] = k + 1 + m + 2;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == 0) && (j1 == n+1)) {
					// левый верхний угол
					if (fabs(A[0][n + 1].ap) > nonzeroEPS) {
						a[ik + id] = A[0][n+1].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[0][n + 1].as) > nonzeroEPS) {
						a[ik + id] = -A[0][n+1].as;
						ja[ik + id] = k + 1 - (m + 2);
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[0][n + 1].ae) > nonzeroEPS) {
						a[ik + id] = -A[0][n+1].ae;
						ja[ik + id] = k + 1 + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == m+1) && (j1 == n+1)) {
					// правый верхний угол.
					if (fabs(A[m + 1][n + 1].ap) > nonzeroEPS) {
						a[ik + id] = A[m+1][n+1].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][n + 1].aw) > nonzeroEPS) {
						a[ik + id] = -A[m+1][n+1].aw;
						ja[ik + id] = k + 1 - 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][n + 1].as) > nonzeroEPS) {
						a[ik + id] = -A[m+1][n+1].as;
						ja[ik + id] = k + 1 - (m + 2);
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((j1 == 0) && ((i1 >= 1) && (i1 < m + 1))) {
					// нижняя граница
					if (fabs(A[i1][0].ap) > nonzeroEPS) {
						a[ik + id] = A[i1][0].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][0].an) > nonzeroEPS) {
						a[ik + id] = -A[i1][0].an;
						ja[ik + id] = k + 1 + m + 2;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((j1 == n+1) && ((i1 >= 1) && (i1 < m + 1))) {
					// верхняя граница
					if (fabs(A[i1][n + 1].ap) > nonzeroEPS) {
						a[ik + id] = A[i1][n + 1].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][n + 1].as) > nonzeroEPS) {
						a[ik + id] = -A[i1][n + 1].as;
						ja[ik + id] = k + 1 - (m + 2);
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == 0) && ((j1>=1)&&(j1<n+1))) {
					// левая граница
					if (fabs(A[0][j1].ap) > nonzeroEPS) {
						a[ik + id] = A[0][j1].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[0][j1].ae) > nonzeroEPS) {
						a[ik + id] = -A[0][j1].ae;
						ja[ik + id] = k + 1 + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == m+1) && ((j1 >= 1) && (j1<n + 1))) {
					// правая граница
					if (fabs(A[m + 1][j1].ap) > nonzeroEPS) {
						a[ik + id] = A[m + 1][j1].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][j1].aw) > nonzeroEPS) {
						a[ik + id] = -A[m + 1][j1].aw;
						ja[ik + id] = k + 1 - 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 >= 1) && (i1 < m + 1) && (j1 >= 1) && (j1 < n + 1)) {
					
					if (fabs(A[i1][j1].ap) > nonzeroEPS) {
						a[ik + id] = A[i1][j1].ap;
						ja[ik + id] = k + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].aw) > nonzeroEPS) {
						a[ik + id] = -A[i1][j1].aw;
						ja[ik + id] = k + 1-1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].as) > nonzeroEPS) {
						a[ik + id] = -A[i1][j1].as;
						ja[ik + id] = k + 1-(m+2);
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].an) > nonzeroEPS) {
						a[ik + id] = -A[i1][j1].an;
						ja[ik + id] = k + 1+m+2;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].ae) > nonzeroEPS) {
						a[ik + id] = -A[i1][j1].ae;
						ja[ik + id] = k + 1 + 1;
						ia[k + id] = min(ik + 1, ia[k + id]);
						ik++;
					}
					k++;
				}
			}
		}

		
		


		


		// TODO : 
		// нужно акуратно прописать выделения и уничтожения памяти с учётом того что было сделано в BiCGStabP.

        // в каждой строке элементы отсортированы по номерам столбцов:
		// Но диагональный элемент всегда на первом месте в строке матрицы.
		int imove=0;
		if (id==0) imove=-1;

		// сортировка ненужна порядок следования любой, но главное чтобы первый в строке был имено диагональный элемент.
       //for (int k=0; k<(maxelm+maxbound); k++) QuickSortCSIR_amg(ja, a, ia[k+1]+1+imove, ia[k+2]-1+imove); // первый элемент всегда диагональный.
		//for (int k=0; k<(maxelm+maxbound); k++) QuickSortCSIR_amg(ja, a, ia[k+1]+imove, ia[k+2]-1+imove); 

		for (int k=1; k<=nnu; k++) ig[k+imove]=ia[k+1+imove]; // инициализация.

		

		//printf("getready ...");
		//getchar();

	    // amg - особенно хорош для поправки давления в SIMPLE алгоритме.
	    // алгоритм 1985 года.
        amg1r5_(a, ia, ja, 
	            u, f, ig, &nda, &ndia,
	            &ndja, &ndu, &ndf, &ndig, 
	            &nnu, &matrix, &iswtch, &iout, 
	            &iprint, &levelx, &ifirst, &ncyc, 
	            &eps, &madapt, &nrd, &nsolco, 
	            &nru, &ecg1, &ecg2, &ewt2, 
	            &nwt, &ntr, &ierr);

		switch (ierr) {
			case 1 : printf("dimension A small\n.");
			//getchar();
				system("pause");
			break;
			case 2 : printf("dimension IA small\n.");
			//getchar();
				system("pause");
			break;
			case 3 : printf("dimension JA small\n.");
			//getchar();
				system("pause");
			break;
			case 4 : printf("dimension U small\n.");
			//getchar();
				system("pause");
			break;
			case 5 : printf("dimension F small\n.");
			//getchar();
				system("pause");
			break;
			case 6 : printf("dimension IG small\n.");
			//getchar();
				system("pause");
			break;
		}

	     // возвращаем решение СЛАУ.
	     for (int i=0; i<maxelm+maxbound; i++) {
	        // обратное копирование.
		    dX0[i]=u[i+1+imove]; 
	     }

		 // инициализация.
		 int k0 = 0;
		
		 for (int j1 = 0; j1 < n + 2; j1++) {
			 for (int i1 = 0; i1 < m + 2; i1++) {
				 u_original[i1][j1] = dX0[k0];
				 k0++;
			 }
		 }
	

	     // освобождение памяти.
		 if (a!=NULL) {
	        delete[] a;
		 }
		 if (ia!=NULL) {
	        delete[] ia;
		 }
		 if (ja!=NULL) {
	        delete[] ja;
		 }
		 if (u!=NULL) {
			 delete[] u;
		 }
		 if (f!=NULL) {
	        delete[] f;
		 }
		 if (ig!=NULL) {
	        delete[] ig;
		 }

		 
		 res_sum=0.0;

		 for (int i = 1; i<m + 1; i++) {
			 for (int j = 1; j<n + 1; j++) {
				 Real buf = 0.0;
				 buf = (A[i][j].aw*u_original[i - 1][j] + A[i][j].as*u_original[i][j - 1] - A[i][j].ap*u_original[i][j] + A[i][j].an*u_original[i][j + 1] + A[i][j].ae*u_original[i + 1][j] + rthdsd[i][j]);
				 buf *= buf;
				 res_sum += buf;
			 }
		 }
		 // Граничные условия :
		 // граничные условия.
		 for (int i = 1; i<m + 1; i++) {
			 Real buf = 0.0;
			 // нижняя граница
			 buf = (A[i][0].an*u_original[i][1] - A[i][0].ap*u_original[i][0] + rthdsd[i][0]);
			 buf *= buf;
			 res_sum += buf;
			 // верхняя граница
			 buf = (A[i][n + 1].as*u_original[i][n] - A[i][n + 1].ap*u_original[i][n + 1] + rthdsd[i][n + 1]);
			 buf *= buf;
			 res_sum += buf;
		 }
		 // граничные условия.
		 for (int j = 1; j<n + 1; j++) {
			 Real buf = 0.0;
			 // левая граница
			 buf = (A[0][j].ae*u_original[1][j] - A[0][j].ap*u_original[0][j] + rthdsd[0][j]);
			 buf *= buf;
			 res_sum += buf;
			 // правая граница
			 buf = (A[m + 1][j].aw*u_original[m][j] - A[m + 1][j].ap*u_original[m + 1][j] + rthdsd[m + 1][j]);
			 buf *= buf;
			 res_sum += buf;
		 }


		 {
			 Real buf = 0.0;
			 // Угловые точки, которые не рассчитываются но визуализируются.
			 // среднее арефметическое
			 // левый нижний угол.
			 buf = (A[0][0].an*u_original[0][1] + A[0][0].ae*u_original[1][0] - A[0][0].ap*u_original[0][0] + rthdsd[0][0]);
			 buf *= buf;
			 res_sum += buf;
			 // правый нижний угол.
			 buf = (A[m + 1][0].aw*u_original[m][0] + A[m + 1][0].an*u_original[m + 1][1] - A[m + 1][0].ap*u_original[m + 1][0] + rthdsd[m + 1][0]);
			 buf *= buf;
			 res_sum += buf;
			 // левый верхний угол.
			 buf = (A[0][n + 1].as*u_original[0][n] + A[0][n + 1].ae*u_original[1][n + 1] - A[0][n + 1].ap*u_original[0][n + 1] + rthdsd[0][n + 1]);
			 buf *= buf;
			 res_sum += buf;
			 // правый верхний угол.
			 buf = (A[m + 1][n + 1].aw*u_original[m][n + 1] + A[m + 1][n + 1].as*u_original[m + 1][n] - A[m + 1][n + 1].ap*u_original[m + 1][n + 1] + rthdsd[m + 1][n + 1]);
			 buf *= buf;
			 res_sum += buf;
		 }

	   res_sum=sqrt(res_sum);
	   //printf("residual finish=%1.4e\n",res_sum);
	   //getchar();
	  // if (bsolid_static_only) {
		   // используется только для теплопередачи в твёрдом теле для ускорения
		   // решения задачи - защита от рестарта.
	    //   finish_residual=res_sum; // значение невязки решённой задачи.
	   //}
	   
		 }

		


         calculation_main_end_time=clock();
		 calculation_search_time +=calculation_main_end_time-calculation_main_start_time;

} // amg_loc_memory
*/

// глобальная память для amg1r5
typedef struct TamgGlobalMemory {
	doublereal *a;
	integer *ia;
	integer *ja;
	doublereal *u;
	doublereal *f;
	integer *ig;
//	integer nda;
//	integer ndia;
//	integer ndja;
//	integer ndu;
//	integer ndf;
//	integer ndig;
} amgGlobalMemory;



amgGlobalMemory amgGM;
bool bfirst_use_amgGM = true;

// адаптивный шаг по времени.
bool bpuasson = true;
Real ebuf = 1.0;
Real en = 1.0;
Real en_1 = 1.0;
Real en_2 = 1.0;
int itol = 0;
bool bfirst_tol = true;
Real start_tolerance = 0.0;

/*
void amg_global_memory(Real** &u_original, Real** rthdsd, MATRNODE** &A,
	int m, int n, int iVar, Real* &dX0)
{
	// dX0 - внешняя память.

	// iVar : 0 potencial, 1 - n, 2 - nu.

	// Замер времени.
	unsigned int calculation_main_start_time = 0; // начало счёта мс.
	unsigned int calculation_main_end_time = 0; // окончание счёта мс.
	unsigned int calculation_search_time = 0; // время нахождения решения.

	calculation_main_start_time = clock(); // момент начала счёта.

	int maxelm = m*n; // число внутренних узлов.
	int maxbound = 2 * (m + n) + 4; // число граничных узлов.

									// На случай если память не была выделена.
	if (dX0 == NULL) {
		dX0 = new Real[maxelm + maxbound];
		for (int i = 0; i<maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}
	else {
		// инициализация.
		int k0 = 0;

		for (int j1 = 0; j1 < n + 2; j1++) {
			for (int i1 = 0; i1 < m + 2; i1++) {
				dX0[k0] = u_original[i1][j1];
				k0++;
			}
		}
	}


	const Real nonzeroEPS = 1e-37; // для отделения вещественного нуля
	Real res_sum = 0.0;
	res_sum = 0.0;
	for (int i = 1; i<m + 1; i++) {
		for (int j = 1; j<n + 1; j++) {
			Real buf = 0.0;
			buf = (A[i][j].aw*u_original[i - 1][j] + A[i][j].as*u_original[i][j - 1] - A[i][j].ap*u_original[i][j] + A[i][j].an*u_original[i][j + 1] + A[i][j].ae*u_original[i + 1][j] + rthdsd[i][j]);
			buf *= buf;
			res_sum += buf;
		}
	}
	// Граничные условия :
	// граничные условия.
	for (int i = 1; i<m + 1; i++) {
		Real buf = 0.0;
		// нижняя граница
		if (bening_condition) {
			buf = (A[i][0].an*u_original[i][1] + A[i][0].ann*u_original[i][2] - A[i][0].ap*u_original[i][0] + rthdsd[i][0]);
		}
		else {
			buf = (A[i][0].an*u_original[i][1] - A[i][0].ap*u_original[i][0] + rthdsd[i][0]);
		}
		buf *= buf;
		res_sum += buf;
		// верхняя граница
		buf = (A[i][n + 1].as*u_original[i][n] - A[i][n + 1].ap*u_original[i][n + 1] + rthdsd[i][n + 1]);
		buf *= buf;
		res_sum += buf;
	}
	// граничные условия.
	for (int j = 1; j<n + 1; j++) {
		Real buf = 0.0;
		// левая граница
		buf = (A[0][j].ae*u_original[1][j] - A[0][j].ap*u_original[0][j] + rthdsd[0][j]);
		buf *= buf;
		res_sum += buf;
		// правая граница
		buf = (A[m + 1][j].aw*u_original[m][j] - A[m + 1][j].ap*u_original[m + 1][j] + rthdsd[m + 1][j]);
		buf *= buf;
		res_sum += buf;
	}


	{
		Real buf = 0.0;
		// Угловые точки, которые не рассчитываются но визуализируются.
		// среднее арефметическое
		// левый нижний угол.
		buf = (A[0][0].an*u_original[0][1] + A[0][0].ae*u_original[1][0] - A[0][0].ap*u_original[0][0] + rthdsd[0][0]);
		buf *= buf;
		res_sum += buf;
		// правый нижний угол.
		buf = (A[m + 1][0].aw*u_original[m][0] + A[m + 1][0].an*u_original[m + 1][1] - A[m + 1][0].ap*u_original[m + 1][0] + rthdsd[m + 1][0]);
		buf *= buf;
		res_sum += buf;
		// левый верхний угол.
		buf = (A[0][n + 1].as*u_original[0][n] + A[0][n + 1].ae*u_original[1][n + 1] - A[0][n + 1].ap*u_original[0][n + 1] + rthdsd[0][n + 1]);
		buf *= buf;
		res_sum += buf;
		// правый верхний угол.
		buf = (A[m + 1][n + 1].aw*u_original[m][n + 1] + A[m + 1][n + 1].as*u_original[m + 1][n] - A[m + 1][n + 1].ap*u_original[m + 1][n + 1] + rthdsd[m + 1][n + 1]);
		buf *= buf;
		res_sum += buf;
	}

	res_sum = sqrt(res_sum);

	if (bpuasson) {
		if (_finite(res_sum) != 0) {
			ebuf = res_sum;
		}
		bpuasson = !bpuasson;
	}
	else if (!bpuasson) {
		if (_finite(res_sum) != 0) {
			if (_finite(ebuf*ebuf + res_sum*res_sum) != 0) {
				en_2 = en_1;
				en_1 = en;
				en = sqrt(ebuf*ebuf + res_sum*res_sum) / (2.0*(m + 2)*(n + 2));
			}
		}
		bpuasson = !bpuasson;
		
		if (bfirst_tol&&(itol>7)) {
			start_tolerance = en;
			bfirst_tol = false;
		}
		else {
			if (bfirst_tol) itol++;
		}
	}

	//printf("residual start=%1.4e\n",res_sum);
	//getchar();

	// результаты тестирования
	// задача, начальная невязка , значение евклидовой нормы невязки при которой решение является полученным.
	// tgf01 5.4357e-1 1.0209e-11
	// CGHV1J с метализацией 3.3667e-1 5.0712e-12
	// tgf02 7.6872e-11 1.434e-11
	// tgf05 1.0871e+0  2.2895e-11
	// резистор на 1мм поликоре 5.0e-2 4.9174e-14
	//Diamond ZUb 4 4.0016e-1  4.64444e-11
	// DiamondZUB 4.0016e-1 1.1443e-8
	// NXP100 4.3399e+0  7.8347e-11 (для решения хватило 8Гб ОЗУ.)


	//if (res_sum>1.0E-10) 
	if (res_sum>1.05*finish_residual) // защита от повторного холостого запуска экономит время конечного пользователя.
	{

		//yes_print_amg=false;
		yes_print_amg = false;



		int id = 0;

		integer ierr = 0;
		doublereal eps = 1.0e-12;

		ierr = 0; // изначальное состояние безошибочное.
				  // Порог точности решения СЛАУ. Значение 1.0E-12 достаточно что проверено в ANSYS icepak.
		eps = 1.0e-12; // рекомендуемое значение которого достаточно. 

					   // Требования к оперативной памяти.
					   //     VECTOR         NEEDED LENGTH (GUESS) 
					   //       A               3*NNA + 5*NNU 
					   //       JA              3*NNA + 5*NNU 
					   //       IA              2.2*NNU 
					   //       U               2.2*NNU 
					   //       F               2.2*NNU 
					   //       IG              5.4*NNU 


		integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.


						 // подсчёт числа ненулевых элементов в матрице.
		nna = 0;
		for (int i = 1; i<m + 1; i++) {
			for (int j = 1; j<n + 1; j++) {
				if (fabs(A[i][j].aw)>nonzeroEPS) (nna)++;
				if (fabs(A[i][j].as)>nonzeroEPS) (nna)++;
				if (fabs(A[i][j].an)>nonzeroEPS) (nna)++;
				if (fabs(A[i][j].ae)>nonzeroEPS) (nna)++;
				if (fabs(A[i][j].ap)>nonzeroEPS) (nna)++;
			}
		}

		for (int i = 1; i<m + 1; i++) {
			// нижняя граница
			if (fabs(A[i][0].an)>nonzeroEPS) (nna)++;
			if (fabs(A[i][0].ap)>nonzeroEPS) (nna)++;
			if (bening_condition) {
				if (fabs(A[i][0].ann)>nonzeroEPS) (nna)++;
			}
			// верхняя граница
			if (fabs(A[i][n + 1].as)>nonzeroEPS) (nna)++;
			if (fabs(A[i][n + 1].ap)>nonzeroEPS) (nna)++;
		}

		// граничные условия.
		for (int j = 1; j<n + 1; j++) {
			// левая граница
			if (fabs(A[0][j].ae)>nonzeroEPS) (nna)++;
			if (fabs(A[0][j].ap)>nonzeroEPS) (nna)++;
			// правая граница
			if (fabs(A[m + 1][j].aw)>nonzeroEPS) (nna)++;
			if (fabs(A[m + 1][j].ap)>nonzeroEPS) (nna)++;
		}
		// Угловые точки, которые не рассчитываются но визуализируются.
		// среднее арефметическое
		// левый нижний угол.
		if (fabs(A[0][0].an)>nonzeroEPS) (nna)++;
		if (fabs(A[0][0].ae)>nonzeroEPS) (nna)++;
		if (fabs(A[0][0].ap)>nonzeroEPS) (nna)++;

		// правый нижний угол.
		if (fabs(A[m + 1][0].aw)>nonzeroEPS) (nna)++;
		if (fabs(A[m + 1][0].an)>nonzeroEPS) (nna)++;
		if (fabs(A[m + 1][0].ap)>nonzeroEPS) (nna)++;

		// левый верхний угол.
		if (fabs(A[0][n + 1].as)>nonzeroEPS) (nna)++;
		if (fabs(A[0][n + 1].ae)>nonzeroEPS) (nna)++;
		if (fabs(A[0][n + 1].ap)>nonzeroEPS) (nna)++;

		// правый верхний угол.
		if (fabs(A[m + 1][n + 1].aw)>nonzeroEPS) (nna)++;
		if (fabs(A[m + 1][n + 1].as)>nonzeroEPS) (nna)++;
		if (fabs(A[m + 1][n + 1].ap)>nonzeroEPS) (nna)++;


		integer nnu = 0; // число неизвестных.
		nnu = maxelm + maxbound;

		
		// Рекомендуемые по умолчанию параметры.
		//integer nda=0; // память под вектор значений матрицы слау.
		//nda=3*(nna)+5*(nnu);
		//integer ndia=0;
		//ndia=(integer)(2.2*(nnu));
		//integer ndja=0;
		//ndja=3*(nna)+5*(nnu);
		//integer ndu=0;
		//ndu=(integer)(2.2*(nnu));
		//integer ndf=0;
		//ndf=(integer)(2.2*(nnu));
		//integer ndig=0;
		//ndig=(integer)(5.4*(nnu));
	

		
		// в двое больше памяти чем рекомендовано.
		//integer nda=0; // память под вектор значений матрицы слау.
		//nda=6*(nna)+10*(nnu);
		//integer ndia=0;
		//ndia=(integer)(4.4*(nnu));
		//integer ndja=0;
		//ndja=6*(nna)+10*(nnu);
		//integer ndu=0;
		//ndu=(integer)(4.4*(nnu));
		//integer ndf=0;
		//ndf=(integer)(4.4*(nnu));
		//integer ndig=0;
		//ndig=(integer)(10.8*(nnu));
	

		// данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
		//doublereal rsize=1.51; // 1048416
		// Вынужденные течения достаточно 2.5. 
		doublereal rsize = 3.0; // на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.

		integer nda = 0; // память под вектор значений матрицы слау.
		nda = (integer)(rsize*(3 * (nna)+5 * (nnu)));
		integer ndia = 0;
		ndia = (integer)(rsize*2.2*(nnu));
		integer ndja = 0;
		ndja = (integer)(rsize*(3 * (nna)+5 * (nnu)));
		integer ndu = 0;
		ndu = (integer)(rsize*2.2*(nnu));
		integer ndf = 0;
		ndf = (integer)(rsize*2.2*(nnu));
		integer ndig = 0;
		ndig = (integer)(rsize*5.4*(nnu));

		//     CLASS 3 - PARAMETERS: 

		//     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). 

		//     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. 

		//                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. 

		//                  2ND DIGIT OF IFIRST  --  ITYPU: 
		//                    =0: NO SETTING OF FIRST APPROXIMATION, 
		//                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, 
		//                    =2: FIRST APPROXIMATION CONSTANT TO ONE, 
		//                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH 
		//                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED 
		//                        BY THE FOLLWING DIGITS. 

		//                  REST OF IFIRST  --  RNDU: 
		//                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN 
		//                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO 
		//                    IFIRST=1372815) 

		//     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE 
		//                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. 

		//                  1ST DIGIT OF NCYC  --  IGAM: 
		//                    =1: V -CYCLE, 
		//                    =2: V*-CYCLE, 
		//                    =3: F -CYCLE, 
		//                    =4: W -CYCLE. 
		//                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE 
		//                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY 
		//                  IGAM V-CYCLES ON THAT PARTICULAR GRID. 

		//                  2ND DIGIT OF NCYC  --  ICGR: 
		//                    =0: NO CONJUGATE GRADIENT, 
		//                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), 
		//                    =2: CONJUGATE GRADIENT (FULL CG). 

		//                  3RD DIGIT OF NCYC  --  ICONV: 
		//                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM 
		//                    (FINEST GRID): 
		//                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY 
		//                        NCYCLE (SEE BELOW) 
		//                    =2: STOP, IF  ||RES|| < EPS 
		//                    =3: STOP, IF  ||RES|| < EPS * |F| 
		//                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| 
		//                    WITH ||RES|| = L2-NORM OF RESIDUAL, 
		//                           EPS     (SEE INPUT PARAMETER EPS) 
		//                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE 
		//                           |U|   = SUPREMUM NORM OF SOLUTION 
		//                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L 
		//                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS 
		//                    AFTER AT MOST NCYCLE CYCLES. 

		//                  REST OF NCYC  --  NCYCLE: 
		//                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR 
		//                    NCYCLE=0: NO CYCLING. 

		//     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE 
		//                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES 
		//                  ARE PERFORMED, REGARDLESS OF EPS. 

		//     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST 
		//                  GRID IN CYCLING: 

		//                  1ST DIGIT OF MADAPT  --  MSEL: 
		//                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP 
		//                        PHASE ARE USED WITHOUT CHECK. 
		//                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED 
		//                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS 
		//                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC 
		//                        (SEE BELOW). 

		//                  REST OF MADAPT  --  FAC 
		//                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART 
		//                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. 
		//                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT 
		//                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 
		//                        BY DEFAULT. 


		//     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): 

		//                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. 

		//                  2ND DIGIT OF NRD  --  NRDX: 
		//                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED 
		//                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS 

		//                  FOLLOWING DIGITS  --  ARRAY NRDTYP: 
		//                    =1: RELAXATION OVER THE F-POINTS ONLY 
		//                    =2: FULL GS SWEEP 
		//                    =3: RELAXATION OVER THE C-POINTS ONLY 
		//                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST 

		//     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: 

		//                  1ST DIGIT  --  NSC: 
		//                    =1: GAUSS-SEIDEL METHOD 
		//                    =2: DIRECT SOLVER (YALE SMP) 

		//                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) 
		//                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). 
		//                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED 
		//                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS 
		//                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) 

		//     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. 

		//         -------------------------------------------------------------- 

		//     CLASS 4 - PARAMETERS: 

		//     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER 
		//     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. 
		//                  THE CHOICE OF THESE PARAMETERS DEPENDS ON 
		//                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) 

		//     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER 
		//                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. 
		//                  THE CHOICE OF THIS PARAMETER DEPENDS ON 
		//                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) 

		//     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION 
		//                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID 
		//                        OPERATORS 
		//                    =1: NO COARSE-GRID OPERATOR TRUNCATION 



		//     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): 

		//          ISWTCH = 4 
		//          IOUT   = 12 
		//          IPRINT = 10606 

		//          LEVELX = 25 
		//          IFIRST = 13 
		//          NCYC   = 10110 
		//          EPS    = 1.D-12 
		//          MADAPT = 27 
		//          NRD    = 1131 
		//          NSOLCO = 110 
		//          NRU    = 1131 

		//          ECG1   = 0. 
		//          ECG2   = 0.25 
		//          EWT2   = 0.35 
		//          NWT    = 2 
		//          NTR    = 0 



		// рекомедуемые параметры по дефолту.

		integer iswtch = 0;
		iswtch = 4;
		integer iout = 0;
		iout = 12;
		integer iprint = 0;
		iprint = 10606;
		integer levelx = 0;
		levelx = 25;
		integer ifirst = 0;
		// начальное приближение :
		// 0 - используется из вне.
		// 1 - нулевое.
		// 2 - единицы.
		// 3 - случайная последовательность.
		ifirst = 10;//13 по умолчанию.
					//ifirst=11; // нулевое начальное приближение.
					//ifirst=10; // вроде как начальное приближение берётся из dX0.
					// но 10 никоим образом не улучшает сходимость.
		integer ncyc = 0;
		ncyc = 10110;
		integer madapt = 0;
		madapt = 27;
		integer nrd = 0;
		nrd = 1131;
		integer nsolco = 0;
		nsolco = 110;
		integer nru = 0;
		nru = 1131;
		doublereal ecg1 = 0.0;
		ecg1 = 0.0;
		doublereal ecg2 = 0.0;
		ecg2 = 0.25;
		doublereal ewt2 = 0.0;
		ewt2 = 0.35;
		integer nwt = 0;
		nwt = 2;
		integer ntr = 0;
		ntr = 0;

		integer matrix = 0;
		//matrix=11; // symmetric SPD.
		matrix = 22;

		if (bfirst_use_amgGM) {

			bfirst_use_amgGM = false;

			// allocate memory.
			//doublereal *amgGM.a = NULL;
			amgGM.a = new doublereal[2*nda + 1];
			if (amgGM.a == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
			//integer *amgGM.ia = NULL;
			amgGM.ia = new integer[2 * ndia + 1];
			if (amgGM.ia == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
			//integer *amgGM.ja = NULL;
			amgGM.ja = new integer[2 * ndja + 1];
			if (amgGM.ja == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
			//doublereal *amgGM.u = NULL;
			amgGM.u = new doublereal[2 * ndu + 1];
			if (amgGM.u == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
			//doublereal *amgGM.f = NULL;
			amgGM.f = new doublereal[2 * ndf + 1];
			if (amgGM.f == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
			//integer *amgGM.ig = NULL;
			amgGM.ig = new integer[2 * ndig + 1];
			if (amgGM.ig == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}

		}

		// Блок инициализации нулём, возможно будет работоспособно и без него.

		for (int k = 0; k <= nda; k++) {
			amgGM.a[k] = 0.0;
		}
		for (int k = 0; k <= ndia; k++) {
			amgGM.ia[k] = 0;
		}
		for (int k = 0; k <= ndja; k++) {
			amgGM.ja[k] = 0;
		}
		for (int k = 0; k <= ndu; k++) {
			amgGM.u[k] = 0.0;
		}
		for (int k = 0; k <= ndf; k++) {
			amgGM.f[k] = 0.0;
		}
		for (int k = 0; k <= ndig; k++) {
			amgGM.ig[k] = 0;
		}


		// обязателная инициализация.
		for (int k = 0; k <= nnu + 1; k++) amgGM.ia[k + id] = nna + 1; // инициализация.
		if (id == 1) amgGM.ia[nnu + 2] = 0;






		// начальное приближение.
		for (int i = 0; i <= ndu; i++) {
			amgGM.u[i] = 0.0;
			if (i<maxelm + maxbound) {
				// обязательно нужно проверить была ли выделена оперативная память. 
				amgGM.u[i + id] = dX0[i];
			}
		}

		// правая часть.
		for (int i = 0; i <= ndf; i++) {
			amgGM.f[i] = 0.0;
		}

		// правая часть .
		// инициализация.
		{
			int k0 = 0;

			for (int j1 = 0; j1 < n + 2; j1++) {
				for (int i1 = 0; i1 < m + 2; i1++) {
					amgGM.f[k0] = rthdsd[i1][j1];
					k0++;
				}
			}
		}

		// см. equation3DtoCRS.

		int ik = 0; // счётчик ненулевых элементов СЛАУ
		int k = 0;

		for (int j1 = 0; j1 < n + 2; j1++) {
			for (int i1 = 0; i1 < m + 2; i1++) {
				if ((i1 == 0) && (j1 == 0)) {
					// левый нижний угол
					if (fabs(A[0][0].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[0][0].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[0][0].an) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[0][0].an;
						amgGM.ja[ik + id] = k + 1 + m + 2;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[0][0].ae) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[0][0].ae;
						amgGM.ja[ik + id] = k + 1 + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == m + 1) && (j1 == 0)) {
					// правый нижний угол
					if (fabs(A[m + 1][0].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[m + 1][0].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][0].aw) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[m + 1][0].aw;
						amgGM.ja[ik + id] = k + 1 - 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][0].an) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[m + 1][0].an;
						amgGM.ja[ik + id] = k + 1 + m + 2;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == 0) && (j1 == n + 1)) {
					// левый верхний угол
					if (fabs(A[0][n + 1].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[0][n + 1].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[0][n + 1].as) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[0][n + 1].as;
						amgGM.ja[ik + id] = k + 1 - (m + 2);
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[0][n + 1].ae) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[0][n + 1].ae;
						amgGM.ja[ik + id] = k + 1 + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == m + 1) && (j1 == n + 1)) {
					// правый верхний угол.
					if (fabs(A[m + 1][n + 1].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[m + 1][n + 1].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][n + 1].aw) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[m + 1][n + 1].aw;
						amgGM.ja[ik + id] = k + 1 - 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][n + 1].as) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[m + 1][n + 1].as;
						amgGM.ja[ik + id] = k + 1 - (m + 2);
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((j1 == 0) && ((i1 >= 1) && (i1 < m + 1))) {
					// нижняя граница
					if (fabs(A[i1][0].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[i1][0].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][0].an) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[i1][0].an;
						amgGM.ja[ik + id] = k + 1 + m + 2;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (bening_condition) {
						if (fabs(A[i1][0].ann) > nonzeroEPS) {
							amgGM.a[ik + id] = -A[i1][0].ann;
							amgGM.ja[ik + id] = k + 1 + m + 2+m+2;
							amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
							ik++;
						}
					}
					k++;
				}
				if ((j1 == n + 1) && ((i1 >= 1) && (i1 < m + 1))) {
					// верхняя граница
					if (fabs(A[i1][n + 1].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[i1][n + 1].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][n + 1].as) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[i1][n + 1].as;
						amgGM.ja[ik + id] = k + 1 - (m + 2);
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == 0) && ((j1 >= 1) && (j1<n + 1))) {
					// левая граница
					if (fabs(A[0][j1].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[0][j1].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[0][j1].ae) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[0][j1].ae;
						amgGM.ja[ik + id] = k + 1 + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 == m + 1) && ((j1 >= 1) && (j1<n + 1))) {
					// правая граница
					if (fabs(A[m + 1][j1].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[m + 1][j1].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[m + 1][j1].aw) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[m + 1][j1].aw;
						amgGM.ja[ik + id] = k + 1 - 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
				if ((i1 >= 1) && (i1 < m + 1) && (j1 >= 1) && (j1 < n + 1)) {

					if (fabs(A[i1][j1].ap) > nonzeroEPS) {
						amgGM.a[ik + id] = A[i1][j1].ap;
						amgGM.ja[ik + id] = k + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].aw) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[i1][j1].aw;
						amgGM.ja[ik + id] = k + 1 - 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].as) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[i1][j1].as;
						amgGM.ja[ik + id] = k + 1 - (m + 2);
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].an) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[i1][j1].an;
						amgGM.ja[ik + id] = k + 1 + m + 2;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					if (fabs(A[i1][j1].ae) > nonzeroEPS) {
						amgGM.a[ik + id] = -A[i1][j1].ae;
						amgGM.ja[ik + id] = k + 1 + 1;
						amgGM.ia[k + id] = min(ik + 1, amgGM.ia[k + id]);
						ik++;
					}
					k++;
				}
			}
		}








		// TODO : 
		// нужно акуратно прописать выделения и уничтожения памяти с учётом того что было сделано в BiCGStabP.

		// в каждой строке элементы отсортированы по номерам столбцов:
		// Но диагональный элемент всегда на первом месте в строке матрицы.
		int imove = 0;
		if (id == 0) imove = -1;

		// сортировка ненужна порядок следования любой, но главное чтобы первый в строке был имено диагональный элемент.
		//for (int k=0; k<(maxelm+maxbound); k++) QuickSortCSIR_amg(ja, a, ia[k+1]+1+imove, ia[k+2]-1+imove); // первый элемент всегда диагональный.
		//for (int k=0; k<(maxelm+maxbound); k++) QuickSortCSIR_amg(ja, a, ia[k+1]+imove, ia[k+2]-1+imove); 

		for (int k = 1; k <= nnu; k++) amgGM.ig[k + imove] = amgGM.ia[k + 1 + imove]; // инициализация.



																		  //printf("getready ...");
																		  //getchar();

																		  // amg - особенно хорош для поправки давления в SIMPLE алгоритме.
																		  // алгоритм 1985 года.
		amg1r5_(amgGM.a, amgGM.ia, amgGM.ja,
			amgGM.u, amgGM.f, amgGM.ig, &nda, &ndia,
			&ndja, &ndu, &ndf, &ndig,
			&nnu, &matrix, &iswtch, &iout,
			&iprint, &levelx, &ifirst, &ncyc,
			&eps, &madapt, &nrd, &nsolco,
			&nru, &ecg1, &ecg2, &ewt2,
			&nwt, &ntr, &ierr);

		switch (ierr) {
		case 1: printf("dimension A small\n.");
			//getchar();
			system("pause");
			break;
		case 2: printf("dimension IA small\n.");
			//getchar();
			system("pause");
			break;
		case 3: printf("dimension JA small\n.");
			//getchar();
			system("pause");
			break;
		case 4: printf("dimension U small\n.");
			//getchar();
			system("pause");
			break;
		case 5: printf("dimension F small\n.");
			//getchar();
			system("pause");
			break;
		case 6: printf("dimension IG small\n.");
			//getchar();
			system("pause");
			break;
		}

		// возвращаем решение СЛАУ.
		for (int i = 0; i<maxelm + maxbound; i++) {
			// обратное копирование.
			dX0[i] = amgGM.u[i + 1 + imove];
		}

		// инициализация.
		int k0 = 0;

		for (int j1 = 0; j1 < n + 2; j1++) {
			for (int i1 = 0; i1 < m + 2; i1++) {
				u_original[i1][j1] = dX0[k0];
				k0++;
			}
		}

		
		// Мы используем глобальную память.
		// освобождение памяти.
		//if (a != NULL) {
			//delete[] a;
	//	}
		//if (ia != NULL) {
			//delete[] ia;
		//}
		//if (ja != NULL) {
			//delete[] ja;
		//}
		//if (u != NULL) {
			//delete[] u;
		//}
		//if (f != NULL) {
			//delete[] f;
		//}
		//if (ig != NULL) {
			//delete[] ig;
		//}
		

		res_sum = 0.0;

		for (int i = 1; i<m + 1; i++) {
			for (int j = 1; j<n + 1; j++) {
				Real buf = 0.0;
				buf = (A[i][j].aw*u_original[i - 1][j] + A[i][j].as*u_original[i][j - 1] - A[i][j].ap*u_original[i][j] + A[i][j].an*u_original[i][j + 1] + A[i][j].ae*u_original[i + 1][j] + rthdsd[i][j]);
				buf *= buf;
				res_sum += buf;
			}
		}
		// Граничные условия :
		// граничные условия.
		for (int i = 1; i<m + 1; i++) {
			Real buf = 0.0;
			// нижняя граница
			if (bening_condition) {
				buf = (A[i][0].an*u_original[i][1] + A[i][0].ann*u_original[i][2] - A[i][0].ap*u_original[i][0] + rthdsd[i][0]);
			}
			else {
				buf = (A[i][0].an*u_original[i][1] - A[i][0].ap*u_original[i][0] + rthdsd[i][0]);
			}
			buf *= buf;
			res_sum += buf;
			// верхняя граница
			buf = (A[i][n + 1].as*u_original[i][n] - A[i][n + 1].ap*u_original[i][n + 1] + rthdsd[i][n + 1]);
			buf *= buf;
			res_sum += buf;
		}
		// граничные условия.
		for (int j = 1; j<n + 1; j++) {
			Real buf = 0.0;
			// левая граница
			buf = (A[0][j].ae*u_original[1][j] - A[0][j].ap*u_original[0][j] + rthdsd[0][j]);
			buf *= buf;
			res_sum += buf;
			// правая граница
			buf = (A[m + 1][j].aw*u_original[m][j] - A[m + 1][j].ap*u_original[m + 1][j] + rthdsd[m + 1][j]);
			buf *= buf;
			res_sum += buf;
		}


		{
			Real buf = 0.0;
			// Угловые точки, которые не рассчитываются но визуализируются.
			// среднее арефметическое
			// левый нижний угол.
			buf = (A[0][0].an*u_original[0][1] + A[0][0].ae*u_original[1][0] - A[0][0].ap*u_original[0][0] + rthdsd[0][0]);
			buf *= buf;
			res_sum += buf;
			// правый нижний угол.
			buf = (A[m + 1][0].aw*u_original[m][0] + A[m + 1][0].an*u_original[m + 1][1] - A[m + 1][0].ap*u_original[m + 1][0] + rthdsd[m + 1][0]);
			buf *= buf;
			res_sum += buf;
			// левый верхний угол.
			buf = (A[0][n + 1].as*u_original[0][n] + A[0][n + 1].ae*u_original[1][n + 1] - A[0][n + 1].ap*u_original[0][n + 1] + rthdsd[0][n + 1]);
			buf *= buf;
			res_sum += buf;
			// правый верхний угол.
			buf = (A[m + 1][n + 1].aw*u_original[m][n + 1] + A[m + 1][n + 1].as*u_original[m + 1][n] - A[m + 1][n + 1].ap*u_original[m + 1][n + 1] + rthdsd[m + 1][n + 1]);
			buf *= buf;
			res_sum += buf;
		}

		res_sum = sqrt(res_sum);

		

		//printf("residual finish=%1.4e\n",res_sum);
		//getchar();
		// if (bsolid_static_only) {
		// используется только для теплопередачи в твёрдом теле для ускорения
		// решения задачи - защита от рестарта.
		//   finish_residual=res_sum; // значение невязки решённой задачи.
		//}

	}




	calculation_main_end_time = clock();
	calculation_search_time += calculation_main_end_time - calculation_main_start_time;

} // amg_global_memory
*/
/*
// Здесь содержится обвязка вызывающая amg1r5.
void amg(equation3D* &sl, equation3D_bon* &slb,
			   int maxelm, int maxbound,
			   Real *dV, Real* &dX0, 
			   Real alpharelax, int iVar, bool bLRfree, QuickMemVorst& m) {
				  // bool bmemory_local=false;

				   //if (bmemory_local) {
					   // локальное выделение памяти , много alloc и free.
					   amg_loc_memory(sl, slb, maxelm,  maxbound, dV, dX0, alpharelax, iVar, bLRfree,m);
				  // }
				//   else {
					   // память выделяется лишь единожды.
					//   amg_global_memory(sl, slb, maxelm,  maxbound, dV, dX0, alpharelax, iVar, bLRfree,m);
				   //}
}
*/