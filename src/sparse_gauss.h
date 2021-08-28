// ���� sparse_gauss.h 
// ���������� ������ ������ ��� ����������� �������.

#define Real double

// ��� �������
typedef struct tagTERM{
	int key;
	Real val;
	// ����������� ���� ���
	// ��� ������, ��� ����������
	// ���������� ����� � �������
	// �� ������� ��� �� ������������.

	//struct tagTERM *left;
	//struct tagTERM *right;
	//int bal;
} TERM;

// ������ ��� ������� � ������� ����
typedef struct tagIRow{
    Real eps0; // ��� ����������� ������������� ����

	// ��������� ��� ���� ������������
	// ������ ��� ���������� ���������� 
	// ����� � �������� �� ������������ �������.
	int POOL_SIZE; // ������ ������� ������� ������� ��������
	int n; // ����� ��������� ���������

	TERM *elm;

} IRow;

// ����������� ������� ����:
// ���������� nxn � ������������ �������������,
// �������� ��������������. ��������� ���������� 
// � ����.
typedef struct tagIMatrix{
    Real eps0; // ��� ����������� ������������� ����

	int n; // ����������� ������� nxn.
    // jp - ������ ������� �����������
    // jm - ������� ������ �����������
    // dd - ������� ���������
	IRow *jp, *jm;
	Real *dd;
} IMatrix;

// ����� ������ � �������� key  
int search_i(TERM* list, int n, int key);

// ������� ������� � ������ ������ key
void deleteAt(TERM* &list, int key, int &n, int &pool);

// �������� ������� �� ���������� : num, val.
void add(TERM* &list, int &n, int &pool, int num, Real val);

// ���������� �������� ������ ����
// ������� ����� key
Real get_val(TERM* list, int n, int key);

// ��������� ����� value � ������
// � ������ key
void modify_add(TERM* &list, int n, int key, Real value);

// ������������� ����� value � ������
// � ������ key
void modify_set(TERM* &list, int n, int key, Real value);

// ������� �� ����������� �������������
// ���������� ��� ��������� ��������
void get_values(TERM *list, int n, int* &indexes, Real* &values);

// ���������� 0 ��� i<=0, i ��� i>0
int sigma(int i);

// �������������� row,col � ���������� ���������� (d,j)
// ��� ������� �����������: d=0..(n-1), j=-(n-d-1)..-1
// ��� ������ �����������:  d=0..(n-1), j=1..(n-d-1)
int getD(int row, int col);
int getJ(int row, int col);

// �������� �������������� ��������� ���������� (d,j) � row,col
int getRow(int d, int j);
int getCol(int d, int j);


// ���������� �������� value � ������ � �������� num
void setValueIRow(IRow *xO, int num, Real value);

// ��������� value � ������������� �������� � ������ num
void addValueIRow(IRow *xO, int num, double value);

// ���������� �������� ������ num
Real getValueIRow(IRow *xO, int num);

// ���������� ��� ��������� ������ ������/�������: 
// ������� ����� - � indexes, �������� � values
int getValuesIRow(IRow *xO, int* &indexes, Real* &values);

// ��������� ������ ��� ����������� �������
void initIMatrix(IMatrix *xO, int n);

// ������������� �������� value � ������ � ������������ [row,col];
// row - ����� ������ �������
// col - ����� ������� �������
void setValueIMatrix(IMatrix *xO, int row, int col, Real value);

// ��������� �������� value � ������ [row,col]
void addValueIMatrix(IMatrix *xO, int row, int col, double value);

// ���������� �������� ������ [row,col]
Real  getValueIMatrix(IMatrix *xO, int  row, int  col);

// ���������� ��������� �������� � ������� ����� ������ d,
// ������� ��������� ������ ������� ���������
int  getJRowIMatrix(IMatrix *xO, int  d, int* &indexes, Real* &values);

// ���������� ��������� �������� � ������� ����� ������� d, 
// ������� ��������� ���� ������� ���������
int  getJColIMatrix(IMatrix *xO, int  d, int* &indexes, Real* &values);

// ������� �����, ������������ ������� x,
// ��������� ������ ��������� ������ b � 
// ���������� ������� xO � ����������� ����������� �������.
void calculateSPARSEgaussArray(IMatrix *xO, Real *x, Real *b);

// ����������� ������� � ������� IMatrix � CSIR ������ 
// ����������� � ����������� ITL ��� ���������� ILU ����������.
void convertIMatrixtoCSIR_ILU_ITL(IMatrix *xO, Real* &U_val, int* &U_ind, int* &U_ptr, Real* &L_val, int* &L_ind, int* &L_ptr);

