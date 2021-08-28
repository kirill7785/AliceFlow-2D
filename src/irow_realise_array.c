// ���� irow_realise_array.c 
// ���������� ���������� ����� � �������� �� ��������.

#include "sparse_gauss.h" // ���������� ���� �������

// ������������ ���������� �� ������� 
// ������������ ������� �������������� �
// ����� ������ ���������� ��� �������.
// 1 ��� ��������� ������ ��������.
const int CLUSTER_SIZE=1; // core 2 quad Q6600 2.4GGH 4cores nodes=2274: i - 10s, 5 - 9s, 10 - 8s. ��� ��� ����� �����

// ����� ������ � ������ key  
int search_i(TERM* list, int n, int key) {
	/*
	// �������� ����� ������ � ������ key
	int i=-1;
	int i1; // ������� ����� for
    // �������� ����� ������ � �������� key
	for (i1=0; i1<n; i1++) if (list[i1].key==key) i=i1;
	return i;
	*/
	// �������� ����� ������ � ������ key
	// ������ list �������������� ������������� �� ���� key.
    int iL=0, iR=n-1;
	int im;
	int i=-1;
	bool found=false;
	while ((!(found)) && (iL<=iR)) {
       im=(int)(iL+iR)/2;
	   if (list[im].key==key) {
		   found=true;
		   i=im;
	   }
	   else if (list[im].key<key) iL=im+1;
	   else iR=im-1;
	}
	return i;

} // search_i

// ������� ������� � ������ ������ key
// ��� �������� �� �������� ��������������� ������� list
void deleteAt(TERM* &list, int key, int &n, int &pool) {

   /* ���������� ��� ��� pool`a:
   // ��������� ������
   TERM *listloc=new TERM[n-1];
   int i1, i2;

   i2=0; // ������� � ��������� �����
   for (i1=0; i1<n; i1++) {
      if (i1!=key) {
		   listloc[i2]=list[i1];
		   i2++;
      }
   }
   
   delete list;
   list = new TERM[n-1];


   // �������� �����������
   for (i1=0; i1<(n-1); i1++) {
	   list[i1]=listloc[i1];
   }
   // ������������ ������
   delete listloc; 
   */

   int ipos=-1,i1; // ������� ����� for
   ipos=search_i( list, n, key);
   //for (i1=0; i1<n; i1++) if (list[i1].key==key) ipos=i1;
   //printf("i=%d",ipos);

   if ((pool-n)<CLUSTER_SIZE) {
        for (i1=ipos; i1<(n-1); i1++) list[i1]=list[i1+1];
		n--;
   }
   else
   {
        // ��������� ������
        TERM *listloc=new TERM[n-1];
        int i2=0; // ������� � ��������� �����

        for (i1=0; i1<n; i1++) {
            if (i1!=ipos) {
		        listloc[i2]=list[i1];
		        i2++;
            }
        }
   
		delete list;
        list = new TERM[n-1];

        // �������� �����������
        for (i1=0; i1<(n-1); i1++) {
	        list[i1]=listloc[i1];
        }
        // ������������ ������
        delete[] listloc;  
		n--;
		pool=n;
   }

} // deleteAt

// �������� ������� �� ���������� : num, val.
// ���������� � pool`�� ������������� ��� ����������
// ������� �������������� �� ��������� ������ � �����������.
// �������� ��������� CLUSTER_SIZE ������� �� ����������� ������� ����
// � ��������� ��� ���� ��� ������������� ������� �� ��� ������ 
// ����������� ���������������� �� ����� ������.
void add(TERM* &list, int &n, int &pool, int num, Real val) {
    TERM *listloc;
    int i1;

	/* ���������� ��� pool`�:
    // ��������� ������
	if (n>0) {
		listloc=new TERM[n];
		// ������ ��������� �����
		for (i1=0; i1<n; i1++) {
			listloc[i1]=list[i1];
		}
        delete list;
	}
	list=new TERM[n+1];
		    
    // �������� �����������
    for (i1=0; i1<n; i1++) {
		list[i1]=listloc[i1];
	}
	list[n].key=num;
	list[n].val=val;
   
	// ������������ ������
	if (n>0) delete listloc;
	n++;
	*/

	/*
	// ������ list �� ����������
	if (n==0) {
		// ���� ������ �� ��� ������
		pool=CLUSTER_SIZE;
		list=new TERM[pool];
		n=1;
        list[n-1].key=num;
	    list[n-1].val=val;  
	}
	 else
	{
		if (n<pool) {
			n++;
            list[n-1].key=num;
	        list[n-1].val=val;
		}
		else
		{
			// n==pool
			pool+=CLUSTER_SIZE;

            listloc=new TERM[n];
		    // ������ ��������� �����
		    for (i1=0; i1<n; i1++) {
			     listloc[i1]=list[i1];
		    }
            delete list; 

            list=new TERM[pool];
			// �������� �����������
            for (i1=0; i1<n; i1++) {
		         list[i1]=listloc[i1];
	        }
	        list[n].key=num;
	        list[n].val=val;
			n++;

			delete listloc;
		}
	}
	*/


	// ������ list ���������� �� ����������� ����� key
    if (n==0) {
		// ���� ������ �� ��� ������
		pool=CLUSTER_SIZE;
		list=new TERM[pool];
		n=1;
        list[n-1].key=num;
	    list[n-1].val=val;  
	}
	 else
	{
		if (n<pool) {
			// ������� �������� � ����������� ���������������
			n++;
            i1=0;
			while ((list[i1].key<num) && (i1<(n-1))) i1++;
			if (i1==(n-1)) {
                 list[n-1].key=num;
	             list[n-1].val=val;
			}
			else
			{
				for (int i2=(n-2); i2>=i1; i2--) list[i2+1]=list[i2];
                list[i1].key=num;
	            list[i1].val=val;
			}
		}
		else
		{
			// n==pool
			pool+=CLUSTER_SIZE;

            listloc=new TERM[n];
		    // ������ ��������� �����
		    for (i1=0; i1<n; i1++) {
			     listloc[i1]=list[i1];
		    }
            delete list; 

            list=new TERM[pool];
			// �������� �����������
            // � ����������� ���������������
			i1=0;
			while ((list[i1].key<num) && (i1<n)) {
				list[i1]=listloc[i1];
				i1++;
			}
			list[i1].key=num;
			list[i1].val=val;
			for (int i2=i1+1; i2<n+1; i2++) list[i2]=listloc[i2-1];
			n++;


			delete[] listloc;
		}
	}

} // add

// ���������� �������� ������ ����
// ������� ����� key
Real get_val(TERM* list, int n, int key) {
	int i=-1;
	//int i1; // ������� ����� for
    // �������� ����� ������ � �������� key
	//for (i1=0; i1<n; i1++) if (list[i1].key==key) i=i1;
	i=search_i( list, n, key);
	if (i!=-1) return (Real) list[i].val;
	return 0.0;
} // get_val

// ��������� ����� value � ������
// � ������ key
void modify_add(TERM* &list, int n, int key, Real value) {
	//int i1; // ������� ����� for
    // �������� ����� ������ � �������� key
	//for (i1=0; i1<n; i1++) if (list[i1].key==key) list[i1].val+=value;
	int i=-1;
	i=search_i( list, n, key); // ����� ������ �������
	list[i].val+=value;

} // modify_add

// ������������� ����� value � ������
// � ������ key
void modify_set(TERM* &list, int n, int key, Real value) {
	//int i1; // ������� ����� for
    // �������� ����� ������ � �������� key
	//for (i1=0; i1<n; i1++) if (list[i1].key==key) list[i1].val=value;
	int i=-1;
	i=search_i( list, n, key); // ����� ������ �������
    list[i].val=value;
} // modify_set


// ������� �� ����������� �������������
void get_values(TERM *list, int n, int* &indexes, Real* &values) {
    int i1; // ������� ����� for
	for (i1=0; i1<n; i1++) {
		indexes[i1]=list[i1].key; 
		values[i1]=list[i1].val;  
	}
} // get_values