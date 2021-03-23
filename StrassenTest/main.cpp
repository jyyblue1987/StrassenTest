#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 


using namespace std; 
typedef int lld; 
  
/* Strassen's Algorithm for matrix multiplication 
   Complexity:    O(n^2.808) */

#define N_0 8
lld** AllocMatrix(int n, int m)
{
	lld *x = new lld[n * m];
	lld **p = new lld*[n];
	memset(x, 0, sizeof(lld) * n * n);

	for(int i = 0; i < n; i++)
		p[i] = x + i * m;

	return p;
}

void FreeMatrix(lld **p)
{
	delete p[0];
	delete p;
}
  
inline lld** MatrixMultiply(lld** a, lld** b, int n, 
                                      int l, int m) 
{ 
    lld** c = AllocMatrix(n, m);

	lld sum = 0;

	lld *aa = NULL;
	lld *bb = NULL;
  
    for (int i = 0; i < n; i++) { 
		for (int j = 0; j < m; j++) { 
            sum = 0; 			
			aa = a[i];
			bb = b[0] + j;
            for (int k = 0; k < l; k++, aa++, bb+=m)
                sum += (*aa) * (*bb); 				
			c[i][j] = sum; 
        } 
    } 
    return c; 
} 

inline void MatrixAdd(lld *a, lld *b, lld *c, int n, int m)
{
	int len = n * m;
	for(int i = 0; i < len; i++, a++, b++, c++)
		*c = *a + *b;
}

inline void MatrixSubstract(lld *a, lld *b, lld *c, int n, int m)
{
	int len = n * m;
	for(int i = 0; i < len; i++, a++, b++, c++)
		*c = *a - *b;
}

  
inline lld** Strassen(lld** a, lld** b, int n,  
                                int l, int m) 
{ 
    if (n <= N_0 || l <= N_0 || m <= N_0 )  
		//if (n == || l == 1 || m == 1)
        return MatrixMultiply(a, b, n, l, m); 

	lld** c = AllocMatrix(n, m);

	int adjN = (n >> 1) + (n & 1); 
	int adjL = (l >> 1) + (l & 1); 
	int adjM = (m >> 1) + (m & 1); 

	lld**** As = new lld***[2]; 
	lld **A = NULL;
	lld *aa = NULL;
	int len = 0;
	
	for (int x = 0; x < 2; x++) { 
		As[x] = new lld**[2]; 
		for (int y = 0; y < 2; y++) { 
			A = AllocMatrix(adjN, adjL);				
			As[x][y] = A;			
			aa = A[0];
			int sx = x * adjN;
			int ex = (x + 1) * adjN;
			if( ex > n )
				ex = n;

			int sy = y * adjL;
			int ey = (y + 1) * adjL;
	
			if( ey > l )
				ey = l;

			len = ey - sy;

			for (int i = sx; i < ex; i++) { 				
				memcpy(aa, a[i] + sy, sizeof(lld) * len);
				aa += adjL;
			} 
		} 
	} 

	lld**** Bs = new lld***[2]; 
	for (int x = 0; x < 2; x++) { 
		Bs[x] = new lld**[2]; 
		for (int y = 0; y < 2; y++) { 			
			A = AllocMatrix(adjL, adjM);				
			Bs[x][y] = A;			
			aa = A[0];
			int sx = x * adjL;
			int ex = (x + 1) * adjL;			
			if( ex > l )
				ex = l;

			int sy = y * adjM;
			int ey = (y + 1) * adjM;
		
			if( ey > m )
				ey = m;

			len = ey - sy;

			for (int i = sx; i < ex; i++) { 
				memcpy(aa, b[i] + sy, sizeof(lld) * len);
				aa += adjM;
			} 
		} 
	} 

	lld*** s = new lld**[10]; 
	for (int i = 0; i < 10; i++) { 
		switch (i) { 
		case 0: 			
			s[i] = AllocMatrix(adjL, adjM);
			MatrixSubstract(Bs[0][1][0], Bs[1][1][0], s[i][0], adjL, adjM);			
			break; 
		case 1: 			
			s[i] = AllocMatrix(adjN, adjL);
			MatrixAdd(As[0][0][0], As[0][1][0], s[i][0], adjN, adjL);			 
			break; 
		case 2: 			
			s[i] = AllocMatrix(adjN, adjL);
			MatrixAdd(As[1][0][0], As[1][1][0], s[i][0], adjN, adjL);			  
			break; 
		case 3: 
			s[i] = AllocMatrix(adjL, adjM);
			MatrixSubstract(Bs[1][0][0], Bs[0][0][0], s[i][0], adjL, adjM);			
			break; 
		case 4: 
			s[i] = AllocMatrix(adjN, adjL);			
			MatrixAdd(As[0][0][0], As[1][1][0], s[i][0], adjN, adjL);			  
			break; 
		case 5: 
			s[i] = AllocMatrix(adjL, adjM);			
			MatrixAdd(Bs[0][0][0], Bs[1][1][0], s[i][0], adjL, adjM);			 
			break; 
		case 6: 
			s[i] = AllocMatrix(adjN, adjL);			
			MatrixSubstract(As[0][1][0], As[1][1][0], s[i][0], adjN, adjL);			  
			break; 
		case 7: 
			s[i] = AllocMatrix(adjL, adjM);			
			MatrixAdd(Bs[1][0][0], Bs[1][1][0], s[i][0], adjL, adjM);			 
			break; 
		case 8: 
			s[i] = AllocMatrix(adjN, adjL);			
			MatrixSubstract(As[0][0][0], As[1][0][0], s[i][0], adjN, adjL);			  
			break; 
		case 9: 			
			s[i] = AllocMatrix(adjL, adjM);		
			MatrixAdd(Bs[0][0][0], Bs[0][1][0], s[i][0], adjL, adjM);			  
			break; 
		} 
	} 

	lld*** p = new lld**[7]; 
	p[0] = Strassen(As[0][0], s[0], adjN, adjL, adjM); 
	p[1] = Strassen(s[1], Bs[1][1], adjN, adjL, adjM); 
	p[2] = Strassen(s[2], Bs[0][0], adjN, adjL, adjM); 
	p[3] = Strassen(As[1][1], s[3], adjN, adjL, adjM); 
	p[4] = Strassen(s[4], s[5], adjN, adjL, adjM); 
	p[5] = Strassen(s[6], s[7], adjN, adjL, adjM); 
	p[6] = Strassen(s[8], s[9], adjN, adjL, adjM); 

	lld *p0 = p[0][0], *p1 = p[1][0], *p2 = p[2][0], *p3 = p[3][0], *p4 = p[4][0], *p5 = p[5][0], *p6 = p[6][0];
	for (int i = 0; i < adjN; i++) { 
		for (int j = 0; j < adjM; j++, p0++, p1++, p2++, p3++, p4++, p5++, p6++) { 
			c[i][j] = *p4 + *p3 - *p1 + *p5; 
			if (j + adjM < m) 
				c[i][j + adjM] = *p0 + *p1; 
			if (i + adjN < n) 
				c[i + adjN][j] = *p2 + *p3;
			if (i + adjN < n && j + adjM < m) 
				c[i + adjN][j + adjM] = *p4 + *p0 - *p2 - *p6;
		} 
	} 

	for (int x = 0; x < 2; x++) { 
		for (int y = 0; y < 2; y++) { 
			FreeMatrix(Bs[x][y]);
		} 
		delete[] Bs[x]; 
	} 
	delete[] Bs;

	for (int x = 0; x < 2; x++) { 
		for (int y = 0; y < 2; y++) { 
			FreeMatrix(As[x][y]);
		} 
		delete[] As[x]; 
	} 
	delete[] As;

	for (int i = 0; i < 10; i++) { 
		FreeMatrix(s[i]);
	} 
	delete[] s; 

	for (int i = 0; i < 7; i++) { 
		FreeMatrix(p[i]);
	} 

	delete[] p; 
   
    return c; 
} 

double T_N(double n, double p, double q);

double T_N(double n, double p, double q) 
{
	if(n <= 2)
		return p * n * n;

	return 7 * T_N(n/2, p, q) + q * n * n;
}

int getCrossOverPointAanlytically(double p, double q) 
{
	int n = N_0;
	double s_n = p * N_0 * N_0 * N_0;
	double t_n = s_n;

	while(true) 
	{
		n += 2;
		t_n = T_N(n, p, q);
		s_n = p * n * n * n;

		if( t_n < s_n )	// cross over point
			break;
	}

	return n;
}
  
int getCrossOverPointExperiemently() 
{
	int n = N_0;
	while(true)
	{
		n *= 2;

		lld **A = AllocMatrix(n, n);
		lld **B = AllocMatrix(n, n);

		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				A[i][j] = rand() % 2;
				B[i][j] = rand() % 2;
			}
		}

		clock_t st = clock(); 
		lld **C_S = MatrixMultiply(A, B, n, n, n);
		clock_t et = clock(); 
		
		clock_t st1 = clock(); 
		lld **C_T = Strassen(A, B, n, n, n);
		clock_t et1 = clock(); 

		// check result is same
		bool flag = false;
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				if( C_S[i][j] != C_T[i][j] )
				{
					flag = true;
					break;
				}

			}
		}
		if( flag == true )
			printf("Result is different\n");

		FreeMatrix(C_S);
		FreeMatrix(C_T);

		FreeMatrix(A);
		FreeMatrix(B);

		printf("n = %d: Strassen = %dms, Standard = %dms\n", n, et1 - st1, et - st);
		if( et1 - st1 + 100 < et - st )	// cross over point
			break;

	}

	return n;
}

int main(int argc, char *argv[]) 
{ 
	if( argc > 3 && atoi(argv[1]) == 0 )
	{
		int n = atoi(argv[2]);

		lld** matA = AllocMatrix(n, n);
		lld** matB = AllocMatrix(n, n);

		// read file
		char * line = NULL;
		size_t len = 0;
		int read = 0;

		FILE *fp = fopen(argv[3], "rt");
		if( fp == NULL )
			return -1;

		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				fscanf(fp, "%d\n", &(matA[i][j]));

		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				fscanf(fp, "%d\n", &(matB[i][j]));
		fclose(fp);

		lld** matC = Strassen(matA, matB, n, n, n); 

		for(int i = 0; i < n; i++)
			printf("%d\n", matC[i][i]);

		FreeMatrix(matA);
		FreeMatrix(matB);
		FreeMatrix(matC);
	}
	else
	{
		// anaylze cross over point
		// standard: S(n) = p * n ^ 3
		// Strassen: T(n) = 7 * T(n/2) + q * n ^2
		int cross_over_point = getCrossOverPointAanlytically(5, 20.0f);
		printf("Therical Cross Over Point = %d\n\n\n", cross_over_point);

		int cross_over_point1 = getCrossOverPointExperiemently();
		printf("Experiement Cross Over Point = %d\n", cross_over_point1);

		lld** matA = AllocMatrix(2, 3); 
		matA[0][0] = 1; 
		matA[0][1] = 2; 
		matA[0][2] = 3; 
		matA[1][0] = 4; 
		matA[1][1] = 5; 
		matA[1][2] = 6; 
  
		lld** matB = AllocMatrix(3, 2); 
		matB[0][0] = 7; 
		matB[0][1] = 8; 
		matB[1][0] = 9; 
		matB[1][1] = 10; 
		matB[2][0] = 11; 
		matB[2][1] = 12; 
  
		lld** matC = Strassen(matA, matB, 2, 3, 2); 
		for (int i = 0; i < 2; i++) { 
			for (int j = 0; j < 2; j++) { 
				printf("%d ", matC[i][j]); 
			} 
			printf("\n"); 
		} 

		// experiement about cross point
		
	}
  
    return 0; 
} 