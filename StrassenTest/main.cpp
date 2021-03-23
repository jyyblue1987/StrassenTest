#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 


using namespace std; 
typedef int lld; 
  
/* Strassen's Algorithm for matrix multiplication 
   Complexity:    O(n^2.808) */


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

  
inline lld** Strassen(lld** a, lld** b, int n,  
                                int l, int m) 
{ 
    if (n == 1 || l == 1 || m == 1)  
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
			for (int j = 0; j < adjL; j++) { 				
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[0][1][j][k] - Bs[1][1][j][k]; 
				} 
			} 
			break; 
		case 1: 			
			s[i] = AllocMatrix(adjN, adjL);
			for (int j = 0; j < adjN; j++) { 				
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][0][j][k] + As[0][1][j][k]; 
				} 
			} 
			break; 
		case 2: 			
			s[i] = AllocMatrix(adjN, adjL);
			for (int j = 0; j < adjN; j++) { 				
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[1][0][j][k] + As[1][1][j][k]; 
				} 
			} 
			break; 
		case 3: 
			s[i] = AllocMatrix(adjL, adjM);
			for (int j = 0; j < adjL; j++) { 				
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[1][0][j][k] - Bs[0][0][j][k]; 
				} 
			} 
			break; 
		case 4: 
			s[i] = AllocMatrix(adjN, adjL);			
			for (int j = 0; j < adjN; j++) { 				
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][0][j][k] + As[1][1][j][k]; 
				} 
			} 
			break; 
		case 5: 
			s[i] = AllocMatrix(adjL, adjM);			
			for (int j = 0; j < adjL; j++) { 				
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[0][0][j][k] + Bs[1][1][j][k]; 
				} 
			} 
			break; 
		case 6: 
			s[i] = AllocMatrix(adjN, adjL);			
			for (int j = 0; j < adjN; j++) { 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][1][j][k] - As[1][1][j][k]; 
				} 
			} 
			break; 
		case 7: 
			s[i] = AllocMatrix(adjL, adjM);			
			for (int j = 0; j < adjL; j++) { 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[1][0][j][k] + Bs[1][1][j][k]; 
				} 
			} 
			break; 
		case 8: 
			s[i] = AllocMatrix(adjN, adjL);			
			for (int j = 0; j < adjN; j++) { 
				for (int k = 0; k < adjL; k++) { 
					s[i][j][k] = As[0][0][j][k] - As[1][0][j][k]; 
				} 
			} 
			break; 
		case 9: 			
			s[i] = AllocMatrix(adjL, adjM);		
			for (int j = 0; j < adjL; j++) { 
				for (int k = 0; k < adjM; k++) { 
					s[i][j][k] = Bs[0][0][j][k] + Bs[0][1][j][k]; 
				} 
			} 
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

	for (int i = 0; i < adjN; i++) { 
		for (int j = 0; j < adjM; j++) { 
			c[i][j] = p[4][i][j] + p[3][i][j] - p[1][i][j] + p[5][i][j]; 
			if (j + adjM < m) 
				c[i][j + adjM] = p[0][i][j] + p[1][i][j]; 
			if (i + adjN < n) 
				c[i + adjN][j] = p[2][i][j] + p[3][i][j]; 
			if (i + adjN < n && j + adjM < m) 
				c[i + adjN][j + adjM] = p[4][i][j] + p[0][i][j] - p[2][i][j] - p[6][i][j]; 
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

int getCrossOverPointAanlytically(double p, double q) 
{
	int n = 1;
	double s_n = p;
	double t_n = s_n;
	while(true) 
	{
		n *= 2;
		t_n = 7 * t_n + q * n * n;
		s_n = p * n * n * n;

		if( t_n < s_n )	// cross over point
			break;
	}

	return n;
}
  
int getCrossOverPointExperiemently() 
{
	int n = 1;
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
		printf("Experiement Cross Over Point = %d\n", cross_over_point);

		lld** matA; 
		matA = new lld*[2]; 
		for (int i = 0; i < 2; i++) 
			matA[i] = new lld[3]; 
		matA[0][0] = 1; 
		matA[0][1] = 2; 
		matA[0][2] = 3; 
		matA[1][0] = 4; 
		matA[1][1] = 5; 
		matA[1][2] = 6; 
  
		lld** matB; 
		matB = new lld*[3]; 
		for (int i = 0; i < 3; i++) 
			matB[i] = new lld[2]; 
		matB[0][0] = 7; 
		matB[0][1] = 8; 
		matB[1][0] = 9; 
		matB[1][1] = 10; 
		matB[2][0] = 11; 
		matB[2][1] = 12; 
  
		lld** matC = Strassen(matA, matB, 2, 3, 2); 
		for (int i = 0; i < 2; i++) { 
			for (int j = 0; j < 2; j++) { 
				printf("%lld ", matC[i][j]); 
			} 
			printf("\n"); 
		} 

		// experiement about cross point
		
	}
  
    return 0; 
} 