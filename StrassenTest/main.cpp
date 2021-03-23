#include <stdio.h>
#include <stdlib.h>
using namespace std; 
typedef int lld; 
  
/* Strassen's Algorithm for matrix multiplication 
   Complexity:    O(n^2.808) */
  
inline lld** MatrixMultiply(lld** a, lld** b, int n, 
                                      int l, int m) 
{ 
    lld** c = new lld*[n]; 
    for (int i = 0; i < n; i++) 
        c[i] = new lld[m]; 

	lld sum = 0;
  
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < m; j++) { 
            sum = 0; 
            for (int k = 0; k < l; k++)
                sum += a[i][k] * b[k][j]; 
			c[i][j] = sum; 
        } 
    } 
    return c; 
} 

lld** AllocMatrix(int n, int m)
{
	lld *x = new lld[n * m];
	lld **p = new lld*[n];

	for(int i = 0; i < n; i++)
		p[i] = x + i * m;

	return p;
}

void FreeMatrix(lld **p)
{
	delete p[0];
	delete p;
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
	for (int x = 0; x < 2; x++) { 
		As[x] = new lld**[2]; 
		for (int y = 0; y < 2; y++) { 
			As[x][y] = AllocMatrix(adjN, adjL);
			for (int i = 0; i < adjN; i++) { 
				for (int j = 0; j < adjL; j++) { 
					int I = i + (x & 1) * adjN; 
					int J = j + (y & 1) * adjL; 
					As[x][y][i][j] = (I < n && J < l) ? a[I][J] : 0; 
				} 
			} 
		} 
	} 

	lld**** Bs = new lld***[2]; 
	for (int x = 0; x < 2; x++) { 
		Bs[x] = new lld**[2]; 
		for (int y = 0; y < 2; y++) { 			
			Bs[x][y] = AllocMatrix(adjL, adjM);
			for (int i = 0; i < adjL; i++) { 				
				for (int j = 0; j < adjM; j++) { 
					int I = i + (x & 1) * adjL; 
					int J = j + (y & 1) * adjM; 
					Bs[x][y][i][j] = (I < l && J < m) ? b[I][J] : 0; 
				} 
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
	}
  
    return 0; 
} 