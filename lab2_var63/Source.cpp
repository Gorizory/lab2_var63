//#define LINEAR
#define CUBIC

#define _USE_MATH_DEFINES

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

#define N_ELEMENTS		40
#ifdef LINEAR
	#define LOCAL_MATR_SIZE	2
#endif
#ifdef CUBIC
	#define LOCAL_MATR_SIZE	4
#endif
#define N_NODES N_ELEMENTS * (LOCAL_MATR_SIZE - 1) + 1
#define X_START			2.0	
#define X_FINISH		10.0	

#define LENGHT (X_FINISH - X_START)
#define ELEMENT_LENGHT  (LENGHT / N_ELEMENTS)

#define LEFT_BOUNDARY_CONDITION -2
#define RIGHT_BOUNDARY_CONDITION 10

using namespace std;

double localMatrix[LOCAL_MATR_SIZE][LOCAL_MATR_SIZE];
double **globalMatrix;
double *vectorB;
double vectorU[N_NODES + 1];

double calcFunc(double x) {
	return (2. * pow(M_E, (-sqrt(31. / 110.) * x - sqrt(62. / 55.))) * (-sqrt(3410.) * pow(M_E, (sqrt(62. / 55.) * x))
		+ 2. * pow(M_E, (sqrt(31. / 110.) * x + sqrt(62. / 55.))) + 2. * pow(M_E, (sqrt(31. / 110.) * x + 9. * sqrt(62. / 55.)))
		+ 153. * pow(M_E, (sqrt(62. / 55.) * x + 4. * sqrt(62. / 55.))) + 153. * pow(M_E, (6. * sqrt(62. / 55.)))
		+ sqrt(3410.) * pow(M_E, (2. * sqrt(310. / 11.)))) / (31. * (1. + pow(M_E, (8. * sqrt(62. / 55.))))));
}

void initLocalMatrix()
{
#ifdef LINEAR
	localMatrix[0][0] = -110.0 / ELEMENT_LENGHT - 31.0 / 3.0 * ELEMENT_LENGHT;
	localMatrix[0][1] = 110.0 / ELEMENT_LENGHT - 31.0 / 6.0 * ELEMENT_LENGHT;
	localMatrix[1][0] = 110.0 / ELEMENT_LENGHT - 31.0 / 6.0 * ELEMENT_LENGHT;
	localMatrix[1][1] = -110.0 / ELEMENT_LENGHT - 31.0 / 3.0 * ELEMENT_LENGHT;
#endif
#ifdef CUBIC
	localMatrix[0][0] = -407.0 / ELEMENT_LENGHT - 248.0 / 105. * ELEMENT_LENGHT;
	localMatrix[0][1] = 2079.0 / 4. / ELEMENT_LENGHT - 1023.0 / 560.0 * ELEMENT_LENGHT;
	localMatrix[0][2] = -297.0 / 2. / ELEMENT_LENGHT + 93.0 / 140.0 * ELEMENT_LENGHT;
	localMatrix[0][3] = 143.0 / 4. / ELEMENT_LENGHT - 589.0 / 1680.0 * ELEMENT_LENGHT;

	localMatrix[1][0] = 2079.0 / 4. / ELEMENT_LENGHT - 1023.0 / 560.0 * ELEMENT_LENGHT;
	localMatrix[1][1] = -1188.0 / ELEMENT_LENGHT - 837.0 / 70.0 * ELEMENT_LENGHT;
	localMatrix[1][2] = 3267.0 / 4. / ELEMENT_LENGHT + 837.0 / 560.0 * ELEMENT_LENGHT;
	localMatrix[1][3] = -297.0 / 2. / ELEMENT_LENGHT + 93.0 / 140.0 * ELEMENT_LENGHT;

	localMatrix[2][0] = -297.0 / 2. / ELEMENT_LENGHT + 93.0 / 140.0 * ELEMENT_LENGHT;
	localMatrix[2][1] = 3267.0 / 4. / ELEMENT_LENGHT + 837.0 / 560.0 * ELEMENT_LENGHT;
	localMatrix[2][2] = -1188.0 / ELEMENT_LENGHT - 837.0 / 70.0 * ELEMENT_LENGHT;
	localMatrix[2][3] = 2079.0 / 4. / ELEMENT_LENGHT - 1023.0 / 560.0 * ELEMENT_LENGHT;

	localMatrix[3][0] = 143.0 / 4. / ELEMENT_LENGHT - 589.0 / 1680.0 * ELEMENT_LENGHT;
	localMatrix[3][1] = -297.0 / 2. / ELEMENT_LENGHT + 93.0 / 140.0 * ELEMENT_LENGHT;
	localMatrix[3][2] = 2079.0 / 4. / ELEMENT_LENGHT - 1023.0 / 560.0 * ELEMENT_LENGHT;
	localMatrix[3][3] = -407.0 / ELEMENT_LENGHT - 248.0 / 105. * ELEMENT_LENGHT;
#endif
}

void initGlobalMatrix()
{
	globalMatrix = (double**)malloc(sizeof(double*)*(N_NODES));
	for (int i = 0; i < N_NODES; i++)
	{
		globalMatrix[i] = (double *)malloc(sizeof(double)*(N_NODES));
	}
	for (int i = 0; i < N_NODES; i++)
	{
		for (int j = 0; j < N_NODES; j++) {
			globalMatrix[i][j] = 0;
		}
	}

#ifdef LINEAR
	for (int i = 0; i < N_NODES - 1; i++)
	{
		globalMatrix[i][i] += localMatrix[0][0];
		globalMatrix[i][i + 1] += localMatrix[0][1];
		globalMatrix[i + 1][i] += localMatrix[1][0];
		globalMatrix[i + 1][i + 1] += localMatrix[1][1];
	}
	for (int i = 0; i < N_NODES - 1; i++) {
		globalMatrix[N_NODES - 1][i] = 0;
	}
	globalMatrix[N_NODES - 1][N_NODES - 1] = 1;
#endif
#ifdef CUBIC
	for (int i = 0; i < N_ELEMENTS; i++)
	{
		globalMatrix[3 * i][3 * i] += localMatrix[0][0];
		globalMatrix[3 * i][3 * i + 1] += localMatrix[0][1];
		globalMatrix[3 * i][3 * i + 2] += localMatrix[0][2];
		globalMatrix[3 * i][3 * i + 3] += localMatrix[0][3];
		
		globalMatrix[3 * i + 1][3 * i] += localMatrix[1][0];
		globalMatrix[3 * i + 1][3 * i + 1] += localMatrix[1][1];
		globalMatrix[3 * i + 1][3 * i + 2] += localMatrix[1][2];
		globalMatrix[3 * i + 1][3 * i + 3] += localMatrix[1][3];

		globalMatrix[3 * i + 2][3 * i] += localMatrix[2][0];
		globalMatrix[3 * i + 2][3 * i + 1] += localMatrix[2][1];
		globalMatrix[3 * i + 2][3 * i + 2] += localMatrix[2][2];
		globalMatrix[3 * i + 2][3 * i + 3] += localMatrix[2][3];

		globalMatrix[3 * i + 3][3 * i] += localMatrix[3][0];
		globalMatrix[3 * i + 3][3 * i + 1] += localMatrix[3][1];
		globalMatrix[3 * i + 3][3 * i + 2] += localMatrix[3][2];
		globalMatrix[3 * i + 3][3 * i + 3] += localMatrix[3][3];
	}
	for (int i = 0; i < N_NODES - 1; i++) {
		globalMatrix[N_NODES - 1][i] = 0;
	}
	globalMatrix[N_NODES - 1][N_NODES - 1] = 1;

	for (int i = 0; i < N_NODES - 1; i++) {
		globalMatrix[(N_NODES - 1) / 2][i] = 0;
	}
	globalMatrix[(N_NODES - 1) / 2][(N_NODES - 1) / 2] = 1;
#endif
}

void initVectorB()
{
	vectorB = (double *)malloc(sizeof(double)*(N_NODES));
	for (int i = 0; i < N_NODES; i++)
	{
		vectorB[i] = 0;
	}
#ifdef LINEAR
	vectorB[0] =  110 * LEFT_BOUNDARY_CONDITION - 2 * ELEMENT_LENGHT;
	for (int i = 1; i < N_NODES - 1; i++) {
		vectorB[i] = -4 * ELEMENT_LENGHT;
	}
	vectorB[N_NODES - 1] = RIGHT_BOUNDARY_CONDITION;
#endif
#ifdef CUBIC
	for (int i = 0; i < N_ELEMENTS; i++) {
		vectorB[3 * i] += -1. / 2. * ELEMENT_LENGHT;
		vectorB[3 * i + 1] += -3. / 2. * ELEMENT_LENGHT;
		vectorB[3 * i + 2] += -3. / 2. * ELEMENT_LENGHT;
		vectorB[3 * i + 3] += -1. / 2. * ELEMENT_LENGHT;
	}
	vectorB[0] += 110 * LEFT_BOUNDARY_CONDITION;
	vectorB[N_NODES - 1] = RIGHT_BOUNDARY_CONDITION;
	vectorB[(N_NODES - 1) / 2] = 12;
#endif
}

void printMatrix(double **A, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			printf("% 4.2f\t", A[i][j]);
		}
		putchar('\n');
	}
}

void printVector(double *A, int size)
{
	for (int j = 0; j < size; j++)
	{
		printf("% 4.2f\n", A[j]);
	}
}

void printResult(double *A, int size)
{
	double maxDelta = 0.;
	ofstream fout("output.csv");
	double x = X_START;
	fout << "x;" << "u counted;" << "u absolute;" << "delta" << endl;
	for (int j = 0; j < size; j++, x += ELEMENT_LENGHT / (LOCAL_MATR_SIZE - 1))
	{
		double u = calcFunc(x);
		double delta = fabs(A[j] - u);
		if (delta > maxDelta) {
			maxDelta = delta;
		}
		fout << x << ";" << A[j] << ";" << u << ";" << delta << ";" << endl;
	}
	fout << ";; max:;" << maxDelta << endl;
	fout.close();
}


int gauss(double **matrA, double *matrB, double *u, int n)
{
	int i, j, row;
	double k;

	for (row = 0; row < n; row++)
	{
		for (i = 0; i < row; i++)
		{
			k = matrA[row][i];
			if (k != 0)
				for (j = 0; j < n; j++)
					matrA[row][j] -= k * matrA[i][j];
			matrB[row] -= k * matrB[i];
		}
		k = matrA[row][row];
		for (j = row; j < n; j++)
			matrA[row][j] /= k;
		matrB[row] /= k;
	}

	for (i = n - 1; i > -1; i--)
	{
		for (j = n - 1; j > i; j--)
		{
			matrB[i] -= matrA[i][j] * u[j];
		}
		u[i] = matrB[i];
	}
	return 0;
}

int main()
{
	initLocalMatrix();
	initGlobalMatrix();
	//printf("Global matrix\n");
	//printMatrix(globalMatrix, N_NODES);
	//printf("Vector b\n");
	initVectorB();
	//printVector(vectorB, N_NODES);
	gauss(globalMatrix, vectorB, vectorU, N_NODES);
	printResult(vectorU, N_NODES);
	
	//getchar();
	return 0;
}

