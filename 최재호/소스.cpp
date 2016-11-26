#include<stdio.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<memory.h>
#include<stdlib.h>
#define nu	0.0000104
#define g	32.2
#define RHO	1.937
#define DENTISET	4608
#define DENTISET2	4608
#define DIM	3
#define DETERMIN_M	-0.000000001
#define DETERMIN_P	0.000000001
#define	Epsilon	0.00000500000004

//전역배열(행렬A,B,dQ)
double A[DIM*DIM] = { 0 };
double B[DIM] = { 0 };
double dQ[DIM] = { 0 };
double C[DIM*DIM] = { 0 };
//파이프
typedef struct _Pipe
{
	double Diameter;
	int length;// 배관의 길이
	double Sigma_k;//부차손실계수의 합
	double Extra_K;//Sigma K + Extra K
	double q;//가정한 배관의 유량 Q*----> 계속 바뀜
	double re;// 레이놀즈수
	double f;//마찰계수
	double le;//등가길이
	double r;//저항계수
	double dq;//DeltaQ
	double Sk_sum_Ek; //Sigma_k+Extra_K
	double result;  //Q
}Pipe;

void init(Pipe * x);
void print_Pipe(Pipe *x);
void Re(Pipe * x);
void F(Pipe * x);
void LE(Pipe * x);
void R(Pipe * x);
void function(Pipe *x);

void PrintMatrix(double *mat, int size); // 행렬 표준출력하기
double Cofactor(double *mat, int size, int i, int j); // i행 j열에 대한 여인수
double Determinant(double *mat, int size); // 행렬식
double* MinorMat(double *majorMat, int size, int row, int col); // i행 j열에 대한 소행렬 구하기
void CofactorMat(double *cfMat, double *oriMat, int size); // 여인수행렬 구하기
void TransposeMat(double *trMat, double *oriMat, int size); // 전치행렬 구하기
void InverseMat(double *inMat, double *oriMat, int size, double det); // 역행렬 구하기
void MultiMat(double *multiMat, double *mat1, double *mat2, int size); // 행렬의 곱

void main() {
	Pipe P1;
	Pipe P2;
	Pipe P3;
	double a, b, c;
	double q1,q2,q3 ; 
	double det, *myMat, *cfMat, *transMat, *invMat, *multiMat;
	int size=3;

	printf("Pipe 초기화............(1/3)\n");
	init(&P1);
	printf("\nPipe 초기화............(2/3)\n");
	init(&P2);
	printf("\nPipe 초기화............(3/3)\n");
	init(&P3);

	printf("\nPipe1의 가정한 유량의 값을 적어주세요 : ");
	scanf("%lf", &q1);		printf("\n");
	P1.q =q1 ;
	printf("Pipe2의 가정한 유량의 값을 적어주세요 : ");
	scanf("%lf", &q2);		printf("\n");
	P2.q = q2;
	printf("Pipe3의 가정한 유량의 값을 적어주세요 : ");
	scanf("%lf", &q3);		printf("\n");
	P3.q = q3;
	
	while (1)
	{
		function(&P1);
		function(&P2);
		function(&P3);
		
		A[0] = 2 * (P1.r)*(P1.q);
		A[1] = 2 * (P2.r)*(P2.q);
		A[2] = 0;  
		A[3] = 2 * (P1.r)*(P1.q);
		A[4] = 0;
		A[5] = 2 * (P3.r)*(P3.q);
		A[6] = 1;
		A[7] = -1;
		A[8] = -1;

		B[1] = (DENTISET/(RHO*g)) - (P1.r)*pow((P1.q),2) - (P2.r)*pow((P2.q),2);
		B[2] = (DENTISET2 / (RHO*g)) - (P1.r)*pow((P1.q), 2) - (P3.r)*pow((P3.q), 2);
		B[3] = 0.0;

		///////////////////
		myMat = A;
		printf("현재 행렬 : \n");
		PrintMatrix(myMat, size); // 입력한 행렬 출력
		printf("\n");
								  // 행렬식
		cfMat = (double*)malloc(sizeof(double)*size*size);
		CofactorMat(cfMat, myMat, size);//여인수 행렬
		transMat = (double*)malloc(sizeof(double)*size*size);
		TransposeMat(transMat, cfMat, size);//전치행렬

		det = Determinant(myMat, size);//행렬수
		printf("\n det(myMat) is %lf\n", det);

		invMat = (double*)malloc(sizeof(double)*size*size);
		InverseMat(C, transMat, size, det);//역행렬
		printf("\n역 행렬 : \n");
		PrintMatrix(C, size);
		printf("\n");
		////////////////////

		P1.dq = C[0] * B[1] + C[1] * B[2] + C[2] * B[3];
		P2.dq = C[3] * B[1] + C[4] * B[2] + C[5] * B[3];
		P3.dq = C[6] * B[1] + C[7] * B[2] + C[8] * B[3];
		
		P1.result = P1.dq + P1.q;
		P2.result = P2.dq + P2.q;
		P3.result = P3.dq + P3.q;

		a=fabs(P1.dq); b = fabs(P2.dq); c = fabs(P3.dq);

		printf("dQ1 = %lf, dQ2 = %lf, dQ3 = %lf\n",P1.dq, P2.dq, P3.dq);
		printf("\n");
	
		if ((a> DETERMIN_M) && (a< DETERMIN_P)&& (b> DETERMIN_M) && (b< DETERMIN_P) && (c> DETERMIN_M) && (c< DETERMIN_P))
		{
			printf("Q1 = %.15lf, Q2 = %.15lf, Q3 = %.15lf\n", P1.result, P2.result, P3.result);
			break;
		}
		else
		{
			printf("Q1 = %.9lf, Q2 = %.9lf, Q3 = %.9lf\n", P1.result, P2.result, P3.result);
			P1.q = P1.result;
			P2.q = P2.result;
			P3.q = P3.result;
		}

	}
	system("pause");

}
//초기화 함수
void init(Pipe * x) {
	//지역변수(함수끝나면 사라짐)
	double Diameter;
	int length;// 배관의 길이
	double Sigma_k;//부차손실계수의 합
	double Extra_K;//Sigma K + Extra K
	//
	printf("Diameter값을 입력해주세요 : ");
	scanf("%lf", &Diameter);
	x->Diameter = Diameter;
	printf("\n");

	printf("배관의 길이값을 입력해주세요 : ");
	scanf("%d", &length);
	x->length = length;
	printf("\n");

	printf("부차손실계수의 합의값(Sigma K)을 입력해주세요 : ");
	scanf("%lf", &Sigma_k);
	x->Sigma_k = Sigma_k;
	printf("\n");

	printf("Extra K값을 입력해주세요 : ");
	scanf("%lf", &Extra_K);
	x->Extra_K = Extra_K;
	printf("\n");

	x->q = 0;
	x->re = 0;
	x->f = 0;
	x->le = 0;
	x->r = 0;
	x->dq = 0;
	x->Sk_sum_Ek = x->Sigma_k + x->Extra_K;
	print_Pipe(x);

}
void print_Pipe(Pipe *x)
{
	printf("Diameter값 : ");
	printf("%f", x->Diameter);
	printf("\n");
	printf("배관의 길이값 : ");
	printf("%9d", x->length);
	printf("\n");
	printf("부차손실계수의 합의값(Sigma K) : ");
	printf("%.9f", x->Sigma_k);
	printf("\n");
	printf("Extra K값 : ");
	printf("%.9f", x->Extra_K);
	printf("\n");
}

//레이놀즈수 구하는 함수
void Re(Pipe * x)
{
	x->re = (4 * (x->q)) / (M_PI*nu*(x->Diameter));
	printf("레이놀즈수 : %.9lf \n", x->re);
}
//마찰계수를 구하는 함수
void F(Pipe * x)
{
	double y = pow(((Epsilon / x->Diameter) / 3.7), (1.11));
	double z = (6.9) / x->re;
	double k = log(y + z);

	x->f = pow((-1.8*k),-2);
	printf("마찰계수 : %.9lf \n", x->f);
}
//등가길이를 구하는 함수
void LE(Pipe * x)
{
	x->le = (x->Diameter)*((x->Sk_sum_Ek) / (x->f));
	printf("등가길이 : %.9lf \n", x->le);
}
//저항계수를 구하는 함수
void R(Pipe * x)
{
	x->r = (8 * (x->f)*((x->length) + (x->le)) / (g * pow(M_PI,2)* pow((x->Diameter), 5) ));
	printf("저항계수 : %.9lf \n", x->r);
}
//함수합
void function(Pipe *x) {
	Re(x);
	F(x);
	LE(x);
	R(x);
}
//행렬

void PrintMatrix(double *mat, int size) // 행렬 표준출력
{
	int i;

	for (i = 0; i<size*size; i++)
	{
		printf("%.10lf ", mat[i]);
		if (i % size == size - 1) { printf("\n"); }
		// size x size 형태로 행렬 출력
	}
}

double Determinant(double *mat, int size) // 행렬식
{
	int j = 0; // 열카운트 변수j
	double det = 0; // 행렬식 저장

	if (size == 1) { // 행렬사이즈가 1이(되)면 원소를 그대로 반환(이경우 원소와 행렬식이 일치한다.)
		return mat[0];
	}
	else {
		for (j = 0; j<size; j++) {
			det += mat[j] * pow(-1, (1 + j + 1))*(Determinant(MinorMat(mat, size, 0, j), size - 1));
			//행렬식 = 1행j열 원소 x 1행j열 여인수 총합.
			//여기서, 1행j열 여인수 = (-1)^(1+j) x 소행렬식의 행렬식이다.
			//(-1)^(1+j)는 pow함수로 구현. 소행렬식의 사이즈가 1이될때까지 재귀호출한다.
		}
	}
	return det; // 구한 행렬식 반환
}

double Cofactor(double *mat, int size, int i, int j) // i행 j열에 대한 여인수
{
	double cf;
	cf = pow(-1, (i + 1) + (j + 1))*(Determinant(MinorMat(mat, size, i, j), size - 1));
	//i행j열 여인수 = (-1)^(i+j) x 소행렬식의 행렬식.
	return cf;
}

double* MinorMat(double *mat, int size, int row, int col) // i행 j열에 대한 소행렬
{
	int i, j, k = 0; // i:행카운트, j:열카운트, k:소행렬 원소 카운트
	double *minorMat = (double*)malloc(sizeof(double)*(size - 1)*(size - 1));
	// 기존행렬보다 한사이즈 작은 소행렬 동적생성.
	for (i = 0; i<size; i++) {
		if (i == row) { continue; }
		for (j = 0; j<size; j++) {
			if (j == col) { continue; }
			minorMat[k] = mat[i*size + j];
			k++;
		}
	}// 기존행렬의 i행원소,j열원소를 제외한 원소를 소행렬에 채운다.

	return minorMat; // 소행렬 반환
}

void TransposeMat(double *trMat, double *oriMat, int size) // 전치행렬 구하기
{
	int i, j;

	for (i = 0; i<size; i++) {
		for (j = 0; j<size; j++) {
			trMat[j*size + i] = oriMat[i*size + j];
			// 행과 열을 바꾼후 전치행렬에 저장.
		}
	}
}

void CofactorMat(double *cfMat, double *oriMat, int size) // 여인수행렬 구하기
{
	int i, j;
	for (i = 0; i<size; i++) {
		for (j = 0; j<size; j++) {
			cfMat[i*size + j] = Cofactor(oriMat, size, i, j);
			// i행j열 원소의 여인수를 구하고 여인수행렬에 넣는다.
		}
	}
}

void InverseMat(double *inMat, double *oriMat, int size, double det) // 역행렬 구하기
{
	int i;
	for (i = 0; i<size*size; i++) {
		inMat[i] = oriMat[i] / det;
		//인자로 받은 수반행렬의 각원소를 행렬식으로 나눈값을 역행렬에 저장.
	}
}

void MultiMat(double *multiMat, double *mat1, double *mat2, int size) // 행렬의 곱
{
	int i, j, q, k = 0; // i,j,q: 두행렬의 행과 열을 카운트할 변수. k:곱결과행렬 원소카운트용
	for (i = 0; i<size*size; i++) {
		multiMat[i] = 0; // 곱결과행렬의 원소를 0으로 초기화한다.
	}
	for (i = 0; i<size; i++) {
		for (j = 0; j<size; j++) {
			for (q = 0; q<size; q++) {
				multiMat[k] += mat1[i*size + q] * mat2[q*size + j];
				// 곱결과행렬=행렬1의 i행q열 x 행렬2의 q행j열 들의 합.
			}
			k++;
		}
	}
}