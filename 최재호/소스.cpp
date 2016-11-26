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

//�����迭(���A,B,dQ)
double A[DIM*DIM] = { 0 };
double B[DIM] = { 0 };
double dQ[DIM] = { 0 };
double C[DIM*DIM] = { 0 };
//������
typedef struct _Pipe
{
	double Diameter;
	int length;// ����� ����
	double Sigma_k;//�����սǰ���� ��
	double Extra_K;//Sigma K + Extra K
	double q;//������ ����� ���� Q*----> ��� �ٲ�
	double re;// ���̳����
	double f;//�������
	double le;//�����
	double r;//���װ��
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

void PrintMatrix(double *mat, int size); // ��� ǥ������ϱ�
double Cofactor(double *mat, int size, int i, int j); // i�� j���� ���� ���μ�
double Determinant(double *mat, int size); // ��Ľ�
double* MinorMat(double *majorMat, int size, int row, int col); // i�� j���� ���� ����� ���ϱ�
void CofactorMat(double *cfMat, double *oriMat, int size); // ���μ���� ���ϱ�
void TransposeMat(double *trMat, double *oriMat, int size); // ��ġ��� ���ϱ�
void InverseMat(double *inMat, double *oriMat, int size, double det); // ����� ���ϱ�
void MultiMat(double *multiMat, double *mat1, double *mat2, int size); // ����� ��

void main() {
	Pipe P1;
	Pipe P2;
	Pipe P3;
	double a, b, c;
	double q1,q2,q3 ; 
	double det, *myMat, *cfMat, *transMat, *invMat, *multiMat;
	int size=3;

	printf("Pipe �ʱ�ȭ............(1/3)\n");
	init(&P1);
	printf("\nPipe �ʱ�ȭ............(2/3)\n");
	init(&P2);
	printf("\nPipe �ʱ�ȭ............(3/3)\n");
	init(&P3);

	printf("\nPipe1�� ������ ������ ���� �����ּ��� : ");
	scanf("%lf", &q1);		printf("\n");
	P1.q =q1 ;
	printf("Pipe2�� ������ ������ ���� �����ּ��� : ");
	scanf("%lf", &q2);		printf("\n");
	P2.q = q2;
	printf("Pipe3�� ������ ������ ���� �����ּ��� : ");
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
		printf("���� ��� : \n");
		PrintMatrix(myMat, size); // �Է��� ��� ���
		printf("\n");
								  // ��Ľ�
		cfMat = (double*)malloc(sizeof(double)*size*size);
		CofactorMat(cfMat, myMat, size);//���μ� ���
		transMat = (double*)malloc(sizeof(double)*size*size);
		TransposeMat(transMat, cfMat, size);//��ġ���

		det = Determinant(myMat, size);//��ļ�
		printf("\n det(myMat) is %lf\n", det);

		invMat = (double*)malloc(sizeof(double)*size*size);
		InverseMat(C, transMat, size, det);//�����
		printf("\n�� ��� : \n");
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
//�ʱ�ȭ �Լ�
void init(Pipe * x) {
	//��������(�Լ������� �����)
	double Diameter;
	int length;// ����� ����
	double Sigma_k;//�����սǰ���� ��
	double Extra_K;//Sigma K + Extra K
	//
	printf("Diameter���� �Է����ּ��� : ");
	scanf("%lf", &Diameter);
	x->Diameter = Diameter;
	printf("\n");

	printf("����� ���̰��� �Է����ּ��� : ");
	scanf("%d", &length);
	x->length = length;
	printf("\n");

	printf("�����սǰ���� ���ǰ�(Sigma K)�� �Է����ּ��� : ");
	scanf("%lf", &Sigma_k);
	x->Sigma_k = Sigma_k;
	printf("\n");

	printf("Extra K���� �Է����ּ��� : ");
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
	printf("Diameter�� : ");
	printf("%f", x->Diameter);
	printf("\n");
	printf("����� ���̰� : ");
	printf("%9d", x->length);
	printf("\n");
	printf("�����սǰ���� ���ǰ�(Sigma K) : ");
	printf("%.9f", x->Sigma_k);
	printf("\n");
	printf("Extra K�� : ");
	printf("%.9f", x->Extra_K);
	printf("\n");
}

//���̳���� ���ϴ� �Լ�
void Re(Pipe * x)
{
	x->re = (4 * (x->q)) / (M_PI*nu*(x->Diameter));
	printf("���̳���� : %.9lf \n", x->re);
}
//��������� ���ϴ� �Լ�
void F(Pipe * x)
{
	double y = pow(((Epsilon / x->Diameter) / 3.7), (1.11));
	double z = (6.9) / x->re;
	double k = log(y + z);

	x->f = pow((-1.8*k),-2);
	printf("������� : %.9lf \n", x->f);
}
//����̸� ���ϴ� �Լ�
void LE(Pipe * x)
{
	x->le = (x->Diameter)*((x->Sk_sum_Ek) / (x->f));
	printf("����� : %.9lf \n", x->le);
}
//���װ���� ���ϴ� �Լ�
void R(Pipe * x)
{
	x->r = (8 * (x->f)*((x->length) + (x->le)) / (g * pow(M_PI,2)* pow((x->Diameter), 5) ));
	printf("���װ�� : %.9lf \n", x->r);
}
//�Լ���
void function(Pipe *x) {
	Re(x);
	F(x);
	LE(x);
	R(x);
}
//���

void PrintMatrix(double *mat, int size) // ��� ǥ�����
{
	int i;

	for (i = 0; i<size*size; i++)
	{
		printf("%.10lf ", mat[i]);
		if (i % size == size - 1) { printf("\n"); }
		// size x size ���·� ��� ���
	}
}

double Determinant(double *mat, int size) // ��Ľ�
{
	int j = 0; // ��ī��Ʈ ����j
	double det = 0; // ��Ľ� ����

	if (size == 1) { // ��Ļ���� 1��(��)�� ���Ҹ� �״�� ��ȯ(�̰�� ���ҿ� ��Ľ��� ��ġ�Ѵ�.)
		return mat[0];
	}
	else {
		for (j = 0; j<size; j++) {
			det += mat[j] * pow(-1, (1 + j + 1))*(Determinant(MinorMat(mat, size, 0, j), size - 1));
			//��Ľ� = 1��j�� ���� x 1��j�� ���μ� ����.
			//���⼭, 1��j�� ���μ� = (-1)^(1+j) x ����Ľ��� ��Ľ��̴�.
			//(-1)^(1+j)�� pow�Լ��� ����. ����Ľ��� ����� 1�̵ɶ����� ���ȣ���Ѵ�.
		}
	}
	return det; // ���� ��Ľ� ��ȯ
}

double Cofactor(double *mat, int size, int i, int j) // i�� j���� ���� ���μ�
{
	double cf;
	cf = pow(-1, (i + 1) + (j + 1))*(Determinant(MinorMat(mat, size, i, j), size - 1));
	//i��j�� ���μ� = (-1)^(i+j) x ����Ľ��� ��Ľ�.
	return cf;
}

double* MinorMat(double *mat, int size, int row, int col) // i�� j���� ���� �����
{
	int i, j, k = 0; // i:��ī��Ʈ, j:��ī��Ʈ, k:����� ���� ī��Ʈ
	double *minorMat = (double*)malloc(sizeof(double)*(size - 1)*(size - 1));
	// ������ĺ��� �ѻ����� ���� ����� ��������.
	for (i = 0; i<size; i++) {
		if (i == row) { continue; }
		for (j = 0; j<size; j++) {
			if (j == col) { continue; }
			minorMat[k] = mat[i*size + j];
			k++;
		}
	}// ��������� i�����,j�����Ҹ� ������ ���Ҹ� ����Ŀ� ä���.

	return minorMat; // ����� ��ȯ
}

void TransposeMat(double *trMat, double *oriMat, int size) // ��ġ��� ���ϱ�
{
	int i, j;

	for (i = 0; i<size; i++) {
		for (j = 0; j<size; j++) {
			trMat[j*size + i] = oriMat[i*size + j];
			// ��� ���� �ٲ��� ��ġ��Ŀ� ����.
		}
	}
}

void CofactorMat(double *cfMat, double *oriMat, int size) // ���μ���� ���ϱ�
{
	int i, j;
	for (i = 0; i<size; i++) {
		for (j = 0; j<size; j++) {
			cfMat[i*size + j] = Cofactor(oriMat, size, i, j);
			// i��j�� ������ ���μ��� ���ϰ� ���μ���Ŀ� �ִ´�.
		}
	}
}

void InverseMat(double *inMat, double *oriMat, int size, double det) // ����� ���ϱ�
{
	int i;
	for (i = 0; i<size*size; i++) {
		inMat[i] = oriMat[i] / det;
		//���ڷ� ���� ��������� �����Ҹ� ��Ľ����� �������� ����Ŀ� ����.
	}
}

void MultiMat(double *multiMat, double *mat1, double *mat2, int size) // ����� ��
{
	int i, j, q, k = 0; // i,j,q: ������� ��� ���� ī��Ʈ�� ����. k:�������� ����ī��Ʈ��
	for (i = 0; i<size*size; i++) {
		multiMat[i] = 0; // ���������� ���Ҹ� 0���� �ʱ�ȭ�Ѵ�.
	}
	for (i = 0; i<size; i++) {
		for (j = 0; j<size; j++) {
			for (q = 0; q<size; q++) {
				multiMat[k] += mat1[i*size + q] * mat2[q*size + j];
				// ��������=���1�� i��q�� x ���2�� q��j�� ���� ��.
			}
			k++;
		}
	}
}