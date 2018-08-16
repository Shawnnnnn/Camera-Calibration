//3D直接线性变换DLT
//
#include "stdafx.h"
#include "DLT.h"
#include <math.h>
#include <vector>

DLT::DLT(void)
{
	x0 = 0;
	y0 = 0;
	fx = 0;
	fy = 0;
	Xs = 0;
	Ys = 0;
	Zs = 0;
	phi = 0;
	omega = 0;
	kappa = 0;
}

DLT::~DLT(void)
{

}
// 方阵求逆函数1,Psr为输入方阵，N为方阵的大小,P为返回的逆矩阵
std::vector<std::vector<double>> inverse( const std::vector<std::vector<double>>& Psr )
{
	double d = 0.0;
	double temp = 0.0;
	int N = Psr.size();//矩阵的行数
	int N1 = Psr[0].size();//矩阵的列数
	std::vector<std::vector<double>> P(N);	
	for (int i = 0; i < N; i++)
	{	P[i].resize(N1);
	for (int j = 0; j < N1; j++)
		P[i][j] = Psr[i][j];
	}
	int *main_row = new int[N];
	int *main_col = new int[N];
	if (N == N1)//方阵才能求逆矩阵
	{
		for (int k = 0; k <= N - 1; k++)
		{
			d = 0.0;
			for (int i = k; i <= N - 1; i++)
			{
				for (int j = k; j <= N - 1; j++)
				{
					temp = abs(P[i][j]);
					if (temp > d)
					{
						d = temp; main_row[k] = i; main_col[k] = j;
					}
				}
			}
			if (main_row[k] != k)
			{
				for (int j = 0; j <= N - 1; j++)
				{
					temp = P[main_row[k]][j]; P[main_row[k]][j] = P[k][j]; P[k][j] = temp;
				}
			}
			if (main_col[k] != k)
			{
				for (int i = 0; i <= N - 1; i++)
				{
					temp = P[i][main_col[k]]; P[i][main_col[k]] = P[i][k]; P[i][k] = temp;
				}
			}
			P[k][k] = 1.0 / P[k][k];
			for (int j = 0; j <= N - 1; j++)
			{
				if (j != k)
				{
					P[k][j] *= P[k][k];
				}
			}

			for (int i = 0; i <= N - 1; i++)
			{
				if (i != k)
				{
					for (int j = 0; j <= N - 1; j++)
					{
						if (j != k)
						{
							P[i][j] -= P[i][k] * P[k][j];
						}
					}
				}
			}
			for (int i = 0; i <= N - 1; i++)
			{
				if (i != k)
				{
					P[i][k] = -P[i][k] * P[k][k];
				}
			}
		}
		for (int k = N - 1; k >= 0; k--)
		{
			if (main_col[k] != k)
			{
				for (int j = 0; j <= N - 1; j++)
				{
					temp = P[k][j]; P[k][j] = P[main_col[k]][j]; P[main_col[k]][j] = temp;
				}
			}
			if (main_row[k] != k)
			{
				for (int i = 0; i <= N - 1; i++)
				{
					temp = P[i][k]; P[i][k] = P[i][main_row[k]]; P[i][main_row[k]] = temp;
				}
			}
		}
	}
	return P;
}

void DLT::DLT_6pts(int dgs, double* imageX, double* imageY, double* spaceX, double* spaceY, double* spaceZ, double* cImageX, double* cImageY)
{
	double *L = new double[11];
	double *initMatrixM  = new double[22*6];
	double *tInitMatrixM = new double[22*6];
	double *initMatrixW  = new double[2*6];

	double initMatrixMTM[121] = {0.0};
	double initMatrixMTW[11]  = {0.0};

	//L系数的初始值计算
	for (int i = 0; i < 6; i++)
	{
		initMatrixM[ 0+22*i] = spaceX[i];
		initMatrixM[ 1+22*i] = spaceY[i];
		initMatrixM[ 2+22*i] = spaceZ[i];
		initMatrixM[ 3+22*i] = 1;
		initMatrixM[ 4+22*i] = 0;
		initMatrixM[ 5+22*i] = 0;
		initMatrixM[ 6+22*i] = 0;
		initMatrixM[ 7+22*i] = 0;
		initMatrixM[ 8+22*i] = imageX[i] * spaceX[i];
		initMatrixM[ 9+22*i] = imageX[i] * spaceY[i];
		initMatrixM[10+22*i] = imageX[i] * spaceZ[i];
		initMatrixM[11+22*i] = 0;
		initMatrixM[12+22*i] = 0;
		initMatrixM[13+22*i] = 0;
		initMatrixM[14+22*i] = 0;
		initMatrixM[15+22*i] = spaceX[i];
		initMatrixM[16+22*i] = spaceY[i];
		initMatrixM[17+22*i] = spaceZ[i];
		initMatrixM[18+22*i] = 1;
		initMatrixM[19+22*i] = imageY[i] * spaceX[i];
		initMatrixM[20+22*i] = imageY[i] * spaceY[i];
		initMatrixM[21+22*i] = imageY[i] * spaceZ[i];

		initMatrixW[0+2*i] = -imageX[i];
		initMatrixW[1+2*i] = -imageY[i];
	}
	
	TransposeMatrix(initMatrixM, tInitMatrixM, 12, 11);
	Multiply(tInitMatrixM, initMatrixM, initMatrixMTM, 11, 12, 11);
	Multiply(tInitMatrixM, initMatrixW, initMatrixMTW, 11, 12, 1);
	Inverse (initMatrixMTM, 11);
	Multiply(initMatrixMTM, initMatrixMTW, L, 11, 11, 1);
	
	//L系数的精确值计算
	double rr3; //r3*r3
	double tempA, tempB, tempC;    //用于解算方便设置的中间变量
	double iterFx, iterFy;    //x向主距，y向主距,上次迭代的fx值,上次迭代的fy的值
	double tempA2, rr;  //误差方程中方便计算的变量 r*r  像点向径r

	double *MatrixM  = new double[24 * dgs];
	double *tMatrixM = new double[24 * dgs];
	double *MatrixW  = new double[2  * dgs];
	double  MatrixMTM[144] = {0.0};
	double  MatrixMTW[12]  = {0.0};

	//计算x0，y0
	rr3    =    1 / ( L[8] * L[8] + L[9] * L[9] + L[10] * L[10] );
	x0     = -rr3 * ( L[0] * L[8] + L[1] * L[9] + L[2] * L[10] );
	y0     = -rr3 * ( L[4] * L[8] + L[5] * L[9] + L[6] * L[10] );
	tempA  =  rr3 * ( L[0] * L[0] + L[1] * L[1] + L[2] * L[2] ) - x0 * x0;
	tempB  =  rr3 * ( L[4] * L[4] + L[5] * L[5] + L[6] * L[6] ) - y0 * y0;
	tempC  =  rr3 * ( L[0] * L[4] + L[1] * L[5] + L[2] * L[6] ) - x0 * y0;

	fx = sqrt(( tempA * tempB - tempC * tempC ) / tempB);
	fy = sqrt(( tempA * tempB - tempC * tempC ) / tempA);
	int iternums = 0;
	do 
	{
		iterFx = fx;
		iterFy = fy;
		for (int i = 0; i < dgs; i++)
		{
			tempA2 = L[8] * spaceX[i] + L[9] * spaceY[i] + L[10] * spaceZ[i] + 1;
			rr = (imageX[i] - x0) * (imageX[i] - x0) + (imageY[i] - y0) * (imageY[i] - y0);

			MatrixM[24*i+0]  = -spaceX[i] / tempA2;
			MatrixM[24*i+1]  = -spaceY[i] / tempA2;
			MatrixM[24*i+2]  = -spaceZ[i] / tempA2;
			MatrixM[24*i+3]  = -1 / tempA2;
			MatrixM[24*i+4]  = 0;
			MatrixM[24*i+5]  = 0;
			MatrixM[24*i+6]  = 0;
			MatrixM[24*i+7]  = 0;
			MatrixM[24*i+8]  = -(imageX[i] * spaceX[i]) / tempA2;
			MatrixM[24*i+9]  = -(imageX[i] * spaceY[i]) / tempA2;
			MatrixM[24*i+10] = -(imageX[i] * spaceZ[i]) / tempA2;
			MatrixM[24*i+11] = -(imageX[i] - x0) * rr;

			MatrixM[24*i+12] = 0;
			MatrixM[24*i+13] = 0;
			MatrixM[24*i+14] = 0;
			MatrixM[24*i+15] = 0;
			MatrixM[24*i+16] = -spaceX[i] / tempA2;
			MatrixM[24*i+17] = -spaceY[i] / tempA2;
			MatrixM[24*i+18] = -spaceZ[i] / tempA2;
			MatrixM[24*i+19] = -1 / tempA2;
			MatrixM[24*i+20] = -(imageY[i] * spaceX[i]) / tempA2;
			MatrixM[24*i+21] = -(imageY[i] * spaceY[i]) / tempA2;
			MatrixM[24*i+22] = -(imageY[i] * spaceZ[i]) / tempA2;
			MatrixM[24*i+23] = -(imageY[i] - y0) * rr;

			MatrixW[2*i+0] = imageX[i] / tempA2;
			MatrixW[2*i+1] = imageY[i] / tempA2;
		}
		TransposeMatrix(MatrixM, tMatrixM, 2 * dgs, 12);
		Multiply(tMatrixM, MatrixM, MatrixMTM, 12, 2 * dgs, 12);
		Multiply(tMatrixM, MatrixW, MatrixMTW, 12, 2 * dgs, 1);

		//std::vector<std::vector<double>> mtm(12), rmtm(12);
		//for (int i = 0; i < 12; i++)
		//{
		//	mtm[i].resize(12);
		//	rmtm[i].resize(12);
		//	for (int j = 0; j < 12; j++)
		//	{
		//		mtm[i][j] = MatrixMTM[i*12+j];
		//	}
		//}
		//rmtm = inverse(mtm);
		//for (int i = 0; i < 12; i++)
		//{
		//	for (int j = 0; j < 12; j++)
		//	{
		//		MatrixMTM[i*12+j] = rmtm[i][j];
		//	}
		//}
		Inverse (MatrixMTM, 12);
		Multiply(MatrixMTM, MatrixMTW, L, 12, 12, 1);

		rr3    =    1 / ( L[8] * L[8] + L[9] * L[9] + L[10] * L[10] );
		x0     = -rr3 * ( L[0] * L[8] + L[1] * L[9] + L[2] * L[10] );
		y0     = -rr3 * ( L[4] * L[8] + L[5] * L[9] + L[6] * L[10] );
		tempA  =  rr3 * ( L[0] * L[0] + L[1] * L[1] + L[2] * L[2] ) - x0 * x0;
		tempB  =  rr3 * ( L[4] * L[4] + L[5] * L[5] + L[6] * L[6] ) - y0 * y0;
		tempC  =  rr3 * ( L[0] * L[4] + L[1] * L[5] + L[2] * L[6] ) - x0 * y0;
		fx = sqrt(( tempA * tempB - tempC * tempC ) / tempB);
		fy = sqrt(( tempA * tempB - tempC * tempC ) / tempA);
		
		iternums++;

	} while (iternums < 1000 && fabs(fx - iterFx) > 0.001 &&  fabs(fy - iterFy) > 0.001);
	
	
	//x0 = -rr3 * ( L[0] * L[8] + L[1] * L[9] + L[2] * L[10] );
	//y0 = -rr3 * ( L[4] * L[8] + L[5] * L[9] + L[6] * L[10] );
	

	//根据计算的11个L系数，反算像方坐标，检验精度
	for (int i = 0; i < dgs; i++)
	{
		double temp1, temp2, temp3;
		temp1 = (L[0] * spaceX[i] + L[1] * spaceY[i] + L[2]  * spaceZ[i] + L[3]);
		temp2 = (L[4] * spaceX[i] + L[5] * spaceY[i] + L[6]  * spaceZ[i] + L[7]);
		temp3 = (L[8] * spaceX[i] + L[9] * spaceY[i] + L[10] * spaceZ[i] + 1);
		cImageX[i] = -( temp1 / temp3 );
		cImageY[i] = -( temp2 / temp3 );
	}

	//求解外方位元素
	double cameraXYZ[3] = {0.0};
	double spaceMatrixA[9] = {0.0};
	double spaceMatrixB[3] = {0.0};

	spaceMatrixA[0] =  L[0];		spaceMatrixA[1] =  L[1];		spaceMatrixA[2] =  L[2];
	spaceMatrixA[3] =  L[4];		spaceMatrixA[4] =  L[5];		spaceMatrixA[5] =  L[6];
	spaceMatrixA[6] =  L[8];		spaceMatrixA[7] =  L[9];		spaceMatrixA[8] =  L[10];
	spaceMatrixB[0] = -L[3];		spaceMatrixB[1] = -L[7];		spaceMatrixB[2] = -1;
	Inverse(spaceMatrixA, 3);
	Multiply(spaceMatrixA, spaceMatrixB, cameraXYZ, 3, 3, 1);

	Xs = cameraXYZ[0];	
	Ys = cameraXYZ[1];	
	Zs = cameraXYZ[2];

	double ds = sqrt(tempA / tempB) - 1;
	double dp = asin(sqrt(tempC * tempC / tempA / tempB));

	/*double RL[3*4]= {X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10],1};
	double RL1[3*3] = {fx,-fx*tan(dp),-x0,0,fx/cos(dp)/(1+ds),-y0,0,0,1};
	double RL2[3*4] = {0};
	Inverse(RL1,3);
	Multiply(RL1,RL,RL2,3,3,4);
	double phi=atan(-RL2[8]/RL2[10]);
	double omega=asin(-RL2[9]);
	double kappa=atan(RL2[1]/RL2[5]);*/

	double b1Exceptb2, temb12;

	a3 = sqrt(rr3) * L[8];
	b3 = sqrt(rr3) * L[9];
	c3 = sqrt(rr3) * L[10];
	temb12 = (L[5] + L[9] * y0) * (1 + ds) * cos(dp);
	b1Exceptb2 = (L[1] + L[9] * x0 + temb12) / temb12;

	phi = atan(a3 / c3);
	omega = asin(-b3);
	kappa = atan(b1Exceptb2);

	//a1 = -0.236567813637011;
	//b1 = 0.97160389366698;
	//c1 = 0.0046414827554309;
	//a2 = 0.0190975007370454;
	//b2 = 0.00942594414114093;
	//c2 = -0.99977319280057;
	//a3 = -0.971427277266087;
	//b3 = -0.236425517633382;
	//c3 = -0.0207850810723557;

	//a1 =  -0.18924753619947;
	//b1 =  0.981928146324741;
	//c1 =  0.00157654612647197;
	//a2 =  0.0192653396526686;
	//b2 =  0.00531826011128786;
	//c2 =  -0.999800261450984;
	//a3 =  -0.981740401903934;
	//b3 =  -0.18917936347458;
	//c3 =  -0.0199236468763905;
	/*phi = atan(a3 / c3);
	omega = asin(-b3);
	kappa = atan(b1 / b2);	*/

	a1 =  cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
	a2 = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
	a3 = -sin(phi) * cos(omega);
	b1 =  cos(omega) * sin(kappa);
	b2 =  cos(omega) * cos(kappa);
	b3 = -sin(omega);
	c1 =  sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
	c2 = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
	c3 = cos(phi) * cos(omega);

}

//矩阵转置,转置后矩阵赋给m2，m为m1的列数，n为m2的列数
void DLT::TransposeMatrix(double *m1,double *m2,int m,int n)
{ 
	int i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			m2[j*m+i]=m1[i*n+j];
}

//矩阵求逆,n阶矩阵a，求出之后逆矩阵依然赋给a
void DLT::Inverse(double *a,int n)
{							
	int i,j,k;
	for(k=0;k<n;k++){
		for(i=0;i<n;i++){
			if(i!=k)
				*(a+i*n+k)=-*(a+i*n+k)/(*(a+k*n+k));
		}
		*(a+k*n+k)=1/(*(a+k*n+k));
		for(i=0;i<n;i++){
			if(i!=k){
				for(j=0;j<n;j++){
					if(j!=k)
						*(a+i*n+j)+=*(a+k*n+j)* *(a+i*n+k);
				}
			}
		}
		for(j=0;j<n;j++){
			if(j!=k)
				*(a+k*n+j)*=*(a+k*n+k);
		}
	}
}

//矩阵相乘，其中result指向输出矩阵
void DLT::Multiply(double *m1,double *m2,double *result, int m,int n,int l)
{						
	int i,j,k;
	for(i=0;i<m;i++){
		for(j=0;j<l;j++){
			result[i*l+j]=0.0;							//输出矩阵初始化
			for(k=0;k<n;k++)
				result[i*l+j]+=m1[i*n+k]*m2[j+k*l];		//输出矩阵赋值
		}
	}
}
