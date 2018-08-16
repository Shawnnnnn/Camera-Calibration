#include "stdafx.h"
#include "DLT3.h"
#include <math.h>
#include <vector>
using namespace std;

DLT3::DLT3(void)
{
	
}

DLT3::~DLT3(void)
{

}

//根据输入的6个控制点的物方xyz坐标、及其像方的像素行列号，进行3维DLT（直接线性变换），求解相机的内外方位元素
//需要输入：控制点的个数dgs，控制点在影像上的行号数组xfhh、列号数组xflh，控制点在物方的X坐标数组wfxzb、物方Y坐标数组wfyzb、物方Z坐标数组wfzzb；
//函数返回：相机的焦距(fhf) fhfx,fhfy，像主点坐标x0,y0;相机摄影中心的XYZ坐标数组fhXYZ（0，1，2三个元素分别是X、Y、Z），相机的外方位角元素phi, omega, kappa；
void DLT3::DLT_cameraPara_6points(int dgs, const std::vector<double>& imageX,const std::vector<double>& imageY, const std::vector<double> &spaceX, const std::vector<double>& spaceY, const std::vector<double>& spaceZ, double &fx, double &fy, double &x0, double &y0, std::vector<double>& cameraXYZ, double &phi, double& omega, double& kappa,std::vector<double>& calculateImageX, std::vector<double>& calculateImageY)
{
	//返回的变量赋初始值
	fx = 0.0;
	fy = 0.0;
	x0 = 0.0;
	y0 = 0.0;
	phi = 0.0; 
	omega = 0.0;
	kappa = 0.0;

	//下面先求解3DDLT的11个系数
	double dltL[11];//3DDLT的11个参数           


	//下面先用前3个和最后3个公共点的两套坐标来求3维DLT的11个参数的近似值
	double temImageX[6], temImageY[6];//临时的像平面坐标x,y数组
	double temSpaceX[6], temSpaceY[6], temSpaceZ[6];//临时的物方三维坐标x,y,z数组

	temImageX[0] = imageX[0]; 
	temImageY[0] = imageY[0];
	temImageX[1] = imageX[1]; 
	temImageY[1] = imageY[1];
	temImageX[2] = imageX[2]; 
	temImageY[2] = imageY[2];
	temImageX[3] = imageX[dgs - 3]; 
	temImageY[3] = imageY[dgs - 3];
	temImageX[4] = imageX[dgs - 2];
	temImageY[4] = imageY[dgs - 2];
	temImageX[5] = imageX[dgs - 1]; 
	temImageY[5] = imageY[dgs - 1];

	temSpaceX[0] = spaceX[0]; 
	temSpaceY[0] = spaceY[0]; 
	temSpaceZ[0] = spaceZ[0];
	temSpaceX[1] = spaceX[1]; 
	temSpaceY[1] = spaceY[1]; 
	temSpaceZ[1] = spaceZ[1];
	temSpaceX[2] = spaceX[2]; 
	temSpaceY[2] = spaceY[2]; 
	temSpaceZ[1] = spaceZ[2];
	temSpaceX[3] = spaceX[dgs - 3]; 
	temSpaceY[3] = spaceY[dgs - 3]; 
	temSpaceZ[3] = spaceZ[dgs - 3];
	temSpaceX[4] = spaceX[dgs - 2];
	temSpaceY[4] = spaceY[dgs - 2];
	temSpaceZ[4] = spaceZ[dgs - 2];
	temSpaceX[5] = spaceX[dgs - 1]; 
	temSpaceY[5] = spaceY[dgs - 1]; 
	temSpaceZ[5] = spaceZ[dgs - 1];

	//11个3维dlt参数的近似值数组
	std::vector<std::vector<double>> initDltL(11); 
	for (int i=0; i<initDltL.size(); i++)
		initDltL[i].resize(1);
	
	//求解参数近似值时的方程的系数值矩阵和常数项矩阵
	std::vector<std::vector<double>> coeffMatrix(11); 
	for (int i=0; i< coeffMatrix.size(); i++)
		coeffMatrix[i].resize(11);
	std::vector<std::vector<double>> constMatrix(11); 
	for (int i=0; i<constMatrix.size(); i++)
		constMatrix[i].resize(1);

	for (int i = 0; i < 5; i++)
	{
		//系数值矩阵
		coeffMatrix[2 * i][ 0] = temSpaceX[i]; 
		coeffMatrix[2 * i][ 1] = temSpaceY[i]; 
		coeffMatrix[2 * i][ 2] = temSpaceZ[i];
		coeffMatrix[2 * i][ 3] = 1; 
		coeffMatrix[2 * i][ 4] = 0; 
		coeffMatrix[2 * i][ 5] = 0; 
		coeffMatrix[2 * i][ 6] = 0; 
		coeffMatrix[2 * i][ 7] = 0;
		coeffMatrix[2 * i][ 8] = temImageX[i] * temSpaceX[i]; 
		coeffMatrix[2 * i][ 9] = temImageX[i] * temSpaceY[i]; 
		coeffMatrix[2 * i][ 10] = temImageX[i] * temSpaceZ[i];

		coeffMatrix[2 * i + 1][ 0] = 0; 
		coeffMatrix[2 * i + 1][ 1] = 0; 
		coeffMatrix[2 * i + 1][ 2] = 0; 
		coeffMatrix[2 * i + 1][ 3] = 0;
		coeffMatrix[2 * i + 1][ 4] = temSpaceX[i]; 
		coeffMatrix[2 * i + 1][ 5] = temSpaceY[i]; 
		coeffMatrix[2 * i + 1][ 6] = temSpaceZ[i];
		coeffMatrix[2 * i + 1][ 7] = 1;
		coeffMatrix[2 * i + 1][ 8] = temImageY[i] * temSpaceX[i]; 
		coeffMatrix[2 * i + 1][ 9] = temImageY[i] * temSpaceY[i]; 
		coeffMatrix[2 * i + 1][ 10] = temImageY[i] * temSpaceZ[i];

		//常数项矩阵
		constMatrix[2 * i][ 0] = -temImageX[i];
		constMatrix[2 * i + 1][ 0] = -temImageY[i];
	}
	coeffMatrix[10][ 0] = temSpaceX[5]; 
	coeffMatrix[10][ 1] = temSpaceY[5]; 
	coeffMatrix[10][ 2] = temSpaceZ[5];
	coeffMatrix[10][ 3] = 1; 
	coeffMatrix[10][ 4] = 0; 
	coeffMatrix[10][ 5] = 0; 
	coeffMatrix[10][ 6] = 0; 
	coeffMatrix[10][ 7] = 0;
	coeffMatrix[10][ 8] = temImageX[5] * temSpaceX[5]; 
	coeffMatrix[10][ 9] = temImageX[5] * temSpaceY[5]; 
	coeffMatrix[10][ 10] = temImageX[5] * temSpaceZ[5];

	constMatrix[10][ 0] = -temImageX[5];

	std::vector<std::vector<double>> nCoeffMatrix;//系数值矩阵的逆阵
	nCoeffMatrix=inverse(coeffMatrix);
	multiply(nCoeffMatrix,constMatrix,initDltL);//求解出的11个参数的近似值

	//多于或等于6个公共点时，在近似值的基础上，再进行平差计算，以多次迭代求解的近似值改正数，从而用近似值加上各次迭代的改正数得到最终的平差值
	if (dgs >= 6)
	{
		std::vector<std::vector<double>> correctMatrix(11);//11个3维dlt参数的近似值的改正数数组
		for (int i=0;i<11;i++)
			correctMatrix[i].resize(1);
		int iterTimes = 0;//迭代次数，用来控制迭代的执行次数
		double maxCorrection = 100000000.0;//改正数中的最大值，用来判断是否需要继续迭代


		//重新定义平差的系数值矩阵和常数项矩阵
		coeffMatrix.clear();
		coeffMatrix.resize(2*dgs);
		for (int i=0;i< coeffMatrix.size();i++)
			coeffMatrix[i].resize(11);

		constMatrix.clear();
		constMatrix.resize(2*dgs);
		for (int i=0;i< constMatrix.size();i++)
			constMatrix[i].resize(1);		


		//下面迭代求解参数近似值的改正数
		double temA = 0, temB, temC, temx, temy;
		std::vector<std::vector<double>> tCoeffMatrix, normalMatrixN, nNormalMatrixN, constNormalMatrixW;//平差时误差方程系数阵的转置矩阵，法方程N矩阵及其逆阵，法方程的常数项矩阵W
		while (iterTimes < 100 && maxCorrection > 0.00001)//迭代次数小于100次且改正值中的最大值大于10e-5时，迭代平差计算
		{
			for (int i = 0; i < dgs; i++)
			{
				temA = initDltL[8][0] * spaceX[i] + initDltL[9][0] * spaceY[i] + initDltL[10][0] * spaceZ[i] + 1;
				temB = initDltL[0][0] * spaceX[i] + initDltL[1][0] * spaceY[i] + initDltL[2][0] * spaceZ[i] + initDltL[3][0];
				temC = initDltL[4][0] * spaceX[i] + initDltL[5][0] * spaceY[i] + initDltL[6][0] * spaceZ[i] + initDltL[7][0];

				temx = temB / temA;
				temy = temC / temA;

				//系数值矩阵
				coeffMatrix[2 * i][0] = -spaceX[i] / temA;
				coeffMatrix[2 * i][1] = -spaceY[i] / temA;
				coeffMatrix[2 * i][2] = -spaceZ[i] / temA;
				coeffMatrix[2 * i][3] = -1 / temA;
				coeffMatrix[2 * i][4] = 0; 
				coeffMatrix[2 * i][5] = 0; 
				coeffMatrix[2 * i][6] = 0; 
				coeffMatrix[2 * i][7] = 0;
				coeffMatrix[2 * i][8] = (temx * spaceX[i]) / temA;
				coeffMatrix[2 * i][9] = (temx * spaceY[i]) / temA;
				coeffMatrix[2 * i][10] = (temx * spaceZ[i]) / temA;

				coeffMatrix[2 * i + 1][0] = 0; 
				coeffMatrix[2 * i + 1][1] = 0; 
				coeffMatrix[2 * i + 1][2] = 0; 
				coeffMatrix[2 * i + 1][3] = 0;
				coeffMatrix[2 * i + 1][4] = -spaceX[i] / temA;
				coeffMatrix[2 * i + 1][5] = -spaceY[i] / temA;
				coeffMatrix[2 * i + 1][6] = -spaceZ[i] / temA;
				coeffMatrix[2 * i + 1][7] = -1 / temA;
				coeffMatrix[2 * i + 1][8] = (temy * spaceX[i]) / temA;
				coeffMatrix[2 * i + 1][9] = (temy * spaceY[i]) / temA;
				coeffMatrix[2 * i + 1][10] = (temy * spaceZ[i]) / temA;

				//常数项矩阵
				constMatrix[2 * i][0] = imageX[i] + temx;
				constMatrix[2 * i + 1][0] = imageY[i] + temy;
			}

			Transpose(coeffMatrix,tCoeffMatrix);//求系数值矩阵的转置矩阵
			multiply(tCoeffMatrix,coeffMatrix,normalMatrixN);//求法方程的N矩阵
			nNormalMatrixN = inverse(normalMatrixN);//求法方程N矩阵的逆阵
			multiply(tCoeffMatrix,constMatrix,constNormalMatrixW);//求法方程的常数项阵W
			multiply(nNormalMatrixN,constNormalMatrixW,correctMatrix);//求解近似值的改正数

			for (int j = 0; j < 11; j++)//利用求解的近似值改正数更新参数的值
				initDltL[j][0] = initDltL[j][0] + correctMatrix[j][0];

			iterTimes = iterTimes + 1;//迭代次数加1

			maxCorrection = abs(correctMatrix[0][0]);
			for (int j = 1; j < 11; j++)//改正数绝对值中的最大值
				if (maxCorrection < abs(correctMatrix[j][0])) 
					maxCorrection = abs(correctMatrix[j][0]);
		}

		for (int i = 0; i < 11; i++)
			dltL[i] = initDltL[i][0];//平差结束后返回的dlt参数数组
	}


	//根据计算的11个L系数，反算像方坐标，检验精度
	calculateImageX.resize(dgs);
	calculateImageY.resize(dgs);

	for (int i = 0; i < dgs; i++)
	{
		double temp1 = (dltL[0] * spaceX[i] + dltL[1] * spaceY[i] + dltL[2]  * spaceZ[i] + dltL[3]);
		double temp2 = (dltL[4] * spaceX[i] + dltL[5] * spaceY[i] + dltL[6]  * spaceZ[i] + dltL[7]);
		double temp3 = (dltL[8] * spaceX[i] + dltL[9] * spaceY[i] + dltL[10] * spaceZ[i] + 1 );
		calculateImageX[i] = -( temp1 / temp3 );
		calculateImageY[i] = -( temp2 / temp3 );
	}

	//下面再利用计算的DLT的11个系数，来求解相机的相关参数
	double temfm,temfz1,temfz2;
	temfm  = dltL[8] * dltL[8] + dltL[9] * dltL[9] + dltL[10] * dltL[10];
	temfz1 = dltL[0] * dltL[8] + dltL[1] * dltL[9] + dltL[2] * dltL[10];
	temfz2 = dltL[4] * dltL[8] + dltL[5] * dltL[9] + dltL[6] * dltL[10];

	x0 = -temfz1 / temfm;//像主点的x坐标x0
	y0 = -temfz2 / temfm;//像主点的y坐标y0

	double rr3 = 1.0 / temfm;//中间值
	double temA2,temB2,temC2;//3个中间值

	temA2 = rr3 * (dltL[0] * dltL[0] + dltL[1] * dltL[1] + dltL[2] * dltL[2]) - x0 * x0;
	temB2 = rr3 * (dltL[4] * dltL[4] + dltL[5] * dltL[5] + dltL[6] * dltL[6]) - y0 * y0;
	temC2 = rr3 * (dltL[0] * dltL[4] + dltL[1] * dltL[5] + dltL[2] * dltL[6]) - x0 * y0;

	fx = sqrt((temA2 * temB2 - temC2 * temC2) / temB2);
	fy = sqrt((temA2 * temB2 - temC2 * temC2) / temA2);

	//fhf = fhfx;// 直接以x方向的主距fx作为相机的焦距f ////Math.Sqrt(fhfx * fhfx + fhfy * fhfy);

	//求解相机中心的三维坐标X、Y、Z
	std::vector<std::vector<double>> coeffMatrix2(3),constMatrix2(3);//求解XYZ的系数矩阵和常数矩阵
	for (int i=0;i<3;i++)
	{
		coeffMatrix2[i].resize(3);
		constMatrix2[i].resize(1);
	}
	coeffMatrix2[0][0] = dltL[0];			coeffMatrix2[0][1] = dltL[1];			coeffMatrix2[0][2] = dltL[2];
	coeffMatrix2[1][0] = dltL[4];			coeffMatrix2[1][1] = dltL[5];			coeffMatrix2[1][2] = dltL[6];
	coeffMatrix2[2][0] = dltL[8];			coeffMatrix2[2][1] = dltL[9];			coeffMatrix2[2][2] = dltL[10];

	constMatrix2[0][0] = -dltL[3];			constMatrix2[1][0] = -dltL[7];			constMatrix2[2][0] = -1;


	std::vector<std::vector<double>> nCoeffMatrix2, temCameraXYZ;
	nCoeffMatrix2 = inverse(coeffMatrix2);
	multiply(nCoeffMatrix2, constMatrix2, temCameraXYZ);

	cameraXYZ.resize(3);
	cameraXYZ[0] = temCameraXYZ[0][0];//相机中心的X坐标
	cameraXYZ[1] = temCameraXYZ[1][0];//相机中心的Y坐标
	cameraXYZ[2] = temCameraXYZ[2][0];//相机中心的Z坐标

	c3 = dltL[10] / sqrt(temfm);
	b3 = dltL[9] / sqrt(temfm);
	a3 = dltL[8] / sqrt(temfm);

	double dp = asin(sqrt(temC2 * temC2 / temA2 / temB2));
	double ds = sqrt(temA2 / temB2) - 1;

	double tempb12 = (dltL[5] + dltL[9] * y0) * (1 + ds) * cos(dp);
	double b1Exceptb2 = (dltL[1] + dltL[9] * x0 + tempb12) / tempb12;

	phi = atan(a3 / c3);
	omega = asin(-b3);
	kappa = atan(b1Exceptb2);

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

// 方阵求逆函数1,Psr为输入方阵，N为方阵的大小,P为返回的逆矩阵
std::vector<std::vector<double>>DLT3::inverse( const std::vector<std::vector<double>>& Psr )
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


/// 矩阵相乘函数
bool DLT3::multiply( const std::vector<std::vector<double>>& jz1, const std::vector<std::vector<double>>& jz2, std::vector<std::vector<double>>& scjz)
{
	scjz.clear();
	int m1, n1, m2, n2;
	double tem;
	m1 = jz1.size(); n1 = jz1[0].size();
	m2 = jz2.size(); n2 = jz2[0].size();
	scjz.resize(m1);
	for (int i = 0; i < m1; i++)
	{
		scjz[i].resize(n2);
		for (int j = 0; j < n2; j++)
		{
			tem = 0;
			for (int k = 0; k < n1; k++)
			{
				tem = tem + jz1[i][k] * jz2[k][j];
			}
			scjz[i][j] = tem;
		}
	}
	return false;
}

bool DLT3::juzhxj( const std::vector<std::vector<double>>& jz1, const std::vector<std::vector<double>>& jz2, std::vector<std::vector<double>>& scjz,int jjbz )
{
	scjz.clear();
	int m1, n1, m2, n2;

	m1 = jz1.size(); n1 = jz1[0].size();
	m2 = jz2.size(); n2 = jz2[0].size();
	scjz.resize(m1);

	if (m1 == m2 && n1 == n2)//两个矩阵的维数相同时，才进行相加或相减的计算，否则，输出矩阵的值全为０
	{
		for (int i = 0; i < m1; i++)
		{
			scjz[i].resize(n2);
			for (int j = 0; j < n2; j++)
			{
				scjz[i][j] = jz1[i][j] + jjbz * jz2[i][j];//当jjbz=1时，为矩阵相加；当jjbz=-1时，为矩阵相减
			}
		}
		return true;

	}
	return false;
}

bool DLT3::Transpose( const std::vector<std::vector<double>>& jz, std::vector<std::vector<double>>& scjz )
{
	scjz.clear();
	int m, n;//jz的行列大小
	m = jz.size(); n = jz[0].size();
	scjz.resize(n);
	for (int i = 0; i < n; i++)
	{
		scjz[i].resize(m);
		for (int j = 0; j < m; j++)
		{
			scjz[i][j] = jz[j][i];//转置
		}
	}
	return false;
}