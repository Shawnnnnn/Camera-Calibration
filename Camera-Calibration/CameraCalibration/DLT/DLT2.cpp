#include "DLT2.h"
#include "Matrix.h"
#include <math.h>
#include <vector>
using namespace std;

DLT2::DLT2(void)
{
	pixelSize = 0.0;	
	imageHeight = imageWidth = 0; 

	weightSpaceCorrd.x = weightSpaceCorrd.Y = weightSpaceCorrd.Z = 0.0;
	weightPixelCorrd.col = weightPixelCorrd.row = 0;

	camereParam.f = 0.0;
	camereParam.x0 = camereParam.y0 = 0.0;
	camereParam.Xs = camereParam.Ys = camereParam.Zs = 0.0;
	camereParam.phi = camereParam.omega = camereParam.kappa = 0.0;
	camereParam.alpha = camereParam.kappaAlpha = 0.0;	

	for (int i = 0; i < 11; i++)
	{
		dltLs[i] = 0.0;
	}

	for (int i = 0; i < 9; i++)
	{
		rotationMatrix[i] = 0.0;
	}

}

DLT2::~DLT2(void)
{

}

void DLT2::InitImageInfo(const int initImageHeight, const int initImageWidth, const double initPixelSize)
{
	imageHeight = initImageHeight;
	imageWidth = initImageWidth;
	pixelSize = initPixelSize;

	weightPixelCorrd.row = (imageHeight - 1) / 2;
	weightPixelCorrd.col = (imageWidth  - 1) / 2;

}

//3DDLT标定算法：根据控制点的图像像素坐标（行列号）和空间坐标，计算相机内外参数
void DLT2::Calibrate(const int numPts, const vector<PixelCoordPoint>& pixelCoordPoints, const vector<CoordPoint3>& spaceCoordPoints)
{
	if (numPts < 6)  //控制点小于6个，无法计算
	{
		return;
	}

	//将图像像素坐标（行列号）转换为图像物理坐标
	vector<CoordPoint2> imageCoordPoints(numPts);

	for (int i = 0; i < numPts; i++)
	{
		imageCoordPoints[i].x = (pixelCoordPoints[i].col - weightPixelCorrd.col) * pixelSize;
		imageCoordPoints[i].y = -(pixelCoordPoints[i].row - weightPixelCorrd.row) * pixelSize;
	}

	//计算图像重心坐标
	for (int i = 0; i < numPts; i++)
	{
		weightSpaceCorrd.x = weightSpaceCorrd.x + spaceCoordPoints[i].x;
		weightSpaceCorrd.Y = weightSpaceCorrd.Y + spaceCoordPoints[i].Y;
		weightSpaceCorrd.Z = weightSpaceCorrd.Z + spaceCoordPoints[i].Z;
	}
	weightSpaceCorrd.x = weightSpaceCorrd.x / numPts;
	weightSpaceCorrd.Y = weightSpaceCorrd.Y / numPts;
	weightSpaceCorrd.Z = weightSpaceCorrd.Z / numPts;
             
	//重心化处理后的空间坐标
	vector<CoordPoint3> weightedSpaceCoordPoints(numPts);

	for (int i = 0; i < numPts; i++)
	{
		weightedSpaceCoordPoints[i].x = spaceCoordPoints[i].x - weightSpaceCorrd.x;
		weightedSpaceCoordPoints[i].Y = spaceCoordPoints[i].Y - weightSpaceCorrd.Y;
		weightedSpaceCoordPoints[i].Z = spaceCoordPoints[i].Z - weightSpaceCorrd.Z;
	}

	//用前3个和最后3个控制点的坐标来求11个L参数的近似值
	vector<CoordPoint2> temImageCoordPoints;  //临时的图像坐标
	vector<CoordPoint3> temSpaceCoordPoints;  //临时的空间坐标

	for (int i = 0; i < 3; i++)
	{
		temImageCoordPoints.push_back(imageCoordPoints[i]);
		temImageCoordPoints.push_back(imageCoordPoints[numPts - i - 1]);

		temSpaceCoordPoints.push_back(weightedSpaceCoordPoints[i]);
		temSpaceCoordPoints.push_back(weightedSpaceCoordPoints[numPts - i - 1]);
	}

	//11个L参数的近似值数组
	mMatrix initDltL(11); 
	for (int i=0; i<initDltL.size(); i++)
		initDltL[i].resize(1);
	
	//求解参数近似值时的方程的系数值矩阵和常数项矩阵
	mMatrix coeffMatrix(11), constMatrix(11); 
	for (int i=0; i< coeffMatrix.size(); i++)
	{
		coeffMatrix[i].resize(11);
		constMatrix[i].resize(1);
	}

	for (int i = 0; i < 5; i++)
	{
		//系数值矩阵
		coeffMatrix[2 * i][ 0] = temSpaceCoordPoints[i].x; 
		coeffMatrix[2 * i][ 1] = temSpaceCoordPoints[i].Y; 
		coeffMatrix[2 * i][ 2] = temSpaceCoordPoints[i].Z;
		coeffMatrix[2 * i][ 3] = 1; 
		coeffMatrix[2 * i][ 4] = 0; 
		coeffMatrix[2 * i][ 5] = 0; 
		coeffMatrix[2 * i][ 6] = 0; 
		coeffMatrix[2 * i][ 7] = 0;
		coeffMatrix[2 * i][ 8] = temImageCoordPoints[i].x * temSpaceCoordPoints[i].x; 
		coeffMatrix[2 * i][ 9] = temImageCoordPoints[i].x * temSpaceCoordPoints[i].Y; 
		coeffMatrix[2 * i][10] = temImageCoordPoints[i].x * temSpaceCoordPoints[i].Z;

		coeffMatrix[2 * i + 1][ 0] = 0; 
		coeffMatrix[2 * i + 1][ 1] = 0; 
		coeffMatrix[2 * i + 1][ 2] = 0; 
		coeffMatrix[2 * i + 1][ 3] = 0;
		coeffMatrix[2 * i + 1][ 4] = temSpaceCoordPoints[i].x; 
		coeffMatrix[2 * i + 1][ 5] = temSpaceCoordPoints[i].Y; 
		coeffMatrix[2 * i + 1][ 6] = temSpaceCoordPoints[i].Z;
		coeffMatrix[2 * i + 1][ 7] = 1;
		coeffMatrix[2 * i + 1][ 8] = temImageCoordPoints[i].y * temSpaceCoordPoints[i].x; 
		coeffMatrix[2 * i + 1][ 9] = temImageCoordPoints[i].y * temSpaceCoordPoints[i].Y; 
		coeffMatrix[2 * i + 1][10] = temImageCoordPoints[i].y * temSpaceCoordPoints[i].Z;

		//常数项矩阵
		constMatrix[2 * i][ 0]     = -temImageCoordPoints[i].x;
		constMatrix[2 * i + 1][ 0] = -temImageCoordPoints[i].y;
	}
	coeffMatrix[10][ 0] = temSpaceCoordPoints[5].x; 
	coeffMatrix[10][ 1] = temSpaceCoordPoints[5].Y; 
	coeffMatrix[10][ 2] = temSpaceCoordPoints[5].Z;
	coeffMatrix[10][ 3] = 1; 
	coeffMatrix[10][ 4] = 0; 
	coeffMatrix[10][ 5] = 0; 
	coeffMatrix[10][ 6] = 0; 
	coeffMatrix[10][ 7] = 0;
	coeffMatrix[10][ 8] = temImageCoordPoints[5].x * temSpaceCoordPoints[5].x;
	coeffMatrix[10][ 9] = temImageCoordPoints[5].x * temSpaceCoordPoints[5].Y;
	coeffMatrix[10][10] = temImageCoordPoints[5].x * temSpaceCoordPoints[5].Z;

	constMatrix[10][ 0] = -temImageCoordPoints[5].x;

	//系数值矩阵的逆阵
	mMatrix nCoeffMatrix;  	
	nCoeffMatrix = Inverse( coeffMatrix );
	Multiply( nCoeffMatrix, constMatrix, initDltL );    //求解出的11个参数的近似值

	//多于或等于控制时，在近似值的基础上，再进行平差计算，以多次迭代求解的近似值改正数，从而用近似值加上各次迭代的改正数得到最终的平差值
	if (numPts >= 6)
	{
		mMatrix correctMatrix(11);    //11个3维dlt参数的近似值的改正数数组
		for (int i = 0; i < 11; i++)
			correctMatrix[i].resize(1);

		int iterTimes = 0;   //迭代次数，用来控制迭代的执行次数
		double maxCorrection = 100000000.0;    //改正数中的最大值，用来判断是否需要继续迭代


		//重新定义平差的系数值矩阵和常数项矩阵
		coeffMatrix.clear();
		coeffMatrix.resize(2 * numPts);
		for (int i = 0; i < coeffMatrix.size(); i++)
			coeffMatrix[i].resize(11);

		constMatrix.clear();
		constMatrix.resize(2 * numPts);
		for (int i=0; i < constMatrix.size(); i++)
			constMatrix[i].resize(1);		


		//下面迭代求解参数近似值的改正数
		double temA = 0, temB, temC, temx, temy;
		mMatrix tCoeffMatrix, normalMatrixN, nNormalMatrixN, constNormalMatrixW;  //平差时误差方程系数阵的转置矩阵，法方程N矩阵及其逆阵，法方程的常数项矩阵W

		while (iterTimes < 100 && maxCorrection > 0.00001)  //迭代次数小于100次且改正值中的最大值大于10e-5时，迭代平差计算
		{
			for (int i = 0; i < numPts; i++)
			{
				temA = initDltL[8][0] * weightedSpaceCoordPoints[i].x + initDltL[9][0] * weightedSpaceCoordPoints[i].Y + initDltL[10][0] * weightedSpaceCoordPoints[i].Z + 1;
				temB = initDltL[0][0] * weightedSpaceCoordPoints[i].x + initDltL[1][0] * weightedSpaceCoordPoints[i].Y + initDltL[2][0]  * weightedSpaceCoordPoints[i].Z + initDltL[3][0];
				temC = initDltL[4][0] * weightedSpaceCoordPoints[i].x + initDltL[5][0] * weightedSpaceCoordPoints[i].Y + initDltL[6][0]  * weightedSpaceCoordPoints[i].Z + initDltL[7][0];

				temx = temB / temA;
				temy = temC / temA;

				//系数值矩阵
				coeffMatrix[2 * i][ 0] = -weightedSpaceCoordPoints[i].x / temA;
				coeffMatrix[2 * i][ 1] = -weightedSpaceCoordPoints[i].Y / temA;
				coeffMatrix[2 * i][ 2] = -weightedSpaceCoordPoints[i].Z / temA;
				coeffMatrix[2 * i][ 3] = -1 / temA;
				coeffMatrix[2 * i][ 4] = 0; 
				coeffMatrix[2 * i][ 5] = 0; 
				coeffMatrix[2 * i][ 6] = 0; 
				coeffMatrix[2 * i][ 7] = 0;
				coeffMatrix[2 * i][ 8] = (temx * weightedSpaceCoordPoints[i].x) / temA;
				coeffMatrix[2 * i][ 9] = (temx * weightedSpaceCoordPoints[i].Y) / temA;
				coeffMatrix[2 * i][10] = (temx * weightedSpaceCoordPoints[i].Z) / temA;

				coeffMatrix[2 * i + 1][ 0] = 0; 
				coeffMatrix[2 * i + 1][ 1] = 0; 
				coeffMatrix[2 * i + 1][ 2] = 0; 
				coeffMatrix[2 * i + 1][ 3] = 0;
				coeffMatrix[2 * i + 1][ 4] = -weightedSpaceCoordPoints[i].x / temA;
				coeffMatrix[2 * i + 1][ 5] = -weightedSpaceCoordPoints[i].Y / temA;
				coeffMatrix[2 * i + 1][ 6] = -weightedSpaceCoordPoints[i].Z / temA;
				coeffMatrix[2 * i + 1][ 7] = -1 / temA;
				coeffMatrix[2 * i + 1][ 8] = (temy * weightedSpaceCoordPoints[i].x) / temA;
				coeffMatrix[2 * i + 1][ 9] = (temy * weightedSpaceCoordPoints[i].Y) / temA;
				coeffMatrix[2 * i + 1][10] = (temy * weightedSpaceCoordPoints[i].Z) / temA;

				//常数项矩阵
				constMatrix[2 * i][0] = imageCoordPoints[i].x + temx;
				constMatrix[2*i+1][0] = imageCoordPoints[i].y + temy;
			}

			Transpose( coeffMatrix, tCoeffMatrix );                         //求系数值矩阵的转置矩阵
			Multiply( tCoeffMatrix, coeffMatrix, normalMatrixN );           //求法方程的N矩阵
			nNormalMatrixN = Inverse( normalMatrixN );                      //求法方程N矩阵的逆阵
			Multiply( tCoeffMatrix, constMatrix, constNormalMatrixW );      //求法方程的常数项阵W
			Multiply( nNormalMatrixN, constNormalMatrixW, correctMatrix );  //求解近似值的改正数

			for (int j = 0; j < 11; j++)                                  //利用求解的近似值改正数更新参数的值
				initDltL[j][0] = initDltL[j][0] + correctMatrix[j][0];

			iterTimes = iterTimes + 1;    //迭代次数加1

			maxCorrection = abs(correctMatrix[0][0]);   //改正数绝对值中的最大值
			for (int j = 1; j < 11; j++)
				if (maxCorrection < abs(correctMatrix[j][0])) 
					maxCorrection = abs(correctMatrix[j][0]);
		}

		for (int i = 0; i < 11; i++)
			dltLs[i] = initDltL[i][0];//平差结束后返回的dlt参数数组
	}
	
	//利用计算的11个L系数，来求解相机的相关参数
	double temfm, temfz1, temfz2;
	temfm  = dltLs[8] * dltLs[8] + dltLs[9] * dltLs[9] + dltLs[10] * dltLs[10];
	temfz1 = dltLs[0] * dltLs[8] + dltLs[1] * dltLs[9] + dltLs[2] * dltLs[10];
	temfz2 = dltLs[4] * dltLs[8] + dltLs[5] * dltLs[9] + dltLs[6] * dltLs[10];

	camereParam.x0 = -temfz1 / temfm;   //像主点的x坐标x0
	camereParam.y0 = -temfz2 / temfm;   //像主点的y坐标y0

	double rr3 = 1.0 / temfm;     //中间值
	double temA2, temB2, temC2;   //3个中间值

	temA2 = rr3 * (dltLs[0] * dltLs[0] + dltLs[1] * dltLs[1] + dltLs[2] * dltLs[2]) - camereParam.x0 * camereParam.x0;
	temB2 = rr3 * (dltLs[4] * dltLs[4] + dltLs[5] * dltLs[5] + dltLs[6] * dltLs[6]) - camereParam.y0 * camereParam.y0;
	temC2 = rr3 * (dltLs[0] * dltLs[4] + dltLs[1] * dltLs[5] + dltLs[2] * dltLs[6]) - camereParam.x0 * camereParam.y0;

	double fx, fy;
	fx = sqrt((temA2 * temB2 - temC2 * temC2) / temB2);
	fy = sqrt((temA2 * temB2 - temC2 * temC2) / temA2);

	camereParam.f = sqrt(fx * fx + fy * fy);

	//求解相机中心的三维坐标X、Y、Z
	mMatrix coeffMatrix2(3), constMatrix2(3);   //求解XYZ的系数矩阵和常数矩阵
	for (int i = 0; i < 3; i++)
	{
		coeffMatrix2[i].resize(3);
		constMatrix2[i].resize(1);
	}
	coeffMatrix2[0][0] = dltLs[0];			coeffMatrix2[0][1] = dltLs[1];			coeffMatrix2[0][2] = dltLs[2];
	coeffMatrix2[1][0] = dltLs[4];			coeffMatrix2[1][1] = dltLs[5];			coeffMatrix2[1][2] = dltLs[6];
	coeffMatrix2[2][0] = dltLs[8];			coeffMatrix2[2][1] = dltLs[9];			coeffMatrix2[2][2] = dltLs[10];

	constMatrix2[0][0] = -dltLs[3];			constMatrix2[1][0] = -dltLs[7];			constMatrix2[2][0] = -1;

	mMatrix nCoeffMatrix2, temCameraCoord;
	nCoeffMatrix2 = Inverse( coeffMatrix2 );
	Multiply( nCoeffMatrix2, constMatrix2, temCameraCoord );

	//因为求3DDLT参数时，物方坐标进行了重心化处理，重心化后的物方坐标相当于是以重心为原点的，所以，用dlt系数计算的相机中心坐标要加上物方重心的坐标，才能和物方坐标的坐标系一致
	camereParam.Xs = temCameraCoord[0][0] + weightSpaceCorrd.x;//相机中心的X坐标
	camereParam.Ys = temCameraCoord[1][0] + weightSpaceCorrd.Y;//相机中心的Y坐标
	camereParam.Zs = temCameraCoord[2][0] + weightSpaceCorrd.Z;//相机中心的Z坐标


	//求解相机的外方位角元素phi, omega, kappa
	//求解相机中心的姿态角：相机方位角fhfwjA，相机倾角fhxjqjaf
	//rotationMatrix旋转矩阵，顺次为a1, b1, c1, a2, b2, c2, a3, b3, c3

	double a1, b1, c1, a2, b2, c2, a3, b3, c3;
		
	a3 = dltLs[8 ] / sqrt(temfm);	
	b3 = dltLs[9 ] / sqrt(temfm);
	c3 = dltLs[10] / sqrt(temfm);

	double dp = asin(sqrt(temC2 * temC2 / temA2 / temB2));
	double ds = sqrt(temA2 / temB2) - 1;

	double temb12 = (dltLs[5] + dltLs[9] * camereParam.x0) * (1 + ds) * cos(dp);
	double b1exceptb2 = (dltLs[1] + dltLs[9] * camereParam.x0 + temb12) / temb12;

	camereParam.phi = atan(a3 / c3);
	camereParam.omega = asin(-b3);
	camereParam.kappa = atan(b1exceptb2);

	camereParam.alpha = acos(c3);             //相机倾角（弧度）
	camereParam.kappaAlpha = atan(a3 / b3);   // 相机方位角（弧度）
	

	//弧度表示的相机姿态角换算成角度
	camereParam.phi *= (180 / PI);
	camereParam.omega *= (180 / PI);
	camereParam.kappa *= (180 / PI);

	camereParam.alpha *= (180 / PI);
	camereParam.kappaAlpha *= (180 / PI);

	//计算旋转矩阵
	rotationMatrix[0] = a1 =  cos(camereParam.phi) * cos(camereParam.kappa) - sin(camereParam.phi) * sin(camereParam.omega) * sin(camereParam.kappa);
	rotationMatrix[3] = a2 = -cos(camereParam.phi) * sin(camereParam.kappa) - sin(camereParam.phi) * sin(camereParam.omega) * cos(camereParam.kappa);
	rotationMatrix[6] = a3 = -sin(camereParam.phi) * cos(camereParam.omega);
	rotationMatrix[1] = b1 =  cos(camereParam.omega) * sin(camereParam.kappa);
	rotationMatrix[4] = b2 =  cos(camereParam.omega) * cos(camereParam.kappa);
	rotationMatrix[7] = b3 = -sin(camereParam.omega);
	rotationMatrix[2] = c1 =  sin(camereParam.phi) * cos(camereParam.kappa) + cos(camereParam.phi) * sin(camereParam.omega) * sin(camereParam.kappa);
	rotationMatrix[5] = c2 = -sin(camereParam.phi) * sin(camereParam.kappa) + cos(camereParam.phi) * sin(camereParam.omega) * cos(camereParam.kappa);
	rotationMatrix[8] = c3 =  cos(camereParam.phi) * cos(camereParam.omega);

}

//根据控制点空间坐标和三维DLT解算的参数，反算图像像素坐标（行列号）
vector<PixelCoordPoint> DLT2::CalculatePixcelCorrd(const vector<CoordPoint3>& spaceCoordPoints)
{
	int numPts = spaceCoordPoints.size();
	vector<CoordPoint3> weightedSpaceCoordPoints(numPts);
	vector<CoordPoint2> imageCoordPoints(numPts);
	vector<PixelCoordPoint>  imagePixelCoordPoints(numPts);

	for(int i = 0; i < numPts; i++)
	{
		weightedSpaceCoordPoints[i].x = spaceCoordPoints[i].x - weightSpaceCorrd.x;
		weightedSpaceCoordPoints[i].Y = spaceCoordPoints[i].Y - weightSpaceCorrd.Y;
		weightedSpaceCoordPoints[i].Z = spaceCoordPoints[i].Z - weightSpaceCorrd.Z;

		double temfz1, temfz2, temfm;

		temfm  = dltLs[8] * weightedSpaceCoordPoints[i].x + dltLs[9] * weightedSpaceCoordPoints[i].Y + dltLs[10] * weightedSpaceCoordPoints[i].Z + 1;
		temfz1 = dltLs[0] * weightedSpaceCoordPoints[i].x + dltLs[1] * weightedSpaceCoordPoints[i].Y + dltLs[ 2] * weightedSpaceCoordPoints[i].Z + dltLs[3];
		temfz2 = dltLs[4] * weightedSpaceCoordPoints[i].x + dltLs[5] * weightedSpaceCoordPoints[i].Y + dltLs[ 6] * weightedSpaceCoordPoints[i].Z + dltLs[7];

		imageCoordPoints[i].x = -temfz1 / temfm;
		imageCoordPoints[i].y = -temfz2 / temfm;

		imagePixelCoordPoints[i].row = weightPixelCorrd.row - imageCoordPoints[i].y / pixelSize;
		imagePixelCoordPoints[i].col = weightPixelCorrd.col + imageCoordPoints[i].x / pixelSize;

	}
	return imagePixelCoordPoints;

}