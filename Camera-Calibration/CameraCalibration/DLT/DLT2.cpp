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

//3DDLT�궨�㷨�����ݿ��Ƶ��ͼ���������꣨���кţ��Ϳռ����꣬��������������
void DLT2::Calibrate(const int numPts, const vector<PixelCoordPoint>& pixelCoordPoints, const vector<CoordPoint3>& spaceCoordPoints)
{
	if (numPts < 6)  //���Ƶ�С��6�����޷�����
	{
		return;
	}

	//��ͼ���������꣨���кţ�ת��Ϊͼ����������
	vector<CoordPoint2> imageCoordPoints(numPts);

	for (int i = 0; i < numPts; i++)
	{
		imageCoordPoints[i].x = (pixelCoordPoints[i].col - weightPixelCorrd.col) * pixelSize;
		imageCoordPoints[i].y = -(pixelCoordPoints[i].row - weightPixelCorrd.row) * pixelSize;
	}

	//����ͼ����������
	for (int i = 0; i < numPts; i++)
	{
		weightSpaceCorrd.x = weightSpaceCorrd.x + spaceCoordPoints[i].x;
		weightSpaceCorrd.Y = weightSpaceCorrd.Y + spaceCoordPoints[i].Y;
		weightSpaceCorrd.Z = weightSpaceCorrd.Z + spaceCoordPoints[i].Z;
	}
	weightSpaceCorrd.x = weightSpaceCorrd.x / numPts;
	weightSpaceCorrd.Y = weightSpaceCorrd.Y / numPts;
	weightSpaceCorrd.Z = weightSpaceCorrd.Z / numPts;
             
	//���Ļ������Ŀռ�����
	vector<CoordPoint3> weightedSpaceCoordPoints(numPts);

	for (int i = 0; i < numPts; i++)
	{
		weightedSpaceCoordPoints[i].x = spaceCoordPoints[i].x - weightSpaceCorrd.x;
		weightedSpaceCoordPoints[i].Y = spaceCoordPoints[i].Y - weightSpaceCorrd.Y;
		weightedSpaceCoordPoints[i].Z = spaceCoordPoints[i].Z - weightSpaceCorrd.Z;
	}

	//��ǰ3�������3�����Ƶ����������11��L�����Ľ���ֵ
	vector<CoordPoint2> temImageCoordPoints;  //��ʱ��ͼ������
	vector<CoordPoint3> temSpaceCoordPoints;  //��ʱ�Ŀռ�����

	for (int i = 0; i < 3; i++)
	{
		temImageCoordPoints.push_back(imageCoordPoints[i]);
		temImageCoordPoints.push_back(imageCoordPoints[numPts - i - 1]);

		temSpaceCoordPoints.push_back(weightedSpaceCoordPoints[i]);
		temSpaceCoordPoints.push_back(weightedSpaceCoordPoints[numPts - i - 1]);
	}

	//11��L�����Ľ���ֵ����
	mMatrix initDltL(11); 
	for (int i=0; i<initDltL.size(); i++)
		initDltL[i].resize(1);
	
	//����������ֵʱ�ķ��̵�ϵ��ֵ����ͳ��������
	mMatrix coeffMatrix(11), constMatrix(11); 
	for (int i=0; i< coeffMatrix.size(); i++)
	{
		coeffMatrix[i].resize(11);
		constMatrix[i].resize(1);
	}

	for (int i = 0; i < 5; i++)
	{
		//ϵ��ֵ����
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

		//���������
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

	//ϵ��ֵ���������
	mMatrix nCoeffMatrix;  	
	nCoeffMatrix = Inverse( coeffMatrix );
	Multiply( nCoeffMatrix, constMatrix, initDltL );    //������11�������Ľ���ֵ

	//���ڻ���ڿ���ʱ���ڽ���ֵ�Ļ����ϣ��ٽ���ƽ����㣬�Զ�ε������Ľ���ֵ���������Ӷ��ý���ֵ���ϸ��ε����ĸ������õ����յ�ƽ��ֵ
	if (numPts >= 6)
	{
		mMatrix correctMatrix(11);    //11��3άdlt�����Ľ���ֵ�ĸ���������
		for (int i = 0; i < 11; i++)
			correctMatrix[i].resize(1);

		int iterTimes = 0;   //�����������������Ƶ�����ִ�д���
		double maxCorrection = 100000000.0;    //�������е����ֵ�������ж��Ƿ���Ҫ��������


		//���¶���ƽ���ϵ��ֵ����ͳ��������
		coeffMatrix.clear();
		coeffMatrix.resize(2 * numPts);
		for (int i = 0; i < coeffMatrix.size(); i++)
			coeffMatrix[i].resize(11);

		constMatrix.clear();
		constMatrix.resize(2 * numPts);
		for (int i=0; i < constMatrix.size(); i++)
			constMatrix[i].resize(1);		


		//�����������������ֵ�ĸ�����
		double temA = 0, temB, temC, temx, temy;
		mMatrix tCoeffMatrix, normalMatrixN, nNormalMatrixN, constNormalMatrixW;  //ƽ��ʱ����ϵ�����ת�þ��󣬷�����N���������󣬷����̵ĳ��������W

		while (iterTimes < 100 && maxCorrection > 0.00001)  //��������С��100���Ҹ���ֵ�е����ֵ����10e-5ʱ������ƽ�����
		{
			for (int i = 0; i < numPts; i++)
			{
				temA = initDltL[8][0] * weightedSpaceCoordPoints[i].x + initDltL[9][0] * weightedSpaceCoordPoints[i].Y + initDltL[10][0] * weightedSpaceCoordPoints[i].Z + 1;
				temB = initDltL[0][0] * weightedSpaceCoordPoints[i].x + initDltL[1][0] * weightedSpaceCoordPoints[i].Y + initDltL[2][0]  * weightedSpaceCoordPoints[i].Z + initDltL[3][0];
				temC = initDltL[4][0] * weightedSpaceCoordPoints[i].x + initDltL[5][0] * weightedSpaceCoordPoints[i].Y + initDltL[6][0]  * weightedSpaceCoordPoints[i].Z + initDltL[7][0];

				temx = temB / temA;
				temy = temC / temA;

				//ϵ��ֵ����
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

				//���������
				constMatrix[2 * i][0] = imageCoordPoints[i].x + temx;
				constMatrix[2*i+1][0] = imageCoordPoints[i].y + temy;
			}

			Transpose( coeffMatrix, tCoeffMatrix );                         //��ϵ��ֵ�����ת�þ���
			Multiply( tCoeffMatrix, coeffMatrix, normalMatrixN );           //�󷨷��̵�N����
			nNormalMatrixN = Inverse( normalMatrixN );                      //�󷨷���N���������
			Multiply( tCoeffMatrix, constMatrix, constNormalMatrixW );      //�󷨷��̵ĳ�������W
			Multiply( nNormalMatrixN, constNormalMatrixW, correctMatrix );  //������ֵ�ĸ�����

			for (int j = 0; j < 11; j++)                                  //�������Ľ���ֵ���������²�����ֵ
				initDltL[j][0] = initDltL[j][0] + correctMatrix[j][0];

			iterTimes = iterTimes + 1;    //����������1

			maxCorrection = abs(correctMatrix[0][0]);   //����������ֵ�е����ֵ
			for (int j = 1; j < 11; j++)
				if (maxCorrection < abs(correctMatrix[j][0])) 
					maxCorrection = abs(correctMatrix[j][0]);
		}

		for (int i = 0; i < 11; i++)
			dltLs[i] = initDltL[i][0];//ƽ������󷵻ص�dlt��������
	}
	
	//���ü����11��Lϵ����������������ز���
	double temfm, temfz1, temfz2;
	temfm  = dltLs[8] * dltLs[8] + dltLs[9] * dltLs[9] + dltLs[10] * dltLs[10];
	temfz1 = dltLs[0] * dltLs[8] + dltLs[1] * dltLs[9] + dltLs[2] * dltLs[10];
	temfz2 = dltLs[4] * dltLs[8] + dltLs[5] * dltLs[9] + dltLs[6] * dltLs[10];

	camereParam.x0 = -temfz1 / temfm;   //�������x����x0
	camereParam.y0 = -temfz2 / temfm;   //�������y����y0

	double rr3 = 1.0 / temfm;     //�м�ֵ
	double temA2, temB2, temC2;   //3���м�ֵ

	temA2 = rr3 * (dltLs[0] * dltLs[0] + dltLs[1] * dltLs[1] + dltLs[2] * dltLs[2]) - camereParam.x0 * camereParam.x0;
	temB2 = rr3 * (dltLs[4] * dltLs[4] + dltLs[5] * dltLs[5] + dltLs[6] * dltLs[6]) - camereParam.y0 * camereParam.y0;
	temC2 = rr3 * (dltLs[0] * dltLs[4] + dltLs[1] * dltLs[5] + dltLs[2] * dltLs[6]) - camereParam.x0 * camereParam.y0;

	double fx, fy;
	fx = sqrt((temA2 * temB2 - temC2 * temC2) / temB2);
	fy = sqrt((temA2 * temB2 - temC2 * temC2) / temA2);

	camereParam.f = sqrt(fx * fx + fy * fy);

	//���������ĵ���ά����X��Y��Z
	mMatrix coeffMatrix2(3), constMatrix2(3);   //���XYZ��ϵ������ͳ�������
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

	//��Ϊ��3DDLT����ʱ���﷽������������Ļ��������Ļ�����﷽�����൱����������Ϊԭ��ģ����ԣ���dltϵ������������������Ҫ�����﷽���ĵ����꣬���ܺ��﷽���������ϵһ��
	camereParam.Xs = temCameraCoord[0][0] + weightSpaceCorrd.x;//������ĵ�X����
	camereParam.Ys = temCameraCoord[1][0] + weightSpaceCorrd.Y;//������ĵ�Y����
	camereParam.Zs = temCameraCoord[2][0] + weightSpaceCorrd.Z;//������ĵ�Z����


	//���������ⷽλ��Ԫ��phi, omega, kappa
	//���������ĵ���̬�ǣ������λ��fhfwjA��������fhxjqjaf
	//rotationMatrix��ת����˳��Ϊa1, b1, c1, a2, b2, c2, a3, b3, c3

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

	camereParam.alpha = acos(c3);             //�����ǣ����ȣ�
	camereParam.kappaAlpha = atan(a3 / b3);   // �����λ�ǣ����ȣ�
	

	//���ȱ�ʾ�������̬�ǻ���ɽǶ�
	camereParam.phi *= (180 / PI);
	camereParam.omega *= (180 / PI);
	camereParam.kappa *= (180 / PI);

	camereParam.alpha *= (180 / PI);
	camereParam.kappaAlpha *= (180 / PI);

	//������ת����
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

//���ݿ��Ƶ�ռ��������άDLT����Ĳ���������ͼ���������꣨���кţ�
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