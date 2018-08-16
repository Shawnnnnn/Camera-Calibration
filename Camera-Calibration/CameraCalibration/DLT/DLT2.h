#pragma once
#include "Define.h"
#include <vector>
using namespace std;


class DLT2
{
public:
	DLT2(void);
	~DLT2(void);
	
	double pixelSize;	//Ӱ������سߴ磬��λΪ����	
	int imageHeight, imageWidth;	//ԭʼͼ��ĸ߶ȡ���ȣ���������	
	CameraParams camereParam;  //������ⷽλԪ��
	double dltLs[11];  //3DDLT��11��Lϵ��
	CoordPoint3 weightSpaceCorrd;  //ͼ�����Ŀռ�����
	PixelCoordPoint weightPixelCorrd;  //ͼ��������������
	double rotationMatrix[9];  //��ת����

	void InitImageInfo(const int initImageHeight, const int initImageWidth, const double initPixelSize);  //���ͼ����Ϣ
	void Calibrate(const int numPts, const vector<PixelCoordPoint>& pixelCoordPoints, const vector<CoordPoint3>& spaceCoordPoints);  //��������궨
	vector<PixelCoordPoint> CalculatePixcelCorrd(const vector<CoordPoint3>& spaceCoordPoints);  //���ݱ궨���������������

};