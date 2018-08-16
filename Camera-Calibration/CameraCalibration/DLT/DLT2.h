#pragma once
#include "Define.h"
#include <vector>
using namespace std;


class DLT2
{
public:
	DLT2(void);
	~DLT2(void);
	
	double pixelSize;	//影像的像素尺寸，单位为毫米	
	int imageHeight, imageWidth;	//原始图像的高度、宽度（行列数）	
	CameraParams camereParam;  //相机内外方位元素
	double dltLs[11];  //3DDLT的11个L系数
	CoordPoint3 weightSpaceCorrd;  //图像重心空间坐标
	PixelCoordPoint weightPixelCorrd;  //图像重心像素坐标
	double rotationMatrix[9];  //旋转矩阵

	void InitImageInfo(const int initImageHeight, const int initImageWidth, const double initPixelSize);  //添加图像信息
	void Calibrate(const int numPts, const vector<PixelCoordPoint>& pixelCoordPoints, const vector<CoordPoint3>& spaceCoordPoints);  //进行相机标定
	vector<PixelCoordPoint> CalculatePixcelCorrd(const vector<CoordPoint3>& spaceCoordPoints);  //根据标定参数反算像点坐标

};