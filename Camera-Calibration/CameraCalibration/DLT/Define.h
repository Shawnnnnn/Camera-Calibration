#pragma once

struct CoordPoint2
{
	double x;
	double y;
};

struct CoordPoint3
{
	double x;
	double Y;
	double Z;
};

struct PixelCoordPoint
{
	int row;
	int col;
};

struct CameraParams
{
	double f;  //主距f（相机摄影中心与相片间的垂直距离）
	double x0, y0;  //像主点坐标x0,y0（像方坐标系）
	double Xs, Ys, Zs;  //相机中心坐标（物方坐标系）
	double phi, omega, kappa;  //外方位角元素ψωκ
	double alpha, kappaAlpha;  //相机倾角α,相机方位角κ

};