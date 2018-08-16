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
	double f;  //����f�������Ӱ��������Ƭ��Ĵ�ֱ���룩
	double x0, y0;  //����������x0,y0��������ϵ��
	double Xs, Ys, Zs;  //����������꣨�﷽����ϵ��
	double phi, omega, kappa;  //�ⷽλ��Ԫ�ئצئ�
	double alpha, kappaAlpha;  //�����Ǧ�,�����λ�Ǧ�

};