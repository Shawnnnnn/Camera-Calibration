//#pragma once
//
//#include <iostream>
//using namespace std;

class DLT
{
public:
	DLT(void);
	~DLT(void);

	double x0, y0;  //像主点
	double fx, fy;  //相片主距
	double Xs, Ys, Zs;  //外方位直线元素，摄影中心S的坐标
	double phi, omega, kappa;  //外方位角元素
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;  //旋转矩阵

	void DLT_6pts(int dgs, double* imageX, double* imageY, double* spaceX, double* spaceY, double* spaceZ, double* cImageX, double* cImageY);
	void TransposeMatrix(double *m1,double *m2,int m,int n);   //矩阵转置,转置后矩阵赋给m2，m为m1的列数，n为m2的列数
	void Inverse(double *a,int n);                             //矩阵求逆,n阶矩阵a，求出之后逆矩阵依然赋给a
	void Multiply(double *m1,double *m2,double *result, int m,int n,int l);  //矩阵相乘，其中result指向输出矩阵

};