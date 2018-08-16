//#pragma once
//
//#include <iostream>
//using namespace std;

class DLT
{
public:
	DLT(void);
	~DLT(void);

	double x0, y0;  //������
	double fx, fy;  //��Ƭ����
	double Xs, Ys, Zs;  //�ⷽλֱ��Ԫ�أ���Ӱ����S������
	double phi, omega, kappa;  //�ⷽλ��Ԫ��
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;  //��ת����

	void DLT_6pts(int dgs, double* imageX, double* imageY, double* spaceX, double* spaceY, double* spaceZ, double* cImageX, double* cImageY);
	void TransposeMatrix(double *m1,double *m2,int m,int n);   //����ת��,ת�ú���󸳸�m2��mΪm1��������nΪm2������
	void Inverse(double *a,int n);                             //��������,n�׾���a�����֮���������Ȼ����a
	void Multiply(double *m1,double *m2,double *result, int m,int n,int l);  //������ˣ�����resultָ���������

};