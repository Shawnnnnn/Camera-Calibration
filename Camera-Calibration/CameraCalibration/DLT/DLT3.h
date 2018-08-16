#pragma once
#include <math.h>
#include <vector>
using namespace std;

class DLT3
{
public:
	DLT3(void);
	~DLT3(void);

	double a1, a2, a3, b1, b2, b3, c1, c2, c3;  //旋转矩阵
	std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& Psr); //方阵求逆矩阵
	bool multiply(const std::vector<std::vector<double>>& jz1, const std::vector<std::vector<double>>& jz2, std::vector<std::vector<double>>& scjz); //矩阵相乘
	bool juzhxj(const std::vector<std::vector<double>>& jz1, const std::vector<std::vector<double>>& jz2, std::vector<std::vector<double>>& scjz,int jjbz);  //两个矩阵jz1,jz2相加或相减的函数(要求，两个输入矩阵的维数相同),当jjbz=1时，为矩阵相加；当jjbz=-1时，为矩阵相减，输出矩阵scjz
	bool Transpose(const std::vector<std::vector<double>>& jz, std::vector<std::vector<double>>& scjz);  //矩阵转置函数

	//根据输入的6个控制点的物方xyz坐标、及其像方的像素行列号，进行3维DLT（直接线性变换），求解相机的内外方位元素
	//需要输入：控制点的个数dgs，控制点在影像上的行号数组xfhh、列号数组xflh，控制点在物方的X坐标数组wfxzb、物方Y坐标数组wfyzb、物方Z坐标数组wfzzb；
	//函数返回：相机的焦距fhf，相机摄影中心的XYZ坐标数组fhXYZ（0，1，2三个元素分别是X、Y、Z），相机的方位角fhfwjA（单位弧度）、相机的倾角fhxjqjaf（单位弧度）
	void DLT_cameraPara_6points(int dgs, const std::vector<double>& xfhh,const std::vector<double>& xflh, const std::vector<double> &wfxzb, const std::vector<double>& wfyzb, const std::vector<double>& wfzzb, double &fhfx, double &fhfy, double &x0, double &y0, std::vector<double>& fhXYZ, double &phi, double& omega, double& kappa, std::vector<double>& calculateImageX, std::vector<double>& calculateImageY);

};