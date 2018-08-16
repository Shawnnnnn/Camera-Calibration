#pragma once
#include <math.h>
#include <vector>
using namespace std;

class DLT3
{
public:
	DLT3(void);
	~DLT3(void);

	double a1, a2, a3, b1, b2, b3, c1, c2, c3;  //��ת����
	std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& Psr); //�����������
	bool multiply(const std::vector<std::vector<double>>& jz1, const std::vector<std::vector<double>>& jz2, std::vector<std::vector<double>>& scjz); //�������
	bool juzhxj(const std::vector<std::vector<double>>& jz1, const std::vector<std::vector<double>>& jz2, std::vector<std::vector<double>>& scjz,int jjbz);  //��������jz1,jz2��ӻ�����ĺ���(Ҫ��������������ά����ͬ),��jjbz=1ʱ��Ϊ������ӣ���jjbz=-1ʱ��Ϊ����������������scjz
	bool Transpose(const std::vector<std::vector<double>>& jz, std::vector<std::vector<double>>& scjz);  //����ת�ú���

	//���������6�����Ƶ���﷽xyz���ꡢ�����񷽵��������кţ�����3άDLT��ֱ�����Ա任���������������ⷽλԪ��
	//��Ҫ���룺���Ƶ�ĸ���dgs�����Ƶ���Ӱ���ϵ��к�����xfhh���к�����xflh�����Ƶ����﷽��X��������wfxzb���﷽Y��������wfyzb���﷽Z��������wfzzb��
	//�������أ�����Ľ���fhf�������Ӱ���ĵ�XYZ��������fhXYZ��0��1��2����Ԫ�طֱ���X��Y��Z��������ķ�λ��fhfwjA����λ���ȣ�����������fhxjqjaf����λ���ȣ�
	void DLT_cameraPara_6points(int dgs, const std::vector<double>& xfhh,const std::vector<double>& xflh, const std::vector<double> &wfxzb, const std::vector<double>& wfyzb, const std::vector<double>& wfzzb, double &fhfx, double &fhfy, double &x0, double &y0, std::vector<double>& fhXYZ, double &phi, double& omega, double& kappa, std::vector<double>& calculateImageX, std::vector<double>& calculateImageY);

};