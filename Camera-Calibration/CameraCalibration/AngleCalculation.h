#pragma once

struct COORDINATION
{
	double x;
	double y;
	double z;
};

class CAngleCalculation
{
public:
	CAngleCalculation();
	CAngleCalculation(COORDINATION camera, COORDINATION irradiationPt);
	~CAngleCalculation();

public:
	static double GetAzimuth(COORDINATION camera, COORDINATION irradiationPt);// ���㷽λ��
	static double GetPitchAngle(COORDINATION camera, COORDINATION irradiationPt);// ���㸩����
};

