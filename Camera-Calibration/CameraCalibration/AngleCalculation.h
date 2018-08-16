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
	static double GetAzimuth(COORDINATION camera, COORDINATION irradiationPt);// 计算方位角
	static double GetPitchAngle(COORDINATION camera, COORDINATION irradiationPt);// 计算俯仰角
};

