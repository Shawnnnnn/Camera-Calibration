#include "stdafx.h"
#include "AngleCalculation.h"

#define M_PI       3.14159265358979323846

CAngleCalculation::CAngleCalculation()
{
}


CAngleCalculation::~CAngleCalculation()
{
}

double CAngleCalculation::GetAzimuth(COORDINATION camera, COORDINATION irradiationPt)
{
	double azimuth = 0.0;
	double tan_a = (irradiationPt.x - camera.x) / (irradiationPt.y - camera.y);
	azimuth = atan(tan_a) * 180 / M_PI;
	if (irradiationPt.y < camera.y)
	{
		return azimuth + 180;
	}
	else if (irradiationPt.x < camera.x)
	{
		return azimuth + 360;
	}
	return azimuth;
}

double CAngleCalculation::GetPitchAngle(COORDINATION camera, COORDINATION irradiationPt)
{
	double pAngle = 0.0;
	double tan_p = (irradiationPt.z - camera.z) / abs(irradiationPt.y - camera.y);
	pAngle = atan(tan_p) * 180 / M_PI;
	if (irradiationPt.z < camera.z && irradiationPt.y < camera.y)
	{
		return pAngle - 180;
	}
	else if (irradiationPt.z > camera.z && irradiationPt.y < camera.y)
	{
		return pAngle + 180;
	}
	return pAngle;
}