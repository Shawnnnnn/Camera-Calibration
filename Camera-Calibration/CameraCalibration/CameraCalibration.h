
// CameraCalibration.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CCameraCalibrationApp: 
// �йش����ʵ�֣������ CameraCalibration.cpp
//

class CCameraCalibrationApp : public CWinApp
{
public:
	CCameraCalibrationApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CCameraCalibrationApp theApp;