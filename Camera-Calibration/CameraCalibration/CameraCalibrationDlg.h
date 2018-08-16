
// CameraCalibrationDlg.h : 头文件
//

#pragma once
#include "afxcmn.h"
#include "DLTCalibrationDlg.h"
#include "CheckerboardCalibrationDlg.h"
#include "ViewAngleDlg.h"


// CCameraCalibrationDlg 对话框
class CCameraCalibrationDlg : public CDialogEx
{
// 构造
public:
	CCameraCalibrationDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_CAMERACALIBRATION_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	CTabCtrl m_tab;
	CDLTCalibrationDlg m_DLTCalibration;
	CCheckerboardCalibrationDlg m_CBoradCalibration;
	afx_msg void OnTcnSelchangeTab(NMHDR *pNMHDR, LRESULT *pResult);
	CViewAngleDlg m_AngleCalculation;
};
