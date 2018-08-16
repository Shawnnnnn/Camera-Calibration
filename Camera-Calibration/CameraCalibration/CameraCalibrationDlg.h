
// CameraCalibrationDlg.h : ͷ�ļ�
//

#pragma once
#include "afxcmn.h"
#include "DLTCalibrationDlg.h"
#include "CheckerboardCalibrationDlg.h"
#include "ViewAngleDlg.h"


// CCameraCalibrationDlg �Ի���
class CCameraCalibrationDlg : public CDialogEx
{
// ����
public:
	CCameraCalibrationDlg(CWnd* pParent = NULL);	// ��׼���캯��

// �Ի�������
	enum { IDD = IDD_CAMERACALIBRATION_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��


// ʵ��
protected:
	HICON m_hIcon;

	// ���ɵ���Ϣӳ�亯��
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
