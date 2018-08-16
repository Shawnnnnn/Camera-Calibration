#pragma once


// CViewAngleDlg 对话框

class CViewAngleDlg : public CDialog
{
	DECLARE_DYNAMIC(CViewAngleDlg)

public:
	CViewAngleDlg(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CViewAngleDlg();

// 对话框数据
	enum { IDD = IDD_ANGLEDLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
private:
	double m_camera_x;
	double m_camera_y;
	double m_camera_z;
	double m_targetPt_x;
	double m_targetPt_y;
	double m_targetPt_z;
	double m_azimuth;
	double m_pitchAngle;
public:
//	afx_msg void OnEnChangeCameraXEdit();
//	afx_msg void OnEnChangeCameraYEdit();
//	afx_msg void OnEnChangeCameraZEdit();
//	afx_msg void OnEnChangeTargetptXEdit();
//	afx_msg void OnEnChangeTargetptYEdit();
//	afx_msg void OnEnChangeTargetptZEdit();
	afx_msg void OnEnKillfocusTargetptXEdit();
	afx_msg void OnEnKillfocusTargetptYEdit();
	afx_msg void OnEnKillfocusTargetptZEdit();
	afx_msg void OnEnKillfocusCameraXEdit();
	afx_msg void OnEnKillfocusCameraYEdit();
	afx_msg void OnEnKillfocusCameraZEdit();
};
