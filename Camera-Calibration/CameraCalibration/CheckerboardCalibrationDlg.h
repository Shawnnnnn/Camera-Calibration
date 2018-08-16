#pragma once


// CCheckerboardCalibrationDlg 对话框

class CCheckerboardCalibrationDlg : public CDialog
{
	DECLARE_DYNAMIC(CCheckerboardCalibrationDlg)

public:
	CCheckerboardCalibrationDlg(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CCheckerboardCalibrationDlg();

// 对话框数据
	enum { IDD = IDD_CKBDLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedCbcbrowserbutton();
	CString m_strCBPicTxt;
	afx_msg void OnEnChangeCbcedit();
	afx_msg void OnBnClickedOk();
	double m_dbx0;
	double m_dby0;
	double m_dbfx;
	double m_dbfy;
};
