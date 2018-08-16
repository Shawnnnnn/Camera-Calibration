#pragma once


// CDLTCalibrationDlg 对话框

class CDLTCalibrationDlg : public CDialog
{
	DECLARE_DYNAMIC(CDLTCalibrationDlg)

public:
	CDLTCalibrationDlg(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CDLTCalibrationDlg();

// 对话框数据
	enum { IDD = IDD_DLTDLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedDltbrowsebutton();
	CString m_strCPointTxt;
	afx_msg void OnBnClickedOk();
	afx_msg void OnEnChangeDltedit();
	double m_dbx0;
	double m_dby0;
	double m_dbfx;
	double m_dbfy;
	double m_dbXs;
	double m_dbYs;
	double m_dbZs;
	double m_dbPhi;
	double m_dbOmega;
	double m_dbKappa;
};
