#pragma once


// CCheckerboardCalibrationDlg �Ի���

class CCheckerboardCalibrationDlg : public CDialog
{
	DECLARE_DYNAMIC(CCheckerboardCalibrationDlg)

public:
	CCheckerboardCalibrationDlg(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~CCheckerboardCalibrationDlg();

// �Ի�������
	enum { IDD = IDD_CKBDLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

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
