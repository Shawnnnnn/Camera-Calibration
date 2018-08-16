
// CameraCalibrationDlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "CameraCalibration.h"
#include "CameraCalibrationDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CCameraCalibrationDlg �Ի���



CCameraCalibrationDlg::CCameraCalibrationDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CCameraCalibrationDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CCameraCalibrationDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TAB, m_tab);
}

BEGIN_MESSAGE_MAP(CCameraCalibrationDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_NOTIFY(TCN_SELCHANGE, IDC_TAB, &CCameraCalibrationDlg::OnTcnSelchangeTab)
END_MESSAGE_MAP()


// CCameraCalibrationDlg ��Ϣ�������

BOOL CCameraCalibrationDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// ���ô˶Ի����ͼ�ꡣ  ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO:  �ڴ���Ӷ���ĳ�ʼ������
	m_tab.InsertItem(0, _T("���Ƶ㶨��"));
	m_tab.InsertItem(1, _T("���̸񶨱�"));
	m_tab.InsertItem(2, _T("���߽Ƕȼ���"));

	m_DLTCalibration.Create(IDD_DLTDLG, GetDlgItem(IDC_TAB));
	m_CBoradCalibration.Create(IDD_CKBDLG, GetDlgItem(IDC_TAB));
	m_AngleCalculation.Create(IDD_ANGLEDLG, GetDlgItem(IDC_TAB));

	//���IDC_TAB�ͻ�����С 
	CRect rs;
	m_tab.GetClientRect(&rs);
	//�����ӶԻ����ڸ������е�λ�� 
	rs.top += 25;
	rs.bottom -= 1;
	rs.left += 1;
	rs.right -= 2;

	//�����ӶԻ���ߴ粢�ƶ���ָ��λ�� 
	m_DLTCalibration.MoveWindow(&rs);
	m_CBoradCalibration.MoveWindow(&rs);
	m_AngleCalculation.MoveWindow(&rs);

	//�ֱ��������غ���ʾ 
	m_DLTCalibration.ShowWindow(true);
	m_CBoradCalibration.ShowWindow(false);
	m_AngleCalculation.ShowWindow(false);

	//����Ĭ�ϵ�ѡ� 
	m_tab.SetCurSel(0);

	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
}

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ  ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void CCameraCalibrationDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR CCameraCalibrationDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


// ��ǩ�ı��¼�
void CCameraCalibrationDlg::OnTcnSelchangeTab(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO:  �ڴ���ӿؼ�֪ͨ����������
	switch (m_tab.GetCurSel())
	{
	case 0:
		m_DLTCalibration.ShowWindow(true);
		m_CBoradCalibration.ShowWindow(false);
		m_AngleCalculation.ShowWindow(false);
		break;
	case 1:
		m_DLTCalibration.ShowWindow(false);
		m_CBoradCalibration.ShowWindow(true);
		m_AngleCalculation.ShowWindow(false);
		break;
	case 2:
		m_DLTCalibration.ShowWindow(false);
		m_CBoradCalibration.ShowWindow(false);
		m_AngleCalculation.ShowWindow(true);
	default:
		break;
	}

	*pResult = 0;
}
