
// CameraCalibrationDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "CameraCalibration.h"
#include "CameraCalibrationDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CCameraCalibrationDlg 对话框



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


// CCameraCalibrationDlg 消息处理程序

BOOL CCameraCalibrationDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO:  在此添加额外的初始化代码
	m_tab.InsertItem(0, _T("控制点定标"));
	m_tab.InsertItem(1, _T("棋盘格定标"));
	m_tab.InsertItem(2, _T("视线角度计算"));

	m_DLTCalibration.Create(IDD_DLTDLG, GetDlgItem(IDC_TAB));
	m_CBoradCalibration.Create(IDD_CKBDLG, GetDlgItem(IDC_TAB));
	m_AngleCalculation.Create(IDD_ANGLEDLG, GetDlgItem(IDC_TAB));

	//获得IDC_TAB客户区大小 
	CRect rs;
	m_tab.GetClientRect(&rs);
	//调整子对话框在父窗口中的位置 
	rs.top += 25;
	rs.bottom -= 1;
	rs.left += 1;
	rs.right -= 2;

	//设置子对话框尺寸并移动到指定位置 
	m_DLTCalibration.MoveWindow(&rs);
	m_CBoradCalibration.MoveWindow(&rs);
	m_AngleCalculation.MoveWindow(&rs);

	//分别设置隐藏和显示 
	m_DLTCalibration.ShowWindow(true);
	m_CBoradCalibration.ShowWindow(false);
	m_AngleCalculation.ShowWindow(false);

	//设置默认的选项卡 
	m_tab.SetCurSel(0);

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CCameraCalibrationDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CCameraCalibrationDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


// 标签改变事件
void CCameraCalibrationDlg::OnTcnSelchangeTab(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO:  在此添加控件通知处理程序代码
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
