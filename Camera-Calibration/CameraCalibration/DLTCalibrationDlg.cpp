// DLTCalibrationDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "CameraCalibration.h"
#include "DLTCalibrationDlg.h"
#include "afxdialogex.h"
#include "DLT/DLT.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;


// CDLTCalibrationDlg 对话框

IMPLEMENT_DYNAMIC(CDLTCalibrationDlg, CDialog)

CDLTCalibrationDlg::CDLTCalibrationDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CDLTCalibrationDlg::IDD, pParent)
	, m_strCPointTxt(_T(""))
	, m_dbx0(0)
	, m_dby0(0)
	, m_dbfx(0)
	, m_dbfy(0)
	, m_dbXs(0)
	, m_dbYs(0)
	, m_dbZs(0)
	, m_dbPhi(0)
	, m_dbOmega(0)
	, m_dbKappa(0)
{

}

CDLTCalibrationDlg::~CDLTCalibrationDlg()
{
}

void CDLTCalibrationDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_DLTEDIT, m_strCPointTxt);
	DDX_Text(pDX, IDC_EDIT_x0, m_dbx0);
	DDX_Text(pDX, IDC_EDIT_y0, m_dby0);
	DDX_Text(pDX, IDC_EDIT_fx, m_dbfx);
	DDX_Text(pDX, IDC_EDIT_fy, m_dbfy);
	DDX_Text(pDX, IDC_EDIT_Xs, m_dbXs);
	DDX_Text(pDX, IDC_EDIT_Ys, m_dbYs);
	DDX_Text(pDX, IDC_EDIT_Zs, m_dbZs);
	DDX_Text(pDX, IDC_EDIT_phi, m_dbPhi);
	DDX_Text(pDX, IDC_EDIT_omega, m_dbOmega);
	DDX_Text(pDX, IDC_EDIT_kappa, m_dbKappa);
}


BEGIN_MESSAGE_MAP(CDLTCalibrationDlg, CDialog)
	ON_BN_CLICKED(IDC_DLTBROWSEBUTTON, &CDLTCalibrationDlg::OnBnClickedDltbrowsebutton)
	ON_BN_CLICKED(IDOK, &CDLTCalibrationDlg::OnBnClickedOk)
	ON_EN_CHANGE(IDC_DLTEDIT, &CDLTCalibrationDlg::OnEnChangeDltedit)
END_MESSAGE_MAP()


// CDLTCalibrationDlg 消息处理程序

// 浏览按钮点击事件
void CDLTCalibrationDlg::OnBnClickedDltbrowsebutton()
{
	// TODO:  在此添加控件通知处理程序代码
	CFileDialog fileDialog(TRUE);//TRUE为打开，FALSE为保存
	fileDialog.m_ofn.lpstrTitle = _T("打开控制点文件");//标题
	fileDialog.m_ofn.lpstrFilter = _T("文本文件(*.txt)\0*.txt\0");//过滤器
	if (IDOK == fileDialog.DoModal())
	{
		m_strCPointTxt = fileDialog.GetPathName();
		UpdateData(FALSE);
	}
	
}

// OK按钮点击事件
void CDLTCalibrationDlg::OnBnClickedOk()
{
	// TODO:  在此添加控件通知处理程序代码
	if (m_strCPointTxt == "")
	{
		MessageBox(_T("请选择文件！"));
	}
	else
	{
		int dgs;
		string s;

		ifstream in(m_strCPointTxt);
		ofstream out("../out_L.txt");
		ofstream out2("../out_L_XYZ.txt");

		if (!in.is_open())
		{
			MessageBox(_T("未能打开文件！"));
			return;
		}
		getline(in, s);
		in >> dgs;

		double* pointID = (double*)malloc(sizeof(double) * dgs);
		double* imageX = (double*)malloc(sizeof(double) * dgs);
		double* imageY = (double*)malloc(sizeof(double) * dgs);
		double* spaceX = (double*)malloc(sizeof(double) * dgs);
		double* spaceY = (double*)malloc(sizeof(double) * dgs);
		double* spaceZ = (double*)malloc(sizeof(double) * dgs);

		double* cImageX = (double*)malloc(sizeof(double) * dgs);
		double* cImageY = (double*)malloc(sizeof(double) * dgs);

		for (int i = 0; i < dgs; i++)
		{
			string s;
			in >> s >> imageX[i] >> imageY[i] >> spaceX[i] >> spaceY[i] >> spaceZ[i];
		}
		in.close();

		if (dgs < 6)
		{
			MessageBox(_T("控制点数量小于6，无法计算！"));
			return;
		}

		DLT mDLT;
		mDLT.DLT_6pts(dgs, imageX, imageY, spaceX, spaceY, spaceZ, cImageX, cImageY);

		if (!out.is_open())
		{
			MessageBox(_T("保存参数文件时出错！"));
			return;
		}

		out << "*****DLT解法*****" << endl;
		out << " x0 = " << mDLT.x0 << endl;
		out << " y0 = " << mDLT.y0 << endl;
		out << " fx = " << mDLT.fx << endl;
		out << " fy = " << mDLT.fy << endl;
		out << " Xs = " << mDLT.Xs << endl;
		out << " Ys = " << mDLT.Ys << endl;
		out << " Zs = " << mDLT.Zs << endl;
		out << " phi = " << mDLT.phi << endl;
		out << " omega = " << mDLT.omega << endl;
		out << " kappa = " << mDLT.kappa << endl;

		out << " a1 = " << mDLT.a1 << endl;
		out << " b1 = " << mDLT.b1 << endl;
		out << " c1 = " << mDLT.c1 << endl;
		out << " a2 = " << mDLT.a2 << endl;
		out << " b2 = " << mDLT.b2 << endl;
		out << " c2 = " << mDLT.c2 << endl;
		out << " a3 = " << mDLT.a3 << endl;
		out << " b3 = " << mDLT.b3 << endl;
		out << " c3 = " << mDLT.c3 << endl;

		out.close();

		m_dbx0 = mDLT.x0;
		m_dby0 = mDLT.y0;
		m_dbfx = mDLT.fx;
		m_dbfy = mDLT.fy;
		m_dbXs = mDLT.Xs;
		m_dbYs = mDLT.Ys;
		m_dbZs = mDLT.Zs;
		m_dbPhi = mDLT.phi;
		m_dbOmega = mDLT.omega;
		m_dbKappa = mDLT.kappa;
		UpdateData(FALSE);

		if (!out2.is_open())
		{
			MessageBox(_T("保存误差文档时出错！"));
			return;
		}
		out2 << "解算X\t解算Y\t原始X\t原始Y\t误差X\t误差Y" << endl;
		for (int i = 0; i < dgs; i++)
		{
			double detaX = imageX[i] - cImageX[i];
			double detaY = imageY[i] - cImageY[i];
			out2 << cImageX[i] << "\t" << cImageY[i] << "\t" << imageX[i] << "\t" << imageY[i] << "\t" << detaX << "\t" << detaY << endl;
		}
		out2.close();

		MessageBox(_T("详细结果文档生成成功！"));
	}

	//CDialog::OnOK();
}

// 编辑框内容发生改变时响应事件
void CDLTCalibrationDlg::OnEnChangeDltedit()
{
	// TODO:  如果该控件是 RICHEDIT 控件，它将不
	// 发送此通知，除非重写 CDialog::OnInitDialog()
	// 函数并调用 CRichEditCtrl().SetEventMask()，
	// 同时将 ENM_CHANGE 标志“或”运算到掩码中。

	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}
