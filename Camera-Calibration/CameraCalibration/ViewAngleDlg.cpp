// ViewAngleDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "CameraCalibration.h"
#include "ViewAngleDlg.h"
#include "afxdialogex.h"
#include "AngleCalculation.h"


// CViewAngleDlg 对话框

IMPLEMENT_DYNAMIC(CViewAngleDlg, CDialog)

CViewAngleDlg::CViewAngleDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CViewAngleDlg::IDD, pParent)
	, m_camera_x(0)
	, m_camera_y(0)
	, m_camera_z(0)
	, m_targetPt_x(0)
	, m_targetPt_y(0)
	, m_targetPt_z(0)
	, m_azimuth(0) 
	, m_pitchAngle(0)
{

}

CViewAngleDlg::~CViewAngleDlg()
{
}

void CViewAngleDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_CAMERA_X_EDIT, m_camera_x);
	DDX_Text(pDX, IDC_CAMERA_Y_EDIT, m_camera_y);
	DDX_Text(pDX, IDC_CAMERA_Z_EDIT, m_camera_z);
	DDX_Text(pDX, IDC_TARGETPT_X_EDIT, m_targetPt_x);
	DDX_Text(pDX, IDC_TARGETPT_Y_EDIT, m_targetPt_y);
	DDX_Text(pDX, IDC_TARGETPT_Z_EDIT, m_targetPt_z);
	DDX_Text(pDX, IDC_AZIMUTH_EDIT, m_azimuth);
	DDX_Text(pDX, IDC_PITCHANGLE_EDIT, m_pitchAngle);
}


BEGIN_MESSAGE_MAP(CViewAngleDlg, CDialog)
	ON_BN_CLICKED(IDOK, &CViewAngleDlg::OnBnClickedOk)
	ON_EN_KILLFOCUS(IDC_TARGETPT_X_EDIT, &CViewAngleDlg::OnEnKillfocusTargetptXEdit)
	ON_EN_KILLFOCUS(IDC_TARGETPT_Y_EDIT, &CViewAngleDlg::OnEnKillfocusTargetptYEdit)
	ON_EN_KILLFOCUS(IDC_TARGETPT_Z_EDIT, &CViewAngleDlg::OnEnKillfocusTargetptZEdit)
	ON_EN_KILLFOCUS(IDC_CAMERA_X_EDIT, &CViewAngleDlg::OnEnKillfocusCameraXEdit)
	ON_EN_KILLFOCUS(IDC_CAMERA_Y_EDIT, &CViewAngleDlg::OnEnKillfocusCameraYEdit)
	ON_EN_KILLFOCUS(IDC_CAMERA_Z_EDIT, &CViewAngleDlg::OnEnKillfocusCameraZEdit)
END_MESSAGE_MAP()


// CViewAngleDlg 消息处理程序


void CViewAngleDlg::OnBnClickedOk()
{
	// TODO:  在此添加控件通知处理程序代码
	COORDINATION camera;
	COORDINATION targetPt;

	camera.x = m_camera_x;
	camera.y = m_camera_y;
	camera.z = m_camera_z;

	targetPt.x = m_targetPt_x;
	targetPt.y = m_targetPt_y;
	targetPt.z = m_targetPt_z;

	CAngleCalculation ac;
	m_azimuth = ac.GetAzimuth(camera, targetPt);
	m_pitchAngle = ac.GetPitchAngle(camera, targetPt);

	UpdateData(FALSE);
	//CDialog::OnOK();
}

void CViewAngleDlg::OnEnKillfocusTargetptXEdit()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}


void CViewAngleDlg::OnEnKillfocusTargetptYEdit()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}


void CViewAngleDlg::OnEnKillfocusTargetptZEdit()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}


void CViewAngleDlg::OnEnKillfocusCameraXEdit()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}


void CViewAngleDlg::OnEnKillfocusCameraYEdit()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}


void CViewAngleDlg::OnEnKillfocusCameraZEdit()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(TRUE);
}