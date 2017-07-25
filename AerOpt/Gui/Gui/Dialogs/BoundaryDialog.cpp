/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		BoundaryDialog.cpp
**
**********************************************/

#include "BoundaryDialog.h"
#include "ui_BoundaryDialog.h"
#include "ProjectData.h"

BoundaryDialog::BoundaryDialog(ProjectData& data, QWidget *parent) :
	QDialog(parent),
	ui(new Ui::BoundaryDialog),
	mData(data)
{
	ui->setupUi(this);
	mMachNo = mData.machNo();
	mReNo = mData.reNo();
	mFreeAlpha = mData.freeAlpha();
	mFreePress = mData.freePress();
	mFreeTemp = mData.freeTemp();

	ui->temp->setValue(mFreeTemp);
	ui->pressure->setValue(mFreePress);
	ui->angle->setValue(mFreeAlpha);
	ui->reynolds->setValue(mReNo);
	ui->mach->setValue(mMachNo);
}

BoundaryDialog::~BoundaryDialog()
{
	delete ui;
}

void BoundaryDialog::accept()
{
	mData.setMachNo(mMachNo);
	mData.setReNo(mReNo);
	mData.setFreeAlpha(mFreeAlpha);
	mData.setFreePress(mFreePress);
	mData.setFreeTemp(mFreeTemp);
	mData.setBoundary(true);

	QDialog::accept();
}

void BoundaryDialog::reject()
{
	QDialog::reject();
	//Do Nothing!
}

void BoundaryDialog::on_temp_valueChanged(int arg1)
{
	mFreeTemp = arg1;
}

void BoundaryDialog::on_pressure_valueChanged(int arg1)
{
	mFreePress = arg1;
}

void BoundaryDialog::on_angle_valueChanged(int arg1)
{
	mFreeAlpha = arg1;
}

void BoundaryDialog::on_reynolds_valueChanged(int arg1)
{
	mReNo = arg1;
}

void BoundaryDialog::on_mach_valueChanged(double arg1)
{
	mMachNo = arg1;
}
