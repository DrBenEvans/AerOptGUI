/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ObjectiveDialog.cpp
**
**********************************************/

#include "ObjectiveDialog.h"
#include "ui_ObjectiveDialog.h"
#include "ProjectData.h"

ObjectiveDialog::ObjectiveDialog(ProjectData& data, QWidget *parent) :
	QDialog(parent),
	ui(new Ui::ObjectiveDialog),
	mData(data)
{
	ui->setupUi(this);
	mObjFunc = mData.objFunc();

	switch (mObjFunc)
	{
		case Enum::ObjFunc::LIFTDRAG :
			ui->liftDrag->setChecked(true);
			break;
		case Enum::ObjFunc::DISTORTION :
//			ui->distortion->setChecked(true);
			break;
		case Enum::ObjFunc::MAXLIFT :
			ui->maxLift->setChecked(true);
			break;
		case Enum::ObjFunc::MINDRAG :
			ui->minDrag->setChecked(true);
			break;
		case Enum::ObjFunc::MAXDOWNFORCE :
			ui->maxDownforce->setChecked(true);
			break;
		case Enum::ObjFunc::MINLIFT :
			ui->minLift->setChecked(true);
			break;
	}
}

ObjectiveDialog::~ObjectiveDialog()
{
	delete ui;
}

void ObjectiveDialog::accept()
{
	mData.setObjFunc(mObjFunc);
	mData.setFunction(true);

	QDialog::accept();
}

void ObjectiveDialog::reject()
{
	QDialog::reject();
	//Do Nothing!
}


void ObjectiveDialog::on_liftDrag_toggled(bool checked)
{
	if (checked) mObjFunc = Enum::ObjFunc::LIFTDRAG;
}

void ObjectiveDialog::on_maxLift_toggled(bool checked)
{
	if (checked) mObjFunc = Enum::ObjFunc::MAXLIFT;
}

void ObjectiveDialog::on_minDrag_toggled(bool checked)
{
	if (checked) mObjFunc = Enum::ObjFunc::MINDRAG;
}

void ObjectiveDialog::on_maxDownforce_toggled(bool checked)
{
	if (checked) mObjFunc = Enum::ObjFunc::MAXDOWNFORCE;
}

void ObjectiveDialog::on_minLift_toggled(bool checked)
{
	if (checked) mObjFunc = Enum::ObjFunc::MINLIFT;
}
