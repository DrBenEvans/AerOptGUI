/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimiserDialog.cpp
**
**********************************************/

#include "OptimiserDialog.h"
#include "ui_OptimiserDialog.h"
#include "ProjectData.h"

OptimiserDialog::OptimiserDialog(ProjectData& data, QWidget *parent) :
	QDialog(parent),
	ui(new Ui::OptimiserDialog),
	mData(data)
{
	ui->setupUi(this);
	mNoAgents = mData.noAgents();
	mNoGens = mData.noGens();
    mNoTop = mData.getNoTop();

	ui->nests->setValue(mNoAgents);
	ui->gens->setValue(mNoGens);
	ui->PercentTop->setValue(mNoTop);
}

OptimiserDialog::~OptimiserDialog()
{
	delete ui;
}

void OptimiserDialog::accept()
{
	mData.setNoAgents(mNoAgents);
	mData.setNoGens(mNoGens);
	mData.setOptimiser(true);
    mData.setNoTop(mNoTop);

	QDialog::accept();
}

void OptimiserDialog::reject()
{
	QDialog::reject();
	//Do Nothing!
}

void OptimiserDialog::on_nests_valueChanged(int arg1)
{
	mNoAgents = arg1;
}

void OptimiserDialog::on_gens_valueChanged(int arg1)
{
	mNoGens = arg1;
}

void OptimiserDialog::on_PercentTop_valueChanged(int arg1)
{
    mNoTop = arg1;
}
