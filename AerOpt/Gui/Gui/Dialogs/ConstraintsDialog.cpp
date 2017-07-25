/*********************************************
**
**	Created on: 	28/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ConstraintsDialog.cpp
**
**********************************************/

#include <QDebug>
#include "ConstraintsDialog.h"
#include "ui_ConstraintsDialog.h"
#include "ProjectData.h"

ConstraintsDialog::ConstraintsDialog(ProjectData& data, QWidget *parent) :
	QDialog(parent),
	ui(new Ui::ConstraintsDialog),
	mData(data)
{
	ui->setupUi(this);
	mIndex = 0;
}

ConstraintsDialog::~ConstraintsDialog()
{
	delete ui;
}

void ConstraintsDialog::setConstraint(const unsigned int index)
{
	mIndex = index;

    BoundaryPoint *point = mData.getControlPoint(0, index);

	qreal x1,y1,x2,y2;
    point->getBoundCoords(&x1,&y1,&x2,&y2);

	ui->xMin->setValue( x1 );
	ui->xMax->setValue( x2 );
	ui->yMin->setValue( y1 );
	ui->yMax->setValue( y2 );
}

void ConstraintsDialog::accept()
{
    BoundaryPoint *point = mData.getControlPoint(0, mIndex);

    point->setBoundCoords(ui->xMin->value(),
						   ui->yMin->value(),
						   ui->xMax->value(),
						   ui->yMax->value());

	QDialog::accept();
}

void ConstraintsDialog::reject()
{
	QDialog::reject();
	//Do Nothing!
}
