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
#include "Mesh.h"

ConstraintsDialog::ConstraintsDialog(Mesh& mesh, QWidget *parent) :
	QDialog(parent),
	ui(new Ui::ConstraintsDialog),
    mMesh(mesh)
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

    BoundaryPoint& point = mMesh.getControlPoint(index);

	qreal x1,y1,x2,y2;
    point.getBoundCoords(&x1,&y1,&x2,&y2);

	ui->xMin->setValue( x1 );
	ui->xMax->setValue( x2 );
	ui->yMin->setValue( y1 );
	ui->yMax->setValue( y2 );
    ui->smoothing->setValue( point.getSmoothing() );
}

void ConstraintsDialog::accept()
{
    BoundaryPoint& point = mMesh.getControlPoint(mIndex);

    point.setBoundCoords(ui->xMin->value(),
						   ui->yMin->value(),
						   ui->xMax->value(),
						   ui->yMax->value());

    point.setSmoothing(ui->smoothing->value());

	QDialog::accept();
}

void ConstraintsDialog::reject()
{
	QDialog::reject();
	//Do Nothing!
}
