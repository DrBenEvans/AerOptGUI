/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		BoundaryDialog.h
**
**********************************************/

#ifndef BOUNDARYDIALOG_H
#define BOUNDARYDIALOG_H

#include <QDialog>

class ProjectData;

namespace Ui {
class BoundaryDialog;
}

class BoundaryDialog : public QDialog
{
	Q_OBJECT

public:
	explicit BoundaryDialog(ProjectData& data, QWidget *parent = 0);
	~BoundaryDialog();

	void accept();
	void reject();

private slots:
	void on_temp_valueChanged(int arg1);

	void on_pressure_valueChanged(int arg1);

	void on_angle_valueChanged(int arg1);

	void on_reynolds_valueChanged(int arg1);

	void on_mach_valueChanged(double arg1);

private:
	Ui::BoundaryDialog *ui;
	float mMachNo;
	float mReNo;
	float mFreeAlpha;
	float mFreePress;
	float mFreeTemp;
	ProjectData& mData;
};

#endif // BOUNDARYDIALOG_H
