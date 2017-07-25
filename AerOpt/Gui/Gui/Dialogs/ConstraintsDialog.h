/*********************************************
**
**	Created on: 	28/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ConstraintsDialog.h
**
**********************************************/

#ifndef CONSTRAINTSDIALOG_H
#define CONSTRAINTSDIALOG_H

#include <QDialog>

namespace Ui
{
	class ConstraintsDialog;
}

class ProjectData;

class ConstraintsDialog : public QDialog
{
	Q_OBJECT

public:
	explicit ConstraintsDialog(ProjectData& data, QWidget *parent = 0);
	~ConstraintsDialog();
	void setConstraint(const unsigned int index);

private:
	Ui::ConstraintsDialog *ui;
	unsigned int mIndex;
	ProjectData& mData;

protected slots:
	void accept();
	void reject();
};

#endif // CONSTRAINTSDIALOG_H
