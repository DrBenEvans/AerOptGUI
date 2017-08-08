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
#include <OptimisationRun.h>

namespace Ui
{
	class ConstraintsDialog;
}

class ConstraintsDialog : public QDialog
{
	Q_OBJECT

public:
    explicit ConstraintsDialog(Mesh& mesh, QWidget *parent = 0);
    ~ConstraintsDialog();
	void setConstraint(const unsigned int index);

private:
	Ui::ConstraintsDialog *ui;
	unsigned int mIndex;
    Mesh& mMesh;

protected slots:
	void accept();
	void reject();
};

#endif // CONSTRAINTSDIALOG_H
