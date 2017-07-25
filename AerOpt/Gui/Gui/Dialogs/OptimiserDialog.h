/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimiserDialog.h
**
**********************************************/

#ifndef OPTIMISERDIALOG_H
#define OPTIMISERDIALOG_H

#include <QDialog>

class ProjectData;

namespace Ui {
class OptimiserDialog;
}

class OptimiserDialog : public QDialog
{
	Q_OBJECT

public:
	explicit OptimiserDialog(ProjectData& data, QWidget *parent = 0);
	~OptimiserDialog();

	void accept();
	void reject();

private slots:
	void on_nests_valueChanged(int arg1);

	void on_gens_valueChanged(int arg1);

    void on_PercentTop_valueChanged(int arg1);

private:
	Ui::OptimiserDialog *ui;
	int mNoAgents;
	int mNoGens;
    int mNoTop;
	ProjectData& mData;
};

#endif // OPTIMISERDIALOG_H
