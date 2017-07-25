/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ObjectiveDialog.h
**
**********************************************/

#ifndef OBJECTIVEDIALOG_H
#define OBJECTIVEDIALOG_H

#include <QDialog>
#include "Enumerations.h"

class ProjectData;

namespace Ui {
class ObjectiveDialog;
}

class ObjectiveDialog : public QDialog
{
	Q_OBJECT

public:
	explicit ObjectiveDialog(ProjectData& data, QWidget *parent = 0);
	~ObjectiveDialog();

	void accept();
	void reject();

private slots:
	void on_liftDrag_toggled(bool checked);

	void on_maxLift_toggled(bool checked);

	void on_minDrag_toggled(bool checked);

	void on_maxDownforce_toggled(bool checked);

	void on_minLift_toggled(bool checked);

private:
	Ui::ObjectiveDialog *ui;
	Enum::ObjFunc mObjFunc;
	ProjectData& mData;
};

#endif // OBJECTIVEDIALOG_H
