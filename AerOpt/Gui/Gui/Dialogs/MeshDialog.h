#ifndef MESHDIALOG_H
#define MESHDIALOG_H

#include <QDialog>
#include "Enumerations.h"

namespace Ui {
class MeshDialog;
}

class OptimisationRun;

class MeshDialog : public QDialog
{
	Q_OBJECT

public:
    explicit MeshDialog(OptimisationRun& data, QWidget *parent = 0);
	~MeshDialog();

	void accept();

private slots:
	void on_course_toggled(bool checked);

	void on_medium_toggled(bool checked);

	void on_fine_toggled(bool checked);

    void on_viscous_toggled(bool checked);

private:
	Ui::MeshDialog *ui;
	Enum::Mesh mMeshDensity;
    OptimisationRun& mData;
};

#endif // MESHDIALOG_H
