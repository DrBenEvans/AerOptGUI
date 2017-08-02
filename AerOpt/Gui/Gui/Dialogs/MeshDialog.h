#ifndef MESHDIALOG_H
#define MESHDIALOG_H

#include <QDialog>
#include "Enumerations.h"
#include "Mesh.h"

namespace Ui {
class MeshDialog;
}

class OptimisationRun;

class MeshDialog : public QDialog
{
	Q_OBJECT

public:
    explicit MeshDialog(QSharedPointer<Mesh> mesh, QWidget *parent = 0);
	~MeshDialog();

	void accept();

private slots:
    void on_viscous_toggled(bool checked);

private:
	Ui::MeshDialog *ui;
    QSharedPointer<Mesh> mMesh;
};

#endif // MESHDIALOG_H
