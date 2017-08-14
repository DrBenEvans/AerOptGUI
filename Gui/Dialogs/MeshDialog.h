#ifndef MESHDIALOG_H
#define MESHDIALOG_H

#include <QtMath>
#include <QDialog>
#include <QGraphicsScene>
#include "ProfileModel.h"
#include "ProfileGraphicsItem.h"
#include "MeshGraphicsItem.h"

#include "Enumerations.h"
#include "Mesh.h"

namespace Ui {
class MeshDialog;
}

class Optimisation;

class MeshDialog : public QDialog
{
	Q_OBJECT

public:
    explicit MeshDialog(std::shared_ptr<Mesh> initMesh, ProfileModel &profileModel, QWidget* parent = 0);
    ~MeshDialog();
    void accept();

public slots:
    void runMesher();
    void setProfile();
    void controlPointStackedWidget();
    void meshStackedWidget();

private slots:
    void on_profile_currentIndexChanged(int index);

    void on_density_currentIndexChanged(int index);

    void on_layers_valueChanged(int arg1);

    void on_thickness_valueChanged(double arg1);

private:
    void on_pushButton_clicked();
    void setMeshActive(bool meshIsActive, bool doToggleProfile = true);

    Ui::MeshDialog *ui;
    std::shared_ptr<Mesh> mMesh;
    QGraphicsScene* mScene;
    ProfileModel& mProfileModel;
    int mScale = 1000;
    ProfileGraphicsItem* mProfileGraphicsItem;
    MeshGraphicsItem* mMeshGraphicsItem;
};

#endif // MESHDIALOG_H
