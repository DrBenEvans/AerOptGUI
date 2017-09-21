#ifndef MESHDIALOG_H
#define MESHDIALOG_H

#include <QtMath>
#include <QDialog>
#include <QGraphicsScene>
#include "ProfileModel.h"
#include "ProfileView.h"
#include "MeshView.h"

#include "Enumerations.h"
#include "MeshDialogModel.h"

namespace Ui {
class MeshDialog;
}

class Optimisation;

class MeshDialog : public QDialog
{
	Q_OBJECT

public:
    explicit MeshDialog(ProfileModel &profileModel, QWidget* parent = 0);
    ~MeshDialog();
    void accept();
    std::vector<BoundaryPoint*> controlPoints();

public slots:
    void runMesher();
    void setProfile();

private slots:
    void on_profile_currentIndexChanged(int index);

    void on_density_currentIndexChanged(int index);

    void on_layers_valueChanged(int arg1);

    void on_thickness_valueChanged(double arg1);

    void on_growthFactor_valueChanged(double arg1);

private:
    void on_pushButton_clicked();
    void setMeshActive(bool meshIsActive, bool doToggleProfile = true);

    Ui::MeshDialog *ui;
    MeshDialogModel* mMeshDialogModel;
    QGraphicsScene* mScene;
    ProfileModel& mProfileModel;
    ViewScaler* mScale;
    ProfileView* mProfileView;
    MeshView* mMeshView;
};

#endif // MESHDIALOG_H
