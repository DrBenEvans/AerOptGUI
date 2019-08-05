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
    explicit MeshDialog(ProfileModel &profileModel, MeshDialogModel *model, QWidget* parent = 0);
    ~MeshDialog();
    void accept();

    /**
     * @brief controlPoints Returns the list of boundary points that are control points.
     * @return List of control points
     */
    std::vector<BoundaryPoint*> controlPoints();

    /**
     * @brief controlPoints Returns the full list of boundary points.
     * @return List of boundary points
     */
    std::vector<BoundaryPoint*> boundaryPoints();

public slots:
    void runMesher();
    void setProfile();

    /**
     * @brief scaleZoomBarValue Scale the value of the zoom bar slider in response to changed in the view transformation.
     * Included validation check which places a limit on hte
     * @param value Scaling value
     */
    void scaleZoomBarValue(double value);

private slots:
    void on_profile_currentIndexChanged(int index);

    void on_density_currentIndexChanged(int index);

    void on_layers_valueChanged(int arg1);

    void on_thickness_valueChanged(double arg1);

    void on_growthFactor_valueChanged(double arg1);

    void on_loadProfileButton_clicked();

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
