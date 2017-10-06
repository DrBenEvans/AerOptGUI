#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Optimisation.h"
#include "ProfileModel.h"
#include "MeshModel.h"
#include "OptimisationModel.h"
#include "ProfileView.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setOptimisationModel(OptimisationModel* optimisationModel);

public slots:
    void newOptimisation();
    void setCurrentOptimisationIndex(int index);
    void optimisationFitnessChanged(int index);
    void optimisationOutputChanged(int index);
    void onSelectedPointChange(int iGen, int agent);

private slots:
    void on_actionShowLog_triggered();
    void revealFilesCurrentOptimisation();
    void clearCurrentSelection();

    void on_loadOptimisationButton_clicked();

private:
    void setMeshViewSimulation(int iGen, int agent);
    void setLogText();
    Optimisation* currentOptimisation();
    Ui::MainWindow *ui;
    ProfileModel mProfileModel;
    OptimisationModel* mOptimisationModel;
    int mCurrentOptimisationIndex = -1;
    MeshModel* mCurrentMeshViewModel;
    ProfileView* mProfileView = 0;
};

#endif // MAINWINDOW_H
