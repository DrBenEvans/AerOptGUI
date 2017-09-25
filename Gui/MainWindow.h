#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Optimisation.h"
#include "ProfileModel.h"
#include "OptimisationModel.h"

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

private slots:
    void on_actionShowLog_triggered();

private:
    void setLogTextFromIndex(int index);
    Ui::MainWindow *ui;
    ProfileModel mProfileModel;
    OptimisationModel* mOptimisationModel;
    int mCurrentOptimisationIndex = -1;
};

#endif // MAINWINDOW_H
