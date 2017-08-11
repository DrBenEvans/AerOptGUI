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

private slots:
    void on_actionShowLog_triggered();

private:
    Ui::MainWindow *ui;
    ProfileModel mProfileModel;
    OptimisationModel* mOptimisationModel;
};

#endif // MAINWINDOW_H
