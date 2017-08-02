#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "OptimisationRun.h"
#include "ProfileLibrary.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_actionNewOptimisation_triggered();
    void on_actionShowLog_triggered();
    void on_actionRunMesher_triggered();

private:
    OptimisationRun& getCurrentOptimisation();

    Ui::MainWindow *ui;
    std::vector<OptimisationRun> mOptimisations;
    ProfileLibrary mProfileLibrary;
};

#endif // MAINWINDOW_H
