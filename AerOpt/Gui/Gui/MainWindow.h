#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "OptimisationRun.h"

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
    Ui::MainWindow *ui;
    std::vector<OptimisationRun> mOptimisations;
    std::vector<Profile> mProfiles;
};

#endif // MAINWINDOW_H
