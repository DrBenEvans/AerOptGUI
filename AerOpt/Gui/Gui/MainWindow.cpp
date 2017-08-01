#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "ConfigSimulationDialog.h"
#include "OptimisationRun.h"
#include "DebugOutput.h"
#include "MeshDialog.h"
#include "Mesh.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    mProfiles.emplace_back("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA0024.prf");
    mProfiles.emplace_back("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA21120.prf");
    mProfiles.emplace_back("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/test.prf");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionNewOptimisation_triggered() {
    mOptimisations.emplace_back();
    OptimisationRun& optimisation = mOptimisations.back();
    ConfigSimulationDialog diag(optimisation,&mProfiles,this);
    diag.exec();
    diag.show();
    ui->canvas->setProfile(*optimisation.getProfileObj());
    ui->canvas->update();
}

void MainWindow::on_actionRunMesher_triggered() {
    QString workDir("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Project");

    const std::list<std::pair<float,float>> profilePoints;
    Mesh mesh(profilePoints);
    MeshDialog diag(mesh, this);
    diag.show();
    diag.exec();
    mesh.runMesher(this,workDir);
}

void MainWindow::on_actionShowLog_triggered() {
    DebugOutput& debugOutput = DebugOutput::Instance();
    debugOutput.show();
}
