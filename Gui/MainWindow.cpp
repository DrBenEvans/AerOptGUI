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

    ProfileModel mProfileModel();
    on_actionNewOptimisation_triggered();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionNewOptimisation_triggered() {
    mOptimisations.emplace_back();
    OptimisationRun& optimisation = mOptimisations.back();
    ConfigSimulationDialog diag(optimisation,mProfileModel,this);
    diag.exec();
    diag.show();
    optimisation.finishConfigure();
    ui->canvas->setProfile(optimisation.getProfileObj());
    ui->canvas->update();
}

void MainWindow::on_actionRunMesher_triggered() {
    QString workDir("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Project");

    ProfilePoints profilePoints;
    OptimisationRun& optimisation = getCurrentOptimisation();
    QSharedPointer<Mesh> mesh = optimisation.getMesh();
    MeshDialog diag(mesh, this);
    diag.show();
    diag.exec();
    ui->canvas->setMesh(mesh);
    connect(mesh.data(),SIGNAL(meshUpdate()),ui->canvas,SLOT(update()));
    mesh->runMesher(workDir);
}

void MainWindow::on_actionShowLog_triggered() {
    DebugOutput& debugOutput = DebugOutput::Instance();
    debugOutput.show();
}

OptimisationRun& MainWindow::getCurrentOptimisation() {
    return mOptimisations.back();
}
