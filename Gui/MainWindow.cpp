#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "ConfigSimulationDialog.h"
#include "Optimisation.h"
#include "DebugOutput.h"
#include "MeshDialog.h"
#include "Mesh.h"
//#include "MeshGraphicsItem.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    mOptimisationModel(nullptr)
{
    ui->setupUi(this);
    connect(ui->newOptimisationButton,&QPushButton::clicked,this,&MainWindow::newOptimisation);
    connect(ui->actionNewOptimisation,&QAction::triggered,this,&MainWindow::newOptimisation);
}

void MainWindow::setOptimisationModel(OptimisationModel* optimisationModel) {
    mOptimisationModel = optimisationModel;
    ui->optimisationComboBox->setModel(mOptimisationModel);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::newOptimisation() {
    auto optimisation = std::make_shared<Optimisation>();
    MeshDialog meshDialog(mProfileModel, this);
    if(meshDialog.exec() == QDialog::Accepted) {
        optimisation->setLabel(QString("Optimisation %1").arg(mOptimisationModel->rowCount()+1));

        ConfigSimulationDialog diag(optimisation, this);
        if(diag.exec() == QDialog::Accepted) {
            optimisation->setControlPoints(meshDialog.controlPoints());
            mOptimisationModel->addOptimisation(optimisation);
            mOptimisationModel->run(optimisation);
        }
    }
}

void MainWindow::on_actionShowLog_triggered() {
    DebugOutput& debugOutput = DebugOutput::Instance();
    debugOutput.show();
}
