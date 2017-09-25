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
    ui->logText->setReadOnly(true);
    connect(ui->newOptimisationButton,&QPushButton::clicked,this,&MainWindow::newOptimisation);
    connect(ui->actionNewOptimisation,&QAction::triggered,this,&MainWindow::newOptimisation);
}

void MainWindow::setOptimisationModel(OptimisationModel* optimisationModel) {
    mOptimisationModel = optimisationModel;
    ui->optimisationComboBox->setModel(mOptimisationModel);
    ui->fitnessPlot->setOptimisationModel(mOptimisationModel);
    connect(ui->optimisationComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::setCurrentOptimisationIndex);
    connect(mOptimisationModel, &OptimisationModel::optimisationDataChanged, this, &MainWindow::optimisationFitnessChanged);
    connect(mOptimisationModel, &OptimisationModel::optimisationOutputChanged, this, &MainWindow::optimisationOutputChanged);
}

void MainWindow::optimisationOutputChanged(int index) {
    if(mCurrentOptimisationIndex == index) {
        setLogTextFromIndex(index);
    }
}

void MainWindow::setLogTextFromIndex(int index) {
    auto opt = mOptimisationModel->optimisation(index);
    ui->logText->setText(opt->outputText());

    // scroll textbox to bottom of text
    QScrollBar *sb = ui->logText->verticalScrollBar();
    sb->setValue(sb->maximum());
}

void MainWindow::optimisationFitnessChanged(int index) {
    if(mCurrentOptimisationIndex == index) {
        ui->fitnessPlot->updatePlot();
    }
}

void MainWindow::setCurrentOptimisationIndex(int index) {
    if(mCurrentOptimisationIndex != index) {
        mCurrentOptimisationIndex = index;
        ui->fitnessPlot->setCurrentOptimisationIndex(index);
        setLogTextFromIndex(index);
    }
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

            int index = mOptimisationModel->rowCount() - 1;
            ui->optimisationComboBox->setCurrentIndex(index);
        }
    }
}

void MainWindow::on_actionShowLog_triggered() {
    DebugOutput& debugOutput = DebugOutput::Instance();
    debugOutput.show();
}
