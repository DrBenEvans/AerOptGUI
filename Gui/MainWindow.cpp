#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "ConfigSimulationDialog.h"
#include "Optimisation.h"
#include "DebugOutput.h"
#include "MeshDialog.h"
#include "MeshView.h"
#include "Mesh.h"
//#include "MeshGraphicsItem.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    mOptimisationModel(nullptr),
    mCurrentMeshViewModel(new MeshModel(this))
{
    ui->setupUi(this);
    ui->logText->setReadOnly(true);
    ui->logText->setLineWrapMode(QTextEdit::NoWrap);
    connect(ui->newOptimisationButton, &QPushButton::clicked, this, &MainWindow::newOptimisation);
    connect(ui->actionNewOptimisation, &QAction::triggered, this, &MainWindow::newOptimisation);
    connect(ui->fitnessPlot, &PlotterWidget::selectedPointChanged, this, &MainWindow::onSelectedPointChange);

    connect(ui->agentSpinBox, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::setSelectedPointFromSpinBox);

    MeshView* meshView = new MeshView(new ViewScaler());
    meshView->setMeshModel(mCurrentMeshViewModel);
    QGraphicsScene* scene = new QGraphicsScene(this);
    scene->addItem(meshView);
    ui->graphicsView->setScene(scene);
}

void MainWindow::onSelectedPointChange(int iGen, int agent) {
    ui->genSpinBox->setValue(iGen+1);
    ui->agentSpinBox->setValue(agent+1);
    setMeshViewSimulation(iGen, agent);
}

void MainWindow::setSelectedPointFromSpinBox() {
    int iGen = ui->genSpinBox->value()-1;
    int agent = ui->agentSpinBox->value()-1;
    ui->fitnessPlot->setCurrentlySelectedPoint(iGen, agent);
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

void MainWindow::setMeshViewSimulation(int iGen, int agent) {
    Optimisation* optimisation = mOptimisationModel->optimisation(mCurrentOptimisationIndex).get();
    Mesh* mesh = optimisation->initMesh().get();
    mCurrentMeshViewModel->setCurrentMesh(mesh);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::newOptimisation() {
    auto optimisation = std::make_shared<Optimisation>();
    MeshDialog meshDialog(mProfileModel, this);
    if(meshDialog.exec() == QDialog::Accepted) {
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
