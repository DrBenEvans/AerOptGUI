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
    connect(ui->actionShowCurrentOptimisationFiles, &QAction::triggered, this, &MainWindow::revealFilesCurrentOptimisation);

    connect(ui->agentSpinBox, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::setSelectedPointFromSpinBox);

    ViewScaler* mScale = new ViewScaler();

    MeshView* meshView = new MeshView(mScale);
    meshView->setMeshModel(mCurrentMeshViewModel);

    mProfileView = new ProfileView(mScale);

    QGraphicsScene* scene = new QGraphicsScene(this);
    scene->addItem(meshView);
    scene->addItem(mProfileView);
    ui->graphicsView->setScene(scene);
}

void MainWindow::revealFilesCurrentOptimisation() {
    mOptimisationModel->revealFiles(mCurrentOptimisationIndex);
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
    connect(mOptimisationModel, &OptimisationModel::optimisationFitnessChanged, this, &MainWindow::optimisationFitnessChanged);
    connect(mOptimisationModel, &OptimisationModel::optimisationOutputChanged, this, &MainWindow::optimisationOutputChanged);
}

void MainWindow::optimisationOutputChanged(int index) {
    if(mCurrentOptimisationIndex == index) {
        setLogText();
    }
}

void MainWindow::setLogText() {
    ui->logText->setText(currentOptimisation()->outputText());

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

        // update fitness plot
        ui->fitnessPlot->setCurrentOptimisationIndex(index);

        // update log text
        setLogText();

        // set profile points
        Mesh* initMesh = currentOptimisation()->initMesh();
        ProfilePoints profilePoints = currentOptimisation()->initMesh()->profilePoints();
        mProfileView->setProfilePoints(profilePoints);
    }
}

Optimisation* MainWindow::currentOptimisation() {
    return mOptimisationModel->optimisation(mCurrentOptimisationIndex);
}

void MainWindow::setMeshViewSimulation(int iGen, int agent) {
    Mesh* mesh = currentOptimisation()->mesh(iGen, agent);
    if(mesh) mCurrentMeshViewModel->setCurrentMesh(mesh);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::newOptimisation() {
    auto optimisation = std::make_shared<Optimisation>();
    MeshDialogModel* meshDialogModel = new MeshDialogModel(this);
    meshDialogModel->setCurrentMesh(optimisation->initMesh());
    MeshDialog meshDialog(mProfileModel, meshDialogModel, this);
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
