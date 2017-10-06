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

    ViewScaler* mScale = new ViewScaler();

    MeshView* meshView = new MeshView(mScale);
    meshView->setMeshModel(mCurrentMeshViewModel);

    mProfileView = new ProfileView(mScale);
    mProfileView->setDrawDots(false);

    QGraphicsScene* scene = new QGraphicsScene(this);
    scene->addItem(meshView);
    scene->addItem(mProfileView);
    ui->graphicsView->setScene(scene);
}

void MainWindow::revealFilesCurrentOptimisation() {
    mOptimisationModel->revealFiles(mCurrentOptimisationIndex);
}

void MainWindow::onSelectedPointChange(int iGen, int agent) {
    QString text = QString("Gen=%1, Index=%2, Fitness=%3");
    text = text.arg(iGen+1).arg(agent+1);
    Optimisation* opt = currentOptimisation();
    if(opt) text = text.arg(opt->fitness(iGen, agent));
    ui->plotText->setText(text);

    setMeshViewSimulation(iGen, agent);
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
        Optimisation* opt = currentOptimisation();
        if(opt) {
            ProfilePoints profilePoints = opt->initProfilePoints();
            mProfileView->setProfilePoints(profilePoints);
        }

        clearCurrentSelection();
    }
}

void MainWindow::clearCurrentSelection() {
    // clear current selection
    ui->fitnessPlot->clearCurrentSelection();
    mCurrentMeshViewModel->setCurrentMesh(nullptr);
    ui->plotText->setText("");
}

Optimisation* MainWindow::currentOptimisation() {
    return mOptimisationModel->optimisation(mCurrentOptimisationIndex);
}

void MainWindow::setMeshViewSimulation(int iGen, int agent) {
    Optimisation* optimisation = currentOptimisation();
    Mesh* mesh = optimisation->mesh(iGen, agent);

    if(!mesh) {
        // try to read again, in case file has been written since last read
        Mesh* mesh = optimisation->mesh(iGen, agent);
        if(!mesh) {
            QMessageBox::warning(this, "Warning", "Data could not be loaded for selected point. Check that file exists and is formatted correctly.", QMessageBox::Ok);
        }
    }
    mCurrentMeshViewModel->setCurrentMesh(mesh);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::newOptimisation() {
    Optimisation *opt = new Optimisation();
    MeshDialogModel *meshDialogModel = new MeshDialogModel(this);
    meshDialogModel->setCurrentMesh(opt->initMesh());
    MeshDialog meshDialog(mProfileModel, meshDialogModel, this);
    if(meshDialog.exec() == QDialog::Accepted) {
        ConfigSimulationDialog diag(opt, this);
        if(diag.exec() == QDialog::Accepted) {
            opt->setControlPoints(meshDialog.controlPoints());
            bool success = mOptimisationModel->run(opt);

            if(!success) {
                QString message("Optimisation %1 run failed! Check log for details.");
                message = message.arg(opt->label());
                QMessageBox::warning(this, "Warning", message, QMessageBox::Ok);
            } else {
                mOptimisationModel->addOptimisation(opt);
            }

            int index = mOptimisationModel->rowCount() - 1;
            ui->optimisationComboBox->setCurrentIndex(index);
        }
    }
}

void MainWindow::on_actionShowLog_triggered() {
    DebugOutput& debugOutput = DebugOutput::Instance();
    debugOutput.show();
}

void MainWindow::on_loadOptimisationButton_clicked()
{
    QString fileName;
    QStringList fileNames;

    QSettings settings;
    QString startLoadPath = settings.value("AerOpt/startLoadPath").toString();

    fileNames.append(QFileDialog::getOpenFileName(this, "Select AerOpt Input File",
                                                  startLoadPath, "AerOpt Input Files (AerOpt_InputParameters.txt)"));
    if (fileNames.size() > 0 && fileNames[0].size() > 0)
    {
        for (const QString &f: fileNames)
        {
            fileName = f;
        }

        qDebug() << "Loading Optimisation: " << fileName;

        QDir fileDir(fileName);
        fileDir.cdUp();
        settings.value("AerOpt/startLoadPath", fileDir.absolutePath());

        QModelIndex index = mOptimisationModel->loadByInputFilePath(fileName);
        if(!mOptimisationModel->isIndexValid(index)) {
            QMessageBox::warning(this, "Warning", "Optimisation load failed.", QMessageBox::Ok);
        } else {
            // setting though ui->optimisationComboBox allows the current index to be maintained
            ui->optimisationComboBox->setCurrentIndex(index.row());
        }
    }
    else
    {
        qWarning() << "Not optimisation selected for load!";
    }
}
