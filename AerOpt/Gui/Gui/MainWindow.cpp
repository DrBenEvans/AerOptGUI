#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "ConfigSimulationDialog.h"
#include "OptimisationRun.h"
#include "DebugOutput.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionNewOptimisation_triggered() {
    OptimisationRun simulation = OptimisationRun();
    ConfigSimulationDialog diag(simulation,this);
    diag.exec();
    diag.show();
    ui->canvas->setData(simulation);
    ui->canvas->update();
}

void MainWindow::on_actionShowLog_triggered() {
    DebugOutput& debugOutput = DebugOutput::Instance();
    debugOutput.show();
}
