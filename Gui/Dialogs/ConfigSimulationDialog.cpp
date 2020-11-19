#include <QDebug>
#include <QMessageBox>
#include "ConfigSimulationDialog.h"
#include "ui_ConfigSimulationDialog.h"
#include "MeshDialog.h"
#include "Enumerations.h"
#include "ClusterLoginDialog.h"

ConfigSimulationDialog::ConfigSimulationDialog(Optimisation *optimisation, QWidget *parent, bool enabled) :
    QDialog(parent),
    ui(new Ui::ConfigSimulationDialog),
    mData(optimisation),
    mEnabled(enabled)
{

    ui->setupUi(this);

    // label setup
    ui->label->setText(mData->label());

    // optimiser setup
    ui->objfunc->setCurrentIndex(mData->objFunc());
    ui->optmethod->setCurrentIndex(mData->getOptimisationMethod());
    ui->nnests->setValue(mData->noAgents());
    ui->ngens->setValue(mData->noGens());
    ui->ndiscard->setValue(mData->getNoTop());

    // flow setup
    ui->temp->setValue(mData->freeTemp());
    ui->pressure->setValue(mData->freePress());
    ui->angle->setValue(mData->freeAlpha());
    ui->reynolds->setValue(mData->reNo());
    ui->mach->setValue(mData->machNo());

    ui->optimisationParameterBox->setEnabled(mEnabled);
    ui->runCondBox->setEnabled(mEnabled);
    ui->projectOptionBox->setEnabled(mEnabled);
    ui->runOnCluster->setEnabled(true);
}

ConfigSimulationDialog::~ConfigSimulationDialog()
{
    delete ui;
}

void ConfigSimulationDialog::validationError(QString message) {
    QMessageBox::warning(this, "Warning", message, QMessageBox::Ok);
}

void ConfigSimulationDialog::accept()
{
    if(mEnabled) {

        // set label
        QString label = ui->label->text();
        if(label.size() == 0) {
            validationError("An optimisation name must be specified");
            return;
        }
        mData->setLabel(label);

        QDir outputDataDirectory = QDir(mData->simulationDirectoryPath());

        if(outputDataDirectory.exists()) {
            validationError("Directory already exists:" + outputDataDirectory.absolutePath());
            return;
        }


        // set optimiser
        mData->setNoAgents(ui->nnests->value());
        mData->setNoGens(ui->ngens->value());
        mData->setNoTop(ui->ndiscard->value());
        mData->setOptimisationMethod(static_cast<Enum::OptMethod>(ui->optmethod->currentIndex()));
        mData->setObjFunc(static_cast<Enum::ObjFunc>(ui->objfunc->currentIndex()));


        // set flow
        mData->setMachNo(ui->mach->value());
        mData->setReNo(ui->reynolds->value());
        mData->setFreeAlpha(ui->angle->value());
        mData->setFreePress(ui->pressure->value());
        mData->setFreeTemp(ui->temp->value());

        mData->runOnCluster = ui->runOnCluster->checkState();

        if(mData->runOnCluster && mData->mClusterPassword.isEmpty()) {

            ClusterLoginDialog loginDiag(mData, this);
            if( loginDiag.exec() == QDialog::Accepted){
                QDialog::accept();
            } else {
                // Return without dismissing the window
                return;
            }

        }

    }

    QDialog::accept();

}



void ConfigSimulationDialog::on_optmethod_currentIndexChanged(int index)
{
    //This line isn't too effective as if the enum indexing changes it breaks
    bool isMCS = (static_cast<Enum::OptMethod>(index) == Enum::OptMethod::MCS);

    if(isMCS) {
        ui->ndiscard->setEnabled(true);
        ui->ndiscard_label->setEnabled(true);
    } else {
        ui->ndiscard->setEnabled(false);
        ui->ndiscard_label->setEnabled(false);
    }

}

void ConfigSimulationDialog::reject()
{
    QDialog::reject();
    //Do Nothing!
}
