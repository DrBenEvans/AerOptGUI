#include <QDebug>
#include "ConfigSimulationDialog.h"
#include "ui_ConfigSimulationDialog.h"
#include <QFileDialog>

ConfigSimulationDialog::ConfigSimulationDialog(OptimisationRun& data, ProfileLibrary &profileLibrary, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ConfigSimulationDialog),
    mData(data),
    mProfileLibrary(profileLibrary)
{
    ui->setupUi(this);

    // label setup
    ui->label->setText(mData.getLabel());

    // profiles QComboBox setup
    ui->profile->setModel(&profileLibrary);

    // optimiser setup
    ui->objfunc->setCurrentIndex(this->objFuncEnumToIndex(mData.objFunc()));
    ui->optmethod->setCurrentIndex(this->optMethodEnumToIndex(mData.getOptimisationMethod()));
    ui->nnests->setValue(mData.noAgents());
    ui->ngens->setValue(mData.noGens());
    ui->ndiscard->setValue(mData.getNoTop());

    // flow setup
    ui->temp->setValue(mData.freeTemp());
    ui->pressure->setValue(mData.freePress());
    ui->angle->setValue(mData.freeAlpha());
    ui->reynolds->setValue(mData.reNo());
    ui->mach->setValue(mData.machNo());

}

ConfigSimulationDialog::~ConfigSimulationDialog()
{
    delete ui;
}

void ConfigSimulationDialog::accept()
{
    // optimiser setup
    mData.setNoAgents(ui->nnests->value());
    mData.setNoGens(ui->ngens->value());
    mData.setNoTop(ui->ndiscard->value());
    mData.setOptimisationMethod(this->indexToOptMethodEnum(ui->optmethod->currentIndex()));
    mData.setObjFunc(this->indexToObjFuncEnum(ui->objfunc->currentIndex()));

    // flow setup
    mData.setMachNo(ui->mach->value());
    mData.setReNo(ui->reynolds->value());
    mData.setFreeAlpha(ui->angle->value());
    mData.setFreePress(ui->pressure->value());
    mData.setFreeTemp(ui->temp->value());

    // profile setup
    QSharedPointer<Profile> profile = ui->profile->currentData(Qt::UserRole).value<QSharedPointer<Profile>>();
    mData.setProfile(profile);

    QDialog::accept();
}

int ConfigSimulationDialog::objFuncEnumToIndex(Enum::ObjFunc enumeration) {
    switch(enumeration) {
    case Enum::ObjFunc::LIFTDRAG:
        return 0;
    case Enum::ObjFunc::MAXLIFT:
        return 1;
    case Enum::ObjFunc::MINDRAG:
        return 2;
    case Enum::ObjFunc::MAXDOWNFORCE:
        return 3;
    case Enum::ObjFunc::MINLIFT:
        return 4;
    case Enum::ObjFunc::FUNCNOTSET:
        qCritical() << "Error: unknown lift function";
        return 998;
    default:
        qCritical() << "Error: unknown lift enumeration";
    }

    return -1;
}

Enum::ObjFunc ConfigSimulationDialog::indexToObjFuncEnum(int index) {
    switch(index) {
    case 0:
        return Enum::ObjFunc::LIFTDRAG;
    case 1:
        return Enum::ObjFunc::MAXLIFT;
    case 2:
        return Enum::ObjFunc::MINDRAG;
    case 3:
        return Enum::ObjFunc::MAXDOWNFORCE;
    case 4:
        return Enum::ObjFunc::MINLIFT;
    default:
        qCritical() << "Error: unknown lift function";
    }

    return Enum::ObjFunc::FUNCNOTSET;

}

int ConfigSimulationDialog::optMethodEnumToIndex(Enum::OptMethod enumeration) {
    switch(enumeration) {
    case Enum::OptMethod::MCS:
        return 0;
    case Enum::OptMethod::DE:
        return 1;
    case Enum::OptMethod::PSO:
        return 2;
    default:
        qCritical() << "Error: unknown optimisation method";
        return 999;
    }

}

Enum::OptMethod ConfigSimulationDialog::indexToOptMethodEnum(int index) {
    switch(index) {
    case 0:
        return Enum::OptMethod::MCS;
    case 1:
        return Enum::OptMethod::DE;
    case 2:
        return Enum::OptMethod::PSO;
    default:
        qCritical() << "Error: unknown optimisation method";
        return Enum::OptMethod::METHODNOTSET;
    }
}


void ConfigSimulationDialog::on_optmethod_currentIndexChanged(int index)
{
    bool isMCS = this->indexToOptMethodEnum(index) == Enum::OptMethod::MCS;

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

void ConfigSimulationDialog::on_pushButton_clicked()
{
    QString fileName;
    QStringList fileNames;

#ifdef Q_OS_UNIX
    // unix load file
    fileNames.append(
                QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath()+"/Documents/Projects/AerOptProject/", "Profile Files (*.prf)")
                );
#endif
#ifdef Q_OS_WIN32
    // do windows stuff here
    fileNames.append(
                QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath(), "Profile Files (*.prf)")
                );
#endif

    if (fileNames.size() > 0)
    {
        for (const QString &f: fileNames)
        {
            fileName = f;
        }

        qDebug() << "File selected: " << fileName;

        mProfileLibrary.addProfileFromFilePath(fileName);
    }
    else
    {
        //Set text Not OK here!
        qWarning() << "Profile data not imported!";
    }
}
