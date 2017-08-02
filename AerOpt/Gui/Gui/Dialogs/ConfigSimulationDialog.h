#ifndef ConfigSimulationDialog_H
#define ConfigSimulationDialog_H

#include <QDialog>
#include "OptimisationRun.h"
#include "ProfileLibrary.h"

namespace Ui {
class ConfigSimulationDialog;
}

class ConfigSimulationDialog : public QDialog
{
    Q_OBJECT


public:
    explicit ConfigSimulationDialog(OptimisationRun& data, ProfileLibrary& profileLibrary, QWidget *parent = 0);
    ~ConfigSimulationDialog();

public slots:
    void accept();
    void reject();

private slots:
    void on_optmethod_currentIndexChanged(int index);

    void on_pushButton_clicked();

private:
    int objFuncEnumToIndex(Enum::ObjFunc enumeration);
    Enum::ObjFunc indexToObjFuncEnum(int index);

    int optMethodEnumToIndex(Enum::OptMethod enumeration);
    Enum::OptMethod indexToOptMethodEnum(int index);

    Ui::ConfigSimulationDialog *ui;
    OptimisationRun& mData;
    ProfileLibrary& mProfileLibrary;
};

#endif // ConfigSimulationDialog_H
