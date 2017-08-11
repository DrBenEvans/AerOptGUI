#ifndef ConfigSimulationDialog_H
#define ConfigSimulationDialog_H

#include <QDialog>
#include "Optimisation.h"
#include "OptimisationModel.h"
#include "ProfileModel.h"

namespace Ui {
class ConfigSimulationDialog;
}

class ConfigSimulationDialog : public QDialog
{
    Q_OBJECT


public:
    explicit ConfigSimulationDialog(std::shared_ptr<Optimisation> optimisation, QWidget *parent = 0);
    ~ConfigSimulationDialog();

public slots:
    void reject();
    void accept();

private slots:
    void on_optmethod_currentIndexChanged(int index);

private:
    int objFuncEnumToIndex(Enum::ObjFunc enumeration);
    Enum::ObjFunc indexToObjFuncEnum(int index);

    int optMethodEnumToIndex(Enum::OptMethod enumeration);
    Enum::OptMethod indexToOptMethodEnum(int index);

    Ui::ConfigSimulationDialog *ui;
    std::shared_ptr<Optimisation> mData;
};

#endif // ConfigSimulationDialog_H
