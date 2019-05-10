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
    explicit ConfigSimulationDialog(Optimisation *optimisation, QWidget *parent = 0, bool enabled = true);
    ~ConfigSimulationDialog();

public slots:
    void reject();
    void accept();

private slots:
    void on_optmethod_currentIndexChanged(int index);

private:
    void validationError(QString message);


    Ui::ConfigSimulationDialog *ui;
    Optimisation *mData;
    bool mEnabled;
};

#endif // ConfigSimulationDialog_H
