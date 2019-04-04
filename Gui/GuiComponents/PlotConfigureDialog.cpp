#include "PlotConfigureDialog.h"
#include "ui_PlotConfigureDialog.h"

PlotConfigureDialog::PlotConfigureDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PlotConfigureDialog)
{
    ui->setupUi(this);
}

PlotConfigureDialog::~PlotConfigureDialog()
{
    delete ui;
}
