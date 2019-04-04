#include "PlotConfigureDialog.h"
#include "ui_PlotConfigureDialog.h"

PlotConfigureDialog::PlotConfigureDialog(int lineSize, int pointSize, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PlotConfigureDialog)
{
    ui->setupUi(this);
    ui->lineSize->setValue(lineSize);
    ui->pointSize->setValue(pointSize);
}

PlotConfigureDialog::~PlotConfigureDialog()
{
    delete ui;
}

int PlotConfigureDialog::getLineSize() {
    return ui->lineSize->value();
}

int PlotConfigureDialog::getPointSize() {
    return ui->pointSize->value();
}
