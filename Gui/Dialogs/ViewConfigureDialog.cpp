#include "ViewConfigureDialog.h"
#include "ui_ViewConfigureDialog.h"

ViewConfigureDialog::ViewConfigureDialog(int lineSize, int pointSize, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ViewConfigureDialog)
{
    ui->setupUi(this);
    ui->lineSize->setValue(lineSize);
    ui->pointSize->setValue(pointSize);
}

ViewConfigureDialog::~ViewConfigureDialog()
{
    delete ui;
}

int ViewConfigureDialog::getLineSize() {
    return ui->lineSize->value();
}

int ViewConfigureDialog::getPointSize() {
    return ui->pointSize->value();
}
