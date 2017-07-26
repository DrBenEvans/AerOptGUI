#include "OptimisationManagerDialog.h"
#include "ui_OptimisationManagerDialog.h"

OptimisationManagerDialog::OptimisationManagerDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OptimisationManagerDialog)
{
    ui->setupUi(this);
}

OptimisationManagerDialog::~OptimisationManagerDialog()
{
    delete ui;
}
