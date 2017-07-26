#include "ConfigureProject.h"
#include "ui_ConfigureProject.h"

ConfigureProject::ConfigureProject(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ConfigureProject)
{
    ui->setupUi(this);
}

ConfigureProject::~ConfigureProject()
{
    delete ui;
}
