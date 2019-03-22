#include <QSettings>
#include "ClusterLoginDialog.h"
#include "ui_ClusterLoginDialog.h"

ClusterLoginDialog::ClusterLoginDialog(Optimisation* opt, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ClusterLoginDialog),
    mData(opt)
{
    QSettings settings;
    QString mUsername = settings.value("Cluster/Username").toString();
    ui->setupUi(this);
    ui->username->setText(mUsername);
    ui->password->setText(mData->mClusterPassword);
}

ClusterLoginDialog::~ClusterLoginDialog()
{
    delete ui;
}

void ClusterLoginDialog::accept()
{
    QSettings settings;

    mData->mClusterPassword = ui->password->text();
    settings.setValue("Cluster/Username", ui->username->text());

    QDialog::accept();
}
