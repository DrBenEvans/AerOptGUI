#include <QSettings>
#include <QMessageBox>
#include "ClusterLoginDialog.h"
#include "ui_ClusterLoginDialog.h"
#include "clusterManager.h"

ClusterLoginDialog::ClusterLoginDialog(Optimisation* opt, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ClusterLoginDialog),
    mData(opt)
{
    QSettings settings;
    QString maddress = settings.value("Cluster/Address").toString();
    QString maccount = settings.value("Cluster/Account").toString();
    QString mUsername = settings.value("Cluster/Username").toString();
    ui->setupUi(this);
    ui->address->setText(maddress);
    ui->account->setText(maccount);
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

    QString password = ui->password->text();
    QString username = ui->username->text();
    QString account = ui->account->text();
    QString address = ui->address->text();

    int failed = sshVerifyPassword(address, username, password);

    if( failed ){
        QMessageBox msgBox;
        msgBox.setText("Unable to login. Please check the password and try again.");
        msgBox.exec();
    } else {
        settings.setValue("Cluster/Username", ui->username->text());
        settings.setValue("Cluster/Account", ui->account->text());
        settings.setValue("Cluster/Address", ui->address->text());
        mData->mClusterPassword = ui->password->text();

        QDialog::accept();
    }
}
