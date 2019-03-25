#include <QSettings>
#include <QMessageBox>
#include "ClusterLoginDialog.h"
#include "ui_ClusterLoginDialog.h"
#include "ClusterFolderChecker.h"

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

    QString password = ui->password->text();
    QString username = ui->username->text();

    int failed = sshVerifyPassword(username, password);

    if( failed ){
        QMessageBox msgBox;
        msgBox.setText("Unable to login. Please check the password and try again.");
        msgBox.exec();
    } else {
        settings.setValue("Cluster/Username", ui->username->text());
        mData->mClusterPassword = ui->password->text();

        QDialog::accept();
    }
}
