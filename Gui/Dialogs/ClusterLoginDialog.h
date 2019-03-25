#ifndef CLUSTERLOGINDIALOG_H
#define CLUSTERLOGINDIALOG_H

#include <QDialog>
#include "Optimisation.h"

namespace Ui {
class ClusterLoginDialog;
}

class ClusterLoginDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ClusterLoginDialog( Optimisation* opt, QWidget *parent = nullptr);
    ~ClusterLoginDialog();

private slots:
    void accept();

private:
    Ui::ClusterLoginDialog *ui;
    Optimisation * mData;
};

#endif // CLUSTERLOGINDIALOG_H
