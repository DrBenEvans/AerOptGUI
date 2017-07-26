#ifndef OPTIMISATIONMANAGERDIALOG_H
#define OPTIMISATIONMANAGERDIALOG_H

#include <QDialog>

namespace Ui {
class OptimisationManagerDialog;
}

class OptimisationManagerDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OptimisationManagerDialog(QWidget *parent = 0);
    ~OptimisationManagerDialog();

private:
    Ui::OptimisationManagerDialog *ui;
};

#endif // OPTIMISATIONMANAGERDIALOG_H
