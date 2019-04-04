#ifndef PLOTCONFIGUREDIALOG_H
#define PLOTCONFIGUREDIALOG_H

#include <QDialog>

namespace Ui {
class PlotConfigureDialog;
}

class PlotConfigureDialog : public QDialog
{
    Q_OBJECT

public:
    explicit PlotConfigureDialog(QWidget *parent = nullptr);
    ~PlotConfigureDialog();

private:
    Ui::PlotConfigureDialog *ui;
};

#endif // PLOTCONFIGUREDIALOG_H
