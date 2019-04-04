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
    explicit PlotConfigureDialog(int lineSize, int pointSize, QWidget *parent = nullptr);
    ~PlotConfigureDialog();
    int getPointSize();
    int getLineSize();

private:
    Ui::PlotConfigureDialog *ui;
};

#endif // PLOTCONFIGUREDIALOG_H
