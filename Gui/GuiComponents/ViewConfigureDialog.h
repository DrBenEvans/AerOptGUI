#ifndef PLOTCONFIGUREDIALOG_H
#define PLOTCONFIGUREDIALOG_H

#include <QDialog>

namespace Ui {
class ViewConfigureDialog;
}

class ViewConfigureDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ViewConfigureDialog(int lineSize, int pointSize, QWidget *parent = nullptr);
    ~ViewConfigureDialog();
    int getPointSize();
    int getLineSize();

private:
    Ui::ViewConfigureDialog *ui;
};

#endif // PLOTCONFIGUREDIALOG_H
