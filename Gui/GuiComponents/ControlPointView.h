#ifndef CONTROLPOINTVIEW_H
#define CONTROLPOINTVIEW_H

#include <QWidget>

namespace Ui {
class ControlPointView;
}

class ControlPointView : public QWidget
{
    Q_OBJECT

public:
    explicit ControlPointView(QWidget *parent = 0);
    ~ControlPointView();

private:
    Ui::ControlPointView *ui;
};

#endif // CONTROLPOINTVIEW_H
