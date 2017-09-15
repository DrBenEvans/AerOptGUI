#ifndef CONTROLPOINTVIEW_H
#define CONTROLPOINTVIEW_H

#include <QWidget>
#include "BoundaryPointModel.h"

namespace Ui {
class ControlPointView;
}

class ControlPointView : public QWidget
{
    Q_OBJECT

public:
    explicit ControlPointView(QWidget *parent = 0);
    ~ControlPointView();

    void updateViewData();
    void setModel(BoundaryPointModel *boundaryPointModel);

public slots:
    void activePointChanged(int index);
    void smoothingValueChanged(double value);
    void controlPointStateChanged(int index, bool isControlPoint);
    void controlPointBoundsChanged(int index);
    void updateModelControlPointState(bool isControlPoint);
    void controlBoundaryChanged();

private:
    void setPointCoords(qreal x, qreal y);
    void controlPointParamsVisible(bool visible);
    Ui::ControlPointView *ui;
    BoundaryPointModel* mBoundaryPointModel;
};

#endif // CONTROLPOINTVIEW_H
