#ifndef CONTROLPOINTBOUNDINGBOXVIEW_H
#define CONTROLPOINTBOUNDINGBOXVIEW_H

#include "ControlPointDragHandle.h"
#include "BoundaryPointModel.h"
#include <QGraphicsItem>
#include "ViewScaler.h"

class ControlPointDragHandle;

class ControlPointBoundingBox : public QGraphicsObject {
    Q_OBJECT
public:
    ControlPointBoundingBox(BoundaryPointModel *model, int index, ViewScaler* scaler, QGraphicsItem *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

    void topLeftMoved(QPointF pos);
    void bottomRightMoved(QPointF pos);
    void setActivated(bool active);

public slots:
    void boundsChanged(int index);

private:
    ControlPointDragHandle* mTopLeft;
    ControlPointDragHandle* mBottomRight;
    bool mActive;
    BoundaryPointModel* mBoundaryPointModel;
    int mBoundaryPointIndex;
    ViewScaler* mScale;
};

#endif // CONTROLPOINTBOUNDINGBOXVIEW_H
