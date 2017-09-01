#ifndef CONTROLPOINTBOUNDINGBOXVIEW_H
#define CONTROLPOINTBOUNDINGBOXVIEW_H

#include "ControlPointDragHandle.h"
#include <QGraphicsItem>

class ControlPointDragHandle;

class ControlPointBoundingBox : public QGraphicsItem {
public:
    ControlPointBoundingBox(QGraphicsItem *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

    void topLeftMoved(QPointF pos);
    void bottomRightMoved(QPointF pos);
private:
    QRectF mControlPointRect;
    ControlPointDragHandle* mTopLeft;
    ControlPointDragHandle* mTopRight;
};

#endif // CONTROLPOINTBOUNDINGBOXVIEW_H
