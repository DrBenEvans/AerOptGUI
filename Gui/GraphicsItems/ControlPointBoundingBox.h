#ifndef CONTROLPOINTBOUNDINGBOXVIEW_H
#define CONTROLPOINTBOUNDINGBOXVIEW_H

#include "ControlPointDragHandle.h"
#include "BoundaryPoint.h"
#include <QGraphicsItem>

class ControlPointDragHandle;

class ControlPointBoundingBox : public QGraphicsItem {
public:
    ControlPointBoundingBox(BoundaryPoint& bp, QGraphicsItem *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

    void topLeftMoved(QPointF pos);
    void bottomRightMoved(QPointF pos);
    void setActivePoint(bool active);
private:
    QRectF controlPointRect() const;
    ControlPointDragHandle* mTopLeft;
    ControlPointDragHandle* mBottomRight;
    bool mActive;
    BoundaryPoint& mBoundaryPoint;
};

#endif // CONTROLPOINTBOUNDINGBOXVIEW_H
