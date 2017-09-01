#ifndef CONTROLPOINTBOUNDINGBOXVIEW_H
#define CONTROLPOINTBOUNDINGBOXVIEW_H

#include "ControlPointDragHandle.h"
#include "BoundaryPoint.h"
#include <QGraphicsItem>

class ControlPointDragHandle;

class ControlPointBoundingBox : public QGraphicsObject {
    Q_OBJECT
public:
    ControlPointBoundingBox(std::shared_ptr<BoundaryPoint> bp, QGraphicsItem *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

    void topLeftMoved(QPointF pos);
    void bottomRightMoved(QPointF pos);
    void setActivePoint(bool active);
public slots:
    void controlRectChanged();

private:
    QRectF controlPointRect() const;
    ControlPointDragHandle* mTopLeft;
    ControlPointDragHandle* mBottomRight;
    bool mActive;
    std::shared_ptr<BoundaryPoint> mBoundaryPoint;
};

#endif // CONTROLPOINTBOUNDINGBOXVIEW_H
