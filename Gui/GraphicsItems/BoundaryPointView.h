#ifndef CONTROLPOINTGRAPHICSITEM_H
#define CONTROLPOINTGRAPHICSITEM_H

#include <QGraphicsItem>
#include "ControlPointBoundingBox.h"
#include "BoundaryPoint.h"

class BoundaryPointView : public QGraphicsObject
{
public:
    BoundaryPointView(BoundaryPoint& boundaryPoint, QGraphicsItem *parent = 0);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

private:
    bool active() const;
    void setActivePoint(bool active);
    bool control();
    void setControl(bool ctl);
    qreal radius() const;
    qreal w(qreal width);
    qreal h(qreal width);

    QRectF mPointRect;
    bool mActive;
    bool mControl = false;
    ControlPointBoundingBox* mControlPointHandles = nullptr;
    BoundaryPoint& mBoundaryPoint;

protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
};

#endif // CONTROLPOINTGRAPHICSITEM_H
