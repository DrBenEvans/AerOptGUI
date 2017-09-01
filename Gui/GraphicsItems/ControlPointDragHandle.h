#ifndef CONTROLPOINTHANDLE_H
#define CONTROLPOINTHANDLE_H

#include "ControlPointBoundingBox.h"
#include <QGraphicsItem>

class ControlPointBoundingBox;

class ControlPointDragHandle : public QGraphicsItem
{
public:
    ControlPointDragHandle(QPointF point, bool top_left, ControlPointBoundingBox *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget */*widget*/) override;

private:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;

    bool mTopLeft;
    ControlPointBoundingBox* mParent;
    QRectF mHandleRect;
    QRectF mHoverRect;
};


#endif // CONTROLPOINTHANDLE_H
