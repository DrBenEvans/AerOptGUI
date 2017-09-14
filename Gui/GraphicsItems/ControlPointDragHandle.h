#ifndef CONTROLPOINTHANDLE_H
#define CONTROLPOINTHANDLE_H

#include "BoundaryPointModel.h"
#include <QGraphicsItem>

class ControlPointBoundingBox;

class ControlPointDragHandle : public QGraphicsItem
{
public:
    ControlPointDragHandle(BoundaryPointModel* model, int index, bool isTopLeft, QGraphicsItem* parent = 0);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget */*widget*/) override;

private:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;

    bool mTopLeft;
    QRectF mHandleRect;
    QRectF mHoverRect;
    BoundaryPointModel* mBoundaryPointModel;
    int mBoundaryPointIndex;
    BoundaryPoint* mBoundaryPoint;
};


#endif // CONTROLPOINTHANDLE_H
