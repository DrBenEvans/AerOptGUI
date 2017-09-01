#include "ControlPointBoundingBox.h"
#include <QPainter>

ControlPointBoundingBox::ControlPointBoundingBox(QGraphicsItem *parent) :
    QGraphicsItem(parent),
    mControlPointRect(-10,-10,20,20)
{
    mTopLeft = new ControlPointDragHandle(mControlPointRect.topLeft(), true, this);
    mTopLeft = new ControlPointDragHandle(mControlPointRect.bottomRight(), false, this);

}

void ControlPointBoundingBox::topLeftMoved(QPointF pos) {
    mControlPointRect.setTopLeft(pos);
    prepareGeometryChange();
}

void ControlPointBoundingBox::bottomRightMoved(QPointF pos) {
    mControlPointRect.setBottomRight(pos);
    prepareGeometryChange();
}

QRectF ControlPointBoundingBox::boundingRect() const {
    return mControlPointRect;
}

void ControlPointBoundingBox::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    painter->setBrush(Qt::NoBrush);
    painter->setPen(QPen(Qt::gray, 1, Qt::DashLine));
    painter->drawRect(mControlPointRect);
}
