#include "ControlPointBoundingBox.h"
#include <QPainter>

ControlPointBoundingBox::ControlPointBoundingBox(QGraphicsItem *parent) :
    QGraphicsItem(parent),
    mControlPointRect(-10,-10,20,20),
    mActive(false)
{
    mTopLeft = new ControlPointDragHandle(mControlPointRect.topLeft(), true, this);
    mBottomRight = new ControlPointDragHandle(mControlPointRect.bottomRight(), false, this);

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
    QRectF r = mControlPointRect;
    qreal margin = 2.0;
    r.adjust(-margin,-margin,margin,margin);
    return r;
}

void ControlPointBoundingBox::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    painter->setBrush(Qt::NoBrush);
    QPen pen;
    if(mActive) {
        pen = QPen(Qt::blue, 2, Qt::DashLine);
    } else {
        pen = QPen(Qt::gray, 2, Qt::DashLine);
    }
    pen.setJoinStyle(Qt::MiterJoin);
    painter->setPen(pen);
    painter->drawRect(mControlPointRect);
}

void ControlPointBoundingBox::setActivePoint(bool active) {
    mActive = active;
    update();
}
