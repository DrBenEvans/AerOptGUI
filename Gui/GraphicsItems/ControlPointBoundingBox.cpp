#include "ControlPointBoundingBox.h"
#include <QPainter>

ControlPointBoundingBox::ControlPointBoundingBox(std::shared_ptr<BoundaryPoint> bp, QGraphicsItem *parent) :
    QGraphicsObject(parent),
    mActive(false),
    mBoundaryPoint(bp)
{
    QRectF rect = controlPointRect();
    mTopLeft = new ControlPointDragHandle(rect.topLeft(), true, this);
    mBottomRight = new ControlPointDragHandle(rect.bottomRight(), false, this);

    // connect boundary point signals
    connect(mBoundaryPoint.get(), &BoundaryPoint::controlRectChanged, this, &ControlPointBoundingBox::controlRectChanged);

}

QRectF ControlPointBoundingBox::controlPointRect() const {
    return mBoundaryPoint->controlPointRect();
}

void ControlPointBoundingBox::controlRectChanged() {
    prepareGeometryChange();
    update();
}

void ControlPointBoundingBox::topLeftMoved(QPointF pos) {
    mBoundaryPoint->setTopLeftBound(pos);
}

void ControlPointBoundingBox::bottomRightMoved(QPointF pos) {
    mBoundaryPoint->setBottomRightBound(pos);
}

QRectF ControlPointBoundingBox::boundingRect() const {
    QRectF r = controlPointRect();
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
    painter->drawRect(controlPointRect());
}

void ControlPointBoundingBox::setActivePoint(bool active) {
    mActive = active;
    update();
}
