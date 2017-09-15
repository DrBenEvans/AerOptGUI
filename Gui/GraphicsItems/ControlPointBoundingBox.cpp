#include "ControlPointBoundingBox.h"
#include <QPainter>

ControlPointBoundingBox::ControlPointBoundingBox(BoundaryPointModel* model, int index, QGraphicsItem *parent) :
    QGraphicsObject(parent),
    mActive(false),
    mBoundaryPointModel(model),
    mBoundaryPointIndex(index)
{
    mTopLeft = new ControlPointDragHandle(mBoundaryPointModel, mBoundaryPointIndex, true, this);
    mBottomRight = new ControlPointDragHandle(mBoundaryPointModel, mBoundaryPointIndex, false, this);
    connect(mBoundaryPointModel, &BoundaryPointModel::controlPointBoundsChanged, this, &ControlPointBoundingBox::boundsChanged);
}

void ControlPointBoundingBox::boundsChanged() {
    prepareGeometryChange();
    update();
}

QRectF ControlPointBoundingBox::boundingRect() const {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->point(mBoundaryPointIndex);
    QRectF r = boundaryPoint->controlPointRect();
    qreal margin = 2.0;
    r.adjust(-margin,-margin,margin,margin);
    return r;
}

void ControlPointBoundingBox::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    QRectF rect = mBoundaryPointModel->point(mBoundaryPointIndex)->controlPointRect();
    painter->setBrush(Qt::NoBrush);
    QPen pen;
    if(mActive) {
        pen = QPen(Qt::blue, 2, Qt::DashLine);
    } else {
        pen = QPen(Qt::gray, 2, Qt::DashLine);
    }
    pen.setJoinStyle(Qt::MiterJoin);
    painter->setPen(pen);
    painter->drawRect(rect);
}

void ControlPointBoundingBox::setActivated(bool active) {
    mActive = active;
    update();
}
