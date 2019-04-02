#include "ControlPointBoundingBox.h"
#include <QPainter>
#include "BoundaryPointModel.h"

ControlPointBoundingBox::ControlPointBoundingBox(BoundaryPointModel* model, int index, ViewScaler* scaler, QGraphicsItem *parent) :
    QGraphicsObject(parent),
    mActive(false),
    mBoundaryPointModel(model),
    mBoundaryPointIndex(index),
    mScale(scaler)
{
    mTopLeft = new ControlPointDragHandle(mBoundaryPointModel, mBoundaryPointIndex, BoundaryPointModel::TOPLEFT, mScale, this);
    mBottomRight = new ControlPointDragHandle(mBoundaryPointModel, mBoundaryPointIndex, BoundaryPointModel::BOTTOMRIGHT, mScale, this);
    connect(mBoundaryPointModel, &BoundaryPointModel::controlPointBoundsChanged, this, &ControlPointBoundingBox::boundsChanged);
}

void ControlPointBoundingBox::boundsChanged(int index) {
    if(index == mBoundaryPointIndex) {
        QRectF rect = mBoundaryPointModel->point(mBoundaryPointIndex)->controlPointRect();
        rect = mScale->toSceneScale(rect);

        mTopLeft->setPos(rect.topLeft());
        mBottomRight->setPos(rect.bottomRight());
        prepareGeometryChange();
        update();
    }
}

QRectF ControlPointBoundingBox::boundingRect() const {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->point(mBoundaryPointIndex);
    QRectF r = boundaryPoint->controlPointRect();
    r = mScale->toSceneScale(r);

    qreal margin = 2.0;
    r.adjust(-margin,-margin,margin,margin);
    return r;
}

void ControlPointBoundingBox::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    QRectF rect = mBoundaryPointModel->point(mBoundaryPointIndex)->controlPointRect();
    rect = mScale->toSceneScale(rect);

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
