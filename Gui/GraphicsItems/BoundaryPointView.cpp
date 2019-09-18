#include <QtWidgets>
#include <QPainter>
#include <QPainter>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include "BoundaryPointView.h"

BoundaryPointView::BoundaryPointView(BoundaryPointModel *model, int index, ViewScaler* scaler, QGraphicsItem* parent) :
    QGraphicsObject(parent),
    mActive(false),
    mBoundaryPointModel(model),
    mBoundaryPointIndex(index),
    mScale(scaler)
{
    setZValue(1);

    setFlags(ItemIsSelectable | ItemIgnoresParentOpacity);
    setAcceptHoverEvents(true);

    qreal r = radius();
    mPointRect = QRectF(-r,-r,2*r,2*r);

    // set position
    BoundaryPoint* bp = mBoundaryPointModel->point(mBoundaryPointIndex);
    setPos(mScale->w(bp->x()), mScale->h(bp->y()));

    connect(mBoundaryPointModel, &BoundaryPointModel::activeIndexChanged, this, &BoundaryPointView::setActivePoint);
    connect(mBoundaryPointModel, &BoundaryPointModel::controlPointStateChanged, this, &BoundaryPointView::refreshControlPointState);

    // Set up click and hold timer
    timer = new QTimer(this);
    timer->setSingleShot(true);
    timer->setInterval(300);
    connect(timer, SIGNAL(timeout()), this, SLOT(timerTimeout()));

}

qreal BoundaryPointView::radius() const {
    return 5.0;
}

QRectF BoundaryPointView::boundingRect() const
{
    qreal margin = 5;
    QRectF r = mPointRect;
    return r.adjusted(-margin, -margin, margin, margin);
}

void BoundaryPointView::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen;
    QBrush brush;
    if(activated()) {
        pen = QPen(Qt::blue, 2, Qt::SolidLine);
    } else {
        pen = QPen(Qt::black, 2, Qt::SolidLine);
    }
    pen.setJoinStyle(Qt::MiterJoin);

    if(mControl) {
        brush = QBrush(Qt::red);
    } else {
        brush = QBrush(Qt::NoBrush);
    }
    painter->setBrush(brush);
    painter->setPen(pen);

    QPainterPath path;
    path.addRect(mPointRect);
    painter->fillPath(path, brush);
    painter->drawPath(path);
}

bool BoundaryPointView::activated() const {
    return mActive;
}

void BoundaryPointView::setActivePoint(int index) {
    bool active = (index == mBoundaryPointIndex);

    // We only need to update a boundary point if the activated value is different
    if (mActive != active){
        mActive = active;
        foreach(auto& view, scene()->views()) {
            if(active) {
                view->setDragMode(QGraphicsView::NoDrag);
            } else {
                view->setDragMode(QGraphicsView::ScrollHandDrag);
            }
        }

        if(mControlPointHandles != nullptr) {
            mControlPointHandles->setActivated(active);
        }

        prepareGeometryChange();
        update();
    }

}


void BoundaryPointView::refreshControlPointState(int index, bool ctl) {
    if (index == mBoundaryPointIndex){

        // Only need to refresh if different
        if(mControl != ctl) {
            //build bounding box, unless it exists
            if(mControlPointHandles == nullptr && ctl)
                mControlPointHandles = new ControlPointBoundingBox(mBoundaryPointModel, mBoundaryPointIndex, mScale, this);

            mControlPointHandles->setVisible(ctl);

            mControl = ctl;
        }
    }

}

bool BoundaryPointView::isControlPoint() {
    return mControl;
}

void BoundaryPointView::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    mBoundaryPointModel->setActiveIndex(mBoundaryPointIndex);
    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::NoDrag);
        view->setCursor(Qt::ArrowCursor);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);
}

void BoundaryPointView::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::ScrollHandDrag);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverLeaveEvent(event);
}

void BoundaryPointView::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) {
    if(mPointRect.contains(event->pos())) {
        mBoundaryPointModel->setControlPointState(mBoundaryPointIndex, !mControl);
    }

    QGraphicsItem::mouseDoubleClickEvent(event);
}

void BoundaryPointView::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    timer->start();
    QGraphicsItem::mousePressEvent(event);
}

void BoundaryPointView::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    timer->stop();
    QGraphicsItem::mouseReleaseEvent(event);
}

void BoundaryPointView::timerTimeout(){
    mBoundaryPointModel->setControlPointState(mBoundaryPointIndex, !mControl);
    if (isControlPoint())
        timer->setInterval(500);
    else {
        timer->setInterval(300);
    }
}
