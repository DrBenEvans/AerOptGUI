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
}

qreal BoundaryPointView::radius() const {
    return 5.0;
}

QRectF BoundaryPointView::boundingRect() const
{
    qreal margin = 5;
    QRectF r = mPointRect;
    r.adjust(-margin, -margin, margin, margin);
    return r;
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

    int r = radius();
    QPainterPath path;
    path.addRect(mPointRect);
    painter->fillPath(path, brush);
    painter->drawPath(path);
}

bool BoundaryPointView::activated() const {
    return mActive;
}

void BoundaryPointView::setActivePoint(int index) {
    setActivated(index == mBoundaryPointIndex);
}

void BoundaryPointView::setActivated(bool active) {
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

void BoundaryPointView::setControl(bool ctl) {
    //build control handles
    if(mControlPointHandles == nullptr && ctl) {
        mControlPointHandles = new ControlPointBoundingBox(mBoundaryPointModel, mBoundaryPointIndex, this);
    }

    mControlPointHandles->setVisible(ctl);

    mControl = ctl;
    mBoundaryPointModel->setActiveIndex(mBoundaryPointIndex);
    mBoundaryPointModel->setControlPointState(mBoundaryPointIndex, ctl);
}

bool BoundaryPointView::control() {
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
        setControl(!mControl);
    }

    QGraphicsItem::mouseDoubleClickEvent(event);
}
