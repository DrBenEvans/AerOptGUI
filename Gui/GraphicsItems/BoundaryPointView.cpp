#include <QtWidgets>
#include <QPainter>
#include <QPainter>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include "BoundaryPointView.h"


BoundaryPointView::BoundaryPointView(int scale, QGraphicsItem* parent) :
    QGraphicsItem(parent),
    mScale(scale),
    mActive(false)
{
    setZValue(1);

    setFlags(ItemIsSelectable | ItemIgnoresParentOpacity);
    setAcceptHoverEvents(true);

    qreal r = radius();
    mPointRect = QRectF(-r,-r,2*r,2*r);
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
    if(active()) {
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

bool BoundaryPointView::active() const {
    return mActive;
}

void BoundaryPointView::setActivePoint(bool active) {
    mActive = active;
    foreach(auto& view, scene()->views()) {
        if(active) {
            view->setDragMode(QGraphicsView::NoDrag);
        } else {
            view->setDragMode(QGraphicsView::ScrollHandDrag);
        }
    }

    if(mControlPointHandles != nullptr) {
        mControlPointHandles->setActivePoint(active);
    }

    prepareGeometryChange();

    update();
}

void BoundaryPointView::setControl(bool ctl) {
    //build control handles
    if(mControlPointHandles == nullptr && ctl) {
        mControlPointHandles = new ControlPointBoundingBox(this);
    }

    mControlPointHandles->setVisible(ctl);

    mControl = ctl;
    setActivePoint(ctl);
}

bool BoundaryPointView::control() {
    return mControl;
}

void BoundaryPointView::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    setActivePoint(true);
    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::NoDrag);
        view->setCursor(Qt::ArrowCursor);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);
}

void BoundaryPointView::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    setActivePoint(false);
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

void BoundaryPointView::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
    printf("%f, %f\n", event->pos());
}
