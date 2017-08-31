#include "BoundaryPointGraphicsItem.h"
#include <QtWidgets>
#include <QPainter>

BoundaryPointGraphicsItem::BoundaryPointGraphicsItem(int scale, QGraphicsItem* parent) :
    QGraphicsItem(parent),
    mScale(scale),
    mActive(false),
    mControlPointRect(-10,-10,20,20)
{
    setZValue(1);

    setFlags(ItemIsSelectable | ItemIgnoresParentOpacity | ItemIgnoresTransformations);
    setAcceptHoverEvents(true);
}

qreal BoundaryPointGraphicsItem::radius() const {
    return 5.0;
}

QRectF BoundaryPointGraphicsItem::controlPointRect() const {
    return mControlPointRect;
}

QRectF BoundaryPointGraphicsItem::boundingRect() const
{
    QRectF r;
    qreal margin = 5;
    if(active() && mControl) {
        r = controlPointRect();
    } else {
        r = pointRect();
    }
    r.adjust(-margin, -margin, margin, margin);
    return r;
}

QRectF BoundaryPointGraphicsItem::pointRect() const {
    qreal r = radius();

    return QRectF(-r,-r,2*r,2*r);
}

void BoundaryPointGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen;
    QBrush brush;
    if(active()) {
        pen = QPen(Qt::green, 2, Qt::SolidLine);
    } else {
        pen = QPen(Qt::blue, 1, Qt::SolidLine);
    }

    if(mControl) {
        brush = QBrush(Qt::red);
    } else {
        brush = QBrush(Qt::NoBrush);
    }
    painter->setBrush(brush);
    painter->setPen(pen);

    int r = radius();
    QPainterPath path;
    //path.addEllipse(QPoint(0,0),r,r);
    path.addRect(pointRect());
    painter->fillPath(path, brush);
    painter->drawPath(path);

    if(mControl) {
        painter->setBrush(Qt::NoBrush);
        painter->setPen(QPen(Qt::gray, 1, Qt::DashLine));
        painter->drawRect(controlPointRect());
    }
}

bool BoundaryPointGraphicsItem::active() const {
    return mActive;
}

void BoundaryPointGraphicsItem::setActive(bool active) {
    mActive = active;
    foreach(auto& view, scene()->views()) {
        if(active) {
            view->setDragMode(QGraphicsView::NoDrag);
        } else {
            view->setDragMode(QGraphicsView::ScrollHandDrag);
        }
    }

    prepareGeometryChange();

    update();
}

void BoundaryPointGraphicsItem::setControl(bool ctl) {
    //build control handles
    if(mTopLeftH == nullptr && ctl) {
        mTopLeftH = new ControlPointHandleGraphicsItem(mControlPointRect.topLeft(), true, this);
        mBottomRightH = new ControlPointHandleGraphicsItem(mControlPointRect.bottomRight(), false, this);
    }

    mTopLeftH->setVisible(ctl);
    mBottomRightH->setVisible(ctl);

    mControl = ctl;
    setActive(ctl);
}
void BoundaryPointGraphicsItem::topLeftMoved(QPointF pos) {
    mControlPointRect.setTopLeft(pos);
    prepareGeometryChange();
}

void BoundaryPointGraphicsItem::bottomRightMoved(QPointF pos) {
    mControlPointRect.setBottomRight(pos);
    prepareGeometryChange();
}

bool BoundaryPointGraphicsItem::control() {
    return mControl;
}

void BoundaryPointGraphicsItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    setActive(true);

    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);
}

void BoundaryPointGraphicsItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    setActive(false);

    // send the event to graphics scene and items
    QGraphicsItem::hoverLeaveEvent(event);
}

void BoundaryPointGraphicsItem::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) {
    setControl(!mControl);

    QGraphicsItem::mouseDoubleClickEvent(event);
}
