#include "ControlPointGraphicsItem.h"
#include <QtWidgets>
#include <QPainter>

ControlPointGraphicsItem::ControlPointGraphicsItem(int scale) :
    mScale(scale),
    mActive(false),
    mRadiusSizeInc(1.0)
{
    setZValue(0);
    setPos(QPointF(0,0));

    setFlags(ItemIsSelectable | ItemIgnoresParentOpacity | ItemIgnoresTransformations);
    setAcceptHoverEvents(true);
}

qreal ControlPointGraphicsItem::radius() const {
    qreal r = 5.0;
    if(mActive) {
        return r*mRadiusSizeInc;
    } else {
        return r;
    }
}


QRectF ControlPointGraphicsItem::boundingRect() const
{
    qreal r = radius()*mRadiusSizeInc;

    return QRectF(-r,-r,2*r,2*r);
}

void ControlPointGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen;
    QBrush brush;
    if(active()) {
        pen = QPen(Qt::red, 2, Qt::SolidLine);
        brush = QBrush(Qt::red);
    } else {
        pen = QPen(Qt::red, 1, Qt::SolidLine);
        brush = QBrush(Qt::NoBrush);
    }
    painter->setBrush(brush);
    painter->setPen(pen);

    int r = radius();
    QPainterPath path;
    path.addEllipse(QPoint(0,0),r,r);
    painter->fillPath(path,brush);
    painter->drawPath(path);
}

bool ControlPointGraphicsItem::active() {
    return mActive;
}

void ControlPointGraphicsItem::setActive(bool active) {
    mActive = active;
    foreach(auto& view, scene()->views()) {
        if(active) {
            view->setDragMode(QGraphicsView::NoDrag);
        } else {
            view->setDragMode(QGraphicsView::ScrollHandDrag);
        }
    }

    update();
}

void ControlPointGraphicsItem::setControl(bool ctl) {
    mControl = ctl;
    setActive(ctl);
}

bool ControlPointGraphicsItem::control() {
    return mControl;
}

void ControlPointGraphicsItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    setActive(true);

    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);
}

void ControlPointGraphicsItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    setActive(false);

    // send the event to graphics scene and items
    QGraphicsItem::hoverLeaveEvent(event);
}

void ControlPointGraphicsItem::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) {
    setControl(!mControl);
}
