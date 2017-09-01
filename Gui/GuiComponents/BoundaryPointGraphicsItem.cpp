#include <QtWidgets>
#include <QPainter>
#include <QPainter>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include "BoundaryPointGraphicsItem.h"


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

ControlPointHandleGraphicsItem::ControlPointHandleGraphicsItem(QPointF point, bool top_left, BoundaryPointGraphicsItem* parent):
    QGraphicsItem(parent),
    mTopLeft(top_left),
    mParent(parent)
{
    setFlags(ItemIsMovable | ItemSendsScenePositionChanges);
    setZValue(1);
    setPos(point);
    setAcceptHoverEvents(true);
}

void ControlPointHandleGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen;
    QBrush brush;

    pen = QPen(Qt::black, 1, Qt::SolidLine);
    brush = QBrush(Qt::red);
    painter->setBrush(brush);
    painter->setPen(pen);

    QPainterPath path;
    path.addRect(boundingRect());
    painter->fillPath(path, brush);
    painter->drawPath(path);
}

QRectF ControlPointHandleGraphicsItem::boundingRect() const
{
    qreal size = 3.0;
    QRectF r(-size,-size,size*2,size*2);
    return r;
}

void ControlPointHandleGraphicsItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);

    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::NoDrag);
        view->setCursor(Qt::SizeFDiagCursor);
    }
}

void ControlPointHandleGraphicsItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    foreach(auto& view, scene()->views()) {
        view->setCursor(Qt::ArrowCursor);
        view->setDragMode(QGraphicsView::ScrollHandDrag);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverLeaveEvent(event);
}

void ControlPointHandleGraphicsItem::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) {
    QGraphicsItem::mouseDoubleClickEvent(event);
}

QVariant ControlPointHandleGraphicsItem::itemChange(GraphicsItemChange change, const QVariant &value)
{
    if (change == ItemPositionChange && scene()) {
        // value is the new position.
        QPointF pos = value.toPointF();
        if(mTopLeft) {
            if(pos.x() > 0) {
                pos.setX(0);
            }

            if(pos.y() > 0) {
                pos.setY(0);
            }

            mParent->topLeftMoved(pos);

            return pos;


        } else {
            if(pos.x() < 0) {
                pos.setX(0);
            }

            if(pos.y() < 0) {
                pos.setY(0);
            }

            mParent->bottomRightMoved(pos);

            return pos;
        }
    }
    return QGraphicsItem::itemChange(change, value);
}
