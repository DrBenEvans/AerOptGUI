#include "ControlPointHandleGraphicsItem.h"
#include <QPainter>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>

ControlPointHandleGraphicsItem::ControlPointHandleGraphicsItem(QPointF point, bool top_left, QGraphicsItem* parent):
    QGraphicsItem(parent),
    mTopLeft(top_left)
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
    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::NoDrag);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);
}

void ControlPointHandleGraphicsItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    foreach(auto& view, scene()->views()) {
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

            return pos;


        } else {
            if(pos.x() < 0) {
                pos.setX(0);
            }

            if(pos.y() < 0) {
                pos.setY(0);
            }

            return pos;
        }
    }
    return QGraphicsItem::itemChange(change, value);
}
