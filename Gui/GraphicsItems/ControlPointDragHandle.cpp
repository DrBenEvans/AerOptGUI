#include "ControlPointDragHandle.h"
#include <QPainter>
#include <QGraphicsView>

ControlPointDragHandle::ControlPointDragHandle(QPointF point, bool top_left, ControlPointBoundingBox *parent):
    QGraphicsItem(parent),
    mTopLeft(top_left),
    mParent(parent)
{
    setFlags(ItemIsMovable | ItemSendsScenePositionChanges);
    setZValue(1);
    setPos(point);
    setAcceptHoverEvents(true);

    // build handle rect
    qreal size = 4.0;
    mHandleRect = QRectF(-size,-size,size*2,size*2);

    // build hover rect
    size = 7.0;
    mHoverRect = QRectF(-size,-size,size*2,size*2);

}

void ControlPointDragHandle::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget */*widget*/)
{
    QPen pen;
    QBrush brush;

    pen = QPen(Qt::black, 1, Qt::SolidLine);
    brush = QBrush(Qt::red);
    painter->setBrush(brush);
    painter->setPen(pen);

    QPainterPath path;
    path.addRect(mHandleRect);
    painter->fillPath(path, brush);
    painter->drawPath(path);
}

QRectF ControlPointDragHandle::boundingRect() const
{
    return mHoverRect;
}

void ControlPointDragHandle::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    // send the event to graphics scene and items
    QGraphicsItem::hoverEnterEvent(event);

    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::NoDrag);
        view->setCursor(Qt::SizeFDiagCursor);
    }
}

void ControlPointDragHandle::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    foreach(auto& view, scene()->views()) {
        view->setDragMode(QGraphicsView::NoDrag);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverLeaveEvent(event);
}

void ControlPointDragHandle::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) {
    QGraphicsItem::mouseDoubleClickEvent(event);
}

QVariant ControlPointDragHandle::itemChange(GraphicsItemChange change, const QVariant &value)
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
