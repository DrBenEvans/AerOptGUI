#include "ControlPointDragHandle.h"
#include <QPainter>
#include <QGraphicsView>
#include "BoundaryPointModel.h"

ControlPointDragHandle::ControlPointDragHandle(BoundaryPointModel* model, int index, BoundaryPointModel::CornerPosition corner, ViewScaler* scaler, QGraphicsItem* parent) :
    QGraphicsItem(parent),
    corner(corner),
    mBoundaryPointModel(model),
    mBoundaryPointIndex(index),
    mBoundaryPoint(mBoundaryPointModel->point(mBoundaryPointIndex)),
    mScale(scaler)
{
    setFlags(ItemIsMovable | ItemSendsScenePositionChanges);
    setZValue(1);
    setAcceptHoverEvents(true);

    QRectF rect = mBoundaryPoint->controlPointRect();
    rect = mScale->toSceneScale(rect);

    switch(corner){

        case BoundaryPointModel::TOPLEFT:
            setPos(rect.topLeft());
            break;

        case BoundaryPointModel::BOTTOMRIGHT:
            setPos(rect.bottomRight());
            break;
    }

    // build handle rect
    qreal size = 4.0;
    mHandleRect = QRectF(-size,-size,size*2,size*2);

    // build hover rect
    size = 7.0;
    mHoverRect = QRectF(-size,-size,size*2,size*2);
}

void ControlPointDragHandle::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget */*widget*/)
{

    QPen pen = QPen(Qt::black, 1, Qt::SolidLine);
    QBrush brush = QBrush(Qt::white);
    pen.setJoinStyle(Qt::MiterJoin);
    painter->setBrush(brush);
    painter->setPen(pen);

    painter->drawRect(mHandleRect);
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
        view->setCursor(Qt::ArrowCursor);
    }

    // send the event to graphics scene and items
    QGraphicsItem::hoverLeaveEvent(event);
}

QVariant ControlPointDragHandle::itemChange(GraphicsItemChange change, const QVariant &value)
{
    if (change == ItemPositionChange && scene()) {
        // value is the new position.
        QPointF pos = value.toPointF();

        switch(corner){
            case BoundaryPointModel::TOPLEFT:
                if(pos.x() > 0) {
                    pos.setX(0);
                }

                if(pos.y() > 0) {
                    pos.setY(0);
                }
                break;

            case BoundaryPointModel::BOTTOMRIGHT:
                if(pos.x() < 0) {
                    pos.setX(0);
                }

                if(pos.y() < 0) {
                    pos.setY(0);
                }
                break;
        }

        mBoundaryPointModel->setControlBoundaryCorner(mBoundaryPointIndex, mScale->fromSceneScale(pos), corner);
        return pos;

    }
    return QGraphicsItem::itemChange(change, value);
}
