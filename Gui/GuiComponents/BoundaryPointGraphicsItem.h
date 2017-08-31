#ifndef CONTROLPOINTGRAPHICSITEM_H
#define CONTROLPOINTGRAPHICSITEM_H

#include <QGraphicsItem>
#include "ControlPointHandleGraphicsItem.h"

class BoundaryPointGraphicsItem : public QGraphicsItem
{
public:
    BoundaryPointGraphicsItem(int scale, QGraphicsItem *parent = 0);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

private:
    bool active() const;
    void setActive(bool active);
    qreal radius() const;
    qreal w(qreal width);
    qreal h(qreal width);
    void setControl(bool ctl);
    bool control();
    QRectF controlPointRect() const;
    QRectF pointRect() const;

    int mScale;
    bool mActive;
    bool mControl = false;
    ControlPointHandleGraphicsItem* mTopLeftH = nullptr;
    ControlPointHandleGraphicsItem* mBottomRightH = nullptr;
    std::vector<ControlPointHandleGraphicsItem> mHandles;
    QRectF mControlPointRect;

protected slots:
    void topLeftMoved(QPointF pos);
    void bottomRightMoved(QPointF pos);

protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;

};

#endif // CONTROLPOINTGRAPHICSITEM_H
