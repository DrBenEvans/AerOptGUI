#ifndef CONTROLPOINTGRAPHICSITEM_H
#define CONTROLPOINTGRAPHICSITEM_H

#include <QGraphicsItem>
#include "ControlPointBoundingBox.h"

class BoundaryPointView : public QGraphicsItem
{
public:
    BoundaryPointView(int scale, QGraphicsItem *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

private:
    bool active() const;
    void setActivePoint(bool active);
    bool control();
    void setControl(bool ctl);
    qreal radius() const;
    qreal w(qreal width);
    qreal h(qreal width);

    QRectF mPointRect;
    int mScale;
    bool mActive;
    bool mControl = false;
    ControlPointBoundingBox* mControlPointHandles = nullptr;

protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
};

#endif // CONTROLPOINTGRAPHICSITEM_H
