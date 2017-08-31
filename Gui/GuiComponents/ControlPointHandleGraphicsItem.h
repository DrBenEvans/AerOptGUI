#ifndef CONTROLPOINTHANDLEGRAPHICSITEM_H
#define CONTROLPOINTHANDLEGRAPHICSITEM_H

#include <QGraphicsItem>

class ControlPointHandleGraphicsItem : public QGraphicsItem
{
public:
    ControlPointHandleGraphicsItem(QPointF point, bool top_left, QGraphicsItem *parent = 0);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

private:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;
    bool mTopLeft;

};

#endif // CONTROLPOINTHANDLEGRAPHICSITEM_H
