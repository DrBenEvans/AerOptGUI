#ifndef CONTROLPOINTGRAPHICSITEM_H
#define CONTROLPOINTGRAPHICSITEM_H

#include <QGraphicsItem>

class ControlPointGraphicsItem : public QGraphicsItem
{
public:
    ControlPointGraphicsItem(int scale);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

private:
    bool active();
    void setActive(bool active);
    qreal radius() const;
    qreal w(qreal width);
    qreal h(qreal width);
    void setControl(bool ctl);
    bool control();

    int mScale;
    bool mActive;
    qreal mRadiusSizeInc;
    bool mControl = false;

protected:
    virtual void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    virtual void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    virtual void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;

};

#endif // CONTROLPOINTGRAPHICSITEM_H
