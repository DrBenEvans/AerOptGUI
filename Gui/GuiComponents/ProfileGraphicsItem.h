#ifndef PROFILEGRAPHICSITEM_H
#define PROFILEGRAPHICSITEM_H

#include "CustomTypes.h"
#include <QGraphicsItem>

class ProfileGraphicsItem : public QGraphicsItem
{
public:
    ProfileGraphicsItem(int scale);

    void setProfilePoints(ProfilePoints profilePoints);
    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

private:
    qreal w(qreal width);
    qreal h(qreal width);

    ProfilePoints mProfile;
    int mScale;
};

#endif // PROFILEGRAPHICSITEM_H
