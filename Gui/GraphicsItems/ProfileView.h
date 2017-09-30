#ifndef PROFILEGRAPHICSITEM_H
#define PROFILEGRAPHICSITEM_H

#include "CustomTypes.h"
#include <QGraphicsItem>
#include "ViewScaler.h"

class ProfileView : public QGraphicsItem
{
public:
    ProfileView(ViewScaler *scale, QGraphicsItem *parent = 0);

    void setProfilePoints(ProfilePoints profilePoints);
    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;
    void setDrawDots(bool drawDots);

private:
    ProfilePoints mProfile;
    ViewScaler* mScale;
    int mDrawDots = true;
};

#endif // PROFILEGRAPHICSITEM_H
