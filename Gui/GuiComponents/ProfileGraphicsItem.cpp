#include "ProfileGraphicsItem.h"
#include <QtWidgets>
#include <QPainter>

ProfileGraphicsItem::ProfileGraphicsItem(int scale) :
    mScale(scale)
{
    setZValue(0);
    setPos(QPointF(-(mScale/2), 0));

    setFlags(ItemIsSelectable);
    setAcceptHoverEvents(true);
}


void ProfileGraphicsItem::setProfilePoints(ProfilePoints profilePoints) {
    mProfile = profilePoints;
    prepareGeometryChange();
    QGraphicsItem::update();
}

QRectF ProfileGraphicsItem::boundingRect() const
{
    float x, y;
    float xmax =  -std::numeric_limits<float>::infinity();
    float ymax =  -std::numeric_limits<float>::infinity();
    float xmin =  std::numeric_limits<float>::infinity();
    float ymin =  std::numeric_limits<float>::infinity();

    for (auto& p : mProfile)
    {
        x = p.first*mScale;
        y = p.second*mScale;

        if(x>xmax)
            xmax = x;
        if(x<xmin)
            xmin = x;
        if(y>ymax)
            ymax = y;
        if(y<ymin)
            ymin = y;

    }
    int margin=3;
    return QRectF(xmin-margin, ymin-margin, (xmax-xmin)+2*margin, (ymax-ymin)+2*margin);
}

void ProfileGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen(Qt::blue, 1, Qt::SolidLine);
    painter->setBrush(Qt::NoBrush);
    painter->setPen(pen);
    QRect rect;

    QPainterPath path;

    auto& firstPoint = mProfile.front();

    path.moveTo( QPoint(w(firstPoint.first), h(firstPoint.second)));
    for (auto& p : mProfile)
    {
        path.lineTo( QPoint(w(p.first), h(p.second)) );
    }
    path.lineTo( QPoint(w(firstPoint.first), h(firstPoint.second)));

    painter->drawPath(path);

    // Draw the dots
    for (auto& p : mProfile)
    {
        path = QPainterPath();
        rect.moveCenter( QPoint(w(p.first), h(p.second)) );
        rect.setWidth( 3 );
        rect.setHeight( 3 );
        path.addRoundedRect(rect,100,100);
        painter->fillPath(path, Qt::blue);
        painter->drawPath(path);
    }
}

qreal ProfileGraphicsItem::w(qreal width) {
    return width*mScale;
}

qreal ProfileGraphicsItem::h(qreal height) {
    return height*mScale;
}
