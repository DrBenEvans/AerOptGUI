#include "ProfileView.h"
#include <QtWidgets>
#include <QPainter>

ProfileView::ProfileView(int scale, QGraphicsItem* parent) :
    QGraphicsItem(parent),
    mScale(scale)
{
    setZValue(0);
    setPos(0,0);
}


void ProfileView::setProfilePoints(ProfilePoints profilePoints) {
    mProfile = profilePoints;
    prepareGeometryChange();
    QGraphicsItem::update();
}

QRectF ProfileView::boundingRect() const
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

void ProfileView::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
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
        QPoint center = QPoint(w(p.first), h(p.second));
        qreal radius = 3.0;
        path.addEllipse(center,radius,radius);
        painter->fillPath(path, Qt::blue);
        painter->drawPath(path);
    }
}

qreal ProfileView::w(qreal width) {
    return width*mScale;
}

qreal ProfileView::h(qreal height) {
    return height*mScale;
}
