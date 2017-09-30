#include "ProfileView.h"
#include <QtWidgets>
#include <QPainter>

ProfileView::ProfileView(ViewScaler* scale, QGraphicsItem* parent) :
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

void ProfileView::setDrawDots(bool drawDots) {
    mDrawDots = drawDots;
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
        x = mScale->w(p.first);
        y = mScale->h(p.second);

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
    QColor color = Qt::red;
    QPen pen(color, 1, Qt::SolidLine);
    painter->setBrush(Qt::NoBrush);
    painter->setPen(pen);
    QRect rect;

    QPainterPath path;

    auto& firstPoint = mProfile.front();

    path.moveTo( QPoint(mScale->w(firstPoint.first), mScale->h(firstPoint.second)));
    for (auto& p : mProfile)
    {
        path.lineTo( QPoint(mScale->w(p.first), mScale->h(p.second)) );
    }
    path.lineTo( QPoint(mScale->w(firstPoint.first), mScale->h(firstPoint.second)));

    painter->drawPath(path);

    // Draw the dots
    if(mDrawDots) {
        for (auto& p : mProfile)
        {
            path = QPainterPath();
            QPoint center = QPoint(mScale->w(p.first), mScale->h(p.second));
            qreal radius = 3.0;
            path.addEllipse(center,radius,radius);
            painter->fillPath(path, color);
            painter->drawPath(path);
        }
    }
}
