#include "Arrow.h"
#include <QtWidgets>
#include <QGraphicsLineItem>
#include <cmath>

const qreal Pi = 3.14;

Arrow::Arrow(QPointF startItem, QPointF endItem, QGraphicsItem *parent, QGraphicsScene *scene) : QGraphicsLineItem(parent)
{
	myStartItem = startItem;
	myEndItem = endItem;
	setFlag(QGraphicsItem::ItemIsSelectable, true);
	myColor = Qt::black;
	setPen(QPen(myColor, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
}

QRectF Arrow::boundingRect() const
{
	qreal extra = ( pen().width() + 20 ) / 2.0;

	return QRectF( line().p1(),
				   QSizeF( line().p2().x()-line().p1().x(), line().p2().y()-line().p1().y() )
				  ).normalized().adjusted(-extra, -extra, extra, extra);
}

QPainterPath Arrow::shape() const
{
	QPainterPath path = QGraphicsLineItem::shape();
	path.addPolygon(arrowHead);
	return path;
}

void Arrow::paint(QPainter *painter)
{
	QPen myPen = pen();
	myPen.setColor(myColor);
	qreal arrowSize = 5;
	painter->setPen(myPen);
	painter->setBrush(myColor);

//	setLine(QLineF(myStartItem, myEndItem));
	setLine(QLineF(myEndItem, myStartItem));

	double angle = ::acos(line().dx() / line().length());
	if (line().dy() >= 0)
		angle = (Pi * 2) - angle;

	QPointF arrowP1 = line().p1() + QPointF(sin(angle + Pi / 3) * arrowSize,
											cos(angle + Pi / 3) * arrowSize);
	QPointF arrowP2 = line().p1() + QPointF(sin(angle + Pi - Pi / 3) * arrowSize,
											cos(angle + Pi - Pi / 3) * arrowSize);

	arrowHead.clear();
	arrowHead << line().p1() << arrowP1 << arrowP2;
	painter->drawLine(line());
	painter->drawPolygon(arrowHead);
	if (isSelected()) {
		painter->setPen(QPen(myColor, 1, Qt::DashLine));
		QLineF myLine = line();
		myLine.translate(0, 4.0);
		painter->drawLine(myLine);
		myLine.translate(0,-8.0);
		painter->drawLine(myLine);
	}
}

