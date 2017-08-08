#ifndef ARROW_H
#define ARROW_H

#include <QGraphicsLineItem>

class QGraphicsPolygonItem;
class QGraphicsLineItem;
class QGraphicsScene;
class QRectF;
class QGraphicsSceneMouseEvent;
class QPainterPath;

class Arrow : public QGraphicsLineItem
{
public:
	enum { Type = UserType + 4 };

	Arrow(QPointF startItem, QPointF endItem, QGraphicsItem *parent = 0, QGraphicsScene *scene = 0);

	int type() const
	{
		return Type;
	}
	void paint(QPainter *painter);
	QRectF boundingRect() const;
	QPainterPath shape() const;
	void setColor(const QColor &color)
	{
		myColor = color;
	}

public slots:

protected:


private:
	QPointF myStartItem;
	QPointF myEndItem;
	QColor myColor;
	QPolygonF arrowHead;
};

#endif // ARROW_H
