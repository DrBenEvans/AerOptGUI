#ifndef VIEWSCALER_H
#define VIEWSCALER_H

#include <QRectF>
#include <QObject>

class ViewScaler : public QObject
{
public:
    ViewScaler(QObject* parent, int scale = 1000);
    qreal w(qreal width);
    qreal h(qreal width);
    QRectF rect(int a1, int a2, int a3, int a4);
    QRectF toSceneScale(QRectF rect);
    QPointF fromSceneScale(QPointF point);

private:
    int mScale;
};

#endif // VIEWSCALER_H
