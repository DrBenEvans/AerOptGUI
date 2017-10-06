#include "ViewScaler.h"

ViewScaler::ViewScaler(QObject* parent, int scale) : QObject(parent),
    mScale(scale)
{

}

qreal ViewScaler::w(qreal width) {
    return width*mScale;
}

qreal ViewScaler::h(qreal height) {
    return height*mScale;
}

QRectF ViewScaler::rect(int a1, int a2, int a3, int a4) {
    return QRectF(a1*mScale, a2*mScale, a3*mScale, a4*mScale);
}

QRectF ViewScaler::toSceneScale(QRectF rect) {
    QRectF scaled = QRectF(rect.topLeft()*mScale, rect.bottomRight()*mScale);
    return scaled;
}

QPointF ViewScaler::fromSceneScale(QPointF point) {
    return point/mScale;
}
