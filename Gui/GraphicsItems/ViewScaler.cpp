#include "ViewScaler.h"

ViewScaler::ViewScaler(int scale) :
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
