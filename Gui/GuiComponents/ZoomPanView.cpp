#include "ZoomPanView.h"
#include <QWheelEvent>
#include "QDebug"
#include <QTouchEvent>

ZoomPanView::ZoomPanView(QWidget* parent) :
    QGraphicsView(parent)
{
    setDragMode(QGraphicsView::ScrollHandDrag);
    setInteractive(true);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setTransformationAnchor(QGraphicsView::AnchorViewCenter);
    setAttribute(Qt::WA_AcceptTouchEvents);
}

void ZoomPanView::touchEvent(QTouchEvent *touchEvent){

    QList<QTouchEvent::TouchPoint> touchPoints = touchEvent->touchPoints();
    if (touchPoints.count() == 2) {
        // determine scale factor
        const QTouchEvent::TouchPoint &touchPoint0 = touchPoints.first();
        const QTouchEvent::TouchPoint &touchPoint1 = touchPoints.last();
        qreal currentScaleFactor =
                QLineF(touchPoint0.pos(), touchPoint1.pos()).length()
                / QLineF(touchPoint0.startPos(), touchPoint1.startPos()).length();
        zoom(currentScaleFactor);
    }

}

void ZoomPanView::wheelEvent(QWheelEvent* event) {

    //Scale the view ie. do the zoom
    qInfo() << event->delta();
    if(event->delta() > 0) {
        // Zoom in
        zoom(ZOOM_SCALE_FACTOR);
    } else {
        //Zoom out
        zoom(1.0 / ZOOM_SCALE_FACTOR);
    }
}

void ZoomPanView::zoom(double scaleFactor){
    scale(scaleFactor, scaleFactor);
}
