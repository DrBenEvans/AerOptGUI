#include "ZoomPanView.h"
#include <QWheelEvent>
#include "QDebug"
#include <QEvent>
#include <QGestureEvent>
#include <QGraphicsItem>

ZoomPanView::ZoomPanView(QWidget* parent) :
    QGraphicsView(parent)
{
    setDragMode(QGraphicsView::ScrollHandDrag);
    setInteractive(true);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setTransformationAnchor(QGraphicsView::AnchorViewCenter);

    grabGestures();
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
    emit scaleChanged(scaleFactor);
}

void ZoomPanView::grabGestures()
{
    QList<Qt::GestureType> gestures;
    gestures << Qt::PinchGesture;
    foreach (Qt::GestureType gesture, gestures)
        grabGesture(gesture);
}

bool ZoomPanView::event(QEvent *event)
{
    if (event->type() == QEvent::Gesture){
        return gestureEvent(static_cast<QGestureEvent*>(event));
    }
    return QWidget::event(event);
}

bool ZoomPanView::gestureEvent(QGestureEvent *event)
{
    if (QGesture *pinch = event->gesture(Qt::PinchGesture))
        pinchTriggered(static_cast<QPinchGesture *>(pinch));
    return true;
}

void ZoomPanView::pinchTriggered(QPinchGesture *gesture)
{
    QPinchGesture::ChangeFlags changeFlags = gesture->changeFlags();

    if (changeFlags & QPinchGesture::ScaleFactorChanged) {
        double pinchScaleFactor = gesture->property("scaleFactor").toDouble();
        zoom(pinchScaleFactor);
        update();
    }

}
