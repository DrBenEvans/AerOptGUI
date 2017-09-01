#include "ZoomPanView.h"
#include <QWheelEvent>
#include "QDebug"

ZoomPanView::ZoomPanView(QWidget* parent) :
    QGraphicsView(parent)
{
    setDragMode(QGraphicsView::ScrollHandDrag);
    setInteractive(true);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setTransformationAnchor(QGraphicsView::AnchorViewCenter);
}

void ZoomPanView::wheelEvent(QWheelEvent* event) {
    //Scale the view ie. do the zoom
    qInfo() << event->delta();
    if(event->delta() > 0) {
        zoomIn();
    } else {
        zoomOut();
    }
}

double ZoomPanView::zoomScaleFactor() {
    return 1.15;
}

void ZoomPanView::zoomIn() {
    double scaleFactor = zoomScaleFactor();
    scale(scaleFactor, scaleFactor);
}

void ZoomPanView::zoomOut() {
    double scaleFactor = zoomScaleFactor();
    scale(1.0 / scaleFactor, 1.0 / scaleFactor);
}
