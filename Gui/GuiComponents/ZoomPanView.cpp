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

void ZoomPanView::zoomIn() {
    scale(ZOOM_SCALE_FACTOR, ZOOM_SCALE_FACTOR);
}

void ZoomPanView::zoomOut() {
    scale(1.0 / ZOOM_SCALE_FACTOR, 1.0 / ZOOM_SCALE_FACTOR);
}
