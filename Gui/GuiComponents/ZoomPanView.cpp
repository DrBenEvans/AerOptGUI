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

/**
  Zoom in or out depending on direction of mouse scroll
 * @brief ZoomPanView::wheelEvent
 * @param event
 */
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
    double scaleFactor = ZOOM_SCALE_FACTOR;
    scale(scaleFactor, scaleFactor);
}

void ZoomPanView::zoomOut() {
    double scaleFactor = ZOOM_SCALE_FACTOR;
    scale(1.0 / scaleFactor, 1.0 / scaleFactor);
}
