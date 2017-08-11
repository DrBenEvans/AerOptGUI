#include "ZoomPanView.h"
#include <QWheelEvent>
#include "QDebug"

ZoomPanView::ZoomPanView(QWidget* parent)
{
    setDragMode(QGraphicsView::ScrollHandDrag);
    setInteractive(false);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setDragMode(QGraphicsView::ScrollHandDrag);
    setTransformationAnchor(QGraphicsView::AnchorViewCenter);
}

void ZoomPanView::wheelEvent(QWheelEvent* event) {

    //Scale the view ie. do the zoom
    double scaleFactor = 1.15; //How fast we zoom
    qInfo() << event->delta();
    if(event->delta() > 0) {
        //Zoom in
        scale(scaleFactor, scaleFactor);
    } else {
        //Zooming out
        scale(1.0 / scaleFactor, 1.0 / scaleFactor);
    }
}
