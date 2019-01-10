#ifndef ZOOMPANVIEW_H
#define ZOOMPANVIEW_H

#include <QGraphicsView>

class ZoomPanView : public QGraphicsView
{
public:
    ZoomPanView(QWidget* parent = 0);

    void wheelEvent(QWheelEvent* event);

private:
    double zoomScaleFactor();
    void zoomIn();
    void zoomOut();
};

#endif // ZOOMPANVIEW_H
