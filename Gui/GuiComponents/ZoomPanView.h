#ifndef ZOOMPANVIEW_H
#define ZOOMPANVIEW_H

#include <QGraphicsView>

class ZoomPanView : public QGraphicsView
{
    public:
        ZoomPanView(QWidget* parent = 0);

        void wheelEvent(QWheelEvent* event);

        const double ZOOM_SCALE_FACTOR = 1.15;

    private:
        void zoomIn();
        void zoomOut();
};

#endif // ZOOMPANVIEW_H
