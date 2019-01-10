#ifndef ZOOMPANVIEW_H
#define ZOOMPANVIEW_H

#include <QGraphicsView>

class ZoomPanView : public QGraphicsView
{
    public:
        ZoomPanView(QWidget* parent = 0);

        /**
         * @brief ZoomPanView::wheelEvent
         * Zoom in or out depending on direction of mouse scroll
         * @param event
         */
        void wheelEvent(QWheelEvent* event);

        /**
         * @brief ZOOM_SCALE_FACTOR
         * Factor by which to rescale view
         */
        const double ZOOM_SCALE_FACTOR = 1.15;

    private:
        /**
         * @brief zoomIn
         * Zoom in to view
         */
        void zoomIn();

        /**
         * @brief zoomOut
         * Zoom out of view
         */
        void zoomOut();
};

#endif // ZOOMPANVIEW_H
