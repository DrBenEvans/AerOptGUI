#ifndef ZOOMPANVIEW_H
#define ZOOMPANVIEW_H

#include <QGraphicsView>
#include <QTouchEvent>

class ZoomPanView : public QGraphicsView
{
    public:
        ZoomPanView(QWidget* parent = 0);

        /**
         * @brief ZoomPanView::wheelEvent Zoom in or out depending on direction of mouse scroll
         * @param event Wheel scroll event
         */
        void wheelEvent(QWheelEvent* event);

        /**
         * @brief touchEvent Pinch-to-Zoom action.
         * @param touchEvent Touch screen event
         */
        void touchEvent(QTouchEvent* touchEvent);

        /**
         * @brief ZOOM_SCALE_FACTOR
         * Factor by which to rescale view
         */
        const double ZOOM_SCALE_FACTOR = 1.15;

    private:

        /**
         * @brief zoom Perform a zoom on the view by the specified scaling factor.
         * @param scaleFactor
         */
        void zoom(double scaleFactor);

};

#endif // ZOOMPANVIEW_H
