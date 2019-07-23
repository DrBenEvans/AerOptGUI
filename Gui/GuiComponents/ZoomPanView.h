#ifndef ZOOMPANVIEW_H
#define ZOOMPANVIEW_H

#include <QGraphicsView>
#include <QtGui>

class QGestureEvent;
class QPinchGesture;

class ZoomPanView : public QGraphicsView
{
    Q_OBJECT
    public:
        ZoomPanView(QWidget* parent = nullptr);

        /**
         * @brief ZoomPanView::wheelEvent Zoom in or out depending on direction of mouse scroll
         * @param event Wheel scroll event
         */
        void wheelEvent(QWheelEvent* event);

        /**
         * @brief ZOOM_SCALE_FACTOR
         * Factor by which to rescale view with mouse wheel gestures
         */
        const double ZOOM_SCALE_FACTOR = 1.15;

        /**
         * @brief grabGestures Determines which gestures to capture.
         * We define it to capture Pinch Gestures only.
         */
        void grabGestures();

    private:

        /**
         * @brief zoom Perform a zoom on the view by the specified scaling factor.
         * @param scaleFactor
         */
        void zoom(double scaleFactor);

        /**
         * @brief gestureEvent Handles different gestures by passing on to relevant subroutine.
         * @param event Gesture event
         * @return true
         */
        bool gestureEvent(QGestureEvent *event);

        /**
         * @brief pinchTriggered Handles pinch gestures to scale the view with zoom according to the size and direction of the pinch.
         */
        void pinchTriggered(QPinchGesture*);

    protected:

        /**
         * @brief event Listens out for events, specifically gestures.
         * @param event
         * @return event
         */
        bool event(QEvent *event) override;


    signals:
        void scaleChanged(double newScale);
};

#endif // ZOOMPANVIEW_H
