#ifndef ZOOMPANVIEW_H
#define ZOOMPANVIEW_H

#include <QGraphicsView>

class ZoomPanView : public QGraphicsView
{
public:
    ZoomPanView(QWidget* parent = 0);

    void wheelEvent(QWheelEvent* event);
};

#endif // ZOOMPANVIEW_H
