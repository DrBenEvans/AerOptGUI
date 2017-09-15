#ifndef CONTROLPOINTBOUNDINGBOXVIEW_H
#define CONTROLPOINTBOUNDINGBOXVIEW_H

#include "ControlPointDragHandle.h"
#include "BoundaryPointModel.h"
#include <QGraphicsItem>

class ControlPointDragHandle;

class ControlPointBoundingBox : public QGraphicsObject {
    Q_OBJECT
public:
    ControlPointBoundingBox(BoundaryPointModel *model, int index, QGraphicsItem *parent);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

    void topLeftMoved(QPointF pos);
    void bottomRightMoved(QPointF pos);
    void setActivated(bool active);

public slots:
    void boundsChanged();

private:
    ControlPointDragHandle* mTopLeft;
    ControlPointDragHandle* mBottomRight;
    bool mActive;
    BoundaryPointModel* mBoundaryPointModel;
    int mBoundaryPointIndex;
};

#endif // CONTROLPOINTBOUNDINGBOXVIEW_H
