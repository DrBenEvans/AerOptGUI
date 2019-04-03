#ifndef CONTROLPOINTHANDLE_H
#define CONTROLPOINTHANDLE_H

#include "BoundaryPointModel.h"
#include <QGraphicsItem>
#include "ViewScaler.h"

class ControlPointBoundingBox;

class ControlPointDragHandle : public QGraphicsItem
{
public:
    ControlPointDragHandle(BoundaryPointModel* model, int index, BoundaryPointModel::CornerPosition corner, ViewScaler* scaler, QGraphicsItem* parent = 0);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget */*widget*/) override;

private:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;

    BoundaryPointModel::CornerPosition corner;
    QRectF mHandleRect;
    QRectF mHoverRect;
    BoundaryPointModel* mBoundaryPointModel;
    int mBoundaryPointIndex;
    BoundaryPoint* mBoundaryPoint;
    ViewScaler* mScale;
};


#endif // CONTROLPOINTHANDLE_H
