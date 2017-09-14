#ifndef CONTROLPOINTGRAPHICSITEM_H
#define CONTROLPOINTGRAPHICSITEM_H

#include <QGraphicsItem>
#include "ControlPointBoundingBox.h"
#include "BoundaryPointModel.h"
#include "ViewScaler.h"

class BoundaryPointView : public QGraphicsObject
{
    Q_OBJECT
public:
    BoundaryPointView(BoundaryPointModel *model, int index, ViewScaler* scaler, QGraphicsItem *parent = 0);

    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

public slots:
    void setActivePoint(int index);

private:
    bool activated() const;
    void setActivated(bool active);
    bool control();
    void setControl(bool ctl);
    qreal radius() const;
    qreal w(qreal width);
    qreal h(qreal width);

    QRectF mPointRect;
    bool mActive;
    bool mControl = false;
    ControlPointBoundingBox* mControlPointHandles = nullptr;
    BoundaryPointModel* mBoundaryPointModel;
    int mBoundaryPointIndex;
    ViewScaler* mScale;

protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
};

#endif // CONTROLPOINTGRAPHICSITEM_H
