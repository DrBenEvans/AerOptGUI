#ifndef MESHGRAPHICSITEM_H
#define MESHGRAPHICSITEM_H

#include <QGraphicsItem>
#include "MeshModel.h"
#include "BoundaryPointModel.h"
#include "ViewScaler.h"
#include "ColorMapper.h"

class MeshView : public QGraphicsObject
{
    Q_OBJECT
public:
    MeshView(ViewScaler *scale, QGraphicsItem *parent = 0);

    void setMeshModel(MeshModel* meshModel);
    void setBoundaryPointModel(BoundaryPointModel* model);
    QRectF boundingRect() const override;
    QPainterPath shape() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;
    void clearBoundaryPoint();

public slots:
    void meshChanged();
    void boundaryPointsReset();

signals:
    void pointActivated(int);

private:
    void calcBoundingBox();
    uint getBrushSize();

    QColor color;
    MeshModel* mMeshModel;
    ViewScaler* mScale;
    ColorMapper mColorMapper;

    float mXmax;
    float mYmax;
    float mXmin;
    float mYmin;
    QRectF mBoundingBox;
    BoundaryPointModel* mBoundaryPointModel;

};

#endif // MESHGRAPHICSITEM_H
