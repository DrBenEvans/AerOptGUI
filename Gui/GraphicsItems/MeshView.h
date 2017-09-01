#ifndef MESHGRAPHICSITEM_H
#define MESHGRAPHICSITEM_H

#include <QGraphicsItem>
#include "Mesh.h"

class MeshView : public QGraphicsItem
{
public:
    MeshView(int scale, QGraphicsItem *parent = 0);

    void setMesh(std::shared_ptr<Mesh> mesh);
    QRectF boundingRect() const override;
    QPainterPath shape() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;
    void showControlPoints(bool visible);
    void setBoundaryPoints(Boundaries& boundaryPoints);

public slots:
    void meshChanged();

private:
    void calcBoundingBox();
    uint getBrushSize();
    qreal w(qreal width);
    qreal h(qreal height);

    QColor color;
    std::shared_ptr<Mesh> mMesh;
    int mScale;

    float mXmax;
    float mYmax;
    float mXmin;
    float mYmin;
    QRectF mBoundingBox;

};

#endif // MESHGRAPHICSITEM_H
