#ifndef BOUNDARYPOINTMODEL_H
#define BOUNDARYPOINTMODEL_H

#include <QObject>
#include "BoundaryPoint.h"

class BoundaryPointModel : public QObject
{
    Q_OBJECT
public:
    explicit BoundaryPointModel(QObject *parent = nullptr);

    void clear();
    void setPoints(std::list<std::pair<float,float>> points);
    void setCurrentBoundaryPoint(int index);
    int count();
    BoundaryPoint* point(int index);
    BoundaryPoint *currentPoint();
    bool currentPointSet();

    enum CornerPosition
    {
        TOPLEFT,
        BOTTOMRIGHT
    };
    void setControlBoundaryCorner(int index, QPointF pos, CornerPosition corner );
    void setActiveIndex(int index);

signals:
    void boundaryPointsReset();
    void controlBoundsChanged(int);
    void activeIndexChanged(int);

public slots:

private:
    std::vector<BoundaryPoint> mBoundaryPoints;
    bool mCurrentPointSet;
    int mCurrentIndex;
    int validIndex(int index);
};

#endif // BOUNDARYPOINTMODEL_H
