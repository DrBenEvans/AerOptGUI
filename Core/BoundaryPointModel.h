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
    int currentIndex();

    enum CornerPosition
    {
        TOPLEFT,
        BOTTOMRIGHT
    };
    void setControlBoundaryCorner(int index, QPointF pos, CornerPosition corner );
    void setActiveIndex(int index);
    void setControlPointState(int index, bool isControlPoint );

signals:
    void boundaryPointsReset();
    void controlPointBoundsChanged(int);
    void activeIndexChanged(int);
    void controlPointStateChanged(int);

public slots:

private:
    std::vector<BoundaryPoint> mBoundaryPoints;
    bool mCurrentPointSet;
    int mCurrentIndex;
    int validIndex(int index);
};

#endif // BOUNDARYPOINTMODEL_H
