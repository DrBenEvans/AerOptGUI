/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimisationRun.h
**
**********************************************/

#ifndef BOUNDARY_POINT_H
#define BOUNDARY_POINT_H

#include <utility>
#include <vector>
#include <QRectF>
#include <QObject>

#include "Enumerations.h"

class BoundaryPoint : public QObject {
    Q_OBJECT

public:
    void setPoint(QPointF ctl_point);
    void setBounds(QRectF ctl_bounds);

    BoundaryPoint(qreal xcoord, qreal ycoord);

    qreal x() const;
    qreal y() const;
    QRectF controlPointRect();
    bool isControlPoint();
    uint getSmoothing();
    uint getSmoothFactor();
    void setSmoothing(uint smoothing);
    void setTopLeftBound(QPointF pos);
    void setBottomRightBound(QPointF pos);

signals:
    void controlRectChanged();

private:
    QPointF mCoord;
    bool mIsControlPoint = false;
    QRectF mControlPointRect;
    uint mSmooth = 0.0;

};

typedef std::vector<std::shared_ptr<BoundaryPoint>> Boundaries;

#endif // BOUNDARY_POINT_H
