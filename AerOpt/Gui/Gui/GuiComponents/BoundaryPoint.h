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
#include <QRectF>

#include "Enumerations.h"

class BoundaryPoint {

public:
    void setPoint(QPointF ctl_point);
    void setBounds(QRectF ctl_bounds);

    BoundaryPoint(QPointF coords);
    BoundaryPoint(qreal xcoord, qreal ycoord);
    BoundaryPoint(QPointF coords, QRectF bounding_box);
    BoundaryPoint(qreal xcoord, qreal ycoord, qreal bb_x, qreal bb_y, qreal bb_width, qreal bb_height);

    qreal x() const;
    qreal y() const;
    void getBoundCoords(qreal *x1,qreal *y1,qreal *x2,qreal *y2);
    void setBoundCoords(qreal x1,qreal y1,qreal x2,qreal y2);
    bool isControlPoint();
    uint getSmoothing();
    uint getSmoothFactor();
    void setSmoothing(uint smoothing);

private:
    QPointF mCoord;
    bool mIsControlPoint = false;
    QRectF mBounds;
    uint mSmooth = 0.0;

};

#endif // BOUNDARY_POINT_H
