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

class BoundaryPoint {

public:
    void setPoint(QPointF ctl_point);

    BoundaryPoint(qreal xcoord, qreal ycoord);

    // point position
    qreal x() const;
    qreal y() const;
    QPointF pos();

    // control points bounds
    void setControlPoint(bool isCtl);
    bool isControlPoint();
    QRectF controlPointRect();
    QRectF setControlPointRect(QRectF ctlPointRect);
    void setTopLeftBound(QPointF pos);
    void setBottomRightBound(QPointF pos);


    uint getSmoothing();
    uint getSmoothFactor();
    void setSmoothing(uint smoothing);

private:
    QPointF mCoord;
    bool mIsControlPoint = false;
    QRectF mControlPointRect;
    uint mSmooth = 0.0;
};

#endif // BOUNDARY_POINT_H
