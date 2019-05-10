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

/**
 * @brief The BoundaryPoint class specifies a sample point where the Mesh meets the profile.
 * A boundary point can be set as a control point which has a bounding box specifying the
 * maximum distance that an optimised point may be from the original BoundaryPoint during
 * the optimisation process.
 */
class BoundaryPoint {

public:
    /**
     * @brief BoundaryPoint Constructor Method.
     * @param xcoord X Coordinate of Boundary Point
     * @param ycoord Y Coordinate of Boundary Point
     */
    BoundaryPoint(qreal xcoord, qreal ycoord);

    // point position
    /**
     * @brief x Returns the x coordinate of the boundary point.
     * @return x coordinate
     */
    qreal x() const;

    /**
     * @brief y Returns the y coordinate of the boundary point.
     * @return y coordinate
     */
    qreal y() const;

    /**
     * @brief x Returns the coordinate point of the boundary point.
     * @return QPointF coordinate point
     */
    QPointF pos();

    // control node bounds
    /**
     * @brief setControlPoint Set whether this boundary point is also a control point.
     * @param isCtl true if point is control point, false if not a control point
     */
    void setControlPoint(bool isCtl);

    /**
     * @brief isControlPoint Returns whether this boundary point is a control point in the current optimisation.
     * @return true iff boundary point is a control point
     */
    bool isControlPoint();

    /**
     * @brief controlPointRect Returns the Control Point bounding box surrounding this boundary point.
     * @return Control Point bounding box
     */
    QRectF controlPointRect();

    /**
     * @brief setControlPointRect Explicitly set the Control Point bounding box for this Boundary Point as a QRect.
     * @param ctlPointRect The Control Point bounding box.
     * @return The newly set Control Point bounding box.
     */
    QRectF setControlPointRect(QRectF ctlPointRect);

    /**
     * @brief setTopLeftBound Set top-left bound for Control Point bounding box which
     * implicitly specifies the minimum x coordinate and maximum y coordinate for the Control Point.
     * @param pos A coordinate point for where the box should extend to.
     */
    void setTopLeftBound(QPointF pos);

    /**
     * @brief setBottomRightBound Set bottom-right bound for Control Point bounding box which
     * implicitly specifies the maximum x coordinate and minimum y coordinate for the Control Point.
     * @param pos A coordinate point for where the box should extend to.
     */
    void setBottomRightBound(QPointF pos);

    /**
     * @brief getSmoothing Returns the value of the smoothing parameter.
     * @return
     */
    uint getSmoothing();

    /**
     * @brief getSmoothFactor Returns the value of the smoothing factor (smoothing parameter / 10).
     * @return
     */
    float getSmoothFactor();

    /**
     * @brief setSmoothing Set the value of the smoothing parameter
     * @param smoothing Value for smoothing parameter.
     */
    void setSmoothing(uint smoothing);

private:
    /**
     * @brief mCoord Coordinate point for this Boundary Point
     */
    QPointF mCoord;

    /**
     * @brief mIsControlPoint States whether this Boundary Point is also a control point in the given optimisation
     */
    bool mIsControlPoint = false;

    /**
     * @brief mControlPointRect Rectangle surrounding Boundary Point that represents the control point parameters.
     */
    QRectF mControlPointRect;

    /**
     * @brief mSmooth Smoothing Parameter
     */
    uint mSmooth = 0.0;
};

#endif // BOUNDARY_POINT_H
