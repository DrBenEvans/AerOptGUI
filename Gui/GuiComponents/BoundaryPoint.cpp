/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimisationRun.cpp
**
**********************************************/

#include <QDebug>
#include <QPointF>
#include <QRectF>
#include <iostream>
#include "BoundaryPoint.h"

BoundaryPoint::BoundaryPoint(QPointF coords, QRectF bounding_box) :
    mCoord(coords),
    mBounds(bounding_box)
{
    mIsControlPoint = false;
}

BoundaryPoint::BoundaryPoint(qreal xcoord, qreal ycoord, qreal bb_x, qreal bb_y, qreal bb_width, qreal bb_height) :
    mCoord(xcoord,ycoord),
    mBounds(bb_x,bb_y,bb_width,bb_height)
{
    mIsControlPoint = false;
}

BoundaryPoint::BoundaryPoint(QPointF coords) :
    mCoord(coords),
    mBounds(0,0,0,0)
{
    mIsControlPoint = false;
}

BoundaryPoint::BoundaryPoint(qreal xcoord, qreal ycoord) :
    mCoord(xcoord,ycoord),
    mBounds(0,0,0,0)
{
}

void BoundaryPoint::getBoundCoords(qreal *x1,qreal *y1,qreal *x2,qreal *y2) {
    mBounds.getCoords(x1,y1,x2,y2);
}

qreal BoundaryPoint::x() const {
    return mCoord.x();
}

qreal BoundaryPoint::y() const {
    return mCoord.y();
}

void BoundaryPoint::setBoundCoords(qreal x1,qreal y1,qreal x2,qreal y2) {
    mIsControlPoint = true;
    mBounds.setCoords(x1,y1,x2,y2);
}

uint BoundaryPoint::getSmoothing() {
    return mSmooth;
}

uint BoundaryPoint::getSmoothFactor() {
    return mSmooth / 10;
}

void BoundaryPoint::setSmoothing(uint smoothing) {
    mSmooth = smoothing;
}

bool BoundaryPoint::isControlPoint() {
    return mIsControlPoint;
}

