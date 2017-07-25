/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ProjectData.cpp
**
**********************************************/

#include <QDebug>
#include <QPointF>
#include <QRectF>
#include <iostream>
#include "BoundaryPoint.h"

BoundaryPoint::BoundaryPoint(QPointF coords, QRectF bounding_box) {
    mCoord = coords;
    mBounds = &bounding_box;
    mIsControlPoint = true;
}

BoundaryPoint::BoundaryPoint(qreal xcoord, qreal ycoord, qreal bb_x, qreal bb_y, qreal bb_width, qreal bb_height) {
    mCoord = QPointF(xcoord,ycoord);
    mBounds = new QRectF(bb_x,bb_y,bb_width,bb_height);
    mIsControlPoint = true;
}

BoundaryPoint::BoundaryPoint(QPointF coords) {
    mCoord = coords;
}

BoundaryPoint::BoundaryPoint(qreal xcoord, qreal ycoord) {
    mCoord = QPointF(xcoord,ycoord);
}

void BoundaryPoint::getBoundCoords(qreal *x1,qreal *y1,qreal *x2,qreal *y2) {
    mBounds->getCoords(x1,y1,x2,y2);
}

qreal BoundaryPoint::x() {
    return mCoord.x();
}

qreal BoundaryPoint::y() {
    return mCoord.y();
}

void BoundaryPoint::setBoundCoords(qreal x1,qreal y1,qreal x2,qreal y2) {
    mIsControlPoint = true;
    if(mBounds==NULL) {
        mBounds = new QRectF(x1,y1,x2,y2);
    } else {
        mBounds->setCoords(x1,y1,x2,y2);
    }
}

uint BoundaryPoint::getSmoothing() {
    return mSmooth;
}

void BoundaryPoint::setSmoothing(uint smoothing) {
    mSmooth = smoothing;
}

bool BoundaryPoint::isControlPoint() {
    if (mIsControlPoint) {
        return false;
    } else {
        return true;
    }

}

