#include "BoundaryPointModel.h"
#include "BoundaryPoint.h"
#include <QSharedPointer>

BoundaryPointModel::BoundaryPointModel(QObject *parent) : QObject(parent),
    mCurrentIndex(-1)
{
}

void BoundaryPointModel::clearPoints() {
    mBoundaryPoints.clear();
    emit boundaryPointsReset();
}

int BoundaryPointModel::controlPointCount() {
    int count = 0;
    for( auto& bp : mBoundaryPoints )
        if(bp.isControlPoint())
            count++;

    return count;
}

std::vector<BoundaryPoint*> BoundaryPointModel::controlPoints() {
    std::vector<BoundaryPoint*> ctlPoints;
    for( auto& bp : mBoundaryPoints )
        if(bp.isControlPoint())
            ctlPoints.push_back(&bp);

    return ctlPoints;
}

std::vector<BoundaryPoint*> BoundaryPointModel::boundaryPoints() {
    std::vector<BoundaryPoint*> boundaryPoints;
    for( auto& bp : mBoundaryPoints) {
        boundaryPoints.push_back(&bp);
    }

    return boundaryPoints;
}

void BoundaryPointModel::setPoints(std::list<std::pair<float, float> > points)
{
    for (auto coord : points)
    {
        mBoundaryPoints.emplace_back(coord.first, coord.second);
    }
    emit boundaryPointsReset();
}

void BoundaryPointModel::setActiveIndex(int index) {
    mCurrentIndex = index;
    emit activeIndexChanged(index);
}

int BoundaryPointModel::count() {
    return mBoundaryPoints.size();
}

int BoundaryPointModel::validIndex(int index) {
    if(index < 0 || index >= count()) {
        return false;
    } else {
        return true;
    }
}

BoundaryPoint* BoundaryPointModel::point(int index) {
    if(validIndex(index)) {
        return &mBoundaryPoints.at(index);
    } else {
        return nullptr;
    }
}

void BoundaryPointModel::setControlPointState(int index, bool isControlPoint ) {
    bool controlPrevState = point(index)->isControlPoint();
    if(controlPrevState != isControlPoint) {
        point(index)->setControlPoint(isControlPoint);
        emit controlPointStateChanged(index, isControlPoint);
    }
}


void BoundaryPointModel::setControlBoundaryCorner(int index, QPointF position, CornerPosition corner) {
    QRectF currentRect = point(index)->controlPointRect();
    if(corner == CornerPosition::TOPLEFT) {
        if(currentRect.topLeft() != position) {
            point(index)->setTopLeftBound(position);
            emit controlPointBoundsChanged(index);
        }
    } else {
        if(currentRect.bottomRight() != position) {
            point(index)->setBottomRightBound(position);
            emit controlPointBoundsChanged(index);
        }
    }
}

int BoundaryPointModel::currentIndex() {
    return mCurrentIndex;
}

void BoundaryPointModel::setControlPointBounds(int index, QRectF ctlBounds) {
    BoundaryPoint* boundaryPoint = point(index);
    if(boundaryPoint) {
        QRectF currentRect = point(index)->controlPointRect();
        if(currentRect != ctlBounds) {
            boundaryPoint->setControlPointRect(ctlBounds);
            emit controlPointBoundsChanged(index);
        }
    }
}
