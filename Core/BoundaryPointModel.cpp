#include "BoundaryPointModel.h"
#include "BoundaryPoint.h"
#include <QSharedPointer>

BoundaryPointModel::BoundaryPointModel(QObject *parent) : QObject(parent),
    mCurrentIndex(-1)
{
}

void BoundaryPointModel::setPoints(std::list<std::pair<float, float> > points)
{
    for (auto coord : points)
    {
        mBoundaryPoints.emplace_back(coord.first, coord.second);
    }
    emit boundaryPointsReset();
}

void BoundaryPointModel::setCurrentBoundaryPoint(int index) {
    mCurrentIndex = index;
}

int BoundaryPointModel::count() {
    return mBoundaryPoints.size();
}

int BoundaryPointModel::validIndex(int index) {
    if(index < 0) {
        return false;
    } else if (index >= count()) {
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

BoundaryPoint* BoundaryPointModel::currentPoint() {
    return point(mCurrentIndex);
}

void BoundaryPointModel::setControlBoundaryCorner(int index, QPointF position, CornerPosition corner ) {
    if(corner == CornerPosition::TOPLEFT) {
        point(index)->setTopLeftBound(position);
    } else {
        point(index)->setBottomRightBound(position);
    }

    emit controlBoundsChanged(index);
}
