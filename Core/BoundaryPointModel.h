#ifndef BOUNDARYPOINTMODEL_H
#define BOUNDARYPOINTMODEL_H

#include <QObject>
#include "BoundaryPoint.h"

/**
 * @brief The BoundaryPointModel class contains a collection of Boundary Points between
 * a profile and enclosing mesh. This class contains operations for interacting these
 * boundary points together or separately and selecting points.
 */
class BoundaryPointModel : public QObject
{
    Q_OBJECT
public:
    explicit BoundaryPointModel(QObject *parent = nullptr);

    /**
     * @brief clearPoints Clears list of Boundary Points (mesh boundary)
     * Emits the signal boundaryPointsReset()
     */
    void clearPoints();

    /**
     * @brief setPoints Set a given list of coordinates to be the new list of boundary points.
     * Changes the list mBoundaryPoints.
     * Emits boundaryPointsReset().
     * @param points List of coordinate points
     */
    void setPoints(std::list<std::pair<float,float>> points);

    /**
     * @brief count Returns the number of boundary points.
     * @return The number of boundary points.
     */
    int count();

    /**
     * @brief point Retrieve the boundary point of given index from mBoundaryPoints
     * @param index Boundary point index
     * @return The boundary point at given index
     */
    BoundaryPoint* point(int index);

    /**
     * @brief currentIndex Returns the index of the point that is currently selected
     * @return Index of currently selected point.
     */
    int currentIndex();

    /**
     * @brief The CornerPosition enum is for representing corners of a
     * Control Point bounding box which have drag handles
     */
    enum CornerPosition
    {
        TOPLEFT,
        BOTTOMRIGHT
    };

    /**
     * @brief setControlBoundaryCorner
     * @param index BoundaryPointIndex
     * @param pos X,Y Position for Control Point relative to boundary point
     * @param corner Corner of the control bounding box we are setting
     */
    void setControlBoundaryCorner(int index, QPointF pos, CornerPosition corner);

    /**
     * @brief setActiveIndex Sets the currently selected index, mCurrentIndex to index
     * Emits activeIndexChanged(index) signal
     * @param index The Index that is selected.
     */
    void setActiveIndex(int index);

    /**
     * @brief setControlPointState Set whether a boundary point is a control point.
     * Emits controlPointStateChanged(index, isControlPoint) on change.
     * @param index The BoundaryPoint index
     * @param isControlPoint Whether this boundary point should be a control point.
     */
    void setControlPointState(int index, bool isControlPoint );

    /**
     * @brief controlPointCount Returns the number of boundary points that are also control points.
     * @return The number of control points
     */
    int controlPointCount();

    /**
     * @brief setControlPointBounds Set the bounding box of a control point.
     * @param index BoundaryPoint index
     * @param ctlBounds Rectangle of Control Point bounding box.
     */
    void setControlPointBounds(int index, QRectF ctlBounds);

    /**
     * @brief controlPoints Returns the list of BoundaryPoints that are control points
     * @return List of ControlPoints
     */
    std::vector<BoundaryPoint *> controlPoints();

    /**
     * @brief boundaryPoints Returns ths list of all boundary points in the model.
     * @return List of boundary points
     */
    std::vector<BoundaryPoint*> boundaryPoints();


signals:
    /**
     * @brief boundaryPointsReset Emits when a new / blank list of boundary points is instantaited.
     */
    void boundaryPointsReset();

    /**
     * @brief controlPointBoundsChanged Emits where the bounding box of a control point changes.
     * Such as when it is created, or changed manually with setControlBoundaryCorner, or setControlPointBounds.
     * @param index The index of the BoundaryPoint that has been changed
     */
    void controlPointBoundsChanged(int index);

    /**
     * @brief activeIndexChanged Emits when the currently selected Boundary Point changes with setActiveIndex.
     * @param index The index of the BoundaryPoint that is selected.
     */
    void activeIndexChanged(int index);

    /**
     * @brief controlPointStateChanged Emits when a BoundaryPoint becomes a control point,
     * or a boundary point is no longer a control point.
     * @param index The boundary point whose control point status is changing.
     * @param isCtl Whether the boundary point is a control point
     */
    void controlPointStateChanged(int index, bool isCtl);

public slots:

private:
    /**
     * @brief mBoundaryPoints List of all Boundary Points for the current profile / mesh
     */
    std::vector<BoundaryPoint> mBoundaryPoints;

    /**
     * @brief mCurrentIndex The index of the currently selected Boundary Point in mBoundaryPoints
     */
    int mCurrentIndex;

    /**
     * @brief validIndex Checks if a given index is valid given the stores list of Boundary Points in mBoundaryPoints.
     * If index is less than zero, or greater than number of points it is invalid.
     * @param index The index to check for validity.
     * @return
     */
    int validIndex(int index);
};

#endif // BOUNDARYPOINTMODEL_H
