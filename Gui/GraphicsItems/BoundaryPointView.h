#ifndef CONTROLPOINTGRAPHICSITEM_H
#define CONTROLPOINTGRAPHICSITEM_H

#include <QGraphicsItem>
#include "ControlPointBoundingBox.h"
#include "BoundaryPointModel.h"
#include "ViewScaler.h"

/**
 * @brief The BoundaryPointView class is means of visually interacting with a single boundary point as a GraphicsItem.
 * Hovering over this item will 'activate' it in the BoundaryPointModel.
 * Double-clicking on this item will convert the boundary point into a control point and provide control point drag-handle functionality.
 */
class BoundaryPointView : public QGraphicsObject
{
    Q_OBJECT
public:
    BoundaryPointView(BoundaryPointModel *model, int index, ViewScaler* scaler, QGraphicsItem *parent = 0);

    /**
     * @brief boundingRect Returns the perimeter rectangle around the boundary point where it can be hovered over.
     * Implemented from QGraphicsItem
     * @return Hoverable Perimeter around Boundary Point
     */
    QRectF boundingRect() const override;

    /**
     * @brief paint Draws the boundary point
     * Implemented from QGraphicsItem
     * @param painter
     * @param item
     * @param widget
     */
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;

public slots:
    /**
     * @brief setActivePoint Executed in response to BoundaryPointModel::activeIndexChanged.
     * If this boundary point has been selected, change mActive to true,
     * disable dragging over the boundary point, and enable control point drag handles.
     * If this boundary point was selected, but another point has been selected, change mActive to false,
     * re-enable dragging over this boundary point, and disable control point drag handles.
     * If the activeIndexChanged value is not this boundary point, then do nothing.
     * @param index The index of the active (selected) boundary point.
     */
    void setActivePoint(int index);

    /**
     * @brief refreshControlPointState Executed in response to BoundaryPointModel::controlPointStateChanged.
     * If boundary point is now a control point, set mControl to true and set controlPointHandles to visible.
     * If boundary point is no longer a control point, set mControl to false and set controlPointHandles to not visible.
     * if control point status has not changed, then do nothing.
     * @param index Index of boundary point that has changed control point status
     * @param ctl Control Point Status: true iff this boundary point is a control point.
     */
    void refreshControlPointState(int index, bool ctl);

private:
    /**
     * @brief activated Returns whether this Boundary Point is currently activated (selected) (AKA mActive).
     * @return true iff activated, otherwise false.
     */
    bool activated() const;

    /**
     * @brief isControlPoint Returns whether this boundary point is a control point (AKA mControl).
     * @return true iff this boundary point is a control point, false otherwise.
     */
    bool isControlPoint();

    /**
     * @brief radius The distance between the precise boundary point location and its bounding box.
     * @return 5.0
     */
    qreal radius() const;

    /**
     * @brief mPointRect This rectangle is drawn by paint to indicate the location of the control point.
     * Its size is determined by radius()
     */
    QRectF mPointRect;

    /**
     * @brief mActive true iff this boundary point is activated (selected), otherwise false.
     */
    bool mActive;

    /**
     * @brief mControl true iff this boundary point is a control point.
     */
    bool mControl = false;
    ControlPointBoundingBox* mControlPointHandles = nullptr;

    /**
     * @brief mBoundaryPointModel Model this view belongs to.
     */
    BoundaryPointModel* mBoundaryPointModel;

    /**
     * @brief mBoundaryPointIndex Relative index of this boundary point.
     */
    int mBoundaryPointIndex;

    /**
     * @brief mScale ViewScaler
     */
    ViewScaler* mScale;

protected:

    /**
     * @brief hoverEnterEvent When hovering over a boundary point, change the cursor to Qt::ArrowCursor
     * and disable dragging i.e. mesh will not move when clicking while hovering.
     * Set the activeIndex of the BoundaryPointModel to the index of the Boundary Point we are hovering over.
     * Overrides hoverEnterEvent method in QGraphicsItem
     * @param event
     */
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event) override;

    /**
     * @brief hoverLeaveEvent After hovering over a boundary point re-enable dragging i.e. mesh will move when clicking on it.
     * Set the activeIndex of the BoundaryPointModel to the index of the Boundary Point we are hovering over.
     * Overrides hoverEnterEvent method in QGraphicsItem.
     * Note: the cursor is changed automatically if a a QGraphicsItem is draggable.
     * @param event
     */
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event) override;

    /**
     * @brief mouseDoubleClickEvent Invert the Control Point Status of this Boundary Point on double-clicking it.
     * @param event QGraphicsSceneMouseEvent Mouse Click event
     */
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
};

#endif // CONTROLPOINTGRAPHICSITEM_H
