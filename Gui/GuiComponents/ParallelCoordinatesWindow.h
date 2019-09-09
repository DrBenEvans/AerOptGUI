#ifndef PARALLELCOORDINATESWINDOW_H
#define PARALLELCOORDINATESWINDOW_H

#include <QDialog>
#include "OptimisationModel.h"
#include "ChartView.h"
#include "chart.h"


namespace Ui {
class ParallelCoordinatesWindow;
}

QT_CHARTS_USE_NAMESPACE
class ParallelCoordinatesWindow : public QDialog
{
    Q_OBJECT

public:
    explicit ParallelCoordinatesWindow(QWidget *parent = nullptr);
    ~ParallelCoordinatesWindow();

    /**
     * @brief setOptimisationModel Sets the optimisation model.
     * @param model The model of loaded optimisations.
     */
    void setOptimisationModel(OptimisationModel* model);

    /**
     * @brief setCurrentOptimisationIndex Set the optimisation to display data for.
     * @param index Index of optimisation in optimisation model.
     */
    void setCurrentOptimisationIndex(int index);

    /**
     * @brief getGraph Returns the Parallel Coordinates graph
     * @return
     */
    QChartView* getGraph();

    /**
     * @brief The GraphType enum Type of graph to plot
     */
    enum GraphType
    {
        BLANK,                 // Blank Graph when an optimisation is not loaded
        CONTROL_POINT,         // Parallel Coordinate Graph of Control Points only
        BOUNDARY_POINT,        // Parallel Coordinate Graph of all boundary points

    };

    /**
     * @brief update Refreshes the graph
     */
    void update();

private slots:

    /**
     * @brief on_switchGraphButton_clicked When switchGraphButton is clicked, alternate between
     * the boundary point graph and the control point graph depending on whether the button is checked.
     */
    void on_switchGraphButton_clicked();

private:
    Ui::ParallelCoordinatesWindow *ui;

    /**
     * @brief mOptimisationModel Model of all loaded optimisations.
     */
    OptimisationModel* mOptimisationModel;

    /**
     * @brief mCurrentOptimisationIndex Index of currently selected optimisation in model.
     */
    int mCurrentOptimisationIndex = -1;

    /**
     * @brief currentOptimisation Returns the currently selected optimisation.
     * @return
     */
    Optimisation* currentOptimisation();

    /**
     * @brief drawGraph Updates chartView to return the specified graph and re-plot.
     * @param type Specifies what kind of graph to draw based on hte GraphType Enum.
     */
    void drawGraph(GraphType type);

    /**
     * @brief plotGraph Draws a parallel coordinates chart of normalised displacement of boundary points.
     * Can either return a chart with all Boundary Points, or just with Boundary Points that are Control Points.
     * @param showAll
     * @return
     */
    Chart* plotGraph(bool showAll);

    /**
     * @brief blankGraph Draws a blank graph for display purposes
     * @return Modified chart
     */
    Chart* blankGraph();

    /**
     * @brief chartView The Parallel Coordinates Graph
     */
    ChartView *chartView = nullptr;

    /**
     * @brief currentlyDisplayedGraph The type of graph we are currently displaying.
     */
    GraphType currentlyDisplayedGraph = BLANK;

    /**
     * @brief normalise Simple feature scaling method of normalisation that transforms x into a value between 0 and 1,
     * based on the minimum and maximum range values given.
     * @param x The value to normalise
     * @param xMin Minimum possible value in range
     * @param xMax Maximum possible value in range
     * @return Normalised value of x
     */
    double normalise(double x, double xMin, double xMax);

    /**
     * @brief rgb Returns a QColor between Red and Blue based on a normalised value
     * @param ratio Normalised value between 0 and 1 where 0 is Red and 1 is Blue
     * @return QColor
     */
    QColor* rgb(double ratio);

};

#endif // PARALLELCOORDINATESWINDOW_H
