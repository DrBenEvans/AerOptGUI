#ifndef PARALLELCOORDINATESPLOTTER_H
#define PARALLELCOORDINATESPLOTTER_H
#include <QWidget>
#include "OptimisationModel.h"
#include "Optimisation.h"
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

QT_CHARTS_USE_NAMESPACE
class ParallelCoordinatesPlotter : public QWidget
{
    Q_OBJECT
public:
    explicit ParallelCoordinatesPlotter(QWidget *parent = nullptr);

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
        CONTROL_POINT,         // Parallel Coordinate Graph of Control Points only
        BOUNDARY_POINT,        // Parallel Coordinate Graph of all boundary points

    };

signals:

public slots:


private:

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
    // Display size
    /**
     * @brief mPointSize The diameter of plotted points.
     */
    int mPointSize = 15;

    /**
     * @brief mLineSize The thickness of plotted lines.
     */
    int mLineSize = 5;

    /**
     * @brief drawGraph
     */
    void drawGraph(GraphType type);

    /**
     * @brief boundaryPointGraph Draws a parallel coordinate graph of all boundary points in the current optimisation.
     */
    QChart* boundaryPointGraph(QChart *chart);

    /**
     * @brief controlPointGraph Draws a parallel coordinate graph of all control points in the current optimisation.
     */
    QChart* controlPointGraph(QChart *chart);

    /**
     * @brief chartView The Parallel Coordinates Graph
     */
    QChartView *chartView = nullptr;

    /**
     * @brief normalise Simple feature scaling method of normalisation that transforms x into a value between 0 and 1,
     * based on the minimum and maximum range values given.
     * @param x The value to normalise
     * @param xMin Minimum possible value in range
     * @param xMax Maximum possible value in range
     * @return Normalised value of x
     */
    double normalise(double x, double xMin, double xMax);

};



#endif // PARALLELCOORDINATESPLOTTER_H
