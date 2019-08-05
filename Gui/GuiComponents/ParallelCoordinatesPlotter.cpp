#include "ParallelCoordinatesPlotter.h"
#include "QWidget"
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>
#include "BoundaryPoint.h"
#include "Mesh.h"
#include <QDebug>

QT_CHARTS_USE_NAMESPACE

ParallelCoordinatesPlotter::ParallelCoordinatesPlotter(QWidget *parent) : QWidget(parent)
{
    if (parent == nullptr){
        this->resize(1200, 500);
    } else{
        this->resize(parent->width(),parent->height());
    }
    //this->resize(1200, 500);
    //QWidget::showFullScreen();
}

void ParallelCoordinatesPlotter::setOptimisationModel(OptimisationModel* model) {
    mOptimisationModel = model;
}

Optimisation* ParallelCoordinatesPlotter::currentOptimisation() {
    return mOptimisationModel->optimisation(mCurrentOptimisationIndex);
}

void ParallelCoordinatesPlotter::setCurrentOptimisationIndex(int index) {
    mCurrentOptimisationIndex = index;
    drawGraph(CONTROL_POINT);
}

QChartView* ParallelCoordinatesPlotter::getGraph(){
    return chartView;
}

void ParallelCoordinatesPlotter::drawGraph(GraphType type){

    // Set up chart object
    QChart *chart = new QChart();
    chart->setTitle("Parallel Coordinates");

    switch(type){
        case CONTROL_POINT:
            chart = controlPointGraph(chart);
            break;

        case BOUNDARY_POINT:
            chart = boundaryPointGraph(chart);
            break;
    }

    // Render widget
    chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setParent(this);
    chartView->resize(this->width(),this->height());
    chartView->show();

}

QChart* ParallelCoordinatesPlotter::boundaryPointGraph(QChart *chart){

    // Initialise axes
    QCategoryAxis *axisX = new QCategoryAxis;
    QValueAxis *axisY = new QValueAxis;


    int bpCount = 0;

    //For each generation
    for(int g = 0; g < numberOfGenerations; g++){

        //Pick best agent (index  = 0)
        Mesh* currentMesh = currentOptimisation()->mesh(g, 0);
        ProfilePoints currentMeshProfile = currentMesh->boundaryPoints();

        // Loop through each boundary point and add a plot to the series
        series = new QLineSeries();
        bpCount = 0;
        for (ProfilePoints::iterator it = currentMeshProfile.begin(); it != currentMeshProfile.end(); ++it) {
            *series << QPointF(2*bpCount, it->first) << QPointF(2*bpCount+1, it->second);
            bpCount++;
        }

        // Attach Last BP series to graph
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        series->setName("Generation " + QString::number(g+1));
        chart->addSeries(series);

    }


QChart* ParallelCoordinatesPlotter::controlPointGraph(QChart *chart) {

    // Initialise axes
    QCategoryAxis *axisX = new QCategoryAxis;
    QValueAxis *axisY = new QValueAxis;


    // Fetch the number of control points
    int numberOfControlPoints  = currentOptimisation()->controlPointCount();
    std::vector<BoundaryPoint*> controlPoints = currentOptimisation()->controlPoints();
    int controlPointIndexes[numberOfControlPoints]; //Index of control points among boundary points

    // Get the indexes of all control points in the mesh
    std::vector<BoundaryPoint*> initialBoundary = currentOptimisation()->initialBoundaryPoints();

    for (int bp = 0; bp < initialBoundary.size(); bp++) {

        for (int cp = 0; cp < numberOfControlPoints; cp++) {

            // If x/y coordinates are the same then this is a control point
            if (initialBoundary.at(bp)->x() == controlPoints[cp]->x()
                    && initialBoundary.at(bp)->y() == controlPoints[cp]->y()) {

                controlPointIndexes[cp] = bp;
                break;
            }
        }
    }

    // Bounding box limit values for normalisation
    double minXValues[numberOfControlPoints];
    double minYValues[numberOfControlPoints];
    double maxXValues[numberOfControlPoints];
    double maxYValues[numberOfControlPoints];

    // Get min/max X/Y bounding box limits
    for(int cp = 0; cp< numberOfControlPoints; cp++){
        minXValues[cp] = controlPoints[cp]->controlPointRect().left();
        minYValues[cp] = controlPoints[cp]->controlPointRect().top();
        maxXValues[cp] = controlPoints[cp]->controlPointRect().right();
        maxYValues[cp] = controlPoints[cp]->controlPointRect().bottom();
    }


    //
    // Set up axes
    //
    // Set up X axis
    for (int i = 0; i < numberOfControlPoints; i++){
        axisX->append("X" + QString::number(i), 2*i);
        axisX->append("Y" + QString::number(i), 2*i+1);
    }
    axisX->setTickCount(numberOfControlPoints);
    //May need to change range back to 0
    axisX->setRange(-1, numberOfControlPoints*2);
    axisX->setTitleText("Control Points");
    axisX->AxisLabelsPositionOnValue;
    chart->addAxis(axisX, Qt::AlignBottom);

    // Set up Y axis
    axisY->setTickCount(11);
    axisY->setRange(0,1); // values are normalised
    axisY->setTitleText("Normalised Values");
    chart->addAxis(axisY, Qt::AlignLeft);



    // Fetch the number of generations
    int numberOfGenerations = currentOptimisation()->noGens();
    QLineSeries *series;

    //For each generation
    for(int g = 0; g < numberOfGenerations; g++){



        //Pick best agent
        //(at current, agents are returned in fitness-descending order from the optimiser so index  = 0)
        Mesh* bestAgent = currentOptimisation()->mesh(g, 0);
        ProfilePoints currentMeshProfile = bestAgent->boundaryPoints();

        // Loop through each control point and add a plot to the series
        series = new QLineSeries();
        for (int cpCount = 0; cpCount < numberOfControlPoints; cpCount++) {

            for (ProfilePoints::iterator it = currentMeshProfile.begin(); it != currentMeshProfile.end(); ++it) {
                if (bpCount == controlPointIndexes[cpCount]) {
                    double xmin = initialBoundary.at(cpCount)->pos().x() + minXValues[cpCount];
                    double xmax = initialBoundary.at(cpCount)->pos().x() + maxXValues[cpCount];
                    double ymin = initialBoundary.at(cpCount)->pos().y() + minYValues[cpCount];
                    double ymax = initialBoundary.at(cpCount)->pos().y()+ maxYValues[cpCount];

                    qInfo() << "BP" + QString::number(bpCount) + "x = " + QString::number(it->first) + ", min = " + QString::number(xmin) + ", max = " + QString::number(xmax) + ", norm = " + QString::number(normalise(it->first, xmin, xmax));
                    //normalise actual X/Y values for this agent w.r.t the bounding box limits
                    //Add normalised values to series

                    *series << QPointF(2*cpCount, normalise(it->first, xmin, xmax));
                    *series << QPointF(2*cpCount+1, normalise(it->second, ymin, ymax));
                    break;
                }

                bpCount++;
            }
        }

        // Add new series to chart
        chart->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        series->setName("Generation " + QString::number(g+1));

    }

    return chart;

}


double ParallelCoordinatesPlotter::normalise(double x, double xMin, double xMax){
    return (x-xMin) / (xMax - xMin);
}

