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
    QChart *chart;

    switch(type){
        case CONTROL_POINT:
            chart->setTitle("Control Point Displacement");
            chart = plotGraph(false);
            break;

        case BOUNDARY_POINT:
            chart->setTitle("Boundary Point Displacement");
            chart = plotGraph(true);
            break;
    }

    // Render widget
    chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setParent(this);
    chartView->resize(this->width(),this->height());
    chartView->show();

}

QChart* ParallelCoordinatesPlotter::plotGraph(bool showAll) {

    QChart *chart = new QChart();

    // Fetch the number of control points
    int numberOfControlPoints  = currentOptimisation()->controlPointCount();
    int numberOfBoundaryPoints = currentOptimisation()->initialBoundaryPoints().size();

    // Fetch the number of generations
    int numberOfGenerations = currentOptimisation()->noGens();

    std::vector<BoundaryPoint*> controlPoints = currentOptimisation()->controlPoints();

    // Get the indexes of all control points in the mesh
    std::vector<BoundaryPoint*> initialBoundary = currentOptimisation()->initialBoundaryPoints();


    //
    // Set up axes
    //
    // Initialise axes
    QCategoryAxis *axisX = new QCategoryAxis;
    QValueAxis *axisY = new QValueAxis;

    // Set up X axis
    int numberOfGraphPoints;
    if (showAll){
        numberOfGraphPoints = numberOfBoundaryPoints;
        axisX->setTitleText("Boundary Points");
    } else {
        numberOfGraphPoints = numberOfControlPoints;
        axisX->setTitleText("Control Points");
    }

    for (int i = 0; i < numberOfGraphPoints; i++){
        axisX->append("X" + QString::number(i), 2*i);
        axisX->append("Y" + QString::number(i), 2*i+1);
    }
    axisX->setTickCount(numberOfGraphPoints);
    //May need to change range back to 0
    axisX->setRange(-1, numberOfGraphPoints*2);
    axisX->setLabelsPosition(axisX->AxisLabelsPositionOnValue);
    chart->addAxis(axisX, Qt::AlignBottom);

    // Set up Y axis
    axisY->setTickCount(21);
    axisY->setRange(0,1); // values are normalised
    axisY->setTitleText("Normalised Displacement Values");
    chart->addAxis(axisY, Qt::AlignLeft);


    // Get minimum and maximum displacement values from original point for normalisation
    double min = 0, max = 0;
    for(int g = 0; g < numberOfGenerations; g++){

        //Get best agent of generation
        //(at current, agents are returned in fitness-descending order from the optimiser so index  = 0)
        Mesh* bestAgent = currentOptimisation()->mesh(g, 0);
        ProfilePoints currentMeshProfile = bestAgent->boundaryPoints();

        // Loop through each profile point and check its XY displacement
        int bpCount = 0;
        for (ProfilePoints::iterator it = currentMeshProfile.begin(); it != currentMeshProfile.end(); ++it) {

            // Check X displacement
            double differenceX = it->first  -  initialBoundary.at(bpCount)->pos().x();
            if(differenceX > max){
                max = differenceX;
            } else if (differenceX < min) {
                min = differenceX;
            }

            // Check Y Displacement
            double differenceY = it->second  -  initialBoundary.at(bpCount)->pos().y();
            if(differenceY > max){
                max = differenceY;
            } else if (differenceY < min) {
                min = differenceY;
            }
            bpCount++;

        }
    }


    //Loop through the best agent of each generation and plot normalised displacement
    for(int g = 0; g < numberOfGenerations; g++){

        //Pick best agent
        //(at current, agents are returned in fitness-descending order from the optimiser so index  = 0)
        Mesh* bestAgent = currentOptimisation()->mesh(g, 0);
        ProfilePoints currentMeshProfile = bestAgent->boundaryPoints();

        // Loop through each point and add a plot to the series
        QLineSeries *series = new QLineSeries();
        int bpCount = 0;
        int cpCount = 0;
        for (ProfilePoints::iterator it = currentMeshProfile.begin(); it != currentMeshProfile.end(); ++it) {

            // Add a point if we are showing all points, or if it is a control point if we are not
            if (showAll || currentOptimisation()->initialBoundaryPoints().at(bpCount)->isControlPoint()) {

                double differenceX = it->first  -  initialBoundary.at(bpCount)->pos().x();
                double differenceY = it->second  -  initialBoundary.at(bpCount)->pos().y();

                //normalise actual X/Y values for this agent w.r.t the maximum displacement values
                //Add normalised values to series
                if (!showAll){
                    *series << QPointF(2*cpCount, normalise(differenceX, min, max));
                    *series << QPointF(2*cpCount+1, normalise(differenceY, min, max));
                    cpCount++;
                } else {
                    *series << QPointF(2*bpCount, normalise(differenceX, min, max));
                    *series << QPointF(2*bpCount+1, normalise(differenceY, min, max));
                }

            }
            bpCount++;
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

