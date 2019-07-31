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
    this->resize(1200, 500);
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
    boundaryPointGraph();
}

QChartView* ParallelCoordinatesPlotter::getGraph(){
    return chartView;
}

void ParallelCoordinatesPlotter::boundaryPointGraph(){

    // Set up chart object
    QChart *chart = new QChart();
    //chart->legend()->hide();
    chart->setTitle("Parallel Coordinates");

    // Fetch the number of generations
    int numberOfGenerations = currentOptimisation()->noGens();
    QLineSeries *series;

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

    // Set up X axis
    chart->addAxis(axisX, Qt::AlignBottom);
    for (int i = 0; i < bpCount; i++){
        axisX->append("BP_" + QString::number(i) + "_X", 2*i);
        axisX->append("BP_" + QString::number(i) + "_Y", 2*i+1);
    }
    axisX->setTickCount(bpCount);
    axisX->setRange(0, bpCount*2 - 1);
    //chart->axes(Qt::Vertical).first()->setRange(0, bpCount*2);

    // Set up Y axis
    axisY->setTickCount(11);
    axisY->setRange(0,1); // values are normalised
    chart->addAxis(axisY, Qt::AlignLeft);

    // Render widget
    chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setParent(this);
    chartView->resize(this->width(),this->height());
    chartView->show();

}

double ParallelCoordinatesPlotter::normalise(double x, double xMin, double xMax){
    return (x-xMin) / (xMax - xMin);
}

