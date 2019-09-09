#include "ParallelCoordinatesWindow.h"
#include "ui_ParallelCoordinatesWindow.h"

#include <QtCharts/QCategoryAxis>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLineSeries>
#include "ChartView.h"

#include "OptimisationModel.h"
#include "Optimisation.h"
#include "chart.h"

QT_CHARTS_USE_NAMESPACE
ParallelCoordinatesWindow::ParallelCoordinatesWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParallelCoordinatesWindow)
{
    ui->setupUi(this);
    //Possibly draw blank graph here
    chartView = new ChartView();
    ui->verticalLayout->addWidget(chartView, 1);
}

ParallelCoordinatesWindow::~ParallelCoordinatesWindow()
{
    delete ui;
}

void ParallelCoordinatesWindow::setOptimisationModel(OptimisationModel* model){
    mOptimisationModel = model;
}

void ParallelCoordinatesWindow::setCurrentOptimisationIndex(int index){
    if(index != mCurrentOptimisationIndex){
        mCurrentOptimisationIndex = index;
        drawGraph(CONTROL_POINT);
    }
}

Optimisation* ParallelCoordinatesWindow::currentOptimisation() {
    return mOptimisationModel->optimisation(mCurrentOptimisationIndex);
}

void ParallelCoordinatesWindow::on_switchGraphButton_clicked()
{
    //ui->verticalLayout->removeItem(ui->verticalLayout->itemAt(1));
    if (ui->switchGraphButton->isChecked()){
        ui->switchGraphButton->setText("Show Control Points Only");
        drawGraph(BOUNDARY_POINT);
    } else {
        ui->switchGraphButton->setText("Show All Boundary Points");
        drawGraph(CONTROL_POINT);
    }
    chartView->repaint();
}

QChartView* ParallelCoordinatesWindow::getGraph(){
    return chartView;
}

void ParallelCoordinatesWindow::drawGraph(GraphType type){

    // Set up chart object
    Chart *chart;

    switch(type){
        case CONTROL_POINT:
            chart = plotGraph(false);
            chart->setTitle("Control Point Displacement");
            break;

        case BOUNDARY_POINT:
            chart = plotGraph(true);
            chart->setTitle("Boundary Point Displacement");
            break;

        default:
            chart = blankGraph();
            break;

    }

    // Render widget
    chartView->setChart(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    //chartView->setParent(this);
}


Chart* ParallelCoordinatesWindow::plotGraph(bool showAll) {

    Chart *chart = new Chart();

    // Fetch the number of control points
    int numberOfControlPoints  = currentOptimisation()->controlPointCount();
    int numberOfBoundaryPoints = currentOptimisation()->initialBoundaryPoints().size();

    // Fetch the number of generations
    int numberOfGenerations = currentOptimisation()->noGens();

    // Fetch min and max fitness values
    std::pair<double,double> fitnessRange = currentOptimisation()->fitnessRangeBestAgents();
    double minFitness = fitnessRange.first;
    double maxFitness  = fitnessRange.second;

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
    axisY->setTickCount(11);
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

            // Add a point if we are showing all points, or only if it is a control point otherwise
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



        // Colour series according to fitness value
        double fitnessValue = currentOptimisation()->fitness(g,0);
        double normalisedFitness = normalise(fitnessValue, minFitness, maxFitness);
        series->setColor(*rgb(normalisedFitness));

        // Add new series to chart
        chart->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        //series->setName("Generation " + QString::number(g+1));

    }
    chart->legend()->setVisible(false);
    chart->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
    return chart;

}


double ParallelCoordinatesWindow::normalise(double x, double xMin, double xMax){
    return (x-xMin) / (xMax - xMin);
}


QColor* ParallelCoordinatesWindow::rgb(double ratio)
{
    //Flip so red is high and blue is low
    ratio = 1-ratio;
    //we want to normalize ratio so that it fits in to 5 regions
    //where each region is 256 units long
    int normalized = int(ratio * 256 * 4);

    //find the distance to the start of the closest region
    int x = normalized % 256;

    int red = 0, grn = 0, blu = 0;
    switch(normalized / 256)
    {
    case 0: red = 255;      grn = x;        blu = 0;       break;//red
    case 1: red = 255 - x;  grn = 255;      blu = 0;       break;//yellow
    case 2: red = 0;        grn = 255;      blu = x;       break;//green
    case 3: red = 0;        grn = 255 - x;  blu = 255;     break;//cyan
    case 4: red = x;        grn = 0;        blu = 255;     break;//blue
    }

    return new QColor(red, grn, blu);
}
