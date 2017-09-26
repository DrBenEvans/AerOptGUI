#include "PlotterWidget.h"
#include "Optimisation.h"

PlotterWidget::PlotterWidget(QWidget *parent) :
    QCustomPlot(parent)
{
    for (uint i = 0; i < 50; ++i) {
        QCPGraph* graph = addGraph();
        graph->setSelectable(QCP::stSingleData);
        mGraphs.push_back(graph);
    }

    connect(this, &PlotterWidget::selectionChangedByUser, this, &PlotterWidget::setCurrentlySelectedPoint);

    clearData();
}

void PlotterWidget::setCurrentlySelectedPoint() {
    int iGen = -1;
    int agent = -1;
    foreach(QCPGraph* graph, selectedGraphs()) {
        for(int i=0; i<mGraphs.size(); i++) {
            if(graph==mGraphs.at(i)) {
                iGen = i;
            }
        }
        QCPDataSelection selection = graph->selection();
        if(selection.dataRangeCount() > 0) {
            QCPDataRange range = selection.dataRange();
            agent = range.begin();
        }
    }

    if(agent > 0 && iGen > 0) {
        emit selectedPointChanged(iGen, agent);
    }
}

void PlotterWidget::setOptimisationModel(OptimisationModel* model) {
    mOptimisationModel = model;
}

void PlotterWidget::setCurrentOptimisationIndex(int index) {
    mCurrentOptimisationIndex = index;
    updatePlot();
}

void PlotterWidget::updatePlot()
{
    std::shared_ptr<Optimisation> opt = mOptimisationModel->optimisation(mCurrentOptimisationIndex);
    bool r = true;

    //Get data for display
    std::vector<std::vector<double>> fitness = opt->fitness();
    uint nGens = fitness.size();
    if(nGens == 0) {
        clearData();
        return;
    }
    int noColumns = fitness.back().size();

    if (r)
    {

        //Compute ranges
        QVector<double> x(nGens);
        for (int i = 0; i < nGens; ++i)
        {
            x[i] = i;
        }


        //for each fitness column do this
        for (int fitCol = 0; fitCol < noColumns; ++fitCol)
        {
            QVector<double> y(nGens);
            for (int i = 0; i < nGens; ++i)
            {
                y[i] = fitness.at(i).at(fitCol);
            }
            //Build graph
            graph(fitCol)->setData(x, y);
        //end do this
        }


        // give the axes some labels:
        xAxis->setLabel("Generation No.");
        yAxis->setLabel("Fitness");
        // set axes ranges, so we see all data:
        rescaleAxes(true);

        if (nGens < 11)
        {
            xAxis->setRange(0, 10);
        }
        else
        {
            xAxis->setRange(0, nGens);
        }

        double min, max;
        max = yAxis->range().upper * 1.05;
        min = yAxis->range().lower * 0.95;

        yAxis->setRange(min, max);
    }

    //Plot graph
    replot();
}

void PlotterWidget::clearData()
{
    // give the axes some labels:
    xAxis->setLabel("Generation No.");
    yAxis->setLabel("Fitness");
    // set axes ranges, so we see all data:
    xAxis->setRange(0, 10);
    yAxis->setRange(0, 1);

    int num = graphCount();
    for (int i = 0; i < num; ++i)
    {
        QVector<double> x(i);
        QVector<double> y(i);
        graph(i)->setData(x, y);
        graph(i)->setPen(QPen(QColor(255, 100, 0)));
        graph(i)->setLineStyle(QCPGraph::lsLine);
        graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
    }

    setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    replot();
}
