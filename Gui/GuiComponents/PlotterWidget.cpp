#include "PlotterWidget.h"
#include "Optimisation.h"

PlotterWidget::PlotterWidget(QWidget *parent) :
    QCustomPlot(parent)
{
    for (uint i = 0; i < 50; ++i) addGraph();

    clearData();

    connect(this, &PlotterWidget::selectionChangedByUser, this, &PlotterWidget::signalCurrentlySelectedPoint);
}

void PlotterWidget::setCurrentlySelectedPoint(int iGen, int agent) {
    std::pair<int, int> agent_gen = getSelection();
    if(agent_gen.first != agent || agent_gen.second != iGen) {
        QCPDataRange dataRange(iGen, iGen+1);
        QCPDataSelection dataSelection;
        dataSelection += dataRange;
        graph(agent)->setSelection(dataSelection);
    }
}

std::pair<int, int> PlotterWidget::getSelection() {
    std::pair<int, int> agent_gen;
    agent_gen.first = -1;
    agent_gen.second = -1;
    foreach(QCPGraph* graphtmp, selectedGraphs()) {
        for(int i=0; i<graphCount(); i++) {
            if(graphtmp==graph(i)) {
                agent_gen.first = i;
            }
        }
        QCPDataSelection selection = graph(agent_gen.first)->selection();
        if(selection.dataRangeCount() > 0) {
            QCPDataRange range = selection.dataRange();
            agent_gen.second = range.begin();
        }
    }

    return agent_gen;
}

void PlotterWidget::signalCurrentlySelectedPoint() {
    std::pair<int, int> agent_gen = getSelection();

    if(agent_gen.first >= 0 && agent_gen.second >= 0) {
        emit selectedPointChanged(agent_gen.second, agent_gen.first);
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
        for (int agent = 0; agent < graphCount(); ++agent)
        {
            if(agent < noColumns) {
                QVector<double> y(nGens);
                for (int i = 0; i < nGens; ++i)
                {
                    y[i] = fitness.at(i).at(agent);
                }
                //Build graph
                graph(agent)->setData(x, y);
                graph(agent)->setVisible(true);
            } else {
                graph(agent)->setVisible(false);
            }
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
        graph(i)->setVisible(false);
        graph(i)->setSelectable(QCP::stSingleData);
    }

    setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    replot();
}
