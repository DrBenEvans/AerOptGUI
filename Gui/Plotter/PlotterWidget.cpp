#include "PlotterWidget.h"

PlotterWidget::PlotterWidget(QCustomPlot *parent) :
    QCustomPlot(parent)
{
    for (uint i = 0; i < 50; ++i) addGraph();

    clearData();
}

void PlotterWidget::setData(const int gen, const QVector<double> nests)
{
    bool r = true;

    //Get data for display
    mGData.push_back(nests);
    r &= mGData.size() == gen;

    if (r)
    {

        //Compute ranges
        QVector<double> x(gen);
        for (int i = 0; i < gen; ++i)
        {
            x[i] = i;
        }


        //for each fitness column do this
        int noColumns = mGData.back().size();
        for (int fitCol = 0; fitCol < noColumns; ++fitCol)
        {
            QVector<double> y(gen);
            for (int i = 0; i < gen; ++i)
            {
                y[i] = mGData.at(i).at(fitCol);
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

        if (gen < 11)
        {
            xAxis->setRange(0, 10);
        }
        else
        {
            xAxis->setRange(0, gen);
        }

        xAxis->setAutoTickStep(false);
        xAxis->setTickStep(1.0);

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
    mGData.clear();

    // give the axes some labels:
    xAxis->setLabel("Generation No.");
    yAxis->setLabel("Fitness");
    // set axes ranges, so we see all data:
    xAxis->setRange(0, 10);
    yAxis->setRange(0, 1);
    xAxis->setAutoTickStep(false);
    xAxis->setTickStep(1.0);

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
