#include "PlotterWidget.h"
#include "Optimisation.h"
#include "ViewConfigureDialog.h"

PlotterWidget::PlotterWidget(QWidget *parent) :
    QCustomPlot(parent)
{
    for (uint i = 0; i < 50; ++i) addGraph();

    clearData();

    // set tick
    QSharedPointer<QCPAxisTickerFixed> fixedTicker(new QCPAxisTickerFixed());
    xAxis->setTicker(fixedTicker);

    fixedTicker->setTickStep(1.0);
    fixedTicker->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

    //connect signals
    connect(this, &PlotterWidget::selectionChangedByUser, this, &PlotterWidget::signalCurrentlySelectedPoint);

    auto signal = QOverload<const QCPRange&, const QCPRange&>::of(&QCPAxis::rangeChanged);
    connect(xAxis, signal, this, &PlotterWidget::validateXRangeChanged);
}

void PlotterWidget::validateXRangeChanged(const QCPRange &range, const QCPRange &oldRange) {
    //fixAxisRange(xAxis, range, mXMin, mXMax, mMinRange);
}

void PlotterWidget::fixAxisRange(QCPAxis *axis, QCPRange range, double lowerBound, double upperBound, double minRange)
{
    QCPRange fixedRange(range);
    bool rangeChanged = false;

    if (fixedRange.lower < lowerBound)
    {
        fixedRange.lower = lowerBound;
        fixedRange.upper = lowerBound + range.size();
        if (fixedRange.upper > upperBound || qFuzzyCompare(range.size(), upperBound-lowerBound))
            fixedRange.upper = upperBound;
        if((fixedRange.upper - fixedRange.lower) < minRange)
            fixedRange.upper = fixedRange.lower + minRange;
        rangeChanged = true;
    } else if (fixedRange.upper > upperBound)
    {
        fixedRange.upper = upperBound;
        fixedRange.lower = upperBound - range.size();
        if (fixedRange.lower < lowerBound || qFuzzyCompare(range.size(), upperBound-lowerBound))
            fixedRange.lower = lowerBound;
        if((fixedRange.upper - fixedRange.lower) < minRange)
            fixedRange.lower = fixedRange.upper - minRange;
        if((fixedRange.upper - fixedRange.lower) < minRange)
            fixedRange.lower = fixedRange.upper - minRange;
        rangeChanged = true;
    }

    if(rangeChanged)
        axis->setRange(fixedRange);
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

void PlotterWidget::clearCurrentSelection() {
    QCPDataSelection dataSelection;
    for(int i=0; i<graphCount(); i++) {
        graph(i)->setSelection(dataSelection);
    }

    replot();
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

    emit selectedPointChanged(agent_gen.second, agent_gen.first);
}

void PlotterWidget::setOptimisationModel(OptimisationModel* model) {
    mOptimisationModel = model;
}

Optimisation* PlotterWidget::currentOptimisation() {
    return mOptimisationModel->optimisation(mCurrentOptimisationIndex);
}

void PlotterWidget::setCurrentOptimisationIndex(int index) {
    mCurrentOptimisationIndex = index;

    updatePlot();
}

void PlotterWidget::updatePlot()
{

    //Get data for display
    Optimisation *opt = currentOptimisation();
    std::vector<std::vector<double>> fitness = opt->allfitness();
    uint nGens = fitness.size();
    if(nGens == 0) {
        clearData();
        return;
    }
    int noColumns = fitness.back().size();


    //Compute ranges
    QVector<double> x(nGens);
    for (int i = 0; i < nGens; ++i)
    {
        // generations are 1-indexed in GUI
        x[i] = i + 1;
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

    std::pair<double, double> range = opt->fitnessRange();
    double rangeNorm = range.second - range.first;
    double max = range.second + rangeNorm*0.1;
    double min = range.first - rangeNorm*0.1;

    yAxis->setRange(min, max);


    if(opt) {
        int genRange = opt->noGens() + 1;
        if (genRange <= mMinRange-mXMin)
        {
            xAxis->setRange(mXMin, mMinRange-mXMin);
        }
        else
        {
            xAxis->setRange(mXMin, genRange-mXMin);
        }
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
    xAxis->setRange(mXMin, 10);
    yAxis->setRange(0, 1);

    QColor color = QColor(255, 100, 0);
    QPen pen = QPen(color, mLineSize);

    int num = graphCount();
    for (int i = 0; i < num; ++i)
    {
        QVector<double> x(i);
        QVector<double> y(i);
        graph(i)->setData(x, y);
        graph(i)->setPen(pen);
        graph(i)->setLineStyle(QCPGraph::lsLine);
        graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, mPointSize));
        graph(i)->setVisible(false);
        graph(i)->setSelectable(QCP::stSingleData);
    }

    setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    replot();
}

int PlotterWidget::getLineSize() {
    return mLineSize;
}

int PlotterWidget::getPointSize() {
    return mPointSize;
}

void PlotterWidget::configureView(int lineSize, int pointSize) {

    mLineSize = lineSize;
    mPointSize = pointSize;

    clearData();
    updatePlot();
}
