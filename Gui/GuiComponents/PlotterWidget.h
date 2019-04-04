#ifndef PLOTTERWIDGET_H
#define PLOTTERWIDGET_H

#include "qcustomplot.h"
#include "OptimisationModel.h"

class PlotterWidget : public QCustomPlot
{
    Q_OBJECT
public:
    PlotterWidget(QWidget* parent = 0);
    void updatePlot();
    void clearData();
    void setOptimisationModel(OptimisationModel* model);
    void setCurrentOptimisationIndex(int index);
    void setCurrentlySelectedPoint(int iGen, int agent);
    void clearCurrentSelection();
    int getLineSize();
    int getPointSize();
    void configureView(int lineSize, int pointSize);

signals:
    void selectedPointChanged(int iGen, int agent);

private slots:
    void validateXRangeChanged(const QCPRange& range, const QCPRange &oldRange);
    void fixAxisRange(QCPAxis *axis, QCPRange range, double lowerBound, double upperBound, double minRange);

private:
    double mXMax = 10;
    double mXMin = 0.7;
    double mYMax = 10000;
    double mMinRange = 3;

    // Display size
    int mPointSize = 15;
    int mLineSize = 5;

    void signalCurrentlySelectedPoint();
    Optimisation* currentOptimisation();
    std::pair<int,int> getSelection();
    OptimisationModel* mOptimisationModel;
    QItemSelectionModel* mSelectionModel;
    int mCurrentOptimisationIndex = -1;
    std::vector<QCPGraph*> mGraphs;
};

#endif // PLOTTERWIDGET_H
