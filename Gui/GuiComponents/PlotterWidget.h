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

signals:
    void selectedPointChanged(int iGen, int agent);

private:
    void signalCurrentlySelectedPoint();
    std::pair<int,int> getSelection();
    OptimisationModel* mOptimisationModel;
    QItemSelectionModel* mSelectionModel;
    uint mCurrentOptimisationIndex;
    std::vector<QCPGraph*> mGraphs;
};

#endif // PLOTTERWIDGET_H
