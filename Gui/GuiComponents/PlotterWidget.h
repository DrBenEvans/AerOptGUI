#ifndef PLOTTERWIDGET_H
#define PLOTTERWIDGET_H

#include "qcustomplot.h"
#include "OptimisationModel.h"

class PlotterWidget : public QCustomPlot
{
public:
    PlotterWidget(QWidget* parent = 0);
    void updatePlot();
    void clearData();
    void setOptimisationModel(OptimisationModel* model);
    void setCurrentOptimisationIndex(int index);

public slots:
    void currentOptimisationChanged(QModelIndex current, QModelIndex prev);

private:
    OptimisationModel* mOptimisationModel;
    QItemSelectionModel* mSelectionModel;
    uint mCurrentOptimisationIndex;
};

#endif // PLOTTERWIDGET_H
