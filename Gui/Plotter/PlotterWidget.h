#ifndef PLOTTERWIDGET_H
#define PLOTTERWIDGET_H

#include "qcustomplot.h"

typedef QVector<QVector<double>> GenData;

class PlotterWidget : public QCustomPlot
{
    Q_OBJECT

public:
    explicit PlotterWidget(QCustomPlot *parent = 0);
    void setData(const int gen, const QVector<double> nests);
    void clearData();

private:
    GenData mGData;
};

#endif // PLOTTERWIDGET_H
