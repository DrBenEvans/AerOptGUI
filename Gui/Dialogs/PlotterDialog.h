#ifndef PLOTTERDIALOG_H
#define PLOTTERDIALOG_H

#include <QDialog>

typedef QVector<QVector<double>> GenData;

namespace Ui {
class PlotterDialog;
}

class PlotterDialog : public QDialog
{
	Q_OBJECT

public:
    explicit PlotterDialog(QWidget *parent = 0);
	~PlotterDialog();
	void setData(const int gen, const QVector<double> nests);
	void clearData();

private:
	Ui::PlotterDialog *ui;

	GenData mGData;
};

#endif // PLOTTERDIALOG_H
