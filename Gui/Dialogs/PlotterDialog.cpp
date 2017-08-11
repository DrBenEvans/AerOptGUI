#include "PlotterDialog.h"
#include "ui_PlotterDialog.h"
#include "Optimisation.h"

PlotterDialog::PlotterDialog(QWidget *parent) :
    QDialog(parent),
	ui(new Ui::PlotterDialog)
{
	ui->setupUi(this);

	for (uint i = 0; i < 50; ++i) ui->widget->addGraph();

	clearData();
}

PlotterDialog::~PlotterDialog()
{
	delete ui;
}

void PlotterDialog::setData(const int gen, const QVector<double> nests)
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
			ui->widget->graph(fitCol)->setData(x, y);
		//end do this
		}


		// give the axes some labels:
		ui->widget->xAxis->setLabel("Generation No.");
		ui->widget->yAxis->setLabel("Fitness");
		// set axes ranges, so we see all data:
		ui->widget->rescaleAxes(true);

		if (gen < 11)
		{
			ui->widget->xAxis->setRange(0, 10);
		}
		else
		{
			ui->widget->xAxis->setRange(0, gen);
		}

		ui->widget->xAxis->setAutoTickStep(false);
		ui->widget->xAxis->setTickStep(1.0);

		double min, max;
		max = ui->widget->yAxis->range().upper * 1.05;
		min = ui->widget->yAxis->range().lower * 0.95;

		ui->widget->yAxis->setRange(min, max);
	}

	//Plot graph
	ui->widget->replot();
}

void PlotterDialog::clearData()
{
	mGData.clear();

	// give the axes some labels:
	ui->widget->xAxis->setLabel("Generation No.");
	ui->widget->yAxis->setLabel("Fitness");
	// set axes ranges, so we see all data:
	ui->widget->xAxis->setRange(0, 10);
	ui->widget->yAxis->setRange(0, 1);
	ui->widget->xAxis->setAutoTickStep(false);
	ui->widget->xAxis->setTickStep(1.0);

	int num = ui->widget->graphCount();
	for (int i = 0; i < num; ++i)
	{
		QVector<double> x(i);
		QVector<double> y(i);
		ui->widget->graph(i)->setData(x, y);
		ui->widget->graph(i)->setPen(QPen(QColor(255, 100, 0)));
		ui->widget->graph(i)->setLineStyle(QCPGraph::lsLine);
		ui->widget->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
	}

	ui->widget->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

	ui->widget->replot();
}
