#include "ControlPointView.h"
#include "ui_ControlPointView.h"

ControlPointView::ControlPointView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ControlPointView)
{
    ui->setupUi(this);
}

ControlPointView::~ControlPointView()
{
    delete ui;
}
