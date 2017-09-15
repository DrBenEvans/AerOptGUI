#include "ControlPointView.h"
#include "ui_ControlPointView.h"

ControlPointView::ControlPointView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ControlPointView),
    mBoundaryPointModel(nullptr)
{
    ui->setupUi(this);

    auto valChangedSignal = QOverload<double>::of(&QDoubleSpinBox::valueChanged);
    connect(ui->smoothingValue, valChangedSignal, this, &ControlPointView::smoothingValueChanged);
    connect(ui->xBoundMin, valChangedSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->xBoundMax, valChangedSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->yBoundMin, valChangedSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->yBoundMax, valChangedSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->controlPointCheckBox, &QCheckBox::toggled, this, &ControlPointView::updateModelControlPointState);

    ui->controlPointCheckBox->setVisible(false);

    controlPointParamsVisible(false);
}

void ControlPointView::setModel(BoundaryPointModel *boundaryPointModel) {
    mBoundaryPointModel = boundaryPointModel;
    connect(mBoundaryPointModel, &BoundaryPointModel::activeIndexChanged, this, &ControlPointView::activePointChanged);
    connect(mBoundaryPointModel, &BoundaryPointModel::controlPointStateChanged, this, &ControlPointView::controlPointStateChanged);
    updateViewData();
}

void ControlPointView::setPointCoords(qreal x, qreal y) {
    QString coordText = QString("(%1, %2)");
    coordText = coordText.arg(x).arg(y);
    ui->coord->setText(coordText);
}

void ControlPointView::updateViewData() {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->currentPoint();
    ui->controlPointCheckBox->setVisible(boundaryPoint);
    adjustSize();
    if(boundaryPoint) {
        ui->controlPointCheckBox->setVisible(true);
        qreal x = boundaryPoint->x();
        qreal y = boundaryPoint->y();
        setPointCoords(x, y);

        QRectF rect = boundaryPoint->controlPointRect();
        ui->yBoundMax->setValue(rect.bottom());
        ui->yBoundMin->setValue(rect.top());
        ui->xBoundMax->setValue(rect.right());
        ui->xBoundMin->setValue(rect.left());

        ui->smoothingValue->setValue(boundaryPoint->getSmoothing());
        ui->controlPointCheckBox->setChecked(boundaryPoint->isControlPoint());
        controlPointParamsVisible(boundaryPoint->isControlPoint());
    }
}

void ControlPointView::controlBoundaryChanged() {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->currentPoint();
    if(boundaryPoint) {
        qreal ymax = ui->yBoundMax->value();
        qreal ymin = ui->yBoundMin->value();
        qreal xmax = ui->xBoundMax->value();
        qreal xmin = ui->xBoundMin->value();

        QRectF ctlBounds = QRectF(xmin, ymin, xmax - xmin, ymax - ymin);
        boundaryPoint->setControlPointRect(ctlBounds);
    }
}

void ControlPointView::smoothingValueChanged(double value) {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->currentPoint();
    if(boundaryPoint) {
        boundaryPoint->setSmoothing(value);
    }
}

void ControlPointView::activePointChanged(int index) {
    setVisible(true);
    updateViewData();
}

void ControlPointView::controlPointParamsVisible(bool visible) {
    ui->groupBox->setVisible(visible);
    ui->groupBox->setAttribute(Qt::WA_TransparentForMouseEvents, !visible);
    adjustSize();
}

void ControlPointView::updateModelControlPointState(bool isControlPoint) {
    if(mBoundaryPointModel) {
        int currentIndex = mBoundaryPointModel->currentIndex();
        mBoundaryPointModel->setControlPointState(currentIndex, isControlPoint);
    }
}

// reflects changes in the model to the view
void ControlPointView::controlPointStateChanged(int index, bool isControlPoint) {
    int currentIndex = mBoundaryPointModel->currentIndex();
    if(index == currentIndex) {
        updateViewData();
    }
}

ControlPointView::~ControlPointView()
{
    delete ui;
}
