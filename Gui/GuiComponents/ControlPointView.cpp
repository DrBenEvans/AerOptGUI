#include "ControlPointView.h"
#include "ui_ControlPointView.h"

ControlPointView::ControlPointView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ControlPointView),
    mBoundaryPointModel(nullptr)
{
    ui->setupUi(this);

    auto dblSpinBoxSignal = QOverload<double>::of(&QDoubleSpinBox::valueChanged);
    auto spinBoxSignal = QOverload<int>::of(&QSpinBox::valueChanged);
    connect(ui->smoothingValue, spinBoxSignal, this, &ControlPointView::smoothingValueChanged);
    connect(ui->xBoundMin, dblSpinBoxSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->xBoundMax, dblSpinBoxSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->yBoundMin, dblSpinBoxSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->yBoundMax, dblSpinBoxSignal, this, &ControlPointView::controlBoundaryChanged);
    connect(ui->controlPointCheckBox, &QCheckBox::toggled, this, &ControlPointView::updateModelControlPointState);

    mInitialControlPointText = ui->coord->text();

    resetView();
}

void ControlPointView::resetView() {
    ui->controlPointCheckBox->setVisible(false);
    ui->coord->setText(mInitialControlPointText);

    controlPointParamsVisible(false);
}

void ControlPointView::setModel(BoundaryPointModel *boundaryPointModel) {
    mBoundaryPointModel = boundaryPointModel;
    connect(mBoundaryPointModel, &BoundaryPointModel::activeIndexChanged, this, &ControlPointView::activePointChanged);
    connect(mBoundaryPointModel, &BoundaryPointModel::controlPointStateChanged, this, &ControlPointView::controlPointStateChanged);
    connect(mBoundaryPointModel, &BoundaryPointModel::controlPointBoundsChanged, this, &ControlPointView::controlPointBoundsChanged);
    connect(mBoundaryPointModel, &BoundaryPointModel::boundaryPointsReset, this, &ControlPointView::resetView);
}

void ControlPointView::setPointCoords(qreal x, qreal y) {
    QString coordText = QString("(%1, %2)");
    coordText = coordText.arg(x).arg(y);
    ui->coord->setText(coordText);
}

void ControlPointView::updateViewData() {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->point(mBoundaryPointIndex);
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

// update changes to the model
void ControlPointView::controlBoundaryChanged() {
    qreal ymax = ui->yBoundMax->value();
    qreal ymin = ui->yBoundMin->value();
    qreal xmax = ui->xBoundMax->value();
    qreal xmin = ui->xBoundMin->value();

    QRectF ctlBounds = QRectF(xmin, ymin, xmax - xmin, ymax - ymin);
    mBoundaryPointModel->setControlPointBounds(mBoundaryPointIndex, ctlBounds);
}

void ControlPointView::smoothingValueChanged(int value) {
    BoundaryPoint* boundaryPoint = mBoundaryPointModel->point(mBoundaryPointIndex);
    if(boundaryPoint) {
        boundaryPoint->setSmoothing(value);
    }
}

void ControlPointView::activePointChanged(int index) {
    mBoundaryPointIndex = index;
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
        mBoundaryPointModel->setControlPointState(mBoundaryPointIndex, isControlPoint);
    }
}

// reflects changes in the model to the view
void ControlPointView::controlPointBoundsChanged(int index) {
    if(index == mBoundaryPointIndex) {
        updateViewData();
    }
}

void ControlPointView::controlPointStateChanged(int index, bool isControlPoint) {
    if(index == mBoundaryPointIndex) {
        updateViewData();
    }
}

ControlPointView::~ControlPointView()
{
    delete ui;
}
