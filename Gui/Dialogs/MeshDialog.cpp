#include "MeshDialog.h"
#include "ui_MeshDialog.h"
#include "ControlPointView.h"
#include <QGraphicsView>
#include <QWheelEvent>
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>
#include <QSettings>
#include "windowsizing.h"

Q_DECLARE_METATYPE(Enum::Mesh)

MeshDialog::MeshDialog(ProfileModel &profileModel, MeshDialogModel* model, QWidget* parent) :
	QDialog(parent),
	ui(new Ui::MeshDialog),
    mMeshDialogModel(model),
    mScene(new QGraphicsScene(this)),
    mProfileModel(profileModel),
    mScale(new ViewScaler(this)),
    mProfileView(new ProfileView(mScale)),
    mMeshView(new MeshView(mScale))
{
	ui->setupUi(this);

    ui->graphicsView->setSceneRect(mScale->rect(-51,-51,102,102));
    ui->graphicsView->setScene(mScene);
    ui->graphicsView->centerOn(mScale->w(0.5),0);

    mMeshView->setMeshModel(mMeshDialogModel);
    mMeshView->setBoundaryPointModel(mMeshDialogModel->boundaryPointModel());

    connect(ui->meshButton,&QPushButton::clicked,this,&MeshDialog::runMesher);

    connect(ui->toggleProfile,&QCheckBox::toggled,[this](bool checked) {
        mProfileView->setVisible(checked);
    });

    connect(ui->zoomSlider,&QSlider::valueChanged,[this](int scale) {
        ui->graphicsView->setTransform(QTransform::fromScale(float(scale)/100, float(scale)/100));;
    });

    connect(ui->graphicsView, &ZoomPanView::scaleChanged,
            this, &MeshDialog::scaleZoomBarValue);

    connect(&mProfileModel, &ProfileModel::newProfileAdded, [this](int index) {
        ui->profile->setCurrentIndex(index);
        ui->profile->update();
        setProfile();
    });

    // profiles QComboBox setup
    ui->profile->setModel(&mProfileModel);

    ui->density->addItem(QString("Coarse"), QVariant::fromValue(Enum::Mesh::COARSE));
    ui->density->addItem(QString("Medium"), QVariant::fromValue(Enum::Mesh::MEDIUM));
    ui->density->addItem(QString("Fine"), QVariant::fromValue(Enum::Mesh::FINE));

    Mesh* mesh = mMeshDialogModel->currentMesh();
    switch (mesh->meshDensity())
	{
        case Enum::Mesh::COARSE :
            ui->density->setCurrentIndex(0);
			break;
		case Enum::Mesh::MEDIUM :
            ui->density->setCurrentIndex(1);
			break;
		case Enum::Mesh::FINE :
            ui->density->setCurrentIndex(2);
            break;
	}

    ui->layers->setValue(mesh->numberBoundaryLayers());
    ui->thickness->setValue(mesh->boundaryLayerThickness());
    ui->growthFactor->setValue(mesh->growthFactor());

    // create control noded view
    ControlPointView* controlPointView = new ControlPointView(ui->graphicsViewContainer);
    controlPointView->setModel(mMeshDialogModel->boundaryPointModel());
    controlPointView->move(10, 10);

    mScene->addItem(mProfileView);
    mScene->addItem(mMeshView);

    setProfile();

    centerAndResizeWindow(this, 0.8, 0.8);
}


void MeshDialog::scaleZoomBarValue(double value){
    ui->zoomSlider->setValue(ui->zoomSlider->value() * value);
}

void MeshDialog::setProfile() {

    Mesh* mesh = mMeshDialogModel->currentMesh();
    if(mesh) {
        ProfilePoints profilePoints = ui->profile->currentData(Qt::UserRole).value<ProfilePoints>();
        mesh->setProfilePoints(profilePoints);
        setMeshActive(false);
        mProfileView->setProfilePoints(mesh->profilePoints());
        mMeshDialogModel->boundaryPointModel()->clearPoints();
        mMeshView->update();
    }
}

MeshDialog::~MeshDialog()
{
	delete ui;
}

std::vector<BoundaryPoint*> MeshDialog::controlPoints() {
    return mMeshDialogModel->boundaryPointModel()->controlPoints();
}

std::vector<BoundaryPoint*> MeshDialog::boundaryPoints() {
    return mMeshDialogModel->boundaryPointModel()->boundaryPoints();
}

void MeshDialog::runMesher() {
    Mesh* mesh = mMeshDialogModel->currentMesh();
    mesh->setMeshDensity(ui->density->currentData().value<Enum::Mesh>());
    mesh->setNumberBoundaryLayers(ui->layers->value());
    mesh->setBoundaryLayerThickness(ui->thickness->value());
    mesh->setGrowthFactor(ui->growthFactor->value());

    // set profile
    bool success = mMeshDialogModel->runMesher();

    if(success)
        setMeshActive(true);
    else {
        QMessageBox::warning(this, "Warning", "Error generating mesh", QMessageBox::Ok);
    }
}

void MeshDialog::setMeshActive(bool meshIsActive, bool doToggleProfile) {
    if(doToggleProfile) {
        ui->toggleProfile->setChecked(!meshIsActive);
    }
    ui->meshButton->setDefault(!meshIsActive);
}

void MeshDialog::accept()
{
    if(mMeshDialogModel->boundaryPointModel()->controlPointCount() == 0) {
        QMessageBox::warning(this, "Message Box",
                                      "No control nodes set. Please select some control nodes before proceeding",
                                      QMessageBox::Ok);
        return;
    }
	QDialog::accept();
}


void MeshDialog::on_profile_currentIndexChanged(int index)
{
   setProfile();
}

void MeshDialog::on_density_currentIndexChanged(int index)
{
    Mesh* mesh = mMeshDialogModel->currentMesh();
    if(mesh->meshDensity() != ui->density->currentData().value<Enum::Mesh>()) {
        setMeshActive(false,false);
    }
}

void MeshDialog::on_layers_valueChanged(int arg1)
{
    setMeshActive(false,false);
}

void MeshDialog::on_thickness_valueChanged(double arg1)
{
    setMeshActive(false,false);
}

void MeshDialog::on_growthFactor_valueChanged(double arg1)
{
    setMeshActive(false,false);
}

void MeshDialog::on_loadProfileButton_clicked()
{
    QSettings settings;
    QString fileName;
    QStringList fileNames;
    QString profilePath = settings.value("AerOpt/profilesDefaultPath").toString();

    fileNames.append(QFileDialog::getOpenFileName(this, "Select Profile File",
                                                  profilePath, "Profile Files (*.prf)"));

    if (fileNames.size() > 0)
    {
        for (const QString &f: fileNames)
        {
            fileName = f;
        }

        qDebug() << "File selected: " << fileName;

        bool success = mProfileModel.addProfileFromFilePath(fileName);
        if(success) {
            QDir dir(fileName);
            dir.cdUp();
            settings.setValue("AerOpt/profilesDefaultPath", dir.absolutePath());
        }

    }
    else
    {
        //Set text Not OK here!
        qWarning() << "Profile data not imported!";
    }
}
