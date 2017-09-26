
#include "MeshDialog.h"
#include "ui_MeshDialog.h"
#include "ControlPointView.h"
#include <QGraphicsView>
#include <QWheelEvent>
#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>

Q_DECLARE_METATYPE(Enum::Mesh)

MeshDialog::MeshDialog(ProfileModel &profileModel, QWidget* parent) :
	QDialog(parent),
	ui(new Ui::MeshDialog),
    mMeshDialogModel(new MeshDialogModel(this)),
    mScene(new QGraphicsScene(this)),
    mScale(new ViewScaler()),
    mProfileModel(profileModel),
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

    connect(mMeshDialogModel,&MeshDialogModel::meshChanged,mMeshView, &MeshView::meshChanged);

    // profiles QComboBox setup
    ui->profile->setModel(&mProfileModel);

    ui->density->addItem(QString("Coarse"), QVariant::fromValue(Enum::Mesh::COURSE));
    ui->density->addItem(QString("Medium"), QVariant::fromValue(Enum::Mesh::MEDIUM));
    ui->density->addItem(QString("Fine"), QVariant::fromValue(Enum::Mesh::FINE));

    Mesh* mesh = mMeshDialogModel->currentMesh();
    switch (mesh->meshDensity())
	{
		case Enum::Mesh::COURSE :
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

    // create control point view
    ControlPointView* controlPointView = new ControlPointView(this);
    controlPointView->setModel(mMeshDialogModel->boundaryPointModel());
    controlPointView->move(30,110);

    mScene->addItem(mProfileView);
    mScene->addItem(mMeshView);

    setProfile();
}

void MeshDialog::setProfile() {

    Mesh* mesh = mMeshDialogModel->currentMesh();
    ProfilePoints profilePoints = ui->profile->currentData(Qt::UserRole).value<ProfilePoints>();
    mesh->setProfilePoints(profilePoints);
    setMeshActive(false);
    mProfileView->setProfilePoints(mesh->profilePoints());
    mMeshDialogModel->boundaryPointModel()->clearPoints();
    mMeshView->update();
}

MeshDialog::~MeshDialog()
{
	delete ui;
}

std::vector<BoundaryPoint*> MeshDialog::controlPoints() {
    return mMeshDialogModel->boundaryPointModel()->controlPoints();
}

void MeshDialog::runMesher() {
    Mesh* mesh = mMeshDialogModel->currentMesh();
    mesh->setMeshDensity(ui->density->currentData().value<Enum::Mesh>());
    mesh->setNumberBoundaryLayers(ui->layers->value());
    mesh->setBoundaryLayerThickness(ui->thickness->value());
    mesh->setGrowthFactor(ui->growthFactor->value());

    // set profile
    mMeshDialogModel->runMesher();

    setMeshActive(true);
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
                                      "No control points set. Please select some control points before proceeding",
                                      QMessageBox::Ok);
        return;
    }
	QDialog::accept();
}


void MeshDialog::on_pushButton_clicked()
{
    QString fileName;
    QStringList fileNames;

#ifdef Q_OS_UNIX
    // unix load file
    fileNames.append(
                QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath()+"/Documents/Projects/AerOptProject/", "Profile Files (*.prf)")
                );
#endif
#ifdef Q_OS_WIN32
    // do windows stuff here
    fileNames.append(
                QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath(), "Profile Files (*.prf)")
                );
#endif

    if (fileNames.size() > 0)
    {
        for (const QString &f: fileNames)
        {
            fileName = f;
        }

        qDebug() << "File selected: " << fileName;

        mProfileModel.addProfileFromFilePath(fileName);
    }
    else
    {
        //Set text Not OK here!
        qWarning() << "Profile data not imported!";
    }
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
