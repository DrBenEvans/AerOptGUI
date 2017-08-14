#include <QDebug>
#include "MeshDialog.h"
#include "ui_MeshDialog.h"
#include <QGraphicsView>
#include <QWheelEvent>
#include <QFileDialog>
#include <ControlPointGraphicsItem.h>

Q_DECLARE_METATYPE(Enum::Mesh)

MeshDialog::MeshDialog(std::shared_ptr<Mesh> initMesh, ProfileModel &profileModel, QWidget* parent) :
	QDialog(parent),
	ui(new Ui::MeshDialog),
    mMesh(initMesh),
    mScene(new QGraphicsScene(this)),
    mProfileModel(profileModel),
    mProfileGraphicsItem(new ProfileGraphicsItem(mScale)),
    mMeshGraphicsItem(new MeshGraphicsItem(mScale))
{
	ui->setupUi(this);

    ui->graphicsView->setSceneRect(QRectF(-51*mScale,-51*mScale,102*mScale,102*mScale));
    ui->graphicsView->setScene(mScene);
    ui->graphicsView->centerOn(0.5*mScale,0);

    connect(ui->meshButton,&QPushButton::clicked,this,&MeshDialog::runMesher);

    connect(ui->toggleProfile,&QCheckBox::toggled,[this](bool checked) {
        mProfileGraphicsItem->setVisible(checked);
    });

    connect(mMesh.get(),&Mesh::meshChanged,[this]() {
        mMeshGraphicsItem->meshChanged();
    });

    // Stacked Widget Controls
    connect(ui->nextButton,&QPushButton::clicked,[this]() {
        controlPointStackedWidget();
    });
    connect(ui->prevButton,&QPushButton::clicked,[this]() {
        meshStackedWidget();
    });
    ui->stackedWidget->setCurrentIndex(0);

    // profiles QComboBox setup
    ui->profile->setModel(&mProfileModel);

    ui->density->addItem(QString("Coarse"), QVariant::fromValue(Enum::Mesh::COURSE));
    ui->density->addItem(QString("Medium"), QVariant::fromValue(Enum::Mesh::MEDIUM));
    ui->density->addItem(QString("Fine"), QVariant::fromValue(Enum::Mesh::FINE));

    switch (mMesh->getMeshDensity())
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

    ui->layers->setValue(mMesh->getNumBoundaryLayers());
    ui->thickness->setValue(mMesh->getBoundaryLayerThickness());

    mScene->addItem(mProfileGraphicsItem);
    mScene->addItem(mMeshGraphicsItem);

    setProfile();
}

void MeshDialog::setProfile() {

    ProfilePoints profilePoints = ui->profile->currentData(Qt::UserRole).value<ProfilePoints>();
    mMesh->setProfilePoints(profilePoints);
    setMeshActive(false);
    mProfileGraphicsItem->setProfilePoints(mMesh->profilePoints());
}

MeshDialog::~MeshDialog()
{
	delete ui;
}

void MeshDialog::runMesher() {
    mMesh->setMeshDensity(ui->density->currentData().value<Enum::Mesh>());
    mMesh->setNumBoundaryLayers(ui->layers->value());
    mMesh->setBoundaryLayerThickness(ui->thickness->value());

    // set profile
    mMesh->runMesher();

    mMeshGraphicsItem->setMesh(mMesh);
    setMeshActive(true);
}

void MeshDialog::setMeshActive(bool meshIsActive, bool doToggleProfile) {
    if(doToggleProfile) {
        ui->toggleProfile->setChecked(!meshIsActive);
    }
    ui->nextButton->setDefault(meshIsActive);
    ui->meshButton->setDefault(!meshIsActive);
    ui->nextButton->setDisabled(!meshIsActive);
}

void MeshDialog::controlPointStackedWidget() {
    ui->stackedWidget->setCurrentIndex(1);
    mMeshGraphicsItem->setOpacity(0.3);
    mProfileGraphicsItem->setOpacity(0.3);
    mMeshGraphicsItem->showControlPoints(true);
}

void MeshDialog::meshStackedWidget() {
    ui->stackedWidget->setCurrentIndex(0);
    mMeshGraphicsItem->setOpacity(1);
    mProfileGraphicsItem->setOpacity(1);
    mMeshGraphicsItem->showControlPoints(false);
}

void MeshDialog::accept()
{
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
    if(mMesh->getMeshDensity() != ui->density->currentData().value<Enum::Mesh>()) {
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
