#include "MeshDialog.h"
#include "ui_MeshDialog.h"

Q_DECLARE_METATYPE(Enum::Mesh)

MeshDialog::MeshDialog(Mesh& mesh, QWidget* parent) :
	QDialog(parent),
	ui(new Ui::MeshDialog),
    mMesh(mesh)
{
	ui->setupUi(this);

    ui->density->addItem(QString("Coarse"), QVariant::fromValue(Enum::Mesh::COURSE));
    ui->density->addItem(QString("Medium"), QVariant::fromValue(Enum::Mesh::MEDIUM));
    ui->density->addItem(QString("Fine"), QVariant::fromValue(Enum::Mesh::FINE));

    switch (mMesh.getMeshDensity())
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

    ui->viscous->setChecked(mMesh.getBoundaryLayerFlag());
    ui->layers->setValue(mMesh.getNumBoundaryLayers());
    ui->thickness->setValue(mMesh.getBoundaryLayerThickness());
}

MeshDialog::~MeshDialog()
{
	delete ui;
}

void MeshDialog::accept()
{
    mMesh.setMeshDensity(ui->density->currentData().value<Enum::Mesh>());

    bool hasBoundaryLayer = ui->viscous->isChecked();
    mMesh.setBoundaryLayerFlag(hasBoundaryLayer);

    if(hasBoundaryLayer) {
        mMesh.setNumBoundaryLayers(ui->layers->value());
        mMesh.setBoundaryLayerThickness(ui->thickness->value());
    }

	QDialog::accept();
}

void MeshDialog::on_viscous_toggled(bool checked)
{
    ui->layers->setEnabled(checked);
    ui->layers_label->setEnabled(checked);
    ui->thickness->setEnabled(checked);
    ui->thickness_label->setEnabled(checked);
}
