#include "MeshDialog.h"
#include "ui_MeshDialog.h"
#include "OptimisationRun.h"

MeshDialog::MeshDialog(OptimisationRun& data, QWidget *parent) :
	QDialog(parent),
	ui(new Ui::MeshDialog),
	mData(data)
{
	ui->setupUi(this);

	mMeshDensity = mData.meshDensity();

	switch (mMeshDensity)
	{
		case Enum::Mesh::COURSE :
			ui->course->setChecked(true);
			break;
		case Enum::Mesh::MEDIUM :
			ui->medium->setChecked(true);
			break;
		case Enum::Mesh::FINE :
			ui->fine->setChecked(true);
			break;
	}

    ui->viscous->setChecked(mData.getBoundaryLayerFlag());
    ui->layers->setValue(mData.getNumBoundaryLayers());
    ui->thickness->setValue(mData.getBoundaryLayerThickness());
}

MeshDialog::~MeshDialog()
{
	delete ui;
}

void MeshDialog::accept()
{
    mData.setMeshDensity(mMeshDensity);

    bool hasBoundaryLayer = ui->viscous->isChecked();
    mData.setBoundaryLayerFlag(hasBoundaryLayer);

    if(hasBoundaryLayer) {
        mData.setNumBoundaryLayers(ui->layers->value());
        mData.setBoundaryLayerThickness(ui->thickness->value());
    }

	QDialog::accept();
}


void MeshDialog::on_course_toggled(bool checked)
{
	if (checked) mMeshDensity = Enum::Mesh::COURSE;
}

void MeshDialog::on_medium_toggled(bool checked)
{
	if (checked) mMeshDensity = Enum::Mesh::MEDIUM;
}

void MeshDialog::on_fine_toggled(bool checked)
{
	if (checked) mMeshDensity = Enum::Mesh::FINE;
}

void MeshDialog::on_viscous_toggled(bool checked)
{
    ui->layers->setEnabled(checked);
    ui->layers_label->setEnabled(checked);
    ui->thickness->setEnabled(checked);
    ui->thickness_label->setEnabled(checked);
}
