#include "MeshDialog.h"
#include "ui_MeshDialog.h"
#include "ProjectData.h"

MeshDialog::MeshDialog(ProjectData& data, QWidget *parent) :
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
}

MeshDialog::~MeshDialog()
{
	delete ui;
}

void MeshDialog::accept()
{
	mData.setMeshDensity(mMeshDensity);

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
