/*********************************************
**
**	Created on: 	08/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		Menu.cpp
**
**********************************************/

#include "Menu.h"
#include "OptimisationRun.h"
#include "Canvas.h"

Menu::Menu(OptimisationRun& data, Canvas& canvas, Enum::TreeType type, QObject *parent) : QObject(parent), mCanvas(canvas), mData(data)
{
	mIndex = 0;
	switch (type)
	{
		case Enum::TreeType::ROOT:
			root(parent);
			break;
		case Enum::TreeType::DATA:
			profile(parent);
			break;
		case Enum::TreeType::MESH:
			mesh(parent);
			break;
		case Enum::TreeType::BOUNDARY:
			boundary(parent);
			break;
		case Enum::TreeType::OPTIMISER:
			optimiser(parent);
			break;
		case Enum::TreeType::RUNTIME:
			runtime(parent);
			break;
		case Enum::TreeType::NODE:
			node(parent);
			break;
		default:
			other(parent);
	}
}

Menu::~Menu()
{
	for (int i = 0; i < mActions.size(); ++i)
	{
		if (mActions.at(i)) delete mActions.at(i);
	}
}

void Menu::exec(const QPoint &pos)
{
	QMenu::exec(mActions, pos);
}

void Menu::exec(const QPoint& pos, const unsigned int& index)
{
	mIndex = index;
	QMenu::exec(mActions, pos);
}

void Menu::root(QObject *parent)
{
	QAction* loadObject = new QAction("Set Project Directory", parent);
	QAction* clear = new QAction("Clear Project", parent);
	//TODO Save project here!
	//To solve end opt thread issue with files.

	mActions.append(loadObject);
	mActions.append(clear);

	QObject::connect(loadObject, SIGNAL( triggered() ), this->parent(), SLOT( loadProject() ));
	QObject::connect(clear, SIGNAL( triggered() ), this->parent(), SLOT( clearProject() ));
}

void Menu::profile(QObject *parent)
{
	QAction* importObject = new QAction("Import Profile", parent);
	QAction* showObject = new QAction("Show profile", parent);
	QAction* hideObject = new QAction("Hide profile", parent);

	mActions.append(importObject);
	mActions.append(showObject);
	mActions.append(hideObject);

    //QObject::connect(importObject, SIGNAL( triggered() ), this->parent(), SLOT( importProfile() ));
	QObject::connect(showObject, SIGNAL( triggered() ), this, SLOT( showProfile() ));
    //QObject::connect(hideObject, SIGNAL( triggered() ), this, SLOT( hideProfile() ));
}

void Menu::mesh(QObject *parent)
{

	QAction* runmeshObject = new QAction("Run Mesher", parent);
	QAction* showObject = new QAction("Show mesh", parent);
	QAction* hideObject = new QAction("Hide mesh", parent);

	mActions.append(runmeshObject);
	mActions.append(showObject);
	mActions.append(hideObject);

    //QObject::connect(runmeshObject, SIGNAL( triggered() ), this->parent(), SLOT( runMesher() ));
    //QObject::connect(showObject, SIGNAL( triggered() ), this, SLOT( showMesh() ));
    //QObject::connect(hideObject, SIGNAL( triggered() ), this, SLOT( hideMesh() ));
}

void Menu::boundary(QObject *parent)
{
    QAction* Object = new QAction("Set Flow Conditions", parent);
	mActions.append(Object);
    //QObject::connect(Object, SIGNAL( triggered() ), this->parent(), SLOT( setBoundary() ));
}

void Menu::optimiser(QObject *parent)
{
	QAction* Object = new QAction("Set Optimiser parameters", parent);
	mActions.append(Object);
    //QObject::connect(Object, SIGNAL( triggered() ), this->parent(), SLOT( setOptimiser() ));
}

void Menu::runtime(QObject *parent)
{
	QAction* runObject = new QAction("Run AerOpt", parent);
	QAction* showObject = new QAction("Show Graph", parent);
	mActions.append(runObject);
	mActions.append(showObject);
	QObject::connect(runObject, SIGNAL( triggered() ), this->parent(), SLOT( runAerOpt() ));
	QObject::connect(showObject, SIGNAL( triggered() ), this->parent(), SLOT( showGraph() ));
}

void Menu::other(QObject *parent)
{

}

void Menu::node(QObject *parent)
{
	QAction* setObject = new QAction("Set Constraints", parent);
	QAction* unsetObject = new QAction("Reset Constraints", parent);

	mActions.append(setObject);
	mActions.append(unsetObject);

	QObject::connect(setObject, SIGNAL( triggered() ), this, SLOT( setConstraints() ));
	QObject::connect(unsetObject, SIGNAL( triggered() ), this, SLOT( resetConstraints() ));
}

void Menu::setConstraints()
{
	mCanvas.setConstraints(mIndex);
}

void Menu::resetConstraints()
{
	mCanvas.resetConstraints(mIndex);
}
