/*********************************************
**
**	Created on: 	08/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		Menu.h
**
**********************************************/

#ifndef MENU_H
#define MENU_H

#include <QObject>
#include <QMenu>
#include <map>
#include <utility>
#include "Enumerations.h"

class Canvas;
class OptimisationRun;

class Menu : public QObject
{
	Q_OBJECT
public:
	explicit Menu(OptimisationRun& profile, Canvas& canvas, Enum::TreeType type, QObject* parent = 0);
	Menu() = delete;
	~Menu();
	void exec(const QPoint& pos);
	void exec(const QPoint& pos, const unsigned int& index);

signals:
	
public slots:
	void setConstraints();
	void resetConstraints();

protected:

private:
	void root(QObject *parent);
	void profile(QObject *parent);
	void mesh(QObject *parent);
	void boundary(QObject *parent);
	void optimiser(QObject *parent);
    void runtime(QObject *);
    void other(QObject *parent);
	void node(QObject *parent);


	QList<QAction*> mActions;
	unsigned int mIndex;
	Canvas& mCanvas;
	OptimisationRun& mData;
};

#endif // MENU_H
