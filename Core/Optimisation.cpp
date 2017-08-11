/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimisationRun.cpp
**
**********************************************/

#include <QDebug>
#include <QPointF>
#include <iostream>
#include "BoundaryPoint.h"
#include "Optimisation.h"

Optimisation::Optimisation() :
    mInitMesh(new Mesh())
{
    mObjFunc = Enum::ObjFunc::LIFTDRAG;

    mMachNo = 0.5;
    mReNo = 0.0;
    mFreeAlpha = 0.0;
    mFreePress = 101325;
    mFreeTemp = 303.15;

    mOptimisationMethod = Enum::OptMethod::MCS;
    mNoAgents = 4;
    mNoGens = 3;
    mNoTop = 75;
}

Optimisation::~Optimisation()
{

}

//PRIVATE FUNCTIONS
Enum::OptMethod Optimisation::getOptimisationMethod() const
{
    return mOptimisationMethod;
}

int Optimisation::getNoTop() const
{
    return mNoTop;
}

void Optimisation::setNoTop(int noTop)
{
    mNoTop = noTop;
}

//Getters and Setters
std::shared_ptr<Mesh> Optimisation::initMesh()
{
    return mInitMesh;
}

void Optimisation::setOptimisationMethod(Enum::OptMethod method)
{
    mOptimisationMethod = method;
}

int Optimisation::noGens() const
{
	return mNoGens;
}

void Optimisation::setNoGens(int noGens)
{
	mNoGens = noGens;
}

int Optimisation::noAgents() const
{
	return mNoAgents;
}

void Optimisation::setNoAgents(int noAgents)
{
	mNoAgents = noAgents;
}

float Optimisation::freeTemp() const
{
	return mFreeTemp;
}

void Optimisation::setFreeTemp(float freeTemp)
{
	mFreeTemp = freeTemp;
}

float Optimisation::freePress() const
{
	return mFreePress;
}

void Optimisation::setFreePress(float freePress)
{
	mFreePress = freePress;
}

float Optimisation::freeAlpha() const
{
	return mFreeAlpha;
}

void Optimisation::setFreeAlpha(float freeAlpha)
{
	mFreeAlpha = freeAlpha;
}

float Optimisation::reNo() const
{
	return mReNo;
}

void Optimisation::setReNo(float reNo)
{
	mReNo = reNo;
}

float Optimisation::machNo() const
{
	return mMachNo;
}

void Optimisation::setMachNo(float machNo)
{
	mMachNo = machNo;
}

Enum::ObjFunc Optimisation::objFunc() const
{
	return mObjFunc;
}

void Optimisation::setObjFunc(const Enum::ObjFunc& objFunc)
{
	mObjFunc = objFunc;
}

QString Optimisation::label() const {
    return mLabel;
}

void Optimisation::setLabel(QString label) {
    mLabel = label;
}
