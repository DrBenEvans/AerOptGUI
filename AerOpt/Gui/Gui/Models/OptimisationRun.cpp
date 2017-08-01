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
#include "OptimisationRun.h"

OptimisationRun::OptimisationRun() : mLabel("Optimisation 1")
{
	clearProject();
}

OptimisationRun::~OptimisationRun()
{

}

void OptimisationRun::clearProject()
{
    mMesh->clear();

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

//Profile Attributes
const std::list<std::pair<float,float>>  OptimisationRun::getProfile() const {
   if(mProfile) {
       mProfile->getProfile();
       auto const& profile = mProfile->getProfile();
       if(profile.size() == 0) {
            return std::list<std::pair<float,float>>();
       } else {
           return profile;
       }
   } else {
       return std::list<std::pair<float,float>>();
   }
}

void OptimisationRun::setProfile(Profile* profile) {
    mProfile = profile;
}

Profile*  OptimisationRun::getProfileObj() const {
    return mProfile;
}


//PRIVATE FUNCTIONS
Enum::OptMethod OptimisationRun::getOptimisationMethod() const
{
    return mOptimisationMethod;
}

int OptimisationRun::getNoTop() const
{
    return mNoTop;
}

void OptimisationRun::setNoTop(int noTop)
{
    mNoTop = noTop;
}

//Getters and Setters
Mesh* OptimisationRun::getMesh()
{
    return mMesh;
}

void OptimisationRun::setOptimisationMethod(Enum::OptMethod method)
{
    mOptimisationMethod = method;
}

int OptimisationRun::noGens() const
{
	return mNoGens;
}

void OptimisationRun::setNoGens(int noGens)
{
	mNoGens = noGens;
}

int OptimisationRun::noAgents() const
{
	return mNoAgents;
}

void OptimisationRun::setNoAgents(int noAgents)
{
	mNoAgents = noAgents;
}

float OptimisationRun::freeTemp() const
{
	return mFreeTemp;
}

void OptimisationRun::setFreeTemp(float freeTemp)
{
	mFreeTemp = freeTemp;
}

float OptimisationRun::freePress() const
{
	return mFreePress;
}

void OptimisationRun::setFreePress(float freePress)
{
	mFreePress = freePress;
}

float OptimisationRun::freeAlpha() const
{
	return mFreeAlpha;
}

void OptimisationRun::setFreeAlpha(float freeAlpha)
{
	mFreeAlpha = freeAlpha;
}

float OptimisationRun::reNo() const
{
	return mReNo;
}

void OptimisationRun::setReNo(float reNo)
{
	mReNo = reNo;
}

float OptimisationRun::machNo() const
{
	return mMachNo;
}

void OptimisationRun::setMachNo(float machNo)
{
	mMachNo = machNo;
}

Enum::ObjFunc OptimisationRun::objFunc() const
{
	return mObjFunc;
}

void OptimisationRun::setObjFunc(const Enum::ObjFunc& objFunc)
{
	mObjFunc = objFunc;
}

QString OptimisationRun::getLabel() {
    return mLabel;
}

void OptimisationRun::setLabel(QString label) {
    mLabel = label;
}
