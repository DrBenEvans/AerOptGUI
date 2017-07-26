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

OptimisationRun::OptimisationRun() : mLabel("Optimisation Label")
{
	clearProject();
}

OptimisationRun::~OptimisationRun()
{

}

void OptimisationRun::clearProject()
{
	mMeshProfile.clear();
	mBoundConnects.clear();
	mControlPoints.clear();
	mMeshPoints.clear();
	mMeshConnectivities.clear();
	mMeshResults.clear();

	mProjectPathSet = false;
	mMeshSet = false;
	mBoundarySet = false;
	mOptimiserSet = false;
	mRunTimeSet = false;
	mRenderProfile = true;
	mRenderMesh = true;

	mMeshProfile = Boundaries(1);
	mBoundConnects = BConnectivities(1);

	mObjFunc = Enum::ObjFunc::LIFTDRAG;
	mMeshDensity = Enum::Mesh::COURSE;

	mMachNo = 0.5;
	mReNo = 0.0;
	mFreeAlpha = 0.0;
	mFreePress = 101325;
	mFreeTemp = 303.15;

    mOptimisationMethod = Enum::OptMethod::MCS;
	mNoAgents = 4;
	mNoGens = 3;
    mNoTop = 75;

    mHasBoundaryLayer = false;
    mNumBoundaryLayers = 10;
    mBoundaryLayerThickness = 0.001;
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

//Mesh Attributes
void OptimisationRun::clearMesh()
{
	mMeshProfile.clear();
	mMeshProfile = Boundaries(1);
	mBoundConnects.clear();
	mBoundConnects = BConnectivities(1);

	mControlPoints.clear();
	mMeshPoints.clear();
	mMeshConnectivities.clear();
	mMeshResults.clear();
}

//mesh boundary points
void OptimisationRun::addBoundaryPoint(const uint& gen, const float& x, const float& y)
{
    BoundaryPoint *point = new BoundaryPoint(x,y,0,0,0,0);

    if ( gen == mMeshProfile.size() )
    {
        mMeshProfile.push_back( std::vector<BoundaryPoint*>() );
    }

    if  ( gen < mMeshProfile.size() )
    {
        mMeshProfile.at(gen).emplace_back(point);
    }
}

void OptimisationRun::addBConnectivity(const uint& gen, const uint& a, const uint& b)
{
	if  ( gen < mBoundConnects.size() )
	{
		mBoundConnects.at(gen).emplace_back( a, b );
	}
	else if ( gen == mBoundConnects.size() )
	{
		mBoundConnects.push_back( std::vector<std::pair<uint, uint>>() );
		mBoundConnects.at(gen).emplace_back( a, b );
	}
}

BoundaryPoint* OptimisationRun::getControlPoint(const uint& genNo, const uint& index)
{
	return mMeshProfile.at(genNo).at(index);
}

bool OptimisationRun::checkBoundaryIntegrity()
{
	bool r = true;



	return r;
}

const Boundaries& OptimisationRun::getBoundary() const
{
	return mMeshProfile;
}

const BConnectivities& OptimisationRun::getBConnects() const
{
	return mBoundConnects;
}

void OptimisationRun::resetBoundary()
{
	//If boundary and bconnects list is greater than 1 then remove all above 1

	if (mMeshProfile.size() > 1)
	{
		mMeshProfile.erase(mMeshProfile.begin()+1, mMeshProfile.end());
	}

	if (mBoundConnects.size() > 1)
	{
		mBoundConnects.erase(mBoundConnects.begin()+1, mBoundConnects.end());
	}
}

//control points
void OptimisationRun::selectControlPoint(const unsigned int& index)
{
	//if control point exists then delete, else add.
	auto it = std::find (mControlPoints.begin(), mControlPoints.end(), index);
	if (it != mControlPoints.end()) mControlPoints.erase(it);
	else mControlPoints.emplace_back(index);
	mControlPoints.unique();
	mControlPoints.sort();
}

bool OptimisationRun::checkControlPointIntegrity()
{
	bool r = true;

	return r;
}

const std::list<unsigned int>& OptimisationRun::getControlPoints() const
{
	return mControlPoints;
}

void OptimisationRun::addMeshPoint(const std::pair<float, float>& pair)
{
	mMeshPoints.emplace_back( pair );
}

void OptimisationRun::clearMeshPoints()
{
	mMeshPoints.clear();
}

const MeshPoints& OptimisationRun::getMeshPoints() const
{
	return mMeshPoints;
}

void OptimisationRun::addMeshConnectivity(const std::tuple<uint, uint, uint>& tuple)
{
	mMeshConnectivities.emplace_back( tuple );
}

void OptimisationRun::clearMeshConnectivities()
{
	mMeshConnectivities.clear();
}

const MeshConnectivities& OptimisationRun::getMeshConnectivities() const
{
	return mMeshConnectivities;
}

void OptimisationRun::addMeshData(const std::tuple<float, float, float, float, float>& tuple)
{
	mMeshResults.emplace_back( tuple );
}

void OptimisationRun::clearMeshData()
{
	mMeshResults.clear();
}

const MeshResults& OptimisationRun::getMeshData() const
{
	return mMeshResults;
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

void OptimisationRun::setBoundaryLayerFlag(bool hasBoundaryLayer) {
    mHasBoundaryLayer = hasBoundaryLayer;
}

void OptimisationRun::setNumBoundaryLayers(int num) {
    mNumBoundaryLayers = num;
}

void OptimisationRun::setBoundaryLayerThickness(qreal thickness) {
    mBoundaryLayerThickness = thickness;
}

bool OptimisationRun::getBoundaryLayerFlag() {
    return mHasBoundaryLayer && mNumBoundaryLayers != 0;
}

int OptimisationRun::getNumBoundaryLayers() {
    if(mHasBoundaryLayer) {
        return mNumBoundaryLayers;
    } else {
        return 0;
    }
}

qreal OptimisationRun::getBoundaryLayerThickness() {
    return mBoundaryLayerThickness;
}



//Getters and Setters
bool OptimisationRun::mesh() const
{
	return mMeshSet;
}

void OptimisationRun::setMesh(bool mesh)
{
	mMeshSet = mesh;
}

bool OptimisationRun::boundary() const
{
	return mBoundarySet;
}

void OptimisationRun::setBoundary(bool boundary)
{
	mBoundarySet = boundary;
}

bool OptimisationRun::optimiser() const
{
	return mOptimiserSet;
}

void OptimisationRun::setOptimiser(bool optimiser)
{
	mOptimiserSet = optimiser;
}

void OptimisationRun::setOptimisationMethod(Enum::OptMethod method)
{
    mOptimisationMethod = method;
}

bool OptimisationRun::runTime() const
{
	return mRunTimeSet;
}

void OptimisationRun::setRunTime(bool runTime)
{
	mRunTimeSet = runTime;
}

bool OptimisationRun::renderProfile() const
{
	return mRenderProfile;
}

void OptimisationRun::setRenderProfile(bool renderProfile)
{
	mRenderProfile = renderProfile;
}

bool OptimisationRun::renderMesh() const
{
	return mRenderMesh;
}

void OptimisationRun::setRenderMesh(bool renderMesh)
{
	mRenderMesh = renderMesh;
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

Enum::Mesh OptimisationRun::meshDensity() const
{
	return mMeshDensity;
}

void OptimisationRun::setMeshDensity(const Enum::Mesh& meshDensity)
{
	mMeshDensity = meshDensity;
}

bool OptimisationRun::projectPathSet() const
{
	return mProjectPathSet;
}

void OptimisationRun::setProjectPathSet(bool projectPathSet)
{
	mProjectPathSet = projectPathSet;
}

QString OptimisationRun::getLabel() {
    return mLabel;
}

void OptimisationRun::setLabel(QString label) {
    mLabel = label;
}
