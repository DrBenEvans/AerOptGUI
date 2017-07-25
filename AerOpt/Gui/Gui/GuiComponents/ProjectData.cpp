/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ProjectData.cpp
**
**********************************************/

#include <QDebug>
#include <QPointF>
#include <iostream>
#include "BoundaryPoint.h"
#include "ProjectData.h"

ProjectData::ProjectData()
{
	clearProject();
}

ProjectData::~ProjectData()
{

}

void ProjectData::clearProject()
{
	mProfile.clear();
	mMeshProfile.clear();
	mBoundConnects.clear();
	mControlPoints.clear();
	mMeshPoints.clear();
	mMeshConnectivities.clear();
	mMeshResults.clear();

	mProfileSet = false;
	mProjectPathSet = false;
	mMeshSet = false;
	mFunctionSet = false;
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

    mOptimisationMethod = 0;
	mNoAgents = 4;
	mNoGens = 3;
    mNoTop = 75;
}

//Profile Attributes
void ProjectData::clearProfile()
{
	mProfile.clear();
}

void ProjectData::addPoint(const float& x, const float& y)
{
	mProfile.emplace_back(x, y);
}

bool ProjectData::checkProfileIntegrity()
{
	bool r = true;
	bool c = false;

	//Points don't duplicate start and finish vertex.
	c = checkDuplicates();
	r &= c;
//	qDebug() << "Passed duplicate point test? " << c;

	c = mProfile.size() >= 6;
	r &= c;
//	qDebug() << "Passed No. points test? " << c;
//	qDebug() << "No. of profile points: " << mProfile.size();

	//Points in sequence and clockwise (check - it could be anticlockwise)
	c = checkClockwise();
	r &= c;
//	qDebug() << "Passed point sequence test? " << c;

	//Points are scaled to x=1 normalisation x[0,...,1] and y is scaled by same qty
	c = checkNormalised();
	r &= c;
//	qDebug() << "Passed normalised range test? " << c;

	return r;
}

const std::list<std::pair<float, float> >& ProjectData::getProfile() const
{
	return mProfile;
}

//Mesh Attributes
void ProjectData::clearMesh()
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
void ProjectData::addBoundaryPoint(const uint& gen, const float& x, const float& y)
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

void ProjectData::addBConnectivity(const uint& gen, const uint& a, const uint& b)
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

BoundaryPoint* ProjectData::getControlPoint(const uint& genNo, const uint& index)
{
	return mMeshProfile.at(genNo).at(index);
}

bool ProjectData::checkBoundaryIntegrity()
{
	bool r = true;



	return r;
}

const Boundaries& ProjectData::getBoundary() const
{
	return mMeshProfile;
}

const BConnectivities& ProjectData::getBConnects() const
{
	return mBoundConnects;
}

void ProjectData::resetBoundary()
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
void ProjectData::selectControlPoint(const unsigned int& index)
{
	//if control point exists then delete, else add.
	auto it = std::find (mControlPoints.begin(), mControlPoints.end(), index);
	if (it != mControlPoints.end()) mControlPoints.erase(it);
	else mControlPoints.emplace_back(index);
	mControlPoints.unique();
	mControlPoints.sort();
}

bool ProjectData::checkControlPointIntegrity()
{
	bool r = true;

	return r;
}

const std::list<unsigned int>& ProjectData::getControlPoints() const
{
	return mControlPoints;
}

void ProjectData::addMeshPoint(const std::pair<float, float>& pair)
{
	mMeshPoints.emplace_back( pair );
}

void ProjectData::clearMeshPoints()
{
	mMeshPoints.clear();
}

const MeshPoints& ProjectData::getMeshPoints() const
{
	return mMeshPoints;
}

void ProjectData::addMeshConnectivity(const std::tuple<uint, uint, uint>& tuple)
{
	mMeshConnectivities.emplace_back( tuple );
}

void ProjectData::clearMeshConnectivities()
{
	mMeshConnectivities.clear();
}

const MeshConnectivities& ProjectData::getMeshConnectivities() const
{
	return mMeshConnectivities;
}

void ProjectData::addMeshData(const std::tuple<float, float, float, float, float>& tuple)
{
	mMeshResults.emplace_back( tuple );
}

void ProjectData::clearMeshData()
{
	mMeshResults.clear();
}

const MeshResults& ProjectData::getMeshData() const
{
	return mMeshResults;
}

//PRIVATE FUNCTIONS
//Profile Attributes
bool ProjectData::checkDuplicates()
{
	bool r = true;

	//Check for duplicates in curve
	//if duplicates = true then return false
	//Actually just remove all duplicate points.
	bool c = true;
	std::list<std::list<std::pair<float,float>>::iterator> iters;
	for (auto i = mProfile.begin(); i != mProfile.end(); ++i)
	{
		for (auto j = i; j != mProfile.end(); ++j)
		{
			if (i == j) continue;

			c = true;
			c &= i->first == j->first;
			c &= i->second == j->second;

			if (c) iters.push_back(j);
		}
	}

//	qDebug() << "Before";
//	for (auto &p : mProfile) qDebug() << p.first << " : " << p.second;

	for (auto i = iters.rbegin(); i != iters.rend(); ++i)
	{
			mProfile.erase(*i);
	}

//	qDebug() << "After";
//	for (auto &p : mProfile) qDebug() << p.first << " : " << p.second;

	return r;
}

bool ProjectData::checkClockwise()
{
	bool r = true;

	//While x is reducing sum1, and while x is increasing sum2

	float sum1 = 0;
	float sum2 = 0;
	for (auto i = std::next(mProfile.begin()); i != mProfile.end(); ++i)
	{
		if (std::prev(i)->first > i->first)//x is reducing?
		{
			sum1 += i->second;
		}
		else
		{
			sum2 += i->second;
		}
	}

	if (sum1 > sum2) mProfile.reverse();

	r &= sum1 != sum2;

	return r;
}

bool ProjectData::checkNormalised()
{
	bool r = true;

	float min =  1000000;
	float max = -1000000;
	float range = 0;
	float scale = 0;

	for (auto i = mProfile.begin(); i != mProfile.end(); ++i)
	{
		if (i->first < min) min = i->first;
		if (i->first > max) max = i->first;
	}

	range = max - min;
	r &= range != 0;

	if (r)
	{
		scale = 1 / range;
		for (auto i = mProfile.begin(); i != mProfile.end(); ++i)
		{
			i->first *= scale;
			i->second *= scale;
		}
	}

	return r;
}

int ProjectData::getOptimisationMethod() const
{
    return mOptimisationMethod;
}

int ProjectData::getNoTop() const
{
    return mNoTop;
}

void ProjectData::setNoTop(int noTop)
{
    mNoTop = noTop;
}



//Getters and Setters
bool ProjectData::profile() const
{
    return mProfileSet;
}

void ProjectData::setProfile(bool profile)
{
	mProfileSet = profile;
}

bool ProjectData::mesh() const
{
	return mMeshSet;
}

void ProjectData::setMesh(bool mesh)
{
	mMeshSet = mesh;
}

bool ProjectData::function() const
{
	return mFunctionSet;
}

void ProjectData::setFunction(bool function)
{
	mFunctionSet = function;
}

bool ProjectData::boundary() const
{
	return mBoundarySet;
}

void ProjectData::setBoundary(bool boundary)
{
	mBoundarySet = boundary;
}

bool ProjectData::optimiser() const
{
	return mOptimiserSet;
}

void ProjectData::setOptimiser(bool optimiser)
{
	mOptimiserSet = optimiser;
}

void ProjectData::setOptimisationMethod(int method_index)
{
    mOptimisationMethod = method_index;
}

bool ProjectData::runTime() const
{
	return mRunTimeSet;
}

void ProjectData::setRunTime(bool runTime)
{
	mRunTimeSet = runTime;
}

bool ProjectData::renderProfile() const
{
	return mRenderProfile;
}

void ProjectData::setRenderProfile(bool renderProfile)
{
	mRenderProfile = renderProfile;
}

bool ProjectData::renderMesh() const
{
	return mRenderMesh;
}

void ProjectData::setRenderMesh(bool renderMesh)
{
	mRenderMesh = renderMesh;
}

int ProjectData::noGens() const
{
	return mNoGens;
}

void ProjectData::setNoGens(int noGens)
{
	mNoGens = noGens;
}

int ProjectData::noAgents() const
{
	return mNoAgents;
}

void ProjectData::setNoAgents(int noAgents)
{
	mNoAgents = noAgents;
}

float ProjectData::freeTemp() const
{
	return mFreeTemp;
}

void ProjectData::setFreeTemp(float freeTemp)
{
	mFreeTemp = freeTemp;
}

float ProjectData::freePress() const
{
	return mFreePress;
}

void ProjectData::setFreePress(float freePress)
{
	mFreePress = freePress;
}

float ProjectData::freeAlpha() const
{
	return mFreeAlpha;
}

void ProjectData::setFreeAlpha(float freeAlpha)
{
	mFreeAlpha = freeAlpha;
}

float ProjectData::reNo() const
{
	return mReNo;
}

void ProjectData::setReNo(float reNo)
{
	mReNo = reNo;
}

float ProjectData::machNo() const
{
	return mMachNo;
}

void ProjectData::setMachNo(float machNo)
{
	mMachNo = machNo;
}

Enum::ObjFunc ProjectData::objFunc() const
{
	return mObjFunc;
}

void ProjectData::setObjFunc(const Enum::ObjFunc& objFunc)
{
	mObjFunc = objFunc;
}

Enum::Mesh ProjectData::meshDensity() const
{
	return mMeshDensity;
}

void ProjectData::setMeshDensity(const Enum::Mesh& meshDensity)
{
	mMeshDensity = meshDensity;
}

bool ProjectData::projectPathSet() const
{
	return mProjectPathSet;
}

void ProjectData::setProjectPathSet(bool projectPathSet)
{
	mProjectPathSet = projectPathSet;
}
