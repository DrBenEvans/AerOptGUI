#include "Profile.h"
#include <fstream>
#include <QDebug>
#include "DebugOutput.h"

Profile::Profile(QObject *parent, QString filePath) :
    QObject(parent),
    mFilePath(filePath)
{
    bool r = true;
    float x, y;

    qInfo() << "Loading profile: " << mFilePath;
    std::ifstream infile(mFilePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();

    if (r)
    {
        clearProfile();
        while (infile >> x >> y)
        {
            addPoint(x, y);
        }
    }
    infile.close();

    r &= checkProfileIntegrity();

    if (!r) clearProfile();
}


QString Profile::getDisplayString() {
    return mFilePath;
}

void Profile::addPoint(const float& x, const float& y)
{
    mProfile.emplace_back(x, y);
}

bool Profile::checkProfileIntegrity()
{
    bool r = true;
    bool c = false;

    //Points don't duplicate start and finish vertex.
    c = checkDuplicates();
    r &= c;
    qDebug() << "Passed duplicate point test? " << c;

    c = mProfile.size() >= 3;
    r &= c;
    qDebug() << "Passed No. points test? " << c;
    qDebug() << "No. of profile points: " << mProfile.size();

    //Points in sequence and clockwise (check - it could be anticlockwise)
    c = checkClockwise();
    r &= c;
    qDebug() << "Passed point sequence test? " << c;

    //Points are scaled to x=1 normalisation x[0,...,1] and y is scaled by same qty
    c = checkNormalised();
    r &= c;
    qDebug() << "Passed normalised range test? " << c;

    return r;
}

void Profile::clearProfile()
{
    mProfile.clear();
}

const std::list<std::pair<float, float> > Profile::getProfile() const
{
    return mProfile;
}

//Profile Attributes
bool Profile::checkDuplicates()
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

bool Profile::checkClockwise()
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

bool Profile::checkNormalised()
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
