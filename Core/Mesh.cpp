#include "Mesh.h"
#include <fstream>
#include <cmath>
#include <QDebug>

Mesh::Mesh(QObject *parent): QObject(parent),
    mMeshDensity(Enum::Mesh::COURSE),
    mProfilePoints()
{
    mMeshProcess.setParent(this);

    // Boundary layer setup
    mNumBoundaryLayers = 10;
    mBoundaryLayerThickness = 0.001;
}

ProfilePoints Mesh::profilePoints() {
    return mProfilePoints;
}

void Mesh::setProfilePoints(ProfilePoints profilePoints) {
    clear();
    mProfilePoints = profilePoints;
    emit meshChanged();
}

void Mesh::runMesher() {
    clear();
    QDir workDir("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/Project");

    bool r = true;

    mMeshPath = workDir;

    QString meshInFile;
    QString meshBacFile;
    QString meshGeoFile;

    QString mesherPath = "/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/FLITE/Mesher/MeshGenerator";

    meshInFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".in").absolutePath();
    meshInFile = QDir::toNativeSeparators(meshInFile);

    meshBacFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".bac").absolutePath();
    meshBacFile = QDir::toNativeSeparators(meshBacFile);

    meshGeoFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".geo").absolutePath();
    meshGeoFile = QDir::toNativeSeparators(meshGeoFile);

    mMeshDatFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".dat").absolutePath();
    mMeshDatFile = QDir::toNativeSeparators(mMeshDatFile);

    r &= createInputFile(meshInFile.toStdString(),
                         meshBacFile.toStdString(),
                         meshGeoFile.toStdString(),
                         mMeshDatFile.toStdString());
    r &= QFile::exists(meshInFile);

    r &= createBacFile(meshBacFile.toStdString());
    r &= QFile::exists(meshBacFile);

    r &= createGeoFile(meshGeoFile.toStdString());
    r &= QFile::exists(meshGeoFile);

    if (r)
    {
        mMeshProcess.setWorkingDirectory( mMeshPath.absolutePath() );
        mMeshProcess.setStandardInputFile( meshInFile );
        mMeshProcess.start( mesherPath );
        mMeshProcess.waitForFinished();
        meshingFinished(mMeshProcess.exitCode(),mMeshProcess.exitStatus());
    }
}

bool Mesh::createBacFile(const std::string& meshBacFile)
{
    bool r = true;

    std::ofstream outfile(meshBacFile, std::ofstream::out);
    r &= outfile.is_open();

    if (r)
    {
        outfile << "title" << std::endl;
        outfile << "   4    2" << std::endl;
        outfile << "   1 -0.600081E+02 -0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "   2  0.600081E+02 -0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "   3  0.600081E+02  0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "   4 -0.600081E+02  0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
        outfile << "  1    1   2   3" << std::endl;
        outfile << "  2    1   3   4" << std::endl;
        outfile << "  sources" << std::endl;
        outfile << "  0  1	 0" << std::endl;
        outfile << "  point" << std::endl;
        outfile << "  line" << std::endl;
        outfile << "InnerRadius" << std::endl;
        switch (getMeshDensity())
        {
            case Enum::Mesh::COURSE :
                outfile << "    0.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
                outfile << "    1.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
                break;
            case Enum::Mesh::MEDIUM :
                outfile << "    0.0  0.0   0.07  0.4  0.6" << std::endl;//<< medium
                outfile << "    1.0  0.0   0.07  0.4  0.6" << std::endl;//<< medium
                break;
            case Enum::Mesh::FINE :
                outfile << "    0.0  0.0   0.02  0.3  0.7" << std::endl;//<< fine
                outfile << "    1.0  0.0   0.02  0.3  0.7" << std::endl;//<< fine
                break;
            default :
                outfile << "    0.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
                outfile << "    1.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
        }
        outfile << " 0 " << std::endl;
        outfile << " 0 " << std::endl;
    }
    outfile.close();

    return r;
}

bool Mesh::createGeoFile(const std::string& meshGeoFile)
{
    bool r = true;
    int noPoints = 0;
    int noSegments = 6;
    int noDomainPoints = 4;
    int nViscSegments = 2;

    std::ofstream outfile(meshGeoFile, std::ofstream::out);
    r &= outfile.is_open();

    noPoints = mProfilePoints.size();
    int nLayers = getNumBoundaryLayers();

    if (r)
    {
        outfile << "npoin nseg nvseg nlayer hmin" << std::endl;
        outfile << int(noPoints + noDomainPoints) << "	  " << noSegments <<  "	  " << nViscSegments << "	  " << nLayers << "   0.000" << std::endl;

        outfile.precision(7);
//		outfile << std::scientific;
        outfile << std::fixed;

        int i = 1;
        for (auto &p : mProfilePoints)
        {
            outfile << i << "	  " << p.first << "	   " << p.second << std::endl;
            ++i;
        }

        outfile << i << " 	 -50.000000	 -50.000000" << std::endl;
        ++i;
        outfile << i << " 	  50.000000	 -50.000000" << std::endl;
        ++i;
        outfile << i << " 	  50.000000	  50.000000" << std::endl;
        ++i;
        outfile << i << " 	 -50.000000	  50.000000" << std::endl;


        int seg = 1;
        int prof1 = std::round(noPoints*0.5);
        int prof2 = noPoints - prof1 + 2;

        outfile << seg << " " << prof1 << " " << 1 << std::endl;
        for (i = 1; i <= prof1; ++i)
        {
            outfile << i << std::endl;
        }
        ++seg;


        outfile << seg << " " << prof2 << " " << 1 << std::endl;
        for (i = prof1; i <= noPoints; ++i)
        {
            outfile << i << std::endl;
        }
        outfile << 1 << std::endl;
        ++seg;


        outfile << seg << " " << 2 << " " << 3 << std::endl;
        outfile << i << " " << i+1 << std::endl;
        ++seg;
        ++i;

        outfile << seg << " " << 2 << " " << 3 << std::endl;
        outfile << i << " " << i+1 << std::endl;
        ++seg;
        ++i;

        outfile << seg << " " << 2 << " " << 3 << std::endl;
        outfile << i << " " << i+1 << std::endl;
        ++seg;
        ++i;

        outfile << seg << " " << 2 << " " << 3 << std::endl;
        outfile << i << " " << i-3 << std::endl;
    }

    qreal thickness = getBoundaryLayerThickness();
    qreal dist = 0;
    for(int i=0; i<nLayers; i++) {
        dist += thickness;
        outfile << dist << std::endl;
    }

    outfile.close();

    return r;
}

bool Mesh::loadMeshProfile(const QString &filePath)
{
    bool r = true;
    int type = 1;
    std::string word1 = "";

    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        while (infile >> word1)
        {
            if (word1 == "clean")
            {
                type = 2;
                break;
            }
        }

        infile.close();

        if (type == 1)
        {
            r &= loadMeshProfileType1(filePath);
        }
        else if (type == 2)
        {
            r &= loadMeshProfileType2(filePath);
        }
    }
    else
    {
        infile.close();
    }

    return r;
}

bool Mesh::loadMeshProfileType1(const QString& filePath)
{
    bool r = true;

    std::list<int> iBounds;
    std::list<std::pair<uint,uint>> bConnectivity;
    std::list<std::pair<float,float>> pPoints;

    resetBConnectivity();

    //Read boundary indicies and connectivity
    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word1("");
        std::string word2("");
        while (infile >> word1)
        {
            if ( word1 == "boundary" )
            {
                infile >> word2;
                break;
            }
        }

        if ( word1 == "boundary" &&  word2 == "sides")
        {
            //Do >> get boundary point indicies
            int i,j,b,t1,t2;
            while ( infile >> i >> j >> t1 >> b >> t2)
            {
                if (b == 1)
                {
                    iBounds.push_back(i);
                    bConnectivity.emplace_back(i,j);
                }
            }
            iBounds.sort();
        }
    }
    infile.close();

    //Read boundary points from indecies
    infile.open(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word3("");
        while (infile >> word3)
        {
            if ( word3 == "coordinates" )
            {
                break;
            }
        }

        if ( word3 == "coordinates" )
        {
            //Do >> get coordinates
            int i;
            float x,y;
            float t1, t2;

            while (infile >> i >> x >> y >> t1 >> t2)
            {
                if (i == iBounds.front())
                {
                    pPoints.emplace_back(x,y);
                    iBounds.pop_front();
                }
            }

            for (auto &pair : pPoints)
            {
                addBoundaryPoint(pair.first, pair.second);
            }

            for (auto& pair : bConnectivity)
            {
                addBConnectivity(pair.first, pair.second);
            }

        }
    }
    infile.close();

    return r;
}

void Mesh::addBConnectivity(const uint& a, const uint& b)
{
    mBoundConnects.emplace_back( a, b );
}

void Mesh::resetBConnectivity()
{
    mBoundConnects.clear();
}

bool Mesh::loadMeshProfileType2(const QString &filePath)
{
    bool r = true;

    std::list<int> iBounds;
    std::list<std::pair<uint,uint>> bConnectivity;
    std::list<std::pair<float,float>> pPoints;

    resetBConnectivity();

    //Read boundary indicies and connectivity
    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word1("");
        std::string word2("");
        while (infile >> word1)
        {
            if ( word1 == "Boundary" )
            {
                infile >> word2;
                break;
            }
        }

        if ( word1 == "Boundary" &&  word2 == "Faces")
        {
            //Do >> get boundary point indicies
            int i,j,b;
            while ( infile >> i >> j >> b )
            {
                if (b == 1)
                {
                    iBounds.push_back(i);
                    bConnectivity.emplace_back(i,j);
                }
            }
            iBounds.sort();
        }
    }
    infile.close();

    //Read boundary points from indecies
    infile.open(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word3("");
        while (infile >> word3)
        {
            if ( word3 == "Coordinates" )
            {
                break;
            }
        }

        if ( word3 == "Coordinates" )
        {
            //Do >> get coordinates
            int i;
            float x,y;

            while (infile >> i >> x >> y)
            {
                if (i == iBounds.front())
                {
                    pPoints.emplace_back(x,y);
                    iBounds.pop_front();
                }
            }

            for (auto &pair : pPoints)
            {
                addBoundaryPoint(pair.first, pair.second);
            }

            for (auto& pair : bConnectivity)
            {
                addBConnectivity(pair.first, pair.second);
            }
        }
    }
    infile.close();

    return r;
}

bool Mesh::loadMesh(const QString& filePath)
{
    bool r = true;
    int type = 1;
    std::string word1 = "";

    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        while (infile >> word1)
        {
            if (word1 == "clean")
            {
                type = 2;
                break;
            }
        }

        infile.close();

        if (type == 1)
        {
            r &= loadMeshType1(filePath );
        }
        else if (type == 2)
        {
            r &= loadMeshType2(filePath);
        }
    }
    else
    {
        infile.close();
    }

    return r;
}

bool Mesh::loadMeshType1(const QString& filePath)
{
    bool r = true;

    //Read mesh connectivity
    std::vector<std::tuple<uint,uint,uint>> mConnectivity;
    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word1("");
        std::string word2("");
        while (infile >> word1)
        {
            if ( word1 == "connectivities" ) break;
        }

        if ( word1 == "connectivities")
        {
            //Do >> get triangle connectivities
            uint i,j,k,t1;
            while ( infile >> word2 >> i >> j >> k >> t1)
            {
                if (word2 != "coordinates")
                {
                    mConnectivity.emplace_back(i,j,k);
                }
                else break;
            }

            clearMeshConnectivities();
            for (auto& tuple : mConnectivity)
            {
                addMeshConnectivity(tuple);
            }
        }
    }
    infile.close();

    //Read mesh points
    std::vector<std::pair<float,float>> mPoints;
    infile.open(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word1("");
        std::string word2("");
        while (infile >> word1)
        {
            if ( word1 == "coordinates" ) break;
        }

        if ( word1 == "coordinates" )
        {
            float x,y,t1,t2;
            while ( infile >> word2 >> x >> y >> t1 >> t2)
            {
                if (word2 != "unknown")
                {
                    mPoints.emplace_back(x,y);
                }
                else break;
            }

            clearMeshPoints();
            for (auto &pair : mPoints)
            {
                addMeshPoint(pair);
            }
        }
    }
    infile.close();

    return r;
}

bool Mesh::loadMeshType2(const QString &filePath)
{
    bool r = true;

    //Read mesh connectivity
    std::vector<std::tuple<uint,uint,uint>> mConnectivity;
    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word1("");
        std::string word2("");
        while (infile >> word1)
        {
            if ( word1 == "Connectivities" ) break;
        }

        if ( word1 == "Connectivities")
        {
            //Do >> get triangle connectivities
            uint i,j,k;
            while ( infile >> word2 >> i >> j >> k)
            {
                if (word2 != "Coordinates")
                {
                    mConnectivity.emplace_back(i,j,k);
                }
                else break;
            }

            clearMeshConnectivities();
            for (auto& tuple : mConnectivity)
            {
                addMeshConnectivity(tuple);
            }
        }
    }
    infile.close();

    //Read mesh points
    std::vector<std::pair<float,float>> mPoints;
    infile.open(filePath.toStdString(), std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        std::string word1("");
        std::string word2("");
        while (infile >> word1)
        {
            if ( word1 == "Coordinates" ) break;
        }

        if ( word1 == "Coordinates" )
        {
            float x,y;
            while ( infile >> word2 >> x >> y)
            {
                if (word2 != "Boundary")
                {
                    mPoints.emplace_back(x,y);
                }
                else break;
            }

            clearMeshPoints();
            for (auto &pair : mPoints)
            {
                addMeshPoint(pair);
            }
        }
    }
    infile.close();

    return r;
}

Enum::Mesh Mesh::getMeshDensity() const
{
    return mMeshDensity;
}

void Mesh::setMeshDensity(const Enum::Mesh& meshDensity)
{
    mMeshDensity = meshDensity;
}

int Mesh::getNumBoundaryLayers() {
    return mNumBoundaryLayers;
}

qreal Mesh::getBoundaryLayerThickness() {
    return mBoundaryLayerThickness;
}


bool Mesh::createInputFile(const std::string& meshInFile,
                               const std::string& meshBacFile,
                               const std::string& meshGeoFile,
                               const std::string& meshDatFile)
{
    bool r = true;

    std::ofstream outfile(meshInFile, std::ofstream::out);
    r &= outfile.is_open();

    if (r)
    {
        outfile << meshGeoFile << std::endl;
        outfile << meshBacFile << std::endl;
        outfile << meshDatFile << std::endl;
        outfile << "0" << std::endl;
        outfile << "1" << std::endl;
        outfile << "N" << std::endl;
        outfile << "1" << std::endl;
        outfile << "1" << std::endl;
        outfile << "1" << std::endl;
        outfile << "1" << std::endl;
        outfile << "3" << std::endl;
        outfile << "5" << std::endl;
    }
    outfile.close();

    return r;
}

void Mesh::stopMesher()
{
    mMeshProcess.kill();
    qInfo() << "Any running mesh jobs have been stopped!";
}

void Mesh::writeStdOutToLog()
{
    QByteArray byteArray = mMeshProcess.readAllStandardOutput();
    QStringList strLines = QString(byteArray).split("\n");

    foreach (QString line, strLines){
        qInfo() << line;
    }
}

void Mesh::writeStdErrToLog()
{
    QByteArray byteArray = mMeshProcess.readAllStandardError();
    QStringList strLines = QString(byteArray).split("\n");

    foreach (QString line, strLines){
        qCritical() << line;
    }
}

void Mesh::meshingFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
    if (exitStatus == QProcess::NormalExit && exitCode == 0)
    {
        bool r = true;

        qInfo() << "Process finished normally";

        //Set and/or check existance of output data file
        r &= mMeshProcess.exitCode() == 0;
        r &= QFile::exists(mMeshDatFile);

        //Now load points into data object
        if (r)
        {
            r &= loadMeshProfile(mMeshDatFile);
            r &= loadMesh(mMeshDatFile);
        }

        //When done set menu to 'Yes'
        if (r)
        {
            qInfo() << "Mesh successfully created.";
        }
        else
        {
            qWarning() << "Mesh not created!";
        }

        emit meshChanged();
    }
    else
    {
        qCritical() << "Process finished abnormally with exit code: " << exitCode;
        qCritical() << "Standard Ouput";
        writeStdOutToLog();
        qCritical() << "Standard Error";
        writeStdErrToLog();
    }
}

void Mesh::clearMeshConnectivities()
{
    mMeshConnectivities.clear();
}

const MeshConnectivities& Mesh::getMeshConnectivities() const
{
    return mMeshConnectivities;
}

const BConnectivities& Mesh::getBConnects() const
{
    return mBoundConnects;
}

Boundaries& Mesh::getMeshBoundary()
{
    return mMeshProfile;
}

void Mesh::selectControlPoint(const unsigned int& index)
{
    //if control point exists then delete, else add.
    auto it = std::find (mControlPoints.begin(), mControlPoints.end(), index);
    if (it != mControlPoints.end()) mControlPoints.erase(it);
    else mControlPoints.emplace_back(index);
    mControlPoints.unique();
    mControlPoints.sort();
}

BoundaryPoint& Mesh::getControlPoint(const uint& index)
{
    return mMeshProfile.at(index);
}

bool Mesh::checkBoundaryIntegrity()
{
    bool r = true;



    return r;
}

void Mesh::resetBoundary()
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
bool Mesh::checkControlPointIntegrity()
{
    bool r = true;

    return r;
}

void Mesh::addMeshPoint(const std::pair<float, float>& pair)
{
    mMeshPoints.emplace_back( pair );
}

void Mesh::clearMeshPoints()
{
    mMeshPoints.clear();
}

const MeshPoints& Mesh::getMeshPoints() const
{
    return mMeshPoints;
}

void Mesh::addMeshConnectivity(const std::tuple<uint, uint, uint>& tuple)
{
    mMeshConnectivities.emplace_back( tuple );
}

void Mesh::addMeshData(const std::tuple<float, float, float, float, float>& tuple)
{
    mMeshResults.emplace_back( tuple );
}

void Mesh::clearMeshData()
{
    mMeshResults.clear();
}

const MeshResults& Mesh::getMeshData() const
{
    return mMeshResults;
}

void Mesh::clear() {
    mMeshProfile.clear();
    mBoundConnects.clear();
    mControlPoints.clear();
    mMeshPoints.clear();
    mMeshConnectivities.clear();
    mMeshResults.clear();
}

void Mesh::addBoundaryPoint(const float& x, const float& y)
{
    mMeshProfile.emplace_back(x,y);
}

bool Mesh::loadResults(const std::string& filePath)
{
    bool r = true;

    float i,rho,u,v,e,p;
    std::ifstream infile(filePath, std::ifstream::in);
    r &= infile.is_open();
    if (r)
    {
        clearMeshData();
        while (infile >> i >> rho >> u >> v >> e >> p)
        {
            addMeshData( std::make_tuple(rho, u, v, e, p) );
        }
    }
    infile.close();

    return r;
}

void Mesh::setNumBoundaryLayers(int num) {
    mNumBoundaryLayers = num;
}

void Mesh::setBoundaryLayerThickness(qreal thickness) {
    mBoundaryLayerThickness = thickness;
}

const std::list<unsigned int>& Mesh::getControlPoints() const
{
    return mControlPoints;
}
