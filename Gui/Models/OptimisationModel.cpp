#include "OptimisationModel.h"

#include "Optimisation.h"
#include <QSettings>
#include <QDebug>
#include <fstream>

using namespace std;

Q_DECLARE_METATYPE(std::shared_ptr<Optimisation>)

OptimisationModel::OptimisationModel(QObject* parent) :
    QAbstractListModel(parent),
    mOptimisations(new vector<std::shared_ptr<Optimisation>>()),
    mSelectionModel(nullptr)
{
    mDirWatcher.setParent(this);
    mOptProcess.setParent(this);

    auto processFinished = QOverload<int>::of(&QProcess::finished);
    connect(&mOptProcess, processFinished, this, &OptimisationModel::aerOptFinished);
}

OptimisationModel::~OptimisationModel()
{
    mOptProcess.kill();
}

QModelIndex OptimisationModel::addOptimisation(std::shared_ptr<Optimisation> optimisation)
{
    int rows = rowCount();
    beginInsertRows(QModelIndex(), rows, rows);
    mOptimisations->push_back(optimisation);
    endInsertRows();
    return index(rows, 0);
}

QModelIndex OptimisationModel::addOptimisation(const Optimisation& optimisation)
{
    int rows = rowCount();
    beginInsertRows(QModelIndex(), rows, rows);
    unique_ptr<Optimisation>newOptimisation(new Optimisation(optimisation));
    mOptimisations->push_back(move(newOptimisation));
    endInsertRows();
    return index(rows, 0);
}

int OptimisationModel::rowCount(const QModelIndex& /*parent*/) const
{
    return mOptimisations->size();
}

QVariant OptimisationModel::data(const QModelIndex& index, int role) const
{
    if (!isIndexValid(index)) {
        return QVariant();
    }

    switch (role) {
        case Qt::DisplayRole:
            return  mOptimisations->at(index.row())->label();
            break;

        case Roles::Object:
            return QVariant::fromValue(mOptimisations->at(index.row()));
            break;

        default:
            return QVariant();
    }
}

bool OptimisationModel::removeRows(int row, int count, const QModelIndex& parent)
{
    if (row < 0
            || row >= rowCount()
            || count < 0
            || (row + count) > rowCount()) {
        return false;
    }

    beginRemoveRows(parent, row, row + count - 1);
    int countLeft = count;
    while(countLeft--) {
        const Optimisation& optimisation = *mOptimisations->at(row + countLeft);
    }
    mOptimisations->erase(mOptimisations->begin() + row,
                    mOptimisations->begin() + row + count);
    endRemoveRows();


    return true;
}

bool OptimisationModel::isIndexValid(const QModelIndex& index) const
{
    if (index.row() < 0
            || index.row() >= rowCount()
            || !index.isValid()) {
        return false;
    }
    return true;
}

std::shared_ptr<Optimisation> OptimisationModel::currentOptimisation() {
    QModelIndex index = selectionModel()->currentIndex();
    return data(index, Roles::Object).value<std::shared_ptr<Optimisation>>();
}

void OptimisationModel::setSelectionModel(QItemSelectionModel* model) {
    mSelectionModel = model;
}

QItemSelectionModel* OptimisationModel::selectionModel() {
    return mSelectionModel;
}

void OptimisationModel::run(std::shared_ptr<Optimisation> optimisation) {
    // load settings
    QSettings settings;
    QString AerOpt = settings.value("AerOpt/executable").toString();
    QString AerOptWorkDir = settings.value("AerOpt/workdir").toString();
    QString AerOptInFile = settings.value("AerOpt/inputFile").toString();
    QString AerOptNodeFile = settings.value("AerOpt/nodeFile").toString();
    QString outputData = settings.value("AerOpt/outputData").toString();
    QString meshDatFile = settings.value("mesher/datFile").toString();
    QString inFolder = settings.value("AerOpt/inFolder").toString();
    QString outFolder = settings.value("AerOpt/outFolder").toString();

    bool r = true;

    //Empty the input and output directories
    r &= FileManipulation::emptyFolder(inFolder);
    r &= FileManipulation::emptyFolder(outFolder);

    //Copy mesh file as input to AerOpt from project directory to input directory
    QString dest = QDir::toNativeSeparators(inFolder + "/Mesh.dat");
    r &= FileManipulation::copyFile(meshDatFile, dest);

    //Set input files in input directory
    r &= createAerOptInFile(AerOptInFile.toStdString(), optimisation);
    r &= createAerOptNodeFile(AerOptNodeFile.toStdString(), optimisation);

    //then run aeropt
    if (r)
    {
        mDirWatcher.removePaths( mDirWatcher.directories() );
        mDirWatcher.addPath(outputData);
        mOptProcess.setWorkingDirectory( AerOptWorkDir );
        mOptProcess.start( AerOpt );
    }
    else
    {
        qWarning() << "Process Aborted!";
    }
}

bool OptimisationModel::createAerOptInFile(const std::string& filePath, std::shared_ptr<Optimisation> optimisation)
{
    bool r = true;
    QSettings settings;
    QString rootDir = settings.value("AerOpt/rootDir").toString();

    std::ofstream outfile(filePath, std::ofstream::out);
    r &= outfile.is_open();

    int controlPointCount = optimisation->controlPointCount();
    r &= controlPointCount > 0;

    if (r)
    {
        std::string xrange;
        std::string mxrange;
        std::string yrange;
        std::string myrange;

        for (auto& point : optimisation->controlPoints())
        {
            QRectF rect = point->controlPointRect();

            xrange += " ";
            xrange += std::to_string(rect.right());

            yrange += " ";
            yrange += std::to_string(rect.bottom());

            mxrange += " ";
            mxrange += std::to_string(rect.left());

            myrange += " ";
            myrange += std::to_string(rect.top());
        }

        outfile << "&inputVariables" << std::endl;
        outfile << "IV%Ma = " << std::to_string( optimisation->machNo() ) << std::endl;
        outfile << "IV%Tamb = " << std::to_string( optimisation->freeTemp() ) << std::endl;// < deg kelvin
        outfile << "IV%Pamb = " << std::to_string( optimisation->freePress() ) << std::endl;// < Pa
        outfile << "IV%R = 287" << std::endl;
        outfile << "IV%gamma = 1.4" << std::endl;
        outfile << "IV%Re = " << std::to_string( optimisation->reNo() ) << std::endl;
        outfile << "IV%xrange =" << xrange << mxrange << std::endl;
        outfile << "IV%yrange =" << yrange << myrange << std::endl;
        outfile << "IV%zrange = 0.00" << std::endl;
        outfile << "IV%angle = 0.0" << std::endl;
        outfile << "IV%Cnconnecttrans = 0" << std::endl;
        outfile << "IV%engFMF = 1.0" << std::endl;
        outfile << "IV%AlphaInflowDirection = " << std::to_string( optimisation->freeAlpha() ) << std::endl;// < angle of attack
        outfile << "IV%turbulencemodel = 0" << std::endl;
        outfile << "IV%Low2Top = " << std::to_string(float(optimisation->getNoTop())/100) << std::endl;
        outfile << "IV%NoSnap = " << std::to_string( optimisation->noAgents() ) << std::endl;
        outfile << "IV%NoCN = " << controlPointCount << std::endl;// < number of control nodes
        outfile << "IV%NoDim = 2" << std::endl;
        outfile << "IV%DoF = 8" << std::endl;// < Degrees freedom
        outfile << "IV%NoG = " << std::to_string( optimisation->noGens() ) << std::endl;// < Generations
        outfile << "IV%NoPOMod = -1" << std::endl;
        outfile << "IV%NoLeviSteps = 100" << std::endl;
        outfile << "IV%NoIter = -3" << std::endl;
        outfile << "IV%delay = 30" << std::endl;
        outfile << "IV%constrain = .True." << std::endl;
        outfile << "IV%AdaptSamp = .FALSE." << std::endl;
        outfile << "IV%waitMax = 48" << std::endl;
        outfile << "IV%maxit = 50000" << std::endl;
        outfile << "IV%Aconst = 0.01" << std::endl;
        outfile << "IV%POD = .false." << std::endl;
        outfile << "IV%NoDelBP = 4" << std::endl;
        outfile << "IV%ObjectiveFunction = " << std::to_string( optimisation->objFunc() ) << std::endl;
        //1 - Lift/Drag
        //2 - Distortion << no point in using this.
        //3 - max Lift
        //4 - min Drag
        //5 - max Downforce
        //6 - min Lift

        outfile << "! Test Parameters for Mesh Deformation" << std::endl;
        outfile << "IV%MeshMovement = 4" << std::endl;
        //1 - 'Linear with Smoothing' << to be tested currently unavailable.
        //2 - 'FDGD with Smoothing' << in validation stage could be good.
        //3 - 'RBF' << produces poor results thus don't use.
        //4 - 'FDGD' << use this for the moment << this is prevous version 1.
        outfile << "IV%Meshtest = .false." << std::endl;

        outfile << "! TEST POD" << std::endl;
        outfile << "IV%multiquadric = .true." << std::endl;
        outfile << "IV%Pol = .true." << std::endl;
        outfile << "IV%meanP = .true." << std::endl;

        outfile << "! Strings" << std::endl;
        outfile << "IV%filepath = '" << rootDir.toStdString() << "'" << std::endl;
        outfile << "IV%filename = 'Geometry'" << std::endl;
        outfile << "IV%runOnCluster = 'N'" << std::endl;

#ifdef Q_OS_UNIX
        // do fancy unix stuff
        outfile << "IV%SystemType = 'L'" << std::endl;
#endif
#ifdef Q_OS_WIN32
        // do windows stuff here
        outfile << "IV%SystemType = 'W'" << std::endl;
#endif

        outfile << "! Login Information" << std::endl;
        outfile << "IV%UserName = 'egnaumann'" << std::endl;
        outfile << "IV%Password = 'Fleur666'" << std::endl;
        outfile << "IV%version = '2.3'" << std::endl;
        outfile << "/" << std::endl;
    }
    else
    {
        qWarning() << "Number of control points is zero or AerOpt Input File could not be created.";
    }
    outfile.close();

    return r;
}

bool OptimisationModel::createAerOptNodeFile(const std::string& filePath, std::shared_ptr<Optimisation> optimisation)
{
    bool r = true;

    std::ofstream outfile(filePath, std::ofstream::out);
    r &= outfile.is_open();

    if (r)
    {
        std::string xrange;
        std::string mxrange;
        std::string yrange;
        std::string myrange;

        for (auto& point : optimisation->controlPoints())
        {
            outfile << std::to_string(point->x())
                    << "	"
                    << std::to_string(point->y())
                    << std::endl;
        }
    }
    outfile.close();

    return r;
}

void OptimisationModel::aerOptFinished(int exitCode) {
    qInfo() << "Optimisation Finished with exit code" << exitCode;
}
