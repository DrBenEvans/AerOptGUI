#include "OptimisationModel.h"

#include "Optimisation.h"
#include <QSettings>
#include <QDebug>
#include <fstream>
#include <QTimer>
#include <QDesktopServices>
#include <QUrl>

using namespace std;

Q_DECLARE_METATYPE(Optimisation*)

OptimisationModel::OptimisationModel(QObject* parent) :
    QAbstractListModel(parent)
{
}

QModelIndex OptimisationModel::addOptimisation(Optimisation *optimisation)
{
    int rows = rowCount();
    beginInsertRows(QModelIndex(), rows, rows);
    mOptimisations.push_back(optimisation);
    optimisation->setModel(this);
    endInsertRows();
    return index(rows, 0);
}

int OptimisationModel::rowCount(const QModelIndex& /*parent*/) const
{
    return mOptimisations.size();
}

QVariant OptimisationModel::data(const QModelIndex& index, int role) const
{
    if (!isIndexValid(index)) {
        return QVariant();
    }

    switch (role) {
        case Qt::DisplayRole:
            return  mOptimisations.at(index.row())->label();

        case Roles::Object:
            return QVariant::fromValue(mOptimisations.at(index.row()));

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
        const Optimisation& optimisation = *mOptimisations.at(row + countLeft);
    }
    mOptimisations.erase(mOptimisations.begin() + row,
                    mOptimisations.begin() + row + count);
    endRemoveRows();


    return true;
}

bool OptimisationModel::isIndexValid(int row) {
    QModelIndex modelIndex = createIndex(row, 0);
    return isIndexValid(modelIndex);
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

Optimisation* OptimisationModel::optimisation(uint index) {
    if(isIndexValid(index)) {
        return mOptimisations.at(index);
    } else {
        return nullptr;
    }
}

bool OptimisationModel::run(Optimisation *optimisation) {
    return optimisation->run();
}

void OptimisationModel::emitOptimisationFitnessChanged(Optimisation* optChanged) {
    for(int i=0; i < mOptimisations.size(); i++) {
        if(optimisation(i) == optChanged) {
            emit optimisationFitnessChanged(i);
            return;
        }
    }
}

void OptimisationModel::emitOptimisationOutputChanged(Optimisation* optChanged) {
    for(int i=0; i < mOptimisations.size(); i++) {
        if(optimisation(i) == optChanged) {
            emit optimisationOutputChanged(i);
            return;
        }
    }
}

void OptimisationModel::revealFiles(int index) {
    if(isIndexValid(index)) {
        QString path = optimisation(index)->simulationDirectoryPath();
        QDesktopServices::openUrl(QUrl::fromLocalFile(path));
    }
}

QModelIndex OptimisationModel::loadByInputFilePath(QString path) {
    Optimisation *optimisation = new Optimisation();
    bool success = optimisation->load(path);
    if(success) {
        return addOptimisation(optimisation);
    } else {
        return QModelIndex();
    }
}
