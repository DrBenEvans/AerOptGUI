#include "OptimisationModel.h"

#include "Optimisation.h"
#include <QSettings>
#include <QDebug>
#include <fstream>
#include <QTimer>

using namespace std;

Q_DECLARE_METATYPE(std::shared_ptr<Optimisation>)

OptimisationModel::OptimisationModel(QObject* parent) :
    QAbstractListModel(parent)
{
}

QModelIndex OptimisationModel::addOptimisation(std::shared_ptr<Optimisation> optimisation)
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
            break;

        case Roles::Object:
            return QVariant::fromValue(mOptimisations.at(index.row()));
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
        const Optimisation& optimisation = *mOptimisations.at(row + countLeft);
    }
    mOptimisations.erase(mOptimisations.begin() + row,
                    mOptimisations.begin() + row + count);
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

std::shared_ptr<Optimisation> OptimisationModel::optimisation(uint index) {
    return mOptimisations.at(index);
}

void OptimisationModel::run(std::shared_ptr<Optimisation> optimisation) {
    optimisation->run();
}

void OptimisationModel::emitOptimisationDataChanged(Optimisation* optimisation) {
    for(int i=0; i<mOptimisations.size(); i++) {
        if(mOptimisations.at(i).get() == optimisation) {
            emit optimisationDataChanged(i);
            return;
        }
    }
}

void OptimisationModel::emitOptimisationOutputChanged(Optimisation* optimisation) {
    for(int i=0; i<mOptimisations.size(); i++) {
        if(mOptimisations.at(i).get() == optimisation) {
            emit optimisationOutputChanged(i);
            return;
        }
    }
}
