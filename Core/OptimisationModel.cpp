#include "OptimisationModel.h"

#include "Optimisation.h"

using namespace std;

Q_DECLARE_METATYPE(std::shared_ptr<Optimisation>)

OptimisationModel::OptimisationModel(QObject* parent) :
    QAbstractListModel(parent),
    mOptimisations(new vector<std::shared_ptr<Optimisation>>()),
    mSelectionModel(nullptr)
{
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
