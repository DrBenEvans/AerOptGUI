#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include <memory>
#include <vector>
#include <QItemSelectionModel>

#include <QAbstractListModel>

#include "Optimisation.h"

class Optimisation;
class DatabaseManager;
class OptimisationModel;

class OptimisationModel : public QAbstractListModel
{
    Q_OBJECT
public:
    enum Roles {
        Object = Qt::UserRole + 1
    };

    OptimisationModel(QObject* parent = 0);

    QModelIndex addOptimisation(const Optimisation &optimisation);
    QModelIndex addOptimisation(std::shared_ptr<Optimisation> optimisation);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role) const override;
    bool removeRows(int row, int count, const QModelIndex& parent) override;

    void setSelectionModel(QItemSelectionModel* model);
    QItemSelectionModel* selectionModel();

    std::shared_ptr<Optimisation> currentOptimisation();

private:
    bool isIndexValid(const QModelIndex& index) const;

private:
    std::unique_ptr<std::vector<std::shared_ptr<Optimisation>>> mOptimisations;
    QItemSelectionModel* mSelectionModel;
};

#endif // SIMULATIONMODEL_H
