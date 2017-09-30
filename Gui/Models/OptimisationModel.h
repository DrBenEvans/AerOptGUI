#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include <memory>
#include <vector>
#include <QItemSelectionModel>

#include <QAbstractListModel>


#include "FileManipulation.h"

class Optimisation;
class DatabaseManager;

class OptimisationModel : public QAbstractListModel
{
    Q_OBJECT
public:
    enum Roles {
        Object = Qt::UserRole + 1
    };

    OptimisationModel(QObject* parent = 0);

    QModelIndex addOptimisation(std::shared_ptr<Optimisation> optimisation);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role) const override;
    bool removeRows(int row, int count, const QModelIndex& parent) override;

    void run(std::shared_ptr<Optimisation> optimisation);

    void emitOptimisationFitnessChanged(Optimisation *optimisation);
    void emitOptimisationOutputChanged(Optimisation *optimisation);

    Optimisation *optimisation(uint index);

    void revealFiles(int index);

signals:
    void optimisationFitnessChanged(int index);
    void optimisationOutputChanged(int index);

private:
    bool isIndexValid(const QModelIndex& index) const;
    bool isIndexValid(int row);
    std::vector<std::shared_ptr<Optimisation>> mOptimisations;
};

#endif // SIMULATIONMODEL_H
