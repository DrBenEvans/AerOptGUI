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

    QModelIndex addOptimisation(Optimisation *optimisation);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role) const override;
    bool removeRows(int row, int count, const QModelIndex& parent) override;

    bool run(Optimisation* optimisation);

    void emitOptimisationFitnessChanged(Optimisation *optChanged);
    void emitOptimisationOutputChanged(Optimisation *optChanged);

    Optimisation *optimisation(uint index);

    void revealFiles(int index);
    QModelIndex loadByInputFilePath(QString path);

    bool isIndexValid(const QModelIndex& index) const;

signals:
    void optimisationFitnessChanged(int index);
    void optimisationOutputChanged(int index);

private:
    bool isIndexValid(int row);
    std::vector<Optimisation*> mOptimisations;
};

#endif // SIMULATIONMODEL_H
