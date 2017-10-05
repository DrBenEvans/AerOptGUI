#ifndef PROFILEMODEL_H
#define PROFILEMODEL_H

#include <QAbstractListModel>
#include <QItemSelectionModel>
#include <memory>
#include "Profile.h"

class ProfileModel : public QAbstractListModel
{
    Q_OBJECT

public:
    ProfileModel(QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
    bool addProfileFromFilePath(QString filePath);
    ProfilePoints getProfileAtIndex(const QModelIndex index) const;
    QString getDisplayStringAtIndex(const QModelIndex index) const;

private:
    std::vector<std::unique_ptr<Profile>> mProfileList;
};

#endif // PROFILEMODEL_H
