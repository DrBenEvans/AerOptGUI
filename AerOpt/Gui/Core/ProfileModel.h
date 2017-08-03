#ifndef PROFILEMODEL_H
#define PROFILEMODEL_H

#include <QAbstractListModel>
#include "Profile.h"

class ProfileModel : public QAbstractListModel
{
    Q_OBJECT

public:
    ProfileModel(QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
    void addProfileFromFilePath(QString filePath);
    ProfileSharedPointer getProfileAtIndex(const QModelIndex index) const;

    std::vector<ProfileSharedPointer> mProfileList;
};

#endif // PROFILEMODEL_H
