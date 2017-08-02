#ifndef PROFILELIBRARY_H
#define PROFILELIBRARY_H

#include <QAbstractListModel>
#include "Profile.h"

class ProfileLibrary : public QAbstractListModel
{
    Q_OBJECT

public:
    ProfileLibrary(QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
    void addProfileFromFilePath(QString filePath);
    ProfileSharedPointer getProfileAtIndex(const QModelIndex index) const;

    std::vector<ProfileSharedPointer> mProfileList;
};

#endif // PROFILELIBRARY_H
