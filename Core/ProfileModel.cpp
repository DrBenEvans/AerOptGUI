#include "ProfileModel.h"
#include <QDebug>

ProfileModel::ProfileModel(QObject *parent) : QAbstractListModel(parent)
{
    addProfileFromFilePath(QString("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA0024.prf"));
    addProfileFromFilePath(QString("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA21120.prf"));
    addProfileFromFilePath(QString("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/test.prf"));
}

void ProfileModel::addProfileFromFilePath(QString filePath) {
    int first = mProfileList.size();
    int last = mProfileList.size();
    beginInsertRows(QModelIndex(),first,last);
    mProfileList.emplace_back(new Profile(this, filePath));
    endInsertRows();
}

int ProfileModel::rowCount(const QModelIndex &parent) const {
    return mProfileList.size();
}

ProfileSharedPointer ProfileModel::getProfileAtIndex(const QModelIndex index) const {
    int i = mProfileList.size() - index.row() - 1;
    return mProfileList.at(i);
}

QVariant ProfileModel::data(const QModelIndex &index, int role) const {
    if(!index.isValid())
        return QVariant();

    if(index.row() >= mProfileList.size() || index.row() < 0)
        return QVariant();

    if(role == Qt::DisplayRole) {
        QString displayString = getProfileAtIndex(index)->getDisplayString();
        qDebug() << displayString;
        return displayString;
    } else if (role == Qt::UserRole) {
        return QVariant::fromValue(getProfileAtIndex(index));
    }

    return QVariant();
}
