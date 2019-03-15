#include "ProfileModel.h"
#include <QDebug>
#include <QSettings>
#include <QFileInfo>

//Container for profile class
ProfileModel::ProfileModel(QObject *parent) : QAbstractListModel(parent)
{
    QSettings settings;

    int size = settings.beginReadArray("profiles");
    for (int i = 0; i < size; ++i) {
        settings.setArrayIndex(i);
        QString profilePath = settings.value("filepath").toString();

        QFileInfo finfo(profilePath);
        if(finfo.exists() && finfo.size() > 0) {
            addProfileFromFilePath(profilePath);
        } else {
            settings.remove("filepath");
        }
    }
    settings.endArray();
}

bool ProfileModel::addProfileFromFilePath(QString filePath) {
    int first = mProfileList.size();
    int last = mProfileList.size();
    std::unique_ptr<Profile> profile(new Profile());
    bool success = profile->setFile(filePath);
    if(success) {
        beginInsertRows(QModelIndex(),first,last);
        mProfileList.push_back(std::move(profile));
        endInsertRows();

        //Emit a signal that a new profile has been added for other processes to see
        emit newProfileAdded(0);
    }

    return success;
}

int ProfileModel::rowCount(const QModelIndex &parent) const {
    return mProfileList.size();
}

ProfilePoints ProfileModel::getProfileAtIndex(const QModelIndex index) const {
    int i = mProfileList.size() - index.row() - 1;
    return mProfileList.at(i)->getProfile();
}

QString ProfileModel::getDisplayStringAtIndex(const QModelIndex index) const {
    int i = mProfileList.size() - index.row() - 1;
    return mProfileList.at(i)->getDisplayString();
}

QVariant ProfileModel::data(const QModelIndex &index, int role) const {
    if(!index.isValid())
        return QVariant();

    if(index.row() >= mProfileList.size() || index.row() < 0)
        return QVariant();

    if(role == Qt::DisplayRole) {
        QString displayString = getDisplayStringAtIndex(index);
        qDebug() << displayString;
        return displayString;

    } else if (role == Qt::UserRole) {
        return QVariant::fromValue(getProfileAtIndex(index));

    } else {
        return QVariant();
    }
}
