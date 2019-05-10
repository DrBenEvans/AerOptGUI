#ifndef PROFILEMODEL_H
#define PROFILEMODEL_H

#include <QAbstractListModel>
#include <QItemSelectionModel>
#include <memory>
#include "Profile.h"

/**
 * @brief The ProfileModel class is a container that provides a means of
 * interacting with several loaded Profiles and adding new profiles from profile files.
 */
class ProfileModel : public QAbstractListModel
{
    Q_OBJECT

public:
    /**
     * @brief ProfileModel Constructor Method.
     * Loads a list of profiles from settings (QSettings).
     * @param parent
     */
    ProfileModel(QObject *parent = nullptr);

    /**
     * @brief rowCount Returns the number of profiles in mProfileList
     * @param parent
     * @return Number of loaded profiles
     */
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;

    /**
     * @brief data Returns either the profile points or the display string for a given profile index.
     * @param index
     * @param role Qt::DisplayRole here will return display string; Qt::UserRole will return profile
     * @return display string OR profile points
     */
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

    /**
     * @brief addProfileFromFilePath Add a new profile from a given filepath.
     * @param filePath Filepath for new profile
     * @return true iff profile is valid and added successfully; otherwise false
     */
    bool addProfileFromFilePath(QString filePath);

    /**
     * @brief getProfileAtIndex Returns the profile at the given index.
     * @param index The profile index
     * @return A list of points comprising the profile
     */
    ProfilePoints getProfileAtIndex(const QModelIndex index) const;

    /**
     * @brief getDisplayStringAtIndex Returns the display string for the profile at the given index.
     * AKA the filepath
     * @param index The profile index
     * @return Profile display string
     */
    QString getDisplayStringAtIndex(const QModelIndex index) const;

signals:
    /**
     * @brief newProfileAdded Emits when a new profile file is added from addProfileFromFilePath().
     * @param index
     */
    void newProfileAdded(int index);

private:
    /**
     * @brief mProfileList List of loaded Profiles.
     */
    std::vector<std::unique_ptr<Profile>> mProfileList;
};

#endif // PROFILEMODEL_H
