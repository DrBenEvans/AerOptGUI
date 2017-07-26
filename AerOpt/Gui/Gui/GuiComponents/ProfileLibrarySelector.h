#ifndef ProfileLibrarySelector_H
#define ProfileLibrarySelector_H

#include <vector>
#include <QString>
#include <QComboBox>
#include "Profile.h"

class ProfileLibrarySelector : public QComboBox
{
public:
    explicit ProfileLibrarySelector(QWidget* parent);
    Profile& getSelectedProfile();
private slots:
    void on_currentIndexChanged(int index);
private:
    std::vector<QString> getProfileNameList();
    QString getIndexedProfileFileName(int);
    std::vector<Profile> mProfiles;
    void addProfile();
    bool addProfileFromFileName(std::string fileName);
};

#endif // ProfileLibrarySelector_H
