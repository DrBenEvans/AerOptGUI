#ifndef ProfileLibrarySelector_H
#define ProfileLibrarySelector_H

#include <vector>
#include <QString>
#include <QComboBox>
#include "Profile.h"

class ProfileLibrarySelector : public QComboBox
{
public:
    explicit ProfileLibrarySelector(std::vector<Profile>* profiles, QWidget* parent);
    explicit ProfileLibrarySelector(QWidget* parent);
    Profile* getSelectedProfile();
private slots:
    void on_currentIndexChanged(int index);
private:
    std::vector<Profile>* mProfiles;
    void addProfile();
};

#endif // ProfileLibrarySelector_H
