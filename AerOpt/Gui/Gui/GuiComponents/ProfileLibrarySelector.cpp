#include "ProfileLibrarySelector.h"
#include "Profile.h"
#include <QFileDialog>
#include <QDebug>

ProfileLibrarySelector::ProfileLibrarySelector(QWidget* parent) : QComboBox(parent)
{

    addProfileFromFileName("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA0024.prf");
    addProfileFromFileName("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA21120.prf");

    addItem(QString("Add New..."));
    for(Profile& profile: mProfiles) {
        this->addItem(QString::fromStdString(profile.getDisplayString()));
    }
    this->setCurrentIndex(1);
}

bool ProfileLibrarySelector::addProfileFromFileName(std::string fileName) {
    Profile profile(fileName);
    mProfiles.push_back(profile);
    return true;
}

std::vector<QString> ProfileLibrarySelector::getProfileNameList() {
    std::vector<QString> string_list;
    string_list.push_back(QString("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA0024.prf"));
    string_list.push_back(QString("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA21120.prf"));
    return string_list;
}

Profile& ProfileLibrarySelector::getSelectedProfile() {
    return mProfiles[this->currentIndex()];
}

void ProfileLibrarySelector::on_currentIndexChanged(int index) {
    if(index==0) {
        this->addProfile();
    }
}

QString ProfileLibrarySelector::getIndexedProfileFileName(int index) {
   return QString('/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/AerOpt/Example_profile_files/NACA21120.prf');
}

void ProfileLibrarySelector::addProfile()
{
    bool success = true;
    QString fileName;
    QStringList fileNames;


#ifdef Q_OS_UNIX
    // do fancy unix stuff
    fileNames.append(
        QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath()+"/Documents/Projects/AerOptProject/", "Profile Files (*.prf)")
                    );
#endif
#ifdef Q_OS_WIN32
    // do windows stuff here
    fileNames.append(
        QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath(), "Profile Files (*.prf)")
                    );
#endif

    if (fileNames.size() > 0)
    {
        for (const QString &f: fileNames)
        {
            fileName = f;
        }

        qInfo() << "File selected: " << fileName;

        success &= addProfileFromFileName(fileName.toStdString());

        if (success)
        {
            //Set text OK here!
            qInfo() << "File successfully loaded.";
        }
        else
        {
            //Set text Not OK here!
            qWarning() << "File failed to load correctly.";
        }
    }
    else
    {
        //Set text Not OK here!
        qWarning() << "Profile data not imported!";
    }
}
