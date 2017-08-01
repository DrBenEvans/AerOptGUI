#include "ProfileLibrarySelector.h"
#include "Profile.h"
#include <QFileDialog>
#include <QDebug>

ProfileLibrarySelector::ProfileLibrarySelector(std::vector<Profile>* profiles, QWidget* parent) : QComboBox(parent), mProfiles(profiles) {
}

ProfileLibrarySelector::ProfileLibrarySelector(QWidget* parent) : QComboBox(parent)
{

    // Add item entries
    addItem(QString("Add New..."));
    for(Profile& profile: *mProfiles) {
        this->addItem(profile.getDisplayString());
    }
    this->setCurrentIndex(1);
}

Profile* ProfileLibrarySelector::getSelectedProfile() {
    Profile& profile = mProfiles->at(this->currentIndex() - 1);
    return &profile;
}

void ProfileLibrarySelector::on_currentIndexChanged(int index) {
    if(index==0) {
        this->addProfile();
    }
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

        mProfiles->emplace_back(fileName);

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
