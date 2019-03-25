#ifndef CLUSTERFOLDERCHECKER_H
#define CLUSTERFOLDERCHECKER_H

#include <iostream>
#include <fstream>
#include <sstream>

#include <QSettings>
#include <QThread>
#include <QTimer>

#define SSH_NO_CPP_EXCEPTIONS
#define LIBSSH_STATIC 1
#include <libssh/libsshpp.hpp>
#include <libssh/sftp.h>
#include <sys/stat.h>
#include <fcntl.h>

class ClusterFolderChecker : public QThread
{
    Q_OBJECT
public:
    ClusterFolderChecker();
    ~ClusterFolderChecker();

    void setWorkingDirectory(QString workDir);
    void setUsername(QString username);
    void setPassword(QString password);

signals:
    void directoryChanged(const QString&);
    void stdOut(const QString);
    void stdErr(const QString);

private:

    std::string workingDirectory = "";
    std::string username = "";
    std::string password = "";

    void run();
    void folderCheckLoop();
    int submitToCluster( std::string AerOptInFile, std::string simulationDirectoryName, std::string username, std::string password );

};


int folderFromCluster( std::string source, std::string destination, std::string username, std::string password );
int sshVerifyPassword( QString username, QString password );
void sshExecute(ssh_session session, std::string command);
int fileToCluster(std::string source, std::string destination, ssh_session session);
ssh_session createSSHSession( std::string username, std::string password );


#endif // CLUSTERFOLDERCHECKER_H
