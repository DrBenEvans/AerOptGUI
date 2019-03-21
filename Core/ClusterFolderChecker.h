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

signals:
    void directoryChanged(const QString&);
    void stdOut(const QString);
    void stdErr(const QString);

private:

    std::string workingDirectory;
    std::string runName;

    void run();

};





ssh_session createSSHSession();

ssh_channel createSSHChannel(ssh_session session);
void sshExecute(std::string command);
sftp_session createSFTPSession(ssh_session session);
int FileToCluster(std::string source, std::string destination);
int fileFromCluster(std::string source, std::string destination);
int folderFromCluster(std::string source, std::string destination);




#endif // CLUSTERFOLDERCHECKER_H
