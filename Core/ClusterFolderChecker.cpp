#include "ClusterFolderChecker.h"
#include <FileManipulation.h>
#include <QDebug>
#include <QDir>
#include <QFileInfo>

ClusterFolderChecker::ClusterFolderChecker()
{
}

ClusterFolderChecker::~ClusterFolderChecker()
{
}


void ClusterFolderChecker::setWorkingDirectory(QString workDir){
    workingDirectory = workDir.toStdString();
}


void ClusterFolderChecker::run(){
    QSettings settings;
    std::string outputFilename;
    std::string localFilename;
    std::string localFolder;
    int line_number=0;
    while (1) {

        outputFilename = workingDirectory + "/Output_Data/output.log";
        QString filePath = QString(("AerOptFiles/"+outputFilename).c_str());
        QDir::toNativeSeparators(filePath);
        localFilename = filePath.toStdString();

        filePath = QString(("AerOptFiles/"+workingDirectory).c_str());
        QDir::toNativeSeparators(filePath);
        localFolder = filePath.toStdString();

        std::cout << localFolder << std::endl;

        folderFromCluster("AerOpt/"+workingDirectory, "AerOptFiles/"+workingDirectory);

        std::ifstream outputfile(localFilename);
        std::string line = "";
        std::string output = "";

        if (!outputfile.is_open())
            perror("File not open");

        for( int l=0; l<line_number; l++ ){
             std::getline(outputfile, line);
        }

        while(std::getline(outputfile, line)){
            output += line + "\n";
            line_number++;
        }

        if (outputfile.bad())
            perror("error while reading file");

        emit stdOut(QString(output.c_str()));
        emit directoryChanged(filePath);

        outputfile.close();

        sleep(10);
    }
}



ssh_session createSSHSession(){
    ssh_session session = ssh_new();
    int rc;
    if (session == NULL)
        exit(-1);

    ssh_options_set(session, SSH_OPTIONS_HOST, "sunbird.swansea.ac.uk");
    std::string username;
    std::string password;
    std::ifstream infile("clusterlogin", std::ifstream::in);
    std::getline(infile, username);
    std::getline(infile, password);

    infile.close();

    ssh_options_set(session, SSH_OPTIONS_USER, username.c_str());

    rc = ssh_connect(session);
    rc &= ssh_userauth_password(session, NULL, password.c_str());
    if (rc != SSH_OK)
    {
      fprintf(stderr, "Error connecting to localhost: %s\n",
              ssh_get_error(session));
      exit(-1);
    }

    return session;
}

ssh_channel createSSHChannel(ssh_session session){
    ssh_channel channel;
    int rc;

    channel = ssh_channel_new(session);
    if (channel == NULL){
        std::cout << "Error creating an ssh channel";
        exit(-1);
    }
    rc = ssh_channel_open_session(channel);
    if (rc != SSH_OK)
    {
      fprintf(stderr, "Error opening ssh channel %s\n",
              ssh_get_error(session));
      ssh_channel_close(channel);
      ssh_channel_free(channel);
      exit(-1);
    }

    return channel;
}

void sshExecute(std::string command){
    int rc;
    ssh_channel channel;
    ssh_session session = createSSHSession();
    channel = createSSHChannel(session);
    rc = ssh_channel_request_exec(channel, command.c_str());
    if (rc != SSH_OK)
    {
      fprintf(stderr, "Failed to execute command 1\n");
      fprintf(stderr, "Error opening ssh channel %s\n",
              ssh_get_error(session));
      ssh_channel_close(channel);
      ssh_channel_free(channel);
      exit(-1);
    }
    ssh_channel_send_eof(channel);
    ssh_disconnect(session);
    ssh_free(session);
}


sftp_session createSFTPSession(ssh_session session){
    sftp_session sftp;
    int rc;

    sftp = sftp_new(session);
    if (sftp == NULL)
    {
      fprintf(stderr, "Error allocating SFTP session: %s\n",
              ssh_get_error(session));
      exit(-1);
    }

    rc = sftp_init(sftp);
    if (rc != SSH_OK)
    {
      fprintf(stderr, "Error initializing SFTP session: %s.\n",
              sftp_get_error(sftp));
      sftp_free(sftp);
      exit(-1);
    }

    return sftp;
}


int FileToCluster(std::string source, std::string destination){
    int rc;
    sftp_session sftp;
    sftp_file file;
    char buffer[16384];
    int nbytes, nwritten;
    int fd;

    ssh_session session = createSSHSession();
    sftp = createSFTPSession(session);

    file = sftp_open(sftp, destination.c_str(), O_WRONLY | O_CREAT, S_IRUSR|S_IWUSR);
    if (file == NULL) {
        fprintf(stderr, "FileToCluster: Can't open file for reading: %s\n",
                ssh_get_error(session));
        return SSH_ERROR;
    }

    fd = open(source.c_str(), O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "FileToCluster: Can't open file for writing: %s\n",
                strerror(errno));
        return SSH_ERROR;
    }

    for (;;) {
          nbytes = read(fd, buffer, sizeof(buffer));
          if (nbytes == 0) {
              break; // EOF
          } else if (nbytes < 0) {
              fprintf(stderr, "FileToCluster: Error while reading file: %s\n",
                      ssh_get_error(session));
              sftp_close(file);
              return SSH_ERROR;
          }
          nwritten = sftp_write(file, buffer, nbytes);
          if (nwritten != nbytes) {
              fprintf(stderr, "FileToCluster: Error writing: %s\n",
                      strerror(errno));
              sftp_close(file);
              return SSH_ERROR;
          }
    }

    rc = sftp_close(file);
    if (rc != SSH_OK) {
        fprintf(stderr, "FileToCluster: Can't close the read file: %s\n",
                ssh_get_error(session));
        return rc;
    }

    ssh_disconnect(session);
    ssh_free(session);
}



int fileFromCluster(std::string source, std::string destination){
    int rc;
    sftp_session sftp;
    sftp_file file;
    char buffer[16384];
    int nbytes, nwritten;
    int fd;

    ssh_session session = createSSHSession();
    sftp = createSFTPSession(session);

    file = sftp_open(sftp, source.c_str(), O_RDONLY, 0);
    if (file == NULL) {
        fprintf(stderr, "fileFromCluster: Can't open file for reading: %s\n",
                ssh_get_error(session));
        fprintf(stderr, "File name: %s\n", source.c_str());
        return SSH_ERROR;
    }

    QString filePath = QString(destination.c_str());
    QFileInfo fileinfo(filePath);
    QDir directory=fileinfo.dir();
    FileManipulation::emptyFolder(QString(directory.path()));

    fd = open(destination.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR|S_IWUSR);
    if (fd < 0) {
        fprintf(stderr, "fileFromCluster: Can't open file for writing: %s\n",
                strerror(errno));
        fprintf(stderr, "File name: %s\n", destination.c_str());
        return SSH_ERROR;
    }


    for (;;) {
          nbytes = sftp_read(file, buffer, sizeof(buffer));
          if (nbytes == 0) {
              break; // EOF
          } else if (nbytes < 0) {
              fprintf(stderr, "fileFromCluster: Error while reading file: %s\n",
                      ssh_get_error(session));
              sftp_close(file);
              return SSH_ERROR;
          }
          nwritten = write(fd, buffer, nbytes);
          if (nwritten != nbytes) {
              fprintf(stderr, "fileFromCluster: Error writing: %s\n",
                      strerror(errno));
              sftp_close(file);
              return SSH_ERROR;
          }
    }

    rc = sftp_close(file);
    if (rc != SSH_OK) {
        fprintf(stderr, "fileFromCluster: Can't close the read file: %s\n",
                ssh_get_error(session));
        return rc;
    }

    ssh_disconnect(session);
    ssh_free(session);

    return 0;
}


int folderFromCluster(std::string source, std::string destination){
    int rc;
    sftp_session sftp;
    sftp_dir directory;

    ssh_session session = createSSHSession();
    sftp = createSFTPSession(session);

    directory = sftp_opendir (sftp, source.c_str());

    while( ! sftp_dir_eof(directory) ){
        sftp_attributes file_attr = sftp_readdir(sftp, directory);
        if(file_attr != NULL){
            std::cout << source << "/" << file_attr->name << " " << file_attr->type << std::endl;
            printf("%d %d\n",file_attr->type,SSH_FILEXFER_TYPE_DIRECTORY);
            std::string subsource = source + "/" + file_attr->name;
            std::string subdestination = destination + "/" + file_attr->name;
            if( file_attr->type == SSH_FILEXFER_TYPE_DIRECTORY){
                if( strcmp(file_attr->name, ".") && strcmp(file_attr->name, "..") ){
                    folderFromCluster(subsource, subdestination);
                }
            } else {
                QString filePath = QDir::toNativeSeparators(subdestination.c_str());
                subdestination = filePath.toStdString();
                std::cout << subsource << " to " << subdestination << std::endl;
                fileFromCluster(subsource, subdestination);
            }
        }
    }

    ssh_disconnect(session);
    ssh_free(session);

    return 0;
}
