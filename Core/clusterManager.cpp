#include "clusterManager.h"
#include <FileManipulation.h>
#include <QDebug>
#include <QSettings>
#include <QDir>
#include <QFile>
#include <QFileInfo>

clusterManager::clusterManager()
{
}

clusterManager::~clusterManager()
{
}


void clusterManager::setWorkingDirectory(QString workDirQString){
    mWorkingDirectory = workDirQString.toStdString();
}

void clusterManager::setClusterAddress(QString addressQString){
    mAddress = addressQString.toStdString();
}

void clusterManager::setUsername(QString usernameQString){
    mUsername = usernameQString.toStdString();
}

void clusterManager::setPassword(QString passwordQString){
    mPassword = passwordQString.toStdString();
}


int clusterManager::submitToCluster(){
    /* Submit a run to the cluster. Copies necessary input files over, creates a full
     * copy of the AerOpt folder on the cluster (currently just assuming its at ~/AerOpt)
     * and starts AerOpt in a screen. */
    QSettings settings;

    ssh_session session = createSSHSession( mAddress, mUsername, mPassword );

    std::string AerOptInFile = settings.value("AerOpt/inputFile").toString().toStdString();
    std::string AerOptNodeFile =  settings.value("AerOpt/nodeFile").toString().toStdString();
    std::string meshDatFile = settings.value("mesher/initMeshFile").toString().toStdString();

    std::string clusterdir = "AerOpt/"+mWorkingDirectory;
    std::string outputDirectory = mWorkingDirectory+"/Output_Data";
    std::string outputfilename = outputDirectory+"/output.log";

    // Delete any conflicting directory and create the new one
    sshExecute(session, "rm -r "+clusterdir);
    sshExecute(session, "mkdir -p "+clusterdir+"/"+outputDirectory);
    sshExecute(session, "mkdir -p "+clusterdir+"/Input_Data/");

    // Copy input files to the cluster
    fileToCluster(AerOptInFile.c_str(),clusterdir+"/Input_Data/AerOpt_InputParameters.txt", session);
    fileToCluster(AerOptNodeFile.c_str(),clusterdir+"/Input_Data/Control_Nodes.txt", session);
    fileToCluster(meshDatFile.c_str(),clusterdir+"/Input_Data/Mesh.dat", session);

    // Copy the input file also to the output directory
    sshExecute(session, "cd "+clusterdir+"; cp Input_Data/AerOpt_InputParameters.txt "+mWorkingDirectory);

    // Copy the AerOpt executable and other necessary binaries
    sshExecute(session, "cd "+clusterdir+"; cp ../AerOpt ./");
    sshExecute(session, "cd "+clusterdir+"; cp -r ../Executables ../executables ./");

    // Create a run script and start it in screen
    sshExecute(session, "cd "+clusterdir+"; echo module load mkl > run.sh");
    sshExecute(session, "cd "+clusterdir+"; echo './AerOpt 2>&1 > "+outputfilename+"' >> run.sh");
    sshExecute(session, "cd "+clusterdir+"; chmod +x run.sh");
    sshExecute(session, "cd "+clusterdir+"; screen -L -d -m ./run.sh ");

    ssh_disconnect(session);
    ssh_free(session);

    return 0;
}


void clusterManager::folderCheckLoop(){
    /* Copy the folder back from the cluster and check for changes in a loop.
     * This essentially replaces folderChecker for local runs. */
    QSettings settings;

    std::string localdirpath = settings.value("AerOpt/workingDirectory").toString().toStdString();
    std::string outputFilename;
    std::string fitnessFilename;
    std::string localFolder;
    int line_number=0;

    // Build local filenames
    std::string clusterdir = "AerOpt/" + mWorkingDirectory+"/";
    outputFilename = mWorkingDirectory + "/Output_Data/output.log";
    QString filePath = QString((localdirpath+outputFilename).c_str());
    QDir::toNativeSeparators(filePath);
    outputFilename = filePath.toStdString();

    filePath = QString((localdirpath+mWorkingDirectory).c_str());
    QDir::toNativeSeparators(filePath);
    localFolder = filePath.toStdString();

    fitnessFilename = mWorkingDirectory + "/FitnessAll.txt";
    filePath = QString((localdirpath+fitnessFilename ).c_str());
    QDir::toNativeSeparators(filePath);
    fitnessFilename = filePath.toStdString();

    // The loop
    while ( 1 ) {

        // Copy the folder over. This will only copy files if the time stamp or
        // file size have changed
        folderFromCluster(clusterdir+mWorkingDirectory, localdirpath+mWorkingDirectory);

        // Check for new lines in the outputfile and send them to the Optimisation
        std::ifstream outputfile(outputFilename);
        std::string line = "";
        std::string output = "";

        if ( outputfile.is_open() ){

            for( int l=0; l<line_number; l++ ){
                std::getline(outputfile, line);
            }

            while(std::getline(outputfile, line)){
                output = line + "\n";
                emit stdOut(QString(output.c_str()));
                line_number++;
            }

            if (outputfile.bad())
                perror("error while reading file");

            outputfile.close();
        }


        // Check for new lines in the fitness file and update plot if they are found
        std::ifstream fitness_file(fitnessFilename);

        if ( fitness_file.is_open() ){
            int lines = 0;
            static int fitness_lines = 0;

            while(std::getline(fitness_file, line)){
                lines++;
            }

            if (fitness_file.bad())
                perror("error while reading file");

            if( lines > fitness_lines ){
                fitness_lines = lines;
                emit directoryChanged(filePath);
            }

            fitness_file.close();
        }


        sleep(5);
    }
}


void clusterManager::run(){
    // This will be run when the Qthread is started. Submit and start checking the folder
    if( submitToCluster()==0 ){
        folderCheckLoop();
    }
}



ssh_session createSSHSession( std::string address, std::string username, std::string password ){
    /* Create a new ssh session and log in. */

    ssh_session session = ssh_new();
    int rc;
    if (session == NULL) {
        printf("Could not create session\n");
        exit(-1);
    }

    ssh_options_set(session, SSH_OPTIONS_HOST, address.c_str());
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
    /* Create an ssh channel for executing a command */

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

void sshExecute(ssh_session session, std::string command){
    /* Execute a command over ssh */
    int rc;
    ssh_channel channel = createSSHChannel(session);

    rc = ssh_channel_request_exec(channel, command.c_str());
    if (rc != SSH_OK)
    {
      fprintf(stderr, "Failed to execute command %s\n",
              ssh_get_error(session));
      std::cout << "Command: " << command << std::endl;
      ssh_channel_close(channel);
      ssh_channel_free(channel);
      exit(-1);
    }
    ssh_channel_close(channel);
    ssh_channel_free(channel);
}


sftp_session createSFTPSession(ssh_session session){
    /* Create an SFTP session for accessing files */

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


int fileToCluster(std::string source, std::string destination, ssh_session session){
    /* Stream the contents of a file to a file in the cluster over sftp */

    int rc;
    sftp_session sftp;
    sftp_file file;
    char buffer[16384];
    int nbytes, nwritten;
    int fd;

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
}



int getClusterFile(std::string source, std::string destination, ssh_session session, sftp_session sftp){
    /* Stream the contents of a file in the cluster to a file over sftp */

    int rc;
    sftp_file file;
    char buffer[16384];
    int nbytes, nwritten;
    FILE *fd;

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
    if(!directory.exists()) {
        QDir().mkpath(directory.path());
    }

    fd = fopen(destination.c_str(), "w+");
    if (fd == NULL) {
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
          nwritten = fwrite(buffer, sizeof(char), nbytes, fd);
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

    fflush(fd);
    fclose(fd);

    return 0;
}


int downloadAndVerifyClusterFile(std::string source, std::string destination, int size, ssh_session session, sftp_session sftp){
    /* Check if a file has matches what should have been downloaded and retry if it doesn't.
     * If everything fails, delete the file */
    struct stat stat_buf;
    int rc = getClusterFile(source, destination, session, sftp);
    rc &= stat(destination.c_str(), &stat_buf);

    int ntry = 0;
    while( !rc && stat_buf.st_size < size/2 ){
        if(ntry>10){
            return 1;
        }
        rc = getClusterFile(source, destination, session, sftp);
        rc &= stat(destination.c_str(), &stat_buf);
        ntry++;
    }
    return rc;
}



int getClusterFolder(std::string source, std::string destination, ssh_session session, sftp_session sftp){
    /* Copy a folder from the cluster recursively, copying only files that have changed. */

    static std::map<std::string, uint64_t> __cluster_file_size_map;
    static std::map<std::string, uint64_t> __cluster_file_time_map;
    sftp_dir directory;

    directory = sftp_opendir (sftp, source.c_str());
    while( ! sftp_dir_eof(directory) ){
        sftp_attributes file_attr = sftp_readdir(sftp, directory);
        if(file_attr != NULL){

            std::string subsource = source + "/" + file_attr->name;
            std::string subdestination = destination + "/" + file_attr->name;

            if( file_attr->type == SSH_FILEXFER_TYPE_DIRECTORY){

                if( strcmp(file_attr->name, ".") && strcmp(file_attr->name, "..") ){
                    getClusterFolder(subsource, subdestination, session, sftp);
                }

            } else {

                QString filePath = QDir::toNativeSeparators(subdestination.c_str());
                subdestination = filePath.toStdString();

                uint64_t old_time = __cluster_file_time_map[subdestination];
                uint64_t old_size = __cluster_file_size_map[subdestination];


                if( old_size != file_attr->size || old_time != file_attr->mtime ){
                    int rc = downloadAndVerifyClusterFile(subsource, subdestination, file_attr->size, session, sftp);

                    if(!rc){
                        __cluster_file_time_map[subdestination] = file_attr->mtime;
                        __cluster_file_size_map[subdestination] = file_attr->size;
                    } else {
                        QString filePath = QString(destination.c_str());
                        QFile file(filePath);
                        file.remove();
                    }
                }
            }
        }
    }
}




int clusterManager::folderFromCluster(std::string source, std::string destination){
    /* Copy a folder from the cluster using the address, username and password
     * provided by the clusterManager */
    sftp_session sftp;

    ssh_session session = createSSHSession(mAddress, mUsername, mPassword);
    sftp = createSFTPSession(session);

    getClusterFolder( source, destination, session, sftp);

    ssh_disconnect(session);
    ssh_free(session);

    return 0;
}

int clusterManager::fileFromCluster(std::string source, std::string destination){
    /* Copy a file from the cluster using the address, username and password
     * provided by the clusterManager */
    sftp_session sftp;

    ssh_session session = createSSHSession(mAddress, mUsername, mPassword);
    sftp = createSFTPSession(session);

    QString filePath = QDir::toNativeSeparators(destination.c_str());
    destination = filePath.toStdString();
    getClusterFile( source, destination, session, sftp);

    ssh_disconnect(session);
    ssh_free(session);

    return 0;
}

int sshVerifyPassword( QString address, QString username, QString password ){
    /* Check the cluster password by attempting to log in. */
    ssh_session session = ssh_new();
    int rc;
    if (session == NULL) {
        fprintf(stderr, "Failed to create ssh session\n");
        exit(-1);
    }

    ssh_options_set(session, SSH_OPTIONS_HOST, address.toStdString().c_str());
    ssh_options_set(session, SSH_OPTIONS_USER, username.toStdString().c_str());

    rc = ssh_connect(session);
    if (rc != SSH_AUTH_SUCCESS) {
      fprintf(stderr, "Failed to connect to the cluster\n");
      exit(-1);
    }

    rc = ssh_userauth_password(session, NULL, password.toStdString().c_str());
    if (rc != SSH_AUTH_SUCCESS)
    {
      fprintf(stderr, "Authentication failed\n");
      return 1;
    }

    return 0;
}

