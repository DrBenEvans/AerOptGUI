#include "ProcessManager.h"
#include <QDebug>
#include <QDir>

ProcessManager::ProcessManager()
{
    mDirWatcher.setParent(this);

    connect(this, SIGNAL(readyReadStandardOutput()), this, SLOT(writeStdOutToLog()));
    connect(this, SIGNAL(readyReadStandardError()), this, SLOT(writeStdErrToLog()));

    connect(&mDirWatcher, &QFileSystemWatcher::directoryChanged, this, &ProcessManager::directoryChanged);

    auto finished = QOverload<int>::of(&QProcess::finished);
    connect(this, finished, this, &ProcessManager::processFinished);
}

ProcessManager::~ProcessManager()
{
    cleanupProcess();
}

void ProcessManager::run(QString process, QString workDir, QString watchPath) {
    mDirWatcher.removePaths( mDirWatcher.directories() );

    // make directory before watching
    QDir outputDir(watchPath);
    if (!outputDir.exists()) {
        outputDir.mkpath(".");
    }

    mDirWatcher.addPath(watchPath);
    setWorkingDirectory(workDir);
    start(process);
}

void ProcessManager::processFinished(int exitStatus) {
    cleanupProcess();
}

void ProcessManager::cleanupProcess() {
    kill();
    mDirWatcher.removePaths( mDirWatcher.directories() );
}

void ProcessManager::writeStdOutToLog()
{
    setReadChannel(QProcess::StandardOutput);
    while(canReadLine()) {
        emit stdOut(QString(readLine()));
    }
}

void ProcessManager::writeStdErrToLog()
{
    setReadChannel(QProcess::StandardError);
    while(canReadLine()) {
        emit stdOut(QString(readLine()));
    }
}
