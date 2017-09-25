#ifndef PROCESSMANAGER_H
#define PROCESSMANAGER_H

#include <QProcess>
#include <QFileSystemWatcher>

class ProcessManager : public QProcess
{
    Q_OBJECT
public:
    ProcessManager();
    ~ProcessManager();
    void run(QString process, QString workDir, QString outputDir);

signals:
    void directoryChanged(const QString&);

public slots:
    void cleanupProcess();

private slots:
    void writeStdOutToLog();
    void writeStdErrToLog();
    void processFinished(int exitCode);

private:
    //Optimiser process
    QFileSystemWatcher mDirWatcher;
};

#endif // PROCESSMANAGER_H
