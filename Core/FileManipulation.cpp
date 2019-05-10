#include "FileManipulation.h"
#include <QDir>
#include <QDebug>

bool FileManipulation::emptyFolder(const QString& path)
{
    bool r = true;

    QDir dir(path);

    // If the directory path doesn't already exist, create it
    if(!dir.exists()) {
        QDir().mkpath(path);
        qInfo() << QString("Creating directory: %1").arg(path);
        if(dir.exists())
            // Empty folder created at directory path - success
            return true;
        else
            return false;
            qInfo() << QString("Failed to create directory: %1").arg(path);

    // If the directory does exist, then delete every file within it
    } else {
        qInfo() << QString("Emptying directory: %1").arg(path);
        dir.setNameFilters(QStringList() << "*.*");
        dir.setFilter(QDir::Files);
        foreach(QString dirFile, dir.entryList())
        {
            r &= dir.remove(dirFile);
        }

        return r;
    }


}

bool FileManipulation::copyFolder(const QString& source, const QString& dest)
{
    bool r = true;

    //Set source folder
    QDir sourcePath(source);

    //Set/create destination folders
    QDir destPath(dest);
    destPath.mkpath(dest);
    qInfo() << QString("Make directory: %1").arg(dest);

    //Copy Input folders
    QStringList filesListSource = sourcePath.entryList(QDir::Files);
    QStringList filesListDest = destPath.entryList(QDir::Files);

    //Delete destination if exists
    foreach (QString filename, filesListDest)
    {
        QFileInfo f(filename);
        QString destname = f.fileName();
        QDir dest = QDir(destPath.absolutePath() + QDir::separator() + destname).absolutePath();
        destname = QDir::toNativeSeparators(dest.path());

        if (QFile::exists(destname))
        {
            QFile::remove(destname);
            qInfo() << QString("Delete existing file: %1").arg(destname);
        }
    }

    //Copy files from source to destination
    foreach (QString filename, filesListSource)
    {
        QFileInfo f(filename);

        QString destname = f.fileName();
        QString sourcename = f.fileName();

        QDir dest = QDir(destPath.absolutePath() + QDir::separator() + destname).absolutePath();
        destname = QDir::toNativeSeparators(dest.path());

        QDir source = QDir(sourcePath.absolutePath() + QDir::separator() + sourcename).absolutePath();
        sourcename = QDir::toNativeSeparators(source.path());

        r &= copyFile( sourcename, destname );
    }

    return r;
}

bool FileManipulation::copyFile(const QString& source, const QString& dest)
{
    bool r = true;

    if (QFile::exists(dest))
    {
        QFile::remove(dest);
        qInfo() << QString("Remove existing file: %1").arg(dest);
    }

    r &= QFile::copy(source, dest);
    if(r) {
        qInfo() << QString("Copy file: %1 -> %2").arg(source).arg(dest);
    } else {
        qWarning() << QString("Copy file FAILED: %1 -> %2").arg(source).arg(dest);
    }

    return r;
}

bool FileManipulation::removeFolder(const QString& path)
{
    bool r = true;
    QDir dir(path);

    if (dir.exists(path))
    {
        // Delete all folders and files at path
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst))
        {

            if (info.isDir()) {
                // If folder then make recursive call
                r = removeFolder(info.absoluteFilePath());
            } else {
                r = QFile::remove(info.absoluteFilePath());
            }

            // If a file could not be deleted, then return error
            if (!r) {
                return r;
            }
        }
        r = dir.rmdir(path);
    }

    if(r) {
        qInfo() << QString("Removed folder: %1").arg(path);
    } else {
        qWarning() << QString("Folder remove failed: %1").arg(path);
    }

    return r;
}
