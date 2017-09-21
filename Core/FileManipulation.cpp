#include "FileManipulation.h"
#include <QDir>

bool FileManipulation::emptyFolder(const QString& path)
{
    bool r = true;

    QDir dir(path);
    dir.setNameFilters(QStringList() << "*.*");
    dir.setFilter(QDir::Files);
    foreach(QString dirFile, dir.entryList())
    {
        dir.remove(dirFile);
    }

    return r;
}

bool FileManipulation::copyFolder(const QString& source, const QString& dest)
{
    bool r = true;

    //Set source folder
    QDir sourcePath(source);

    //Set/create destination folders
    QDir destPath(dest);
    destPath.mkpath(dest);

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

        r &= QFile::copy( sourcename, destname );
    }

    return r;
}

bool FileManipulation::copyFile(const QString& source, const QString& dest)
{
    bool r = true;

    if (QFile::exists(dest))
    {
        QFile::remove(dest);
    }

    r &= QFile::copy(source, dest);

    return r;
}

bool FileManipulation::removeFolder(const QString& path)
{
    bool r = true;
    QDir dir(path);

    if (dir.exists(path))
    {
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst))
        {
            if (info.isDir())
            {
                r = removeFolder(info.absoluteFilePath());
            }
            else
            {
                r = QFile::remove(info.absoluteFilePath());
            }

            if (!r)
            {
                return r;
            }
        }
        r = dir.rmdir(path);
    }
    return r;
}
