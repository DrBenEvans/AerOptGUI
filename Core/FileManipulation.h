#ifndef FILEMANIPULATION_H
#define FILEMANIPULATION_H

#include <QString>


class FileManipulation
{
public:
    static bool emptyFolder(const QString& path);
    static bool copyFolder(const QString& source, const QString& dest);
    static bool copyFile(const QString& source, const QString& dest);
    static bool removeFolder(const QString& path);
};

#endif // FILEMANIPULATION_H
