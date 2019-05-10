#ifndef FILEMANIPULATION_H
#define FILEMANIPULATION_H

#include <QString>

/**
 * @brief The FileManipulation class provides methods for generic file manipulation tasks including:
 * * Empty the contents of a folder
 * * Copy a folder from a source to a destination
 * * Copy a file from a souce to a destination
 * * Delete a folder
 */
class FileManipulation
{
public:
    /**
     * @brief emptyFolder Creates an empty folder at the specified filepath.
     * If the directory already exists then delete all files in it, otherwise try to create the directory.
     * @param path The filepath for the empty directory.
     * @return true iff empty directory now exists at the given path
     */
    static bool emptyFolder(const QString& path);

    /**
     * @brief copyFolder copy a folder and its contents from its source to a new destination.
     * If a folder already exists at the destination with the same name as the source, it will be overwritten.
     * @param source The filepath of the folder to be copied
     * @param dest The filepath where the folder should be copied to.
     * @return true iff folder copied successfully
     */
    static bool copyFolder(const QString& source, const QString& dest);

    /**
     * @brief copyFile Copies a file from its source to a specified destination.
     * If a file of the same name already exists at the destination, then this is replaced with the source file.
     * @param source The filepath of the file to be copied.
     * @param dest Destination filepath where file should be copied to
     * @return true iff file copied successfullt.
     */
    static bool copyFile(const QString& source, const QString& dest);

    /**
     * @brief removeFolder Recursively delete all files and folders at the specified directory path
     * @warning This is a dangerous function, use with caution.
     * @param path The filepath to the
     * @return true iff folder and all files contained within have been deleted.
     */
    static bool removeFolder(const QString& path);
};

#endif // FILEMANIPULATION_H
