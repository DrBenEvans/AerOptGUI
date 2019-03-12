#ifndef PROFILE_H
#define PROFILE_H

#include <list>
#include <QString>
#include "CustomTypes.h"

/**
 * @brief The Profile class containts a list of coordinate points for an optimisation profile.
 * The list of points can be set using setFile() and the profile points can be retrieved with getProfile().
 */
class Profile {

public:
    /**
     * @brief setFile Loads new profile data from the specified file and updates the list of
     * profile points to those of the new profile if it is valid.
     * @param filePath Filepath of the profile
     * @return true iff new profile file is valid and has been set
     */
    bool setFile(QString filePath);

    /**
     * @brief Profile::getProfile Returns the list of profile points
     * @return list of profile points
     */
    const ProfilePoints getProfile() const;

    /**
     * @brief getDisplayString Returns the file path for the profile.
     * @return Profile file filepath
     */
    QString getDisplayString();

private:
    /**
     * @brief clearProfile Clears the current profile data.
     */
    void clearProfile();

    /**
     * @brief addPoint Adds a point to the current profile.
     * @param x Dimension x location of point.
     * @param y Dimension y location of point.
     */
    void addPoint(const float& x, const float& y);

    /**
     * @brief checkProfileIntegrity Carries out several tests to confirm validity of profile.
     * @return True if integrity of profile is OK.
     * Checks the integrity of the current profile by carrying out a series of tests
     * including checkDuplicates, checkClockwise, checkNormalised,
     * checking the number of points exceeds the minimum value of 3
     */
    bool checkProfileIntegrity();
    /**
     * @brief getProfile
     * @return A reference to the profile as a series of x y pairs.
     */

    /**
     * @brief checkDuplicates Integrity check: Checks for duplicate points in the profile and removes them.
     * @return true
     */
    bool checkDuplicates();

    /**
     * @brief checkClockwise Integrity check: Check points are in correct clockwise sequence
     * @return true iff points are in correct sequence
     */
    bool checkClockwise();

    /**
     * @brief checkNormalised Integrity check: Scales points to x=1 normalisation x[0,...,1], y is scaled by same quantity
     * @return true iff points are scalable and have been scaled
     */
    bool checkNormalised();

    /**
     * @brief mProfile List of profile points with X and Y coordinates
     */
    ProfilePoints mProfile;

    /**
     * @brief mFilePath Filepath of raw profile file
     */
    QString mFilePath;
};

#endif // PROFILE_H
