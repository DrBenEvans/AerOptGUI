#ifndef PROFILE_H
#define PROFILE_H

#include <list>
#include <QString>

class Profile
{
public:
    Profile(QString filePath);
    const std::list<std::pair<float,float>> getProfile() const;
    QString getDisplayString();

private:
    /**
     * @brief clearProfile
     * Clears the current profile data.
     */
    void clearProfile();
    /**
     * @brief addPoint
     * @param x Dimension x location of point.
     * @param y Dimension y location of point.
     * Adds points to the current profile.
     */
    void addPoint(const float& x, const float& y);
    /**
     * @brief checkProfileIntegrity
     * @return True if integrity of profile is OK.
     * Checks the integrity of the current profile.
     */
    bool checkProfileIntegrity();
    /**
     * @brief getProfile
     * @return A reference to the profile as a series of x y pairs.
     */

    bool checkDuplicates();
    bool checkClockwise();
    bool checkNormalised();
    std::list<std::pair<float,float>> mProfile;
    QString mFilePath;
};

#endif // PROFILE_H
