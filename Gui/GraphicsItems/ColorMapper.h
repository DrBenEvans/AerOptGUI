#ifndef COLORMAPPER_H
#define COLORMAPPER_H

#include <QColor>
#include <vector>
#include <cmath>

/**
 * @brief RGBA Tuple for storing a QColor
 */
typedef std::tuple<int,int,int,int> rgba;

/**
 * @brief The ColorMapper class is used to convert scaled pressure readings to a corresponding colour where colouring a mesh.
 */
class ColorMapper
{
public:
    /**
     * @brief ColorMapper Constructor Method.
     * Establishes colourmap values.
     */
    ColorMapper();

    /**
     * @brief color Return a QColor corresponding to a given pressure value.
     * If the scale has not been set, then a default colour is returned.
     * @param value Pressure
     * @return QColor corresponding to pressure value.
     */
    const QColor color(const float& value) const;

    /**
     * @brief setMin Set the minimum pressure value for colour scaling.
     * @param min Minimum pressure value.
     */
    void setMin(const float& min);

    /**
     * @brief setMax Set the maximum pressure value for colour scaling.
     * @param max Maximum pressure value.
     */
    void setMax(const float& max);

private:
    /**
     * @brief transferFunction Normalise RGB parameters of a QColor by a normalised pressure value.
     * @param value Normalised pressure value.
     * @param colour RGBA tuple to be Normalised
     */
    void transferFunction(const float& value, rgba& colour) const;
    /**
     * @brief mMin Minimum pressure value.
     */
    float mMin;
    /**
     * @brief mMax Maximum pressure value.
     */
    float mMax;

    /**
     * @brief colourmap Reference vector of different mesh colours.
     */
    std::vector<std::tuple<float,float,float,float>> colourmap;
};

#endif // COLORMAPPER_H
