#ifndef COLORMAPPER_H
#define COLORMAPPER_H

#include <QColor>
#include <vector>
#include <cmath>

typedef std::tuple<int,int,int,int> rgba;

class ColorMapper
{
public:
    ColorMapper();
    const QColor color(const float& value) const;
    void setMin(const float& min);
    void setMax(const float& max);

private:
    void transferFunction(const float& value, rgba& colour) const;
    float mMin;
    float mMax;
    std::vector<std::tuple<float,float,float,float>> colourmap;
};

#endif // COLORMAPPER_H
