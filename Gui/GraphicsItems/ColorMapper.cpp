#include "ColorMapper.h"
#include <tuple>

ColorMapper::ColorMapper()
{
    colourmap.emplace_back( 0,         0,        0,        0.423499   );
    colourmap.emplace_back( 0.0668895, 0,        0.119341, 0.529244   );
    colourmap.emplace_back( 0.133779,  0,        0.238697, 0.634974   );
    colourmap.emplace_back( 0.200669,  0,        0.346853, 0.687877   );
    colourmap.emplace_back( 0.267559,  0,        0.450217, 0.718135   );
    colourmap.emplace_back( 0.334448,  0,        0.553552, 0.664836   );
    colourmap.emplace_back( 0.401338,  0,        0.651087, 0.51931    );
    colourmap.emplace_back( 0.468227,  0.115846, 0.724788, 0.35285    );
    colourmap.emplace_back( 0.535117,  0.326772, 0.781201, 0.140185   );
    colourmap.emplace_back( 0.602006,  0.522759, 0.79852,  0.0284581  );
    colourmap.emplace_back( 0.668897,  0.703166, 0.788678, 0.00885023 );
    colourmap.emplace_back( 0.735786,  0.845121, 0.751141, 0          );
    colourmap.emplace_back( 0.802675,  0.955734, 0.690822, 0          );
    colourmap.emplace_back( 0.869565,  0.995407, 0.56791,  0.0618448  );
    colourmap.emplace_back( 0.936455,  0.987716, 0.403403, 0.164858   );
    colourmap.emplace_back( 1,         0.980407, 0.247105, 0.262699   );
}

void ColorMapper::setMin(const float& min) {
    mMin = min;
}

void ColorMapper::setMax(const float& max) {
    mMax = max;
}

const QColor ColorMapper::color(const float& value) const
{
    rgba&& colour = std::make_tuple(255, 255, 255, 255);

    if (mMin <  mMax)
    {
        //Normalise pressure value
        float norm = (value - mMin) / (mMax - mMin);

        //compute rgba from normalised value
        transferFunction(norm, colour);
    }

    QColor qcolour(std::get<0>(colour),
                   std::get<1>(colour),
                   std::get<2>(colour),
                   std::get<3>(colour));

    return qcolour;
}

void ColorMapper::transferFunction(const float& value, rgba& colour) const
{
    uint i;
    for (i = 0; i < colourmap.size()-1; i++)
    {
        if (value < std::get<0>(colourmap.at(i+1))) break;
    }

    if (i < colourmap.size()-1)
    {

        float diff = std::get<0>(colourmap.at(i+1)) - std::get<0>(colourmap.at(i+0));
        float v = value - std::get<0>(colourmap.at(i+0));
        float linear = v / diff;

        float diffr = std::get<1>(colourmap.at(i+1)) - std::get<1>(colourmap.at(i+0));
        float r = (diffr * linear) + std::get<1>(colourmap.at(i+0));

        float diffg = std::get<2>(colourmap.at(i+1)) - std::get<2>(colourmap.at(i+0));
        float g = (diffg * linear) + std::get<2>(colourmap.at(i+0));

        float diffb = std::get<3>(colourmap.at(i+1)) - std::get<3>(colourmap.at(i+0));
        float b = (diffb * linear) + std::get<3>(colourmap.at(i+0));

        std::get<0>(colour) = std::round(255.0 * r);
        std::get<1>(colour) = std::round(255.0 * g);
        std::get<2>(colour) = std::round(255.0 * b);
    }
}
