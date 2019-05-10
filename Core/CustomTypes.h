#ifndef CUSTOMTYPES_H
#define CUSTOMTYPES_H

#include <QMetaType>
#include <list>

//Abbreviate long type names

typedef std::vector<std::tuple<uint,uint,uint>> MeshConnectivities;
typedef std::vector<std::pair<uint,uint>> BConnectivities;
typedef std::vector<std::pair<float,float>> MeshPoints;

// 5 elements
// 0 = rho, 1 = u, 2 = v, 3 = energy, 4 = pressure,
typedef std::vector<std::tuple<float,float,float,float,float>> MeshResults;

typedef std::list<std::pair<float,float>> ProfilePoints;
Q_DECLARE_METATYPE(ProfilePoints)

#endif // CUSTOMTYPE_H
