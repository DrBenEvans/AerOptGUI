#ifndef CUSTOMTYPES_H
#define CUSTOMTYPES_H

#include <QMetaType>
#include "BoundaryPoint.h"
#include <list>

//Abbreviate long type names

typedef std::vector<std::tuple<uint,uint,uint>> MeshConnectivities;
typedef std::vector<std::pair<uint,uint>> BConnectivities;
typedef std::vector<std::pair<float,float>> MeshPoints;
typedef std::vector<std::tuple<float,float,float,float,float>> MeshResults;

typedef std::list<std::pair<float,float>> ProfilePoints;
Q_DECLARE_METATYPE(ProfilePoints)

#endif // CUSTOMTYPE_H
