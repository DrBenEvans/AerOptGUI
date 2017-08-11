#ifndef SIMULATION_H
#define SIMULATION_H

#include "CustomTypes.h"
#include "Mesh.h"

class Simulation
{
public:
    Simulation();

    int generation() const;
    void setGeneration(int gen);
    int id() const;
    void setId(int id);
    std::shared_ptr<Mesh> mesh();
    void setOptimisationId(int id);
    ProfilePoints profilePoints() const;
    void setProfilePoints(ProfilePoints profile);

private:
    int mGeneration;
    int mId;
    int mOptimisationId;
    ProfilePoints mProfilePoints;
    std::shared_ptr<Mesh> mMesh;
};

#endif // SIMULATION_H
