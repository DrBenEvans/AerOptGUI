#include "Simulation.h"

Simulation::Simulation()
{

}

int Simulation::generation() const {
    return mGeneration;
}

void Simulation::setGeneration(int gen) {
    mGeneration = gen;
}

int Simulation::id() const {
    return mId;
}

void Simulation::setId(int id) {
    mId = id;
}

std::shared_ptr<Mesh> Simulation::mesh() {
    return mMesh;
}

void Simulation::setOptimisationId(int id) {
    mOptimisationId = id;
}

ProfilePoints Simulation::profilePoints() const {
    return mProfilePoints;
}

void Simulation::setProfilePoints(ProfilePoints profile) {
    mProfilePoints = profile;
}
