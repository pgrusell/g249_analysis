#pragma once
#include <string>

class FairRunSim;

struct SimPaths
{
    std::string vmcworkdir;
    std::string geomdir;
    std::string pardir;
    std::string configdir;
};

class ISimGeometry
{
public:
    virtual ~ISimGeometry() = default;
    virtual void configure(FairRunSim &run, const SimPaths &paths) = 0;
};