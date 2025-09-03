#pragma once

class FairPrimaryGenerator;

class ISimGenerator
{
public:
    virtual ~ISimGenerator() = default;
    virtual FairPrimaryGenerator *build() = 0;
};