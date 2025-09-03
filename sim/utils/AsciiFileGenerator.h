#include "ISimGenerator.h"

#include <FairPrimaryGenerator.h>
#include <R3BAsciiGenerator.h>
#include <TSystem.h>
#include <FairLogger.h>
#include <string>
#include <optional>

class AsciiFileGenerator : public ISimGenerator
{
public:
    explicit AsciiFileGenerator(std::string eventFile,
                                double dx = 0., double dy = 0., double dz = 0.)
        : fEventFile(std::move(eventFile)), fDx(dx), fDy(dy), fDz(dz)
    {
    }

    AsciiFileGenerator &withXYZ(double x, double y, double z)
    {
        fXYZ = XYZ{x, y, z};
        return *this;
    }

    FairPrimaryGenerator *build() override
    {
        if (gSystem->AccessPathName(fEventFile.c_str()))
        {
            FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
                                             "AsciiFileGenerator: event file '%s' not found (continuing anyway).",
                                             fEventFile.c_str());
        }

        auto *prim = new FairPrimaryGenerator();

        auto *gen = new R3BAsciiGenerator(fEventFile.c_str());
        gen->SetDxDyDz(fDx, fDy, fDz);

        if (fXYZ.has_value())
        {
            gen->SetXYZ(fXYZ->x, fXYZ->y, fXYZ->z);
        }

        prim->AddGenerator(gen);
        return prim;
    }

private:
    struct XYZ
    {
        double x{}, y{}, z{};
    };

    std::string fEventFile;
    double fDx{0.}, fDy{0.}, fDz{0.};
    std::optional<XYZ> fXYZ{};
};