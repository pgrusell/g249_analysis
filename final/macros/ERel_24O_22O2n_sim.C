#include "../utils/ERelCalculator.h"
#include <memory>
#include <vector>

void ERel_24O_22O2n_sim()
{
    std::unique_ptr<IERelCalculator> calc;
    calc = std::make_unique<ERelSimCalculator>();
    calc->getData({"sim_O23_O22n.root", "sim_O24_O23n.root"}, {1000080230, 1000080220, 2112});
    calc->Calculate();
    calc->plotErel();
}