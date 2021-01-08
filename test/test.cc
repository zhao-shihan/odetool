#include <fstream>
#include <cmath>
#include "odetool.hh"
using namespace odetool;

inline odefloat f(list_odefloat y, odefloat x) {
    return -sin(y[0]);
}

int main(int argc, char const* argv[]) {
    list_odefunction eqn{ f };
    list_odefloat x = DiscretizeX(0, 10, 0.01);
    list_list_odefloat y = DiscretizeY(1, x.size());
    y[0][0] = 1;
    RK45(eqn, x, y);

    std::ofstream fout("result.csv");
    if (fout.is_open()) {
        fout << "x,y\n";
        for (size_t i = 0; i < x.size(); ++i) {
            fout << x[i] << ',' << y[0][i] << std::endl;
        }
        fout.close();
    }

    return 0;
}

