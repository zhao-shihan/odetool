#include <fstream>
#include <cmath>
#include "odetool.hh"
using namespace odetool;

const size_t object_num = 1000;

const size_t dof = 3 * object_num;

inline odefloat Velocity(list_odefloat y, odefloat, size_t i) {
    return y[dof + i];
}

inline odefloat GravityX(list_odefloat y, odefloat, size_t eqn_index) {
    odefloat ax = 0;
    odefloat d = 0;
    size_t i = eqn_index - dof;
    for (size_t j = 0; j < i; j += 3) {
        d = sqrt((y[i] - y[j]) * (y[i] - y[j]) +
            (y[i + 1] - y[j + 1]) * (y[i + 1] - y[j + 1]) +
            (y[i + 2] - y[j + 2]) * (y[i + 2] - y[j + 2]));
        ax += (y[j] - y[i]) / (d * d * d);
    }
    for (size_t j = i + 3; j < dof; j += 3) {
        d = sqrt((y[i] - y[j]) * (y[i] - y[j]) +
            (y[i + 1] - y[j + 1]) * (y[i + 1] - y[j + 1]) +
            (y[i + 2] - y[j + 2]) * (y[i + 2] - y[j + 2]));
        ax += (y[j] - y[i]) / (d * d * d);
    }
    return ax;
}

inline odefloat GravityY(list_odefloat y, odefloat, size_t eqn_index) {
    odefloat ay = 0;
    odefloat d = 0;
    size_t i = eqn_index - dof;
    for (size_t j = 1; j < i; j += 3) {
        d = sqrt((y[i - 1] - y[j - 1]) * (y[i - 1] - y[j - 1]) +
            (y[i] - y[j]) * (y[i] - y[j]) +
            (y[i + 1] - y[j + 1]) * (y[i + 1] - y[j + 1]));
        ay += (y[j] - y[i]) / (d * d * d);
    }
    for (size_t j = i + 3; j < dof; j += 3) {
        d = sqrt((y[i - 1] - y[j - 1]) * (y[i - 1] - y[j - 1]) +
            (y[i] - y[j]) * (y[i] - y[j]) +
            (y[i + 1] - y[j + 1]) * (y[i + 1] - y[j + 1]));
        ay += (y[j] - y[i]) / (d * d * d);
    }
    return ay;
}

inline odefloat GravityZ(list_odefloat y, odefloat, size_t eqn_index) {
    odefloat az = 0;
    odefloat d = 0;
    size_t i = eqn_index - dof;
    for (size_t j = 2; j < i; j += 3) {
        d = sqrt((y[i - 2] - y[j - 2]) * (y[i - 2] - y[j - 2]) +
            (y[i - 1] - y[j - 1]) * (y[i - 1] - y[j - 1]) +
            (y[i] - y[j]) * (y[i] - y[j]));
        az += (y[j] - y[i]) / (d * d * d);
    }
    for (size_t j = i + 3; j < dof; j += 3) {
        d = sqrt((y[i - 2] - y[j - 2]) * (y[i - 2] - y[j - 2]) +
            (y[i - 1] - y[j - 1]) * (y[i - 1] - y[j - 1]) +
            (y[i] - y[j]) * (y[i] - y[j]));
        az += (y[j] - y[i]) / (d * d * d);
    }
    return az;
}

int main(int argc, char const* argv[]) {
    list_odefunction eqn(2 * dof);
    for (size_t i = 0; i < dof; ++i) {
        eqn[i] = Velocity;
    }
    for (size_t i = dof; i < 2 * dof; ++i) {
        eqn[i] = GravityX;
        ++i;
        eqn[i] = GravityY;
        ++i;
        eqn[i] = GravityZ;
    }

    list_odefloat t = DiscretizeX(0, 10, 0.0001);

    list_list_odefloat y = DiscretizeY(eqn.size(), t.size());
    for (size_t i = 0; i < dof; ++i) {
        if (i % 3 == 0) {
            *(y[i].begin()) = (odefloat)i;
        } else {
            *(y[i].begin()) = 0.0;
        }
    }
    for (size_t i = dof; i < 2 * dof; ++i) {
        if (i % 3 == 1) {
            if ((i / 3) % 2 == 0) {
                *(y[i].begin()) = 0.2;
            } else {
                *(y[i].begin()) = -0.2;
            }
        } else {
            *(y[i].begin()) = 0.0;
        }
    }

    RK45(eqn, t, y);

    std::ofstream fout("result.csv");
    if (fout.is_open()) {
        fout << 't';
        for (size_t i = 0; i < dof; ++i) {
            fout << ",x" << i / 3 + 1;
            ++i;
            fout << ",y" << i / 3 + 1;
            ++i;
            fout << ",z" << i / 3 + 1;
        }
        for (size_t i = 0; i < dof; ++i) {
            fout << ",vx" << i / 3 + 1;
            ++i;
            fout << ",vy" << i / 3 + 1;
            ++i;
            fout << ",vz" << i / 3 + 1;
        }
        fout << std::endl;
        for (size_t n = 0; n < t.size(); n += 10) {
            fout << t[n];
            for (size_t i = 0; i < 2 * dof; ++i) {
                fout << ',' << y[i][n];
            }
            fout << std::endl;
        }
        fout.close();
    }

    return 0;
}

