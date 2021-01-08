#ifndef ODETOOL_H
#define ODETOOL_H

#include <vector>
#include <omp.h>

namespace odetool {

    typedef double                                            odefloat;
    typedef std::vector<odefloat>                             list_odefloat;
    typedef std::vector<std::vector<odefloat>>                list_list_odefloat;
    typedef std::vector<odefloat(*)(list_odefloat, odefloat)> list_odefunction;

    inline list_odefloat DiscretizeX(const odefloat& x_begin, const odefloat& x_end, const odefloat& h) {
        size_t n = (x_end - x_begin) / h + 1;
        list_odefloat x(n);
#pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            x[i] = x_begin + i * h;
        }
        return x;
    }

    inline list_list_odefloat DiscretizeY(const size_t& eqn_size, const size_t& x_size) {
        list_list_odefloat y(eqn_size);
#pragma omp parallel for
        for (list_list_odefloat::iterator it = y.begin(); it != y.end(); ++it) {
            it->resize(x_size);
        }
        return y;
    }

    void RK45(const list_odefunction& eqn, const list_odefloat& x, list_list_odefloat& y) {
        list_list_odefloat k(4);
#pragma omp parallel for
        for (list_list_odefloat::iterator it = k.begin(); it != k.end(); ++it) {
            it->resize(eqn.size());
        }
        list_list_odefloat temp_y(4);
#pragma omp parallel for
        for (list_list_odefloat::iterator it = temp_y.begin(); it != temp_y.end(); ++it) {
            it->resize(eqn.size());
        }

        for (size_t n = 0; n < x.size() - 1; ++n) {
            odefloat h = x[n + 1] - x[n];

#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                temp_y[0][i] = y[i][n];
            }
#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                k[0][i] = eqn[i](temp_y[0], x[n]) * h;
            }

#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                temp_y[1][i] = temp_y[0][i] + 0.5 * k[0][i];
            }
#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                k[1][i] = eqn[i](temp_y[1], x[n] + 0.5 * h) * h;
            }

#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                temp_y[2][i] = temp_y[0][i] + 0.5 * k[1][i];
            }
#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                k[2][i] = eqn[i](temp_y[2], x[n] + 0.5 * h) * h;
            }

#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                temp_y[3][i] = temp_y[0][i] + k[2][i];
            }
#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                k[3][i] = eqn[i](temp_y[3], x[n] + h) * h;
            }

#pragma omp parallel for
            for (size_t i = 0; i < eqn.size(); ++i) {
                y[i][n + 1] = y[i][n] + 0.166666666666667 * (k[0][i] + 2.0 * k[1][i] + 2.0 * k[2][i] + k[3][i]);
            }
        }
    }

} // namespace odetool

#endif // ODETOOL_H

