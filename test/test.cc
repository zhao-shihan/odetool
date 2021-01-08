#include "odetool.hh"

int main(int argc, char const* argv[]) {
    list_odeeqn eqn(3);
    ODEtool ode45(eqn, (odefloat)0.1, (odefloat)1.0, (odefloat)0.01);
    return 0;
}

