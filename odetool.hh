#ifndef ODETOOL_H
#define ODETOOL_H

#include <vector>

#define odelist vector
using std::odelist;

typedef double                                         odefloat;

typedef odelist<odefloat(*)(list_odefloat, odefloat)>  list_odeeqn;
typedef odelist<odefloat>                              list_odefloat;
typedef odelist<list_odefloat>                         list_list_odefloat;

typedef list_odeeqn::iterator                          it_list_odeeqn;
typedef list_list_odefloat::iterator                   it_list_list_odefloat;
typedef list_odefloat::iterator                        it_list_odefloat;
typedef odelist<it_list_odefloat>                      list_it_list_odefloat;
typedef odelist<it_list_odefloat>::iterator            it_list_it_list_odefloat;

class ODEtool {
private:
    list_odeeqn eqn_;
    it_list_odeeqn it_eqn_;

    size_t n;

    list_odefloat list_x_;
    it_list_odefloat it_list_x_;

    list_list_odefloat list_list_y_;
    it_list_list_odefloat it_list_list_y_;
    list_it_list_odefloat list_it_list_y_;
    it_list_it_list_odefloat it_list_it_list_y_;

    list_list_odefloat list_list_k_;
    it_list_list_odefloat it_list_list_k_;
    list_it_list_odefloat list_it_list_k_;
    it_list_it_list_odefloat it_list_it_list_k_;

public:
    ODEtool(const list_odeeqn& eqn, odefloat x_begin, odefloat x_end, const odefloat& h);
    ~ODEtool();

    void RK45(const list_odefloat& init_cond);
};

ODEtool::ODEtool(const list_odeeqn& eqn, odefloat x_begin, odefloat x_end, const odefloat& h) :
    eqn_(eqn),
    it_eqn_(eqn_.begin()),
    list_list_y_(eqn.size()),
    it_list_list_y_(list_list_y_.begin()),
    list_it_list_y_(eqn.size()),
    it_list_it_list_y_(list_it_list_y_.begin()),
    list_list_k_(eqn.size()),
    it_list_list_k_(list_list_k_.begin()) {
    if (x_begin > x_end) {
        odefloat temp = x_begin;
        x_begin = x_end;
        x_end = temp;
    }
    n = (x_end - x_begin) / h + 1;

    list_x_.resize(n);
    it_list_x_ = list_x_.begin();
    *it_list_x_ = x_begin;
    for (++it_list_x_; it_list_x_ != list_x_.end(); ++it_list_x_) {
        *it_list_x_ = *(it_list_x_ - 1) + h;
    }
    for (; it_list_list_y_ != list_list_y_.end(); ++it_list_list_y_) {
        it_list_list_y_->resize(n);
        *it_list_it_list_y_ = it_list_list_y_->begin();
        ++it_list_it_list_y_;
    }
}

ODEtool::~ODEtool() {}

void ODEtool::RK45(const list_odefloat& init_cond) {
    list_it_list_k_.resize(4);
    it_list_it_list_k_ = list_it_list_k_.begin();
    for (it_list_list_k_ = list_list_k_.begin(); it_list_list_k_ != list_list_k_.end(); ++it_list_list_k_) {
        it_list_list_k_->resize(4);
    }
    list_odefloat::const_iterator it_init_cond = init_cond.begin();
    for (it_list_it_list_y_ = list_it_list_y_.begin(); it_list_it_list_y_ != list_it_list_y_.end();
        ++it_list_it_list_y_) {
        **it_list_it_list_y_ = *it_init_cond;
        ++it_init_cond;
    }

    for (it_list_x_ = list_x_.begin(); it_list_x_ != list_x_.end(); ++it_list_x_) {

    }
}

#endif // ODETOOL_H
