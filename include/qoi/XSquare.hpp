//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_X_SQUARE_HPP
#define SZ_QOI_X_SQUARE_HPP

#include <algorithm>
#include "QoI.hpp"
#include "def.hpp"

namespace SZ {
    template<class T>
    class QoI_X_Square : public concepts::QoIInterface<T, 1> {

    public:
        QoI_X_Square(T tolerance, size_t num_elements, T global_eb) : 
            tolerance(tolerance),
            num_rest(num_elements),
            global_eb(global_eb) {
                printf("global_eb = %.4f\n", global_eb);
            }

        T interpret_eb(T data) const {
            if (data == 0) return std::min(sqrt(tolerance), global_eb);
            T eb = - data + sqrt(data * data + tolerance);
            if ((data < 0) && (data * data - tolerance >= 0)){
                if (eb > - data - sqrt(data * data - tolerance))
                    eb = - data - sqrt(data * data - tolerance);
            }
            return std::min(eb, global_eb);
        }

        void update_tolerance(T data, T dec_data){}

        void print(){}

    private:
        T tolerance;
        T global_eb;
        int num_rest;
    };
}
#endif 