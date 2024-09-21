//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_X_CUBIC_HPP
#define SZ_QOI_X_CUBIC_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_X_Cubic : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Cubic(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 9;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            
            if (data == 0)
                return global_eb;

            double b = fabs(6*data);
            double a = fabs(0.5*b*data);//datatype may be T

            
            T eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            if (eb < global_eb)
                count +=1
            //T eb = data >= 0 ? std::cbrt(data*data*data+tolerance)-data : data - std::cbrt(data*data*data-tolerance);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(data*data - dec_data*dec_data) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){std::cout<<"Count: "<<count<<std::endl;}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

    private:
        T tolerance;
        T global_eb;
    };
}
#endif 
