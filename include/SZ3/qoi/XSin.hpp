//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_X_SIN_HPP
#define SZ_QOI_X_SIN_HPP

#include <algorithm>
#include <cmath>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_X_Sin : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Sin(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 10;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data)  {
            

            double a = fabs(cos(data));//datatype may be T
            double b = fabs(sin(data));
            T eb;
            if( b !=0)
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (a!=0)
                eb = tolerance/a;
            else 
                eb = global_eb;
            if (eb < global_eb)
                count ++;

            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) {
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
        int count =0;
    };
}
#endif 
