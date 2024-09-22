//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_FX_HPP
#define SZ_QOI_FX_HPP

#include <algorithm>
#include <cmath>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"
#include <symengine/expression.h>
#include <symengine/parser.h>
#include <symengine/symbol.h>
#include <symengine/derivative.h>
#include <symengine/eval.h> 
using SymEngine::Expression;
using SymEngine::Symbol;
using SymEngine::parse;
using SymEngine::diff;
using SymEngine::RealDouble;
using SymEngine::evalf;
using SymEngine::map_basic_num;
using SymEngine::down_cast;
using SymEngine::RCP;
using SymEngine::Basic;

namespace SZ {
    template<class T, uint N>
    class QoI_FX : public concepts::QoIInterface<T, N> {

    public:
        QoI_FX(T tolerance, T global_eb, std::string ff = "x^2") : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 14;
            std::cout<<"init 1 "<< std::endl;
            RCP<const Basic>  x = Symbol("x").rcp_from_this();
    
            f = Expression(ff);
             std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            df = f.diff(x);
             std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            ddf = df.diff(x);
             std::cout<<"init 4 "<< std::endl;
              std::cout<<"f: "<< f<<std::endl;
            std::cout<<"df: "<< df<<std::endl;
            std::cout<<"ddf: "<< ddf<<std::endl;
           // RCP<const Basic> result = evalf(df.subs(map_basic_basic({{x,RealDouble(2).rcp_from_this()}})),53, SymEngine::EvalfDomain::Real);
            RCP<const Realdouble> value(2);
            map_basic_num mbb=  {{x,value}};
            std::cout<<"init 5 "<< std::endl;
             double result = (double)df.subs(mbb); 
           
            std::cout<<"Eval res: "<<result<<std::endl;
            //SymEngine::RCP<const Basic> result = evalf(df,53, SymEngine::EvalfDomain::Real);
            //std::cout<< (down_cast<const RealDouble &>(*result)).as_double()<<std::endl;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            

            double a = fabs(cos(data));//datatype may be T
            double b = fabs(sin(data));
            T eb;
            if( b !=0)
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (a!=0)
                eb = tolerance/a;
            else 
                eb = global_eb;
           
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(sin(data) - sin(dec_data)) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

    private:
        T tolerance;
        T global_eb;
        //Symbol x;
        Expression f;
        Expression df;
        Expression ddf;
     
    };
}
#endif 
