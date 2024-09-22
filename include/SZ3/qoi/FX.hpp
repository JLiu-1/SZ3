//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_FX_HPP
#define SZ_QOI_FX_HPP

#include <algorithm>
#include <cmath>
#include <functional>
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
using SymEngine::symbol;
using SymEngine::parse;
using SymEngine::diff;
using SymEngine::RealDouble;
using SymEngine::Integer;
using SymEngine::evalf;
using SymEngine::map_basic_basic;
using SymEngine::down_cast;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::real_double;
using SymEngine::eval_double;

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
           // std::cout<<"init 1 "<< std::endl;
            
            Expression f;
            Expression df;
            Expression ddf;
             x = symbol("x");
    
            f = Expression(ff);
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            df = f.diff(x);
            // std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            ddf = df.diff(x);
            std::cout<<"f: "<< f<<std::endl;
            std::cout<<"df: "<< df<<std::endl;
            std::cout<<"ddf: "<< ddf<<std::endl;
  
            func = convert_expression_to_function(f, x);
            deri_1 = convert_expression_to_function(df, x);
            deri_2 = convert_expression_to_function(ddf, x);
            // std::cout<<"init 4 "<< std::endl;
              
           // RCP<const Basic> result = evalf(df.subs(map_basic_basic({{x,RealDouble(2).rcp_from_this()}})),53, SymEngine::EvalfDomain::Real);
           // RCP<const Symbol> value = symbol("2");
           // map_basic_basic mbb=  {{x,value}};
            //std::cout<<"init 5 "<< std::endl;
             //double result = (double)df.subs({{x,real_double(2)}}); 
           
           // std::cout<<"Eval res: "<<result<<std::endl;
            //SymEngine::RCP<const Basic> result = evalf(df,53, SymEngine::EvalfDomain::Real);
            //std::cout<< (down_cast<const RealDouble &>(*result)).as_double()<<std::endl;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            

            double a = fabs(deri_1(data));//datatype may be T
            double b = fabs(deri_2(data));
           // 
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and b !=0 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;
           // std::cout<<data<<" "<<a<<" "<<b<<" "<<eb<<" "<<global_eb<<std::endl; 
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(func(data) - func(dec_data)) < tolerance);
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

        inline double evaluate(const Expression & func, T val) const{
            
            return (double)func.subs({{x,real_double(val)}}); 

        } 
        std::function<double(T)> convert_expression_to_function(const Basic &expr, const RCP<const Symbol> &x) {
            //std::cout<<SymEngine::type_code_name(expr.get_type_code())<<std::endl;
            // x
            if (SymEngine::is_a<const SymEngine::Symbol>(expr)) {
                return [](T x_value) { return x_value; };
            }
            // c
            else if (SymEngine::is_a<const RealDouble>(expr) or SymEngine::is_a<const Integer>(expr)) {
                double constant_value = eval_double(expr);
                return [constant_value](T) { return constant_value; };
            }
            // +
            else if (SymEngine::is_a<SymEngine::Add>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) + right(x_value);
                };
            }
            // -
            /*
            else if (SymEngine::is_a<SymEngine::sub>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) - right(x_value);
                };
            }*/
            // *
            else if (SymEngine::is_a<SymEngine::Mul>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) * right(x_value);
                };
            }
            // /
            /*
            else if (SymEngine::is_a<SymEngine::div>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) / right(x_value);
                };
            }*/
            // pow
            else if (SymEngine::is_a<SymEngine::Pow>(expr)) {
                auto args = expr.get_args();
                auto base = convert_expression_to_function(Expression(args[0]), x);
                auto exponent = convert_expression_to_function(Expression(args[1]), x);
                return [base, exponent](T x_value) {
                    return std::pow(base(x_value), exponent(x_value));
                };
            }
            // sin
            else if (SymEngine::is_a<SymEngine::Sin>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::sin(arg(x_value));
                };
            }
            // cos
            else if (SymEngine::is_a<SymEngine::Cos>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::cos(arg(x_value));
                };
            }

            else if (SymEngine::is_a<SymEngine::Tan>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::tan(arg(x_value));
                };
            }

            else if (SymEngine::is_a<SymEngine::Sinh>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::sinh(arg(x_value));
                };
            }
            // cos
            else if (SymEngine::is_a<SymEngine::Cosh>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::cosh(arg(x_value));
                };
            }

            else if (SymEngine::is_a<SymEngine::Tanh>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::tanh(arg(x_value));
                };
            }

            else if (SymEngine::is_a<SymEngine::Sign>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return (x_value > 0) - (0 > x_value);
                };
            }
            //  log
            else if (SymEngine::is_a<SymEngine::Log>(expr)) {
                auto args = expr.get_args();
                auto arg = convert_expression_to_function(Expression(args[0]), x);

                if (args.size() == 2) { // base log
                    auto base = convert_expression_to_function(Expression(args[1]), x);
                    return [arg, base](T x_value) {
                        return std::log(arg(x_value)) / std::log(base(x_value));
                    };
                } else { // ln
                    return [arg](T x_value) {
                        return std::log(arg(x_value));
                    };
                }
            }

            throw std::runtime_error("Unsupported expression type");
        }
        RCP<const Symbol>  x;
        T tolerance;
        T global_eb;
        std::function<double(T)> func;
        std::function<double(T)> deri_1;
        std::function<double(T)> deri_2;
     
    };
}
#endif 
