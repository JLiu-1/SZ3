

#ifndef SZ3_LINEAR_QUANTIZER_OMP_HPP
#define SZ3_LINEAR_QUANTIZER_OMP_HPP
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

template <class T>
class LinearQuantizerOMP : public concepts::QuantizerOMPInterface<T, int> {
   public:
    LinearQuantizerOMP() : error_bound(1), double_error_bound(2), double_error_bound_reciprocal(0.5), radius(32768) {}

    LinearQuantizerOMP(double eb, int r = 32768) : error_bound(eb),double_error_bound(2*eb), double_error_bound_reciprocal(0.5 / eb), radius(r) {
        assert(eb != 0);
    }

    double get_eb() const { return error_bound; }

    void set_eb(double eb) {
        error_bound = eb;
        double_error_bound = 2 * eb;
        double_error_bound_reciprocal = 1.0 / double_error_bound;
    }

    std::pair<int, int> get_out_range() const override { return std::make_pair(0, radius * 2); }

    // quantize the data with a prediction value, and returns the quantization index and the decompressed data
    // int quantize(T data, T pred, T& dec_data);
    ALWAYS_INLINE int quantize_and_overwrite(T &data, T pred, size_t data_idx) override {

        T diff = data - pred;
        int quant_index = std::llrint(diff * this->double_error_bound_reciprocal);
        if (std::abs(quant_index) < this->radius ) {
            //if (diff < 0) 
            //    quant_index = -quant_index;
            T decompressed_data = pred + quant_index * this->double_error_bound;

   
            // if data is NaN, the error is NaN, and NaN <= error_bound is false
            if (fabs(decompressed_data - data) <= this->error_bound) {
                data = decompressed_data;
                
               // auto quant_index_shifted = ;
                
                return this->radius + quant_index;
            } else {
                if(data_idx == 22448644)
              std::cout<<"bb"<<std::endl;
                save_unpred(data, data_idx);
                return 0;
            }
        } else {
            if(data_idx == 22448644)
              std::cout<<"pp"<<std::endl;
            save_unpred(data, data_idx);
            return 0;
        }
    }

    // recover the data using the quantization index
    ALWAYS_INLINE T recover(T pred, int quant_index) override {
        if (quant_index ) {
            return recover_pred(pred, quant_index);
        } else {
            //return recover_unpred();
            return T(0);
        }
    }

    ALWAYS_INLINE T recover_pred(T pred, int quant_index) {
        return pred + (quant_index - this->radius) * this->double_error_bound;
    }

    //ALWAYS_INLINE T recover_unpred() { return unpred[index++]; }

    ALWAYS_INLINE int save_unpred(T ori, size_t data_idx){
        #ifdef _OPENMP
        #pragma omp critical 
        {
            unpred.push_back(ori);
            unpred_idx.push_back(data_idx);
        }
        #else
            unpred.push_back(ori);
            unpred_idx.push_back(data_idx);
        #endif
        return 0;
    }

    void unpack_unpred(T *data) const{

        #ifdef _OPENMP
           #pragma omp parallel for
        #endif
        for(size_t i = 0; i < unpred.size(); i++)
            data[unpred_idx[i]] = unpred[i];

    }



    size_t size_est() { return unpred.size() * sizeof(T); }

    void save(unsigned char *&c) const override {
        write(uid, c);
        write(this->error_bound, c);
        write(this->radius, c);
        size_t unpred_size = unpred.size();
        write(unpred_size, c);
        if (unpred_size > 0) {
            assert (unpred_size == unpred_idx.size());
            write(unpred.data(), unpred_size, c);
            write(unpred_idx.data(), unpred_size, c);       
        }
    }

    void load(const unsigned char *&c, size_t &remaining_length) override {
        uchar uid_read;
        read(uid_read, c, remaining_length);
        if (uid_read != uid) {
            throw std::invalid_argument("LinearQuantizer uid mismatch");
        }
        double eb;
        read(eb, c, remaining_length);
        set_eb(eb);
        read(this->radius, c, remaining_length);
        size_t unpred_size = 0;
        read(unpred_size, c, remaining_length);
        if (unpred_size > 0) {
            unpred.resize(unpred_size);
            read(unpred.data(), unpred_size, c, remaining_length);
            unpred_idx.resize(unpred_size);
            read(unpred_idx.data(), unpred_size, c, remaining_length);
        }
        index = 0;
    }

    void print() override {
        printf("[LinearQuantizer] error_bound = %.8G, radius = %d, unpred = %zu\n", error_bound, radius, unpred.size());
    }

   private:
    std::vector<T> unpred;
    std::vector<size_t> unpred_idx;
    size_t index = 0;  // used in decompression only
    uchar uid = 0b11;

    double error_bound;
    double double_error_bound;
    double double_error_bound_reciprocal;
    int radius;  // quantization interval radius
};

}  // namespace SZ3
#endif