#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

#include "SZ3/def.hpp"
#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/api/impl/SZImplOMP.hpp"
#include <cmath>

namespace SZ3 {
    template<class T, uint N>
    void SZ_compress_impl(Config &conf, const T *data, uchar *dst, size_t &dstLen) {
#ifndef _OPENMP
        conf.openmp=false;
#endif
        if (conf.openmp) {
            //dataCopy for openMP is handled by each thread
            SZ_compress_OMP<T, N>(conf, data, dst, dstLen);
        } else {
            std::vector<T> dataCopy(data, data + conf.num);
            SZ_compress_dispatcher<T, N>(conf, dataCopy.data(), dst, dstLen);
        }
    }


    template<class T, uint N>
    void SZ_decompress_impl(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {


#ifndef _OPENMP
        conf.openmp=false;
#endif
        if (conf.openmp) {
            SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
        } else {
            SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
        }
    }
}
#endif