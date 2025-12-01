#ifdef _OPENMP

#ifndef SZ3_INTERPOLATION_DECOMPOSITION_OMP_HPP
#define SZ3_INTERPOLATION_DECOMPOSITION_OMP_HPP
#include <omp.h>
#include <cmath>
#include <cstring>
#include <immintrin.h>
#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/quantizer/Quantizer_Omp.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/BlockwiseIterator.hpp"


namespace SZ3 {
template <class T, uint N, class QuantizerOMP>
class InterpolationDecomposition_OMP : public concepts::DecompositionInterface<T, int, N> {
   public:
    InterpolationDecomposition_OMP(const Config &conf, QuantizerOMP quantizer) : quantizer(quantizer) {
        static_assert(std::is_base_of<concepts::QuantizerOMPInterface<T, int>, QuantizerOMP>::value,
                      "must implement the quantizer interface");
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        init();
        /*
        if(N==3){
            std::vector<int> quant_inds_vec_reordered(num_elements);
            auto* __restrict dst = quant_inds_vec_reordered.data();
            auto const* __restrict src = quant_inds_vec.data();

            #pragma omp parallel for collapse(3) schedule(static)
            for (size_t x0 = 0; x0 < original_dimensions[0]; ++x0) {
                for (size_t y0 = 0; y0 < original_dimensions[1]; ++y0) {
                    for (size_t z0 = 0; z0 < original_dimensions[2]; ++z0) {
                        const size_t idx = x0 * original_dim_offsets[0] + y0 * original_dim_offsets[1] + z0;

                        size_t x = x0, y = y0, z = z0;

                        const int max_level = interp_level - 1;
                        unsigned level = 0;
                        if (max_level > 0) {
                            unsigned tzx = x ? __builtin_ctzll(x) : 32;
                            unsigned tzy = y ? __builtin_ctzll(y) : 32;
                            unsigned tzz = z ? __builtin_ctzll(z) : 32;
                            unsigned l = std::min<unsigned>(max_level,
                                                            std::min(tzx, std::min(tzy, tzz)));
                            level = l;
                            x >>= l;
                            y >>= l;
                            z >>= l;
                        }

                        size_t reordered_idx =
                            x * reduced_dim_offsets[level][0] +
                            y * reduced_dim_offsets[level][1] +
                            z;

                        if (level < (unsigned)max_level) {

                            size_t t0 = ((x + 1) >> 1) * reduced_dim_offsets[level + 1][0];
                            size_t t1 = ((y + 1) >> 1) * reduced_dim_offsets[level + 1][1];
                            size_t t2 = ((z + 1) >> 1);

                            reordered_idx += level_prefix[level]
                                           - t0
                                           - ((x % 2 == 0) ? t1 : 0)
                                           - ((x % 2 == 0 && y % 2 == 0) ? t2 : 0);
                        }

                        dst[idx] = src[reordered_idx];
                    }
                }
            }
           
            quant_inds.clear();
            quant_inds.shrink_to_fit();
            quant_inds= std::move( quant_inds_vec_reordered);



        }*/


        auto default_nThreads = omp_get_max_threads();
        //std::cout<<"max threads: "<<default_nThreads<<std::endl;

        size_t max_usable_threads = default_nThreads;
        for (uint i = 1; i < N; ++i) 
            max_usable_threads = std::min(max_usable_threads, original_dimensions[i]);
        omp_set_num_threads(max_usable_threads);

        nThreads = omp_get_max_threads(); // for safety
        //std::cout<<"used threads: "<<nThreads<<std::endl;


        buffer_len =  (max_dim + 2 * AVX_256_parallelism - max_dim % AVX_256_parallelism) ;
        size_t total_buffer_len = buffer_len * nThreads;
        interp_buffer_1 = new T[total_buffer_len];

        interp_buffer_2 = new T[total_buffer_len];
        interp_buffer_3 = new T[total_buffer_len];
        interp_buffer_4 = new T[total_buffer_len];

        pred_buffer = new T[total_buffer_len];
        #pragma omp parallel for
        for(size_t i =0;i<total_buffer_len;++i)
            pred_buffer[i] = interp_buffer_1[i] = interp_buffer_2[i] = interp_buffer_3[i] = interp_buffer_4[i] = T(0);

        #pragma omp parallel for
        for(size_t i =0; i<num_elements;++i)
            dec_data[i] = T(0);
        quantizer.unpack_unpred(dec_data);

        this->quant_inds = quant_inds.data();
        double eb = quantizer.get_eb();
        //visited.resize(num_elements);
        auto start_level = interp_level;
        if (anchor_stride == 0) {                                               // check whether used anchor points
            *dec_data += quantizer.recover(0, this->quant_inds[0]);  // no anchor points
        } else {
           // recover_anchor_grid(dec_data);  // recover anchor points, not needed because because all outliers were previously unpacked.
            start_level--;
        }
        
        for (int level = start_level; level > 0; level--) {
            // set level-wise error bound
            if (eb_alpha < 0) {
                if (level >= 3) {
                    quantizer.set_eb(eb * eb_ratio);
                } else {
                    quantizer.set_eb(eb);
                }
            } else if (eb_alpha >= 1) {
                double cur_ratio = pow(eb_alpha, level - 1);
                if (cur_ratio > eb_beta) {
                    cur_ratio = eb_beta;
                }
                quantizer.set_eb(eb / cur_ratio);
            }
            size_t stride = 1U << (level - 1);
            auto interp_block_size = blocksize * stride;
            auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                dec_data, std::begin(original_dimensions), std::end(original_dimensions), interp_block_size, 0);
            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();
            for (auto block = inter_begin; block != inter_end; ++block) {
                auto end_idx = block.get_global_index();
                for (uint i = 0; i < N; ++i) {
                    end_idx[i] += interp_block_size;
                    if (end_idx[i] > original_dimensions[i] - 1) {
                        end_idx[i] = original_dimensions[i] - 1;
                    }
                }
                interpolation(
                    dec_data, block.get_global_index(), end_idx, interpolators[interp_id],
                    [&](const std::array<size_t,N>&idx_array,size_t idx, T &d, T pred, int level) { d += quantizer.recover(pred, quant_inds[index_mapping(idx_array,idx,level)]);},// no need to use idx. the outliers will be unpacked separately (todo).
                    direction_sequence_id, stride);
            }
        }
        quantizer.postdecompress_data();
       

      


        delete [] interp_buffer_1;
        delete [] interp_buffer_2;
        delete [] interp_buffer_3;
        delete [] interp_buffer_4;
        delete [] pred_buffer;

        omp_set_num_threads(default_nThreads);
        return dec_data;
    }

    // compress given the error bound
    std::vector<int> compress(const Config &conf, T *data) override {




        std::copy_n(conf.dims.begin(), N, original_dimensions.begin());
        
        interp_id = conf.interpAlgo;
        direction_sequence_id = conf.interpDirection;
        anchor_stride = conf.interpAnchorStride;
        blocksize = 256;  // a empirical value. Can be very large but not helpful
        eb_alpha = conf.interpAlpha;
        eb_beta = conf.interpBeta;

        init();



        auto default_nThreads = omp_get_max_threads();
        //std::cout<<"max threads: "<<default_nThreads<<std::endl;

        size_t max_usable_threads = default_nThreads;
        for (uint i = 1; i < N; ++i) 
            max_usable_threads = std::min(max_usable_threads, original_dimensions[i]);
        omp_set_num_threads(max_usable_threads);

        nThreads = omp_get_max_threads(); // for safety
        //std::cout<<"used threads: "<<nThreads<<std::endl;
   

        buffer_len =  (max_dim + 2 * AVX_256_parallelism - max_dim % AVX_256_parallelism) ;
        size_t total_buffer_len = buffer_len * nThreads;
        interp_buffer_1 = new T[total_buffer_len];
        interp_buffer_2 = new T[total_buffer_len];
        interp_buffer_3 = new T[total_buffer_len];
        interp_buffer_4 = new T[total_buffer_len];

        pred_buffer = new T[total_buffer_len];
        #pragma omp parallel for
        for(size_t i =0;i < total_buffer_len;++i)
            pred_buffer[i] = interp_buffer_1[i] = interp_buffer_2[i] = interp_buffer_3[i] = interp_buffer_4[i] = T(0);
       

        std::vector<int> quant_inds_vec(num_elements);
        //visited.resize(num_elements);
        quant_inds = quant_inds_vec.data();
        double eb = quantizer.get_eb();
        auto start_level = interp_level;
        if (anchor_stride == 0) {  // check whether to use anchor points
            quant_inds[0] = quantizer.quantize_and_overwrite(*data, 0, 0);  // no
        } else {
            build_anchor_grid(data);  // losslessly saving anchor points
            start_level--;
        }

        for (int level = start_level; level > 0; level--) {
            double cur_eb = eb;
            // set level-wise error bound
            if (eb_alpha < 0) {
                if (level >= 3) {
                    cur_eb = eb * eb_ratio;
                } else {
                    cur_eb = eb;
                }
            } else if (eb_alpha >= 1) {
                double cur_ratio = pow(eb_alpha, level - 1);
                if (cur_ratio > eb_beta) {
                    cur_ratio = eb_beta;
                }
                cur_eb = eb / cur_ratio;
            }
            quantizer.set_eb(cur_eb);
            size_t stride = 1U << (level - 1);

            auto interp_block_size = blocksize * stride;

            auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                data, std::begin(original_dimensions), std::end(original_dimensions), interp_block_size, 0);

            auto inter_begin = inter_block_range->begin();
            auto inter_end = inter_block_range->end();

            for (auto block = inter_begin; block != inter_end; ++block) {
                auto end_idx = block.get_global_index();
                for (uint i = 0; i < N; ++i) {
                    end_idx[i] += interp_block_size;
                    if (end_idx[i] > original_dimensions[i] - 1) {
                        end_idx[i] = original_dimensions[i] - 1;
                    }
                }

                interpolation(
                    data, block.get_global_index(), end_idx, interpolators[interp_id],
                    [&](const std::array<size_t,N>&idx_array, size_t idx, T &d, T pred, int level) {
                        
                        quant_inds[index_mapping(idx_array, idx,level)] = (quantizer.quantize_and_overwrite(d, pred, idx));
                       
                    },
                    direction_sequence_id, stride);
            }
        }
        quantizer.set_eb(eb);
        quantizer.postcompress_data();


        




        delete [] interp_buffer_1;
        delete [] interp_buffer_2;
        delete [] interp_buffer_3;
        delete [] interp_buffer_4;
        delete [] pred_buffer;

        omp_set_num_threads(default_nThreads);
        /*
        if(N==3){
            std::vector<int> quant_inds_vec_reordered(num_elements);
            auto* __restrict dst = quant_inds_vec_reordered.data();
            auto const* __restrict src = quant_inds_vec.data();

            #pragma omp parallel for collapse(3) schedule(static)
            for (size_t x0 = 0; x0 < original_dimensions[0]; ++x0) {
                for (size_t y0 = 0; y0 < original_dimensions[1]; ++y0) {
                    for (size_t z0 = 0; z0 < original_dimensions[2]; ++z0) {
                        const size_t idx = x0 * original_dim_offsets[0] + y0 * original_dim_offsets[1] + z0;

                        size_t x = x0, y = y0, z = z0;

                        const int max_level = interp_level - 1;
                        unsigned level = 0;
                        if (max_level > 0) {
                            unsigned tzx = x ? __builtin_ctzll(x) : 32;
                            unsigned tzy = y ? __builtin_ctzll(y) : 32;
                            unsigned tzz = z ? __builtin_ctzll(z) : 32;
                            unsigned l = std::min<unsigned>(max_level,
                                                            std::min(tzx, std::min(tzy, tzz)));
                            level = l;
                            x >>= l;
                            y >>= l;
                            z >>= l;
                        }

                        size_t reordered_idx =
                            x * reduced_dim_offsets[level][0] +
                            y * reduced_dim_offsets[level][1] +
                            z;

                        if (level < (unsigned)max_level) {

                            size_t t0 = ((x + 1) >> 1) * reduced_dim_offsets[level + 1][0];
                            size_t t1 = ((y + 1) >> 1) * reduced_dim_offsets[level + 1][1];
                            size_t t2 = ((z + 1) >> 1);

                            reordered_idx += level_prefix[level]
                                           - t0
                                           - ((x % 2 == 0) ? t1 : 0)
                                           - ((x % 2 == 0 && y % 2 == 0) ? t2 : 0);
                        }

                        dst[reordered_idx] = src[idx];
                    }
                }
            }

         
           
            return quant_inds_vec_reordered;

        }*/


        return quant_inds_vec;
    }

    void save(uchar *&c) override {
        write(original_dimensions.data(), N, c);
        write(blocksize, c);
        write(interp_id, c);
        write(direction_sequence_id, c);
        write(anchor_stride, c);
        write(eb_alpha, c);
        write(eb_beta, c);
        

        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(original_dimensions.data(), N, c, remaining_length);
        read(blocksize, c, remaining_length);
        read(interp_id, c, remaining_length);
        read(direction_sequence_id, c, remaining_length);
        read(anchor_stride, c, remaining_length);
        read(eb_alpha, c, remaining_length);
        read(eb_beta, c, remaining_length);
        

        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    void init() {
       
        quant_index = 0;
        assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
        assert((anchor_stride & anchor_stride - 1) == 0 && "Anchor stride should be 0 or 2's exponentials");
        num_elements = 1;
        interp_level = -1;
	    bool use_anchor = false;
        max_dim = 1;
        for (uint i = 0; i < N; ++i) {
            if (interp_level < ceil(log2(original_dimensions[i]))) {
                interp_level = static_cast<int>(ceil(log2(original_dimensions[i])));
            }
    	    if (original_dimensions[i] > anchor_stride)
    	        use_anchor = true;
            num_elements *= original_dimensions[i];
            max_dim = std::max(max_dim,original_dimensions[i]);
        }

        if (!use_anchor)
            anchor_stride = 0;
        if (anchor_stride > 0) {
            int max_interpolation_level = static_cast<int>(log2(anchor_stride)) + 1;
            if (max_interpolation_level <= interp_level) {
                interp_level = max_interpolation_level;
            }
        }

        original_dim_offsets[N - 1] = 1;
        for (int i = N - 2; i >= 0; --i) {
            original_dim_offsets[i] = original_dim_offsets[i + 1] * original_dimensions[i + 1];
        }

        dim_sequences = std::vector<std::array<int, N>>();
        auto sequence = std::array<int, N>();
        for (uint i = 0; i < N; ++i) {
            sequence[i] = i;
        }
        do {
            dim_sequences.push_back(sequence);
        } while (std::next_permutation(sequence.begin(), sequence.end()));

        if constexpr (N==3){
       
            auto d_size = original_dimensions;
            reduced_dim_offsets.resize(interp_level );
            level_prefix.resize(interp_level , 0);
            
            int level = 0;
            while(level < interp_level){
                //grid_leaps[level][0] = 1;
                reduced_dim_offsets[level][2] = 1;
                reduced_dim_offsets[level][1] = d_size[2];
                reduced_dim_offsets[level][0] = d_size[1] * d_size[2];
              
               
                
                
                if(level + 1 < interp_level ){
                    d_size[0] = (d_size[0] + 1) >> 1;
                    d_size[1] = (d_size[1] + 1) >> 1;
                    d_size[2] = (d_size[2] + 1) >> 1;
                    level_prefix[level] = d_size[0] *  d_size[1] * d_size[2];
                }
                ++level;
            }  
        }
    }

    ALWAYS_INLINE size_t index_mapping(const std::array<size_t,N>&idx_array, const size_t &idx,const int &level){
        if constexpr (N!=3){
            return idx;
        }
        const size_t dim0 = original_dim_offsets[0];
        const size_t dim1 = original_dim_offsets[1];

        // 拆 idx -> (x0, y0, z0)，每个 dim 只做一次除法
        //size_t x = idx / dim0;
        //size_t r  = idx - x * dim0;  // r = idx % dim0

       // size_t y = r / dim1;
       // size_t z = r - y * dim1;    // z0 = r % dim1
       // return x * original_dim_offsets[0] + y * original_dim_offsets[1] + z; 
       // size_t x = x0, y = y0, z = z0;

        size_t x = idx_array[0],y=idx_array[1],z=idx_array[2];
        const int max_level = interp_level - 1;

        x >>= level;
        y >>= level;
        z >>= level;
        /*
        unsigned level = 0;
        if (max_level > 0) {
            size_t v = x | y | z;

            unsigned tz;
            //if (v != 0) {
                // 注意：__builtin_ctzll 参数不能为 0, but we will not process idx in this function when it is 0 
                tz = static_cast<unsigned>(__builtin_ctzll(v));
           // } else {
            //    return 0;
            //}

            // level = min(max_level, tz)
            level = (tz < max_level) ? tz : max_level;

           
        }*/

        size_t reordered_idx =
            x * reduced_dim_offsets[level][0] +
            y * reduced_dim_offsets[level][1] +
            z;
            
        if (level < (unsigned)max_level) {

            size_t t0 = ((x + 1) >> 1) * reduced_dim_offsets[level + 1][0];

            reordered_idx += level_prefix[level]- t0;

            if( (x & 1) == 0 ){
                reordered_idx -= ((y + 1) >> 1) * reduced_dim_offsets[level + 1][1];
                if( (y & 1) == 0){
                    reordered_idx -= ((z + 1) >> 1);
                }
            }

            
         
        }

         /*
            reordered_idx += level_prefix[level]
                           - t0
                           - ((x % 2 == 0) ? t1 : 0)
                           - ((x % 2 == 0 && y % 2 == 0) ? t2 : 0);*/
        return reordered_idx;
    }

    void build_anchor_grid(T *data) {  // store anchor points. steplength: anchor_stride on each dimension 
        std::array<size_t, N> strides;
        std::array<size_t, N> begins{0};
        std::fill(strides.begin(), strides.end(), anchor_stride);
        foreach_omp
            <T, N>(data, 0, begins, original_dimensions, strides, original_dim_offsets,
                   [&](T *d,const std::array<size_t, N> idx_array) { auto idx = d - data; quant_inds[idx] = quantizer.save_unpred( *d, idx);});
    }

    void recover_anchor_grid(T *data) {  // recover anchor points. steplength: anchor_stride on each dimension
        std::array<size_t, N> strides;
        std::array<size_t, N> begins{0};
        std::fill(strides.begin(), strides.end(), anchor_stride);
        foreach_omp
            <T, N>(data, 0, begins, original_dimensions, strides, original_dim_offsets, [&](T *d,const std::array<size_t, N> idx_array) {
                *d = quantizer.recover_unpred(d - data);
                //quant_index++;
            });
    }

    /**
     * Do interpolations along a certain dimension, and move through that dimension only.
     * This is the original API, described in the ICDE'21 paper.
     * @tparam QuantizeFunc
     * @param data
     * @param begin
     * @param end
     * @param stride
     * @param interp_func
     * @param quantize_func
     * @return
     */
    template <class QuantizeFunc>
    double interpolation_1d(T *data, size_t begin, size_t end, size_t stride, const std::string &interp_func,
                            QuantizeFunc &&quantize_func) {
        size_t n = (end - begin) / stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0;

        std::array<size_t,N> idx;

        size_t stride3x = 3 * stride;
        size_t stride5x = 5 * stride;
        if (interp_func == "linear" || n < 5) {
            // if (pb == PB_predict_overwrite) {
            #pragma omp parallel for
            for (size_t i = 1; i < n - 1; i += 2) {
                T *d = data + begin + i * stride;
                quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),0);
            }
            if (n % 2 == 0) {
                T *d = data + begin + (n - 1) * stride;
                if (n < 4) {
                    quantize_func(idx, d-data, *d, *(d - stride),0);
                } else {
                    quantize_func(idx, d-data, *d, interp_linear1(*(d - stride3x), *(d - stride)),0);
                }
            }
            // }
        } else {
            T *d;
            size_t i;
            #pragma omp parallel for
            for (i = 3; i < n - 3; i += 2) {
                d = data + begin + i * stride;
                quantize_func(idx, d-data,*d,
                              interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),0);
            }
            d = data + begin + stride;
            quantize_func(idx, d-data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),0);

            d = data + begin + i * stride;
            quantize_func(idx, d-data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),0);
            if (n % 2 == 0) {
                d = data + begin + (n - 1) * stride;
                quantize_func(idx, d-data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)),0);
            }
        }

        return predict_error;
    }

    /**
     * Do all interpolations along a certain dimension on the full data grid. Moving on the fastest-dim.
     * This is the new API, described in the SIGMOD'24 paper.
     * @tparam QuantizeFunc
     * @param data
     * @param begin_idx
     * @param end_idx
     * @param direction
     * @param strides
     * @param math_stride
     * @param interp_func
     * @param quantize_func
     * @return
     */
    template <class QuantizeFunc>
    double interpolation_1d_fastest_dim_first(T *data, const std::array<size_t, N> &begin_idx,
                                              const std::array<size_t, N> &end_idx, const size_t &direction,
                                              std::array<size_t, N> &strides, const size_t &math_stride,
                                              const std::string &interp_func, QuantizeFunc &&quantize_func) {
        for (size_t i = 0; i < N; ++i) {
            if (end_idx[i] < begin_idx[i]) return 0;
        }
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        int level = __builtin_ctzll(math_stride);
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0.0;
        size_t offset = 0;
        size_t stride = math_stride * original_dim_offsets[direction];
        std::array<size_t, N> begins = begin_idx, ends = end_idx, dim_offsets;
        /*
        for (size_t i = 0; i < N; ++i) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            dim_offsets[i] = original_dim_offsets[i];
            offset += original_dim_offsets[i] * begin_idx[i];
        }*/
        dim_offsets[direction] = stride;
        size_t stride2x = 2 * stride;
        if (interp_func == "linear" ) {
            begins[direction] = math_begin_idx + math_stride;
            ends[direction] = math_begin_idx + math_stride * (n - 1);
            //strides[direction] = 2;
            foreach_omp
                <T, N>(data, 0, begins, ends, strides, original_dim_offsets,
                       [&](T *d, const std::array<size_t,N> &idx) { quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level); });
            if (n % 2 == 0) {
                begins[direction] = ends[direction];
                ends[direction] += math_stride;
                foreach_omp //todo: this is infficient when direction = 0
                    <T, N>(data, 0, begins, ends, strides, original_dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (n < 3)
                            quantize_func(idx, d-data, *d, *(d - stride),level);
                        else
                            quantize_func(idx, d-data, *d, interp_linear1(*(d - stride2x), *(d - stride)),level);
                    });
            }
        } else {
            size_t stride3x = 3 * stride;
            size_t i_start = 3;
            begins[direction] = math_begin_idx + i_start * math_stride;
            ends[direction] = (n >= 3) ?  math_begin_idx + (n - 3) * math_stride : math_begin_idx;
            //strides[direction] = 2;
            foreach_omp 
                <T, N>(data, 0, begins, ends, strides, original_dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                    quantize_func(idx, d-data, *d,
                                  interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),level);
            });
            std::vector<size_t> boundaries;
            boundaries.push_back(1);
            if (n % 2 == 1 && n > 3) {
                boundaries.push_back(n - 2);
            }
            if (n % 2 == 0 && n > 4) {
                boundaries.push_back(n - 3);
            }
            if (n % 2 == 0 && n > 2) {
                boundaries.push_back(n - 1);
            }
            for (auto boundary : boundaries) {
                begins[direction] = math_begin_idx + boundary * math_stride;
                ends[direction] = begins[direction] + math_stride;
                
                foreach_omp //todo: this is infficient when direction = 0
                    <T, N>(data, 0, begins, ends, strides, original_dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (boundary >= 3) {
                            if (boundary + 3 < n)
                                quantize_func(
                                    idx, d-data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),level);
                            else if (boundary + 1 < n)
                                quantize_func(idx, d-data, *d,
                                              interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),level);
                            else
                                quantize_func(idx, d-data, *d, interp_linear1(*(d - stride3x), *(d - stride)),level);
                        } else {
                            if (boundary + 3 < n)
                                quantize_func(idx, d-data, *d,
                                              interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),level);
                            
                            else if (boundary + 1 < n)
                               
                                quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level);
                            
                            else
                                quantize_func(idx, d-data, *d, *(d - stride),level);
                            
                        }
                    });
            }
        }
        return predict_error;
    }
    
    void avx_interp_cubic(const T * a,const T * b,const T * c,const T * d,T * p, const size_t &len){
       // assert(len <= max_dim);
         constexpr bool is_float  = std::is_same_v<T, float>;
        constexpr bool is_double = std::is_same_v<T, double>;

        size_t i = 0;

        if constexpr (is_float) {
            const size_t step = AVX_256_parallelism;
            const __m256 nine  = _mm256_set1_ps(9.0f);
            const __m256 factor = _mm256_set1_ps(1.0f / 16.0f);

            for (; i  < len; i += step) {
                __m256 va = _mm256_loadu_ps(a + i);
                __m256 vb = _mm256_loadu_ps(b + i);
                __m256 vc = _mm256_loadu_ps(c + i);
                __m256 vd = _mm256_loadu_ps(d + i);

                 __m256 sum = _mm256_add_ps(vb, vc); 
                 sum = _mm256_mul_ps(sum, nine); 
                 sum = _mm256_sub_ps(sum, va); 
                sum = _mm256_sub_ps(sum, vd); 

               
                      
                sum = _mm256_mul_ps(sum, factor);        

                _mm256_storeu_ps(p + i, sum);
            }
        }
        else if constexpr (is_double) {
            const size_t step = AVX_256_parallelism;
            const __m256d nine  = _mm256_set1_pd(9.0);
            const __m256d factor = _mm256_set1_pd(1.0 / 16.0);

            for (; i  < len; i += step) {
                __m256d va = _mm256_loadu_pd(a + i);
                __m256d vb = _mm256_loadu_pd(b + i);
                __m256d vc = _mm256_loadu_pd(c + i);
                __m256d vd = _mm256_loadu_pd(d + i);

                __m256d sum = _mm256_add_pd(vb, vc); 
                 sum = _mm256_mul_pd(sum, nine); 
                 sum = _mm256_sub_pd(sum, va); 
                sum = _mm256_sub_pd(sum, vd); 

               
                      
                sum = _mm256_mul_pd(sum, factor);    
                _mm256_storeu_pd(p + i, sum);
            }
        }
        /*
        for (; i < len; i++) {
            p[i] = (-a[i] + T(9) * b[i] + T(9) * c[i] - d[i]) / T(16);
        }*/


    }

    void avx_interp_cubic_1D(const T * buf, T * p, const size_t &len){
       // assert(len <= max_dim);
         constexpr bool is_float  = std::is_same_v<T, float>;
        constexpr bool is_double = std::is_same_v<T, double>;
        if(len == 1)
            return;

        auto odd_len = len / 2;
        auto even_len = len - odd_len;

        if(even_len < 2)
            p[0] = (buf[0]);

        else if(even_len < 3)
            p[0] = interp_linear(buf[0], buf[1]);
        else
            p[0] = interp_quad_1(buf[0], buf[1], buf[2]) ;
        size_t i = 0;

        if constexpr (is_float) {
            const size_t step = AVX_256_parallelism;
            const __m256 nine  = _mm256_set1_ps(9.0f);
            const __m256 factor = _mm256_set1_ps(1.0f / 16.0f);

            for (; i + 3  <= even_len; i += step) { // 3 is not AVX_256_parallelism - 1 !!
                __m256 va = _mm256_loadu_ps(buf + i);
                __m256 vb = _mm256_loadu_ps(buf + i + 1);
                __m256 vc = _mm256_loadu_ps(buf + i + 2);
                __m256 vd = _mm256_loadu_ps(buf + i + 3);

                 __m256 sum = _mm256_add_ps(vb, vc); 
                 sum = _mm256_mul_ps(sum, nine); 
                 sum = _mm256_sub_ps(sum, va); 
                sum = _mm256_sub_ps(sum, vd);                       
                sum = _mm256_mul_ps(sum, factor);        

                _mm256_storeu_ps(p + i + 1, sum);
            }
        }
        else if constexpr (is_double) {
            const size_t step = AVX_256_parallelism;
            const __m256d nine  = _mm256_set1_pd(9.0);
            const __m256d factor = _mm256_set1_pd(1.0 / 16.0);

            for (; i + 3 <= even_len; i += step) { // 3 is not AVX_256_parallelism - 1 !!
                __m256d va = _mm256_loadu_pd(buf + i);
                __m256d vb = _mm256_loadu_pd(buf + i + 1);
                __m256d vc = _mm256_loadu_pd(buf + i + 2);
                __m256d vd = _mm256_loadu_pd(buf + i + 3);

                __m256d sum = _mm256_add_pd(vb, vc); 
                 sum = _mm256_mul_pd(sum, nine); 
                 sum = _mm256_sub_pd(sum, va); 
                sum = _mm256_sub_pd(sum, vd); 
                sum = _mm256_mul_pd(sum, factor);    
                _mm256_storeu_pd(p + i + 1, sum);
            }
        }
        if(odd_len > 1){
            if(odd_len < even_len){//the only boundary is p[len- 1] 
                //odd_len < even_len so even_len > 2
                p[odd_len - 1] = interp_quad_2(buf[even_len - 3], buf[even_len - 2], buf[even_len - 1]);

            }
            else{//the boundary points are is p[len -2 ] and p[len -1 ]
                if(odd_len > 2){ //len - 2
                 //odd_len = even_len so even_len > 2
                    p[odd_len - 2] = interp_quad_2(buf[even_len - 3],  buf[even_len - 2], buf[even_len - 1]);
                }
                //len -1
                //odd_len = even_len so even_len > 1
                    p[odd_len - 1] = interp_linear1(buf[even_len - 2], buf[even_len - 1]);
                

            }
        }
        /*
        for (; i < len; i++) {
            p[i] = (-a[i] + T(9) * b[i] + T(9) * c[i] - d[i]) / T(16);
        }*/


    }

    template <class QuantizeFunc>
    double interpolation_1d_simd_3d_x(T *data, const std::array<size_t, N> &begin_idx,
                                              const std::array<size_t, N> &end_idx, const size_t &direction,
                                              std::array<size_t, N> &strides, const size_t &math_stride,
                                              const std::string &interp_func, QuantizeFunc &&quantize_func) {
        assert(direction==0  && N==3);
        for (size_t i = 0; i < N; ++i) {
            if (end_idx[i] < begin_idx[i]) return 0;
        }
        int level = __builtin_ctzll(math_stride);
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        if (n <= 1) {
            return 0;
        }
        std::array<size_t,N> idx;
        double predict_error = 0.0;
        size_t offset = 0;
        size_t stride = math_stride * original_dim_offsets[direction];
        std::array<size_t, N> begins, ends, dim_offsets;
        for (size_t i = 0; i < N; i++) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            dim_offsets[i] = original_dim_offsets[i];
            offset += original_dim_offsets[i] * begin_idx[i];
        }
        dim_offsets[direction] = stride;
        size_t stride2x = 2 * stride;
        if (interp_func == "linear") {
            begins[direction] = 1;
            ends[direction] = n - 1;
            strides[direction] = 2;
            foreach_omp
                <T, N>(data, offset, begins, ends, strides, dim_offsets,
                       [&](T *d, const std::array<size_t,N> &idx) { quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level); });
            if (n % 2 == 0) {
                begins[direction] = n - 1;
                ends[direction] = n;
                foreach_omp //todo: this is infficient when direction = 0
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (n < 3)
                            quantize_func(idx, d-data, *d, *(d - stride),level);
                        else
                            quantize_func(idx, d-data, *d, interp_linear1(*(d - stride2x), *(d - stride)),level);
                    });
            }
        } else {
            size_t stride3x = 3 * stride;
            size_t i_start = 3;
            begins[direction] = i_start;
            ends[direction] = (n >= 3) ? (n - 3) : 0;
            strides[direction] = 2;
            size_t vector_len = ends[2] > begins[2] ? (ends[2]-begins[2]-1)/strides[2] + 1 : 0;

            #pragma omp parallel for
            for (size_t j = begins[1]; j < ends[1]; j += strides[1]) {
                auto tid = omp_get_thread_num();
                auto buffer_offset = buffer_len * tid;
                auto cur_buffer_1 = interp_buffer_1 + buffer_offset;
                auto cur_buffer_2 = interp_buffer_2 + buffer_offset;
                auto cur_buffer_3 = interp_buffer_3 + buffer_offset;
                auto cur_buffer_4 = interp_buffer_4 + buffer_offset; 
                auto cur_pred_buffer = pred_buffer + buffer_offset;
                

               
                
                for(size_t i = begins[0]; i < ends[0]; i += strides[0]){

                    auto cur_ij_offset = offset + i * dim_offsets[0] + j * dim_offsets[1];
                    size_t buffer_idx = 0;
                    if( i == begins[0]){
                        for (size_t k = begins[2]; k < ends[2]; k += strides[2]) {
                            auto cur_offset =  cur_ij_offset + k;
                            //if (visited[cur_offset- stride3x]==0 or visited[cur_offset+ stride3x]==0)
                            //   std::cout<<"e1 "<<i<<" "<<j<<" "<<k<<" "<<stride3x<<std::endl;
                           cur_buffer_1[buffer_idx] = data[cur_offset -  stride3x];
                           //cur_buffer_1[buffer_idx] = data[0];
                            cur_buffer_2[buffer_idx] = data[cur_offset - stride];
                           // cur_buffer_2[buffer_idx] = data[0];
                          cur_buffer_3[buffer_idx] = data[cur_offset + stride];
                           //cur_buffer_3[buffer_idx] = data[0];
                            cur_buffer_4[buffer_idx] = data[cur_offset +  stride3x];
                          //  cur_buffer_4[buffer_idx] = data[0];
                            ++buffer_idx;

                        }
                    }
                    
                    else{
                        auto temp_buffer = cur_buffer_1;
                        cur_buffer_1 = cur_buffer_2;
                        cur_buffer_2 = cur_buffer_3;
                        cur_buffer_3 = cur_buffer_4;
                        cur_buffer_4 = temp_buffer;

                        buffer_idx = 0;
                        for (size_t k = begins[2]; k < ends[2]; k += strides[2]) {
                            auto cur_offset =  cur_ij_offset + stride3x + k;
                           // if (visited[cur_offset]==0 )
                            //   std::cout<<"e2 "<<i<<" "<<j<<" "<<k<<" "<<cur_offset<<std::endl;
                           
                           cur_buffer_4[buffer_idx++] = data[cur_offset];

                            //cur_buffer_4[buffer_idx++] = data[0];

                        }
                    }
                    
                    avx_interp_cubic(cur_buffer_1,cur_buffer_2,cur_buffer_3,cur_buffer_4, cur_pred_buffer, vector_len);
                    buffer_idx = 0;
                    for (size_t k = begins[2]; k < ends[2]; k += strides[2]){
                        auto pred = cur_pred_buffer[buffer_idx++];
                        auto d = data + cur_ij_offset + k;
                      // if (d-data < 0 || d-data>=num_elements)
                      //      std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                        quantize_func(idx, d-data, *d,pred,level);

                    }
                    
                }
            }
            std::vector<size_t> boundaries;
            boundaries.push_back(1);
            if (n % 2 == 1 && n > 3) {
                boundaries.push_back(n - 2);
            }
            if (n % 2 == 0 && n > 4) {
                boundaries.push_back(n - 3);
            }
            if (n % 2 == 0 && n > 2) {
                boundaries.push_back(n - 1);
            }
            for (auto boundary : boundaries) {
                begins[direction] = boundary;
                ends[direction] = boundary + 1;
                foreach_omp
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (boundary >= 3) {
                            if (boundary + 3 < n)
                                quantize_func(
                                    idx, d-data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),level);
                            else if (boundary + 1 < n)
                                quantize_func(idx, d-data, *d,
                                              interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),level);
                            else
                                quantize_func(idx, d-data, *d, interp_linear1(*(d - stride3x), *(d - stride)),level);
                        } else {
                            if (boundary + 3 < n)
                                quantize_func(idx, d-data, *d,
                                              interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),level);
                            else if (boundary + 1 < n)
                                quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level);
                            else
                                quantize_func(idx, d-data, *d, *(d - stride),level);
                        }
                    });
            }
        }
        return predict_error;
    }

    template <class QuantizeFunc>
    double interpolation_1d_simd_3d_y(T *data, const std::array<size_t, N> &begin_idx,
                                              const std::array<size_t, N> &end_idx, const size_t &direction,
                                              std::array<size_t, N> &strides, const size_t &math_stride,
                                              const std::string &interp_func, QuantizeFunc &&quantize_func) {
        assert(direction==1  && N==3);
        std::array<size_t,N> idx;
        int level = __builtin_ctzll(math_stride);
        for (size_t i = 0; i < N; ++i) {
            if (end_idx[i] < begin_idx[i]) return 0;
        }
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0.0;
        size_t offset = 0;
        size_t stride = math_stride * original_dim_offsets[direction];
        std::array<size_t, N> begins, ends, dim_offsets;
        for (size_t i = 0; i < N; ++i) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            dim_offsets[i] = original_dim_offsets[i];
            offset += original_dim_offsets[i] * begin_idx[i];
        }
        dim_offsets[direction] = stride;
        size_t stride2x = 2 * stride;
        if (interp_func == "linear") {
            begins[direction] = 1;
            ends[direction] = n - 1;
            strides[direction] = 2;
            foreach_omp
                <T, N>(data, offset, begins, ends, strides, dim_offsets,
                       [&](T *d, const std::array<size_t,N> &idx) { quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level); });
            if (n % 2 == 0) {
                begins[direction] = n - 1;
                ends[direction] = n;
                foreach_omp
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (n < 3)
                            quantize_func(idx, d-data, *d, *(d - stride),level);
                        else
                            quantize_func(idx, d-data, *d, interp_linear1(*(d - stride2x), *(d - stride)),level);
                    });
            }
        } else {
            size_t stride3x = 3 * stride;
            size_t i_start = 3;
            begins[direction] = i_start;
            ends[direction] = (n >= 3) ? (n - 3) : 0;
            strides[direction] = 2;
            size_t vector_len = ends[2] > begins[2] ? (ends[2]-begins[2]-1)/strides[2] + 1 : 0;

            #pragma omp parallel for
            for (size_t i = begins[0]; i < ends[0]; i += strides[0]) {
                auto tid = omp_get_thread_num();
                auto buffer_offset = buffer_len * tid;
                auto cur_buffer_1 = interp_buffer_1 + buffer_offset;
                auto cur_buffer_2 = interp_buffer_2 + buffer_offset;
                auto cur_buffer_3 = interp_buffer_3 + buffer_offset;
                auto cur_buffer_4 = interp_buffer_4 + buffer_offset; 
                auto cur_pred_buffer = pred_buffer + buffer_offset;
                

               
                
                for(size_t j = begins[1]; j < ends[1]; j += strides[1]){

                    auto cur_ij_offset = offset + i * dim_offsets[0] + j * dim_offsets[1];
                    size_t buffer_idx = 0;
                    if( j == begins[1]){
                        for (size_t k = begins[2]; k < ends[2]; k += strides[2]) {
                            auto cur_offset =  cur_ij_offset + k;
                            //if (visited[cur_offset- stride3x]==0 or visited[cur_offset+ stride3x]==0)
                            //   std::cout<<"e1 "<<i<<" "<<j<<" "<<k<<" "<<stride3x<<std::endl;
                           cur_buffer_1[buffer_idx] = data[cur_offset -  stride3x];
                           //cur_buffer_1[buffer_idx] = data[0];
                            cur_buffer_2[buffer_idx] = data[cur_offset - stride];
                           // cur_buffer_2[buffer_idx] = data[0];
                          cur_buffer_3[buffer_idx] = data[cur_offset + stride];
                           //cur_buffer_3[buffer_idx] = data[0];
                            cur_buffer_4[buffer_idx] = data[cur_offset +  stride3x];
                          //  cur_buffer_4[buffer_idx] = data[0];
                            buffer_idx++;

                        }
                    }
                    
                    else{
                        auto temp_buffer = cur_buffer_1;
                        cur_buffer_1 = cur_buffer_2;
                        cur_buffer_2 = cur_buffer_3;
                        cur_buffer_3 = cur_buffer_4;
                        cur_buffer_4 = temp_buffer;

                        buffer_idx = 0;
                        for (size_t k = begins[2]; k < ends[2]; k += strides[2]) {
                            auto cur_offset =  cur_ij_offset + stride3x + k;
                           // if (visited[cur_offset]==0 )
                            //   std::cout<<"e2 "<<i<<" "<<j<<" "<<k<<" "<<cur_offset<<std::endl;
                           
                           cur_buffer_4[buffer_idx++] = data[cur_offset];

                            //cur_buffer_4[buffer_idx++] = data[0];

                        }
                    }
                    
                    avx_interp_cubic(cur_buffer_1,cur_buffer_2,cur_buffer_3,cur_buffer_4, cur_pred_buffer, vector_len);
                    buffer_idx = 0;
                    for (size_t k = begins[2]; k < ends[2]; k += strides[2]){
                        auto pred = cur_pred_buffer[buffer_idx++];
                        auto d = data + cur_ij_offset + k;
                      // if (d-data < 0 || d-data>=num_elements)
                      //      std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                        quantize_func(idx, d-data, *d, pred,level);

                    }
                    
                }
            }
            std::vector<size_t> boundaries;
            boundaries.push_back(1);
            if (n % 2 == 1 && n > 3) {
                boundaries.push_back(n - 2);
            }
            if (n % 2 == 0 && n > 4) {
                boundaries.push_back(n - 3);
            }
            if (n % 2 == 0 && n > 2) {
                boundaries.push_back(n - 1);
            }
            for (auto boundary : boundaries) {
                begins[direction] = boundary;
                ends[direction] = boundary + 1;
                foreach_omp
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (boundary >= 3) {
                            if (boundary + 3 < n)
                                quantize_func(
                                    idx, d-data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),level);
                            else if (boundary + 1 < n)
                                quantize_func(idx, d-data, *d,
                                              interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),level);
                            else
                                quantize_func(idx, d-data, *d, interp_linear1(*(d - stride3x), *(d - stride)),level);
                        } else {
                            if (boundary + 3 < n)
                                quantize_func(idx, d-data, *d,
                                              interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),level);
                            else if (boundary + 1 < n)
                                quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level);
                            else
                                quantize_func(idx, d-data, *d, *(d - stride),level);
                        }
                    });
            }
        }
        return predict_error;
    }


    template <class QuantizeFunc>
    double interpolation_1d_simd_3d_z(T *data, const std::array<size_t, N> &begin_idx,
                                              const std::array<size_t, N> &end_idx, const size_t &direction,
                                              std::array<size_t, N> &strides, const size_t &math_stride,
                                              const std::string &interp_func, QuantizeFunc &&quantize_func) {
        assert(direction==2 && N==3);
        std::array<size_t,N> idx;
        int level = __builtin_ctzll(math_stride);
        for (size_t i = 0; i < N; ++i) {
            if (end_idx[i] < begin_idx[i]) return 0;
        }
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0.0;
        size_t offset = 0;
        size_t stride = math_stride * original_dim_offsets[direction];
        std::array<size_t, N> begins, ends, dim_offsets;
        for (size_t i = 0; i < N; ++i) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            dim_offsets[i] = original_dim_offsets[i];
            offset += original_dim_offsets[i] * begin_idx[i];
        }
        dim_offsets[direction] = stride;
        size_t stride2x = 2 * stride;
        if (interp_func == "linear") {
            begins[direction] = 1;
            ends[direction] = n - 1;
            strides[direction] = 2;
            foreach_omp
                <T, N>(data, offset, begins, ends, strides, dim_offsets,
                       [&](T *d, const std::array<size_t,N> &idx) { quantize_func(idx, d-data, *d, interp_linear(*(d - stride), *(d + stride)),level); });
            if (n % 2 == 0) {
                begins[direction] = n - 1;
                ends[direction] = n;
                foreach_omp
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d, const std::array<size_t,N> &idx) {
                        if (n < 3)
                            quantize_func(idx, d-data, *d, *(d - stride),level);
                        else
                            quantize_func(idx, d-data, *d, interp_linear1(*(d - stride2x), *(d - stride)),level);
                    });
            }
        } else {
            //size_t stride3x = 3 * stride;
            //size_t i_start = 3;
            begins[direction] = 1;
            ends[direction] = n;
            strides[direction] = 2;

            #pragma omp parallel for
            for (size_t i = begins[0]; i < ends[0]; i += strides[0]) {
                
                for(size_t j = begins[1]; j < ends[1]; j += strides[1]){
                    auto tid = omp_get_thread_num();
                    auto buffer_offset = buffer_len * tid;
                    auto cur_buffer = interp_buffer_1 + buffer_offset;
                    auto cur_pred_buffer = pred_buffer + buffer_offset;

                    auto cur_ij_offset = offset + i * dim_offsets[0] + j * dim_offsets[1];
                    size_t odd_len = n/2;//, even_len = n - odd_len;
                        
                    for (size_t k = 0; k < n; k += 2) {
                        auto cur_offset = cur_ij_offset + k * dim_offsets[2];
                        cur_buffer[k/2] = data[cur_offset];
                    }
                    
                    avx_interp_cubic_1D(cur_buffer,cur_pred_buffer, n);
                    for (size_t k = 0; k < odd_len; ++k ){
                        auto pred = cur_pred_buffer[k];
                        auto d = data + cur_ij_offset + (2 * k + 1) * dim_offsets[2];
                      // if (d-data < 0 || d-data>=num_elements)
                      //      std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                        quantize_func(idx, d-data, *d,pred,level);

                    }
                    
                }
            }
        }
        return predict_error;
    }




    template <class QuantizeFunc>
    double interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end,
                         const std::string &interp_func, QuantizeFunc &&quantize_func, const int direction,
                         size_t stride = 1) {
        if constexpr (N == 1) {  // old API
            return interpolation_1d(data, begin[0], end[0], stride, interp_func, quantize_func);
        } else if constexpr (N == 2) {  // old API
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dim_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                size_t begin_offset =
                    begin[dims[0]] * original_dim_offsets[dims[0]] + j * original_dim_offsets[dims[1]];
                predict_error += interpolation_1d(
                    data, begin_offset, begin_offset + (end[dims[0]] - begin[dims[0]]) * original_dim_offsets[dims[0]],
                    stride * original_dim_offsets[dims[0]], interp_func, quantize_func);
            }
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                size_t begin_offset =
                    i * original_dim_offsets[dims[0]] + begin[dims[1]] * original_dim_offsets[dims[1]];
                predict_error += interpolation_1d(
                    data, begin_offset, begin_offset + (end[dims[1]] - begin[dims[1]]) * original_dim_offsets[dims[1]],
                    stride * original_dim_offsets[dims[1]], interp_func, quantize_func);
            }
            return predict_error;
        } else if constexpr (N == 3 || N == 4) {  // new API (for faster speed)
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dim_sequences[direction];
            std::array<size_t, N> strides;
            std::array<size_t, N> begin_idx = begin, end_idx = end;
            strides[dims[0]] = 1;
            size_t max_interp_seq_length = 0;
             for (uint i = 0; i < N; ++i) 
                max_interp_seq_length = std::max(max_interp_seq_length, (end[i]-begin[i])/stride );
            for (uint i = 1; i < N; ++i) {
                begin_idx[dims[i]] = (begin[dims[i]] ? begin[dims[i]] + stride2x : 0);
                strides[dims[i]] = stride2x;
            }
            if(0){//if(N==3  && max_interp_seq_length >= 2 * AVX_256_parallelism * nThreads){//avx
                if(direction ==0 ){//xyz
                    predict_error += interpolation_1d_simd_3d_x(data, begin_idx, end_idx, dims[0], strides, stride, interp_func, quantize_func);
                    //predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[0], strides, stride, interp_func, quantize_func);
                    begin_idx[1] = begin[1];
                    begin_idx[0] = (begin[0] ? begin[0] + stride : 0);
                    strides[0] = stride;
                    predict_error += interpolation_1d_simd_3d_y(data, begin_idx, end_idx, dims[1], strides, stride, interp_func, quantize_func);
                   // predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[1], strides, stride, interp_func, quantize_func);
                    begin_idx[2] = begin[2];
                    begin_idx[1] = (begin[1] ? begin[1] + stride : 0);
                    strides[1] = stride;
                    //predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[2], strides, stride, interp_func, quantize_func);
                    predict_error += interpolation_1d_simd_3d_z(data, begin_idx, end_idx, dims[2], strides, stride, interp_func, quantize_func);
                }
                else{//zyx
                    predict_error += interpolation_1d_simd_3d_z(data, begin_idx, end_idx, dims[0], strides, stride, interp_func, quantize_func);
                    //predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[0], strides, stride, interp_func, quantize_func);
                    begin_idx[1] = begin[1];
                    begin_idx[2] = (begin[2] ? begin[2] + stride : 0);
                    strides[2] = stride;
                    predict_error += interpolation_1d_simd_3d_y(data, begin_idx, end_idx, dims[1], strides, stride, interp_func, quantize_func);
                   // predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[1], strides, stride, interp_func, quantize_func);
                    begin_idx[0] = begin[0];
                    begin_idx[1] = (begin[1] ? begin[1] + stride : 0);
                    strides[1] = stride;
                    predict_error += interpolation_1d_simd_3d_x(data, begin_idx, end_idx, dims[2], strides, stride, interp_func, quantize_func);
                    //predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[2], strides, stride, interp_func, quantize_func);
                }
            }

            else{
                predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[0], strides, stride, interp_func, quantize_func);
                for (uint i = 1; i < N; ++i) {
                begin_idx[dims[i]] = begin[dims[i]];
                begin_idx[dims[i - 1]] = (begin[dims[i - 1]] ? begin[dims[i - 1]] + stride : 0);
                strides[dims[i - 1]] = stride;
               
                predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[i], strides, stride, interp_func, quantize_func);
                }
            }

            
            return predict_error;
        } else {
            throw std::runtime_error("Unsupported dimension in InterpolationDecomposition");
        }
    }

    int interp_level = -1;
    int interp_id;
    uint blocksize;
    std::vector<std::string> interpolators = {"linear", "cubic"};
    int *quant_inds;
    size_t quant_index = 0;
    double max_error;
    QuantizerOMP quantizer;
    size_t num_elements;
    std::array<size_t, N> original_dimensions;
    std::array<size_t, N> original_dim_offsets;
    std::vector<std::array<int, N>> dim_sequences;
    int direction_sequence_id;
    size_t anchor_stride = 0;
    double eb_alpha = -1;
    double eb_beta = -1;
    double eb_ratio = 0.5;  // To be deprecated
    const size_t AVX_256_parallelism = 32 / sizeof(T);
    size_t max_dim = 1;
    size_t nThreads = 1;
    size_t buffer_len = 1024;
    T *interp_buffer_1,*interp_buffer_2,*interp_buffer_3,*interp_buffer_4,*pred_buffer;
    //std::vector<int> visited;

    std::vector<size_t> level_prefix;
    std::vector<std::array<size_t,N> >reduced_dim_offsets;

    
};

template <class T, uint N, class QuantizerOMP>
InterpolationDecomposition_OMP<T, N, QuantizerOMP> make_decomposition_interpolation_omp(const Config &conf, QuantizerOMP quantizer) {
    return InterpolationDecomposition_OMP<T, N, QuantizerOMP>(conf, quantizer);
}

}  // namespace SZ3

#endif
#endif