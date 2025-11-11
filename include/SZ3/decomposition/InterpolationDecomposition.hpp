#ifndef SZ3_INTERPOLATION_DECOMPOSITION_HPP
#define SZ3_INTERPOLATION_DECOMPOSITION_HPP

#include <cmath>
#include <cstring>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {
template <class T, uint N, class Quantizer>
class InterpolationDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    InterpolationDecomposition(const Config &conf, Quantizer quantizer) : quantizer(quantizer) {
        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quantizer interface");
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        init();

        this->quant_inds = quant_inds.data();
        double eb = quantizer.get_eb();

        if (anchor_stride == 0) {                                               // check whether used anchor points
            *dec_data = quantizer.recover(0, this->quant_inds[quant_index++]);  // no anchor points
        } else {
            recover_anchor_grid(dec_data);  // recover anchor points
            interp_level--;
        }

        for (int level = interp_level; level > 0 && level <= interp_level; level--) {
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
            bool use_gather_scatter = (stride == 1) && (N <= 3) && (max_dim >=35); //todo: try different conditions
            if(use_gather_scatter){
                gather(dec_data, stride);
                interpolation_gathered(
                    dec_data, interpolators[interp_id],
                    [&](size_t idx, T &d, T pred) { d = quantizer.recover(pred, quant_inds[quant_index++]); },
                    direction_sequence_id, stride);
                scatter(dec_data, stride);
                
            }
            else{
                auto interp_block_size = blocksize * stride;
                auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                    dec_data, std::begin(original_dimensions), std::end(original_dimensions), interp_block_size, 0);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto end_idx = block.get_global_index();
                    for (uint i = 0; i < N; i++) {
                        end_idx[i] += interp_block_size;
                        if (end_idx[i] > original_dimensions[i] - 1) {
                            end_idx[i] = original_dimensions[i] - 1;
                        }
                    }
                    interpolation(
                        dec_data, block.get_global_index(), end_idx, interpolators[interp_id],
                        [&](size_t idx, T &d, T pred) { d = quantizer.recover(pred, quant_inds[quant_index++]); },
                        direction_sequence_id, stride);
                }
            }
            
                
        }
        quantizer.postdecompress_data();
        delete []buffer;
        delete []aligned_buffer;
        return dec_data;
    }

    // compress given the error bound
    std::vector<int> compress(const Config &conf, T *data) override {
        std::copy_n(conf.dims.begin(), N, original_dimensions.begin());

        interp_id = conf.interpAlgo;
        direction_sequence_id = conf.interpDirection;
        anchor_stride = conf.interpAnchorStride;
        blocksize = 32;  // a empirical value. Can be very large but not helpful
        eb_alpha = conf.interpAlpha;
        eb_beta = conf.interpBeta;

        init();


        std::vector<int> quant_inds_vec(num_elements);
        quant_inds = quant_inds_vec.data();
        double eb = quantizer.get_eb();
        if (anchor_stride == 0) {  // check whether to use anchor points
            quant_inds[quant_index++] = quantizer.quantize_and_overwrite(*data, 0);  // no
        } else {
            build_anchor_grid(data);  // losslessly saving anchor points
            interp_level--;
        }

        for (int level = interp_level; level > 0 && level <= interp_level; level--) {

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

            bool use_gather_scatter = (stride == 1) && (N <= 3) && (max_dim >=35); //todo: try different conditions
            
            if(use_gather_scatter){//cannot use blocked interp here
                Timer timer(true);
                gather(data, stride);
                timer.stop("Gather");

                timer.start();
                T * test = new T[conf.num];
                if(N==3){
                    
                    auto even_len_x =  (original_dimensions[0] - 1)/2 + 1;
                    auto even_len_y =  (original_dimensions[1] - 1)/2 + 1;
                    auto even_len_z =  (original_dimensions[2] - 1)/2 + 1;
                    for (size_t i =0 ;i < original_dimensions[0];i++){
                        auto ii = i % 2 == 0 ? i /2 : even_len_x + i/2;
                        for (size_t j =0 ;j < original_dimensions[1];j++){
                            auto jj = j % 2 == 0 ? j /2 : even_len_y + j/2;
                            for (size_t k =0 ;k < original_dimensions[2];k++){
                                
                                auto kk = k % 2 == 0 ? k /2 : even_len_z + k/2;
                                test[ii*original_dim_offsets[0]+jj*original_dim_offsets[1]+kk] = data[i*original_dim_offsets[0]+j*original_dim_offsets[1]+k];

                            }
                        }
                    }
                }
                timer.stop("Gather2");
                
                interpolation_gathered(
                        tedt, interpolators[interp_id],
                        [&](size_t idx, T &d, T pred) {
                            quant_inds[quant_index++] = (quantizer.quantize_and_overwrite(d, pred));
                        },
                        direction_sequence_id, stride);
                timer.start();
                scatter(data, stride);
                timer.stop("Scatter");
                timer.start();
                if(N==3){
                    
                    auto even_len_x =  (original_dimensions[0] - 1)/2 + 1;
                    auto even_len_y =  (original_dimensions[1] - 1)/2 + 1;
                    auto even_len_z =  (original_dimensions[2] - 1)/2 + 1;
                    for (size_t i =0 ;i < original_dimensions[0];i++){
                        auto ii = i >= even_len_x ? (i-even_len_x) *2 + 1:  i * 2;
                        for (size_t j =0 ;j < original_dimensions[1];j++){
                            auto jj =  j >= even_len_y ? (j-even_len_y) *2 + 1:  j * 2;
                            for (size_t k =0 ;k < original_dimensions[2];k++){
                                
                                auto kk =  k >= even_len_z ? (k-even_len_z) *2 + 1:  k * 2;
                                data[ii*original_dim_offsets[0]+jj*original_dim_offsets[1]+kk] = test[i*original_dim_offsets[0]+j*original_dim_offsets[1]+k];

                            }
                        }
                    }
                }
                 delete []test;
                 timer.stop("Scatter2");


            }
            else{
                auto interp_block_size = blocksize * stride;

                auto inter_block_range = std::make_shared<multi_dimensional_range<T, N>>(
                    data, std::begin(original_dimensions), std::end(original_dimensions), interp_block_size, 0);

                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();

                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto end_idx = block.get_global_index();
                    for (uint i = 0; i < N; i++) {
                        end_idx[i] += interp_block_size;
                        if (end_idx[i] > original_dimensions[i] - 1) {
                            end_idx[i] = original_dimensions[i] - 1;
                        }
                    }

                    interpolation(
                        data, block.get_global_index(), end_idx, interpolators[interp_id],
                        [&](size_t idx, T &d, T pred) {
                            quant_inds[quant_index++] = (quantizer.quantize_and_overwrite(d, pred));
                        },
                        direction_sequence_id, stride);
                }
            }

        }
        quantizer.set_eb(eb);
        quantizer.postcompress_data();
        delete []buffer;
        delete []aligned_buffer;
        std::cout<<quant_inds_vec.size()<<std::endl;
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
        for (uint i = 0; i < N; i++) {
            max_dim = std::max(max_dim, original_dimensions[i]);
            if (interp_level < ceil(log2(original_dimensions[i]))) {
                interp_level = static_cast<int>(ceil(log2(original_dimensions[i])));
            }
    	    if (original_dimensions[i] > anchor_stride)
    	        use_anchor = true;
                num_elements *= original_dimensions[i];
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
        for (int i = N - 2; i >= 0; i--) {
            original_dim_offsets[i] = original_dim_offsets[i + 1] * original_dimensions[i + 1];
        }

        dim_sequences = std::vector<std::array<int, N>>();
        auto sequence = std::array<int, N>();
        for (uint i = 0; i < N; i++) {
            sequence[i] = i;
        }
        do {
            dim_sequences.push_back(sequence);
        } while (std::next_permutation(sequence.begin(), sequence.end()));

   
        buffer = new T [max_dim * column_num];
        size_t alignment = 32;  // 256 bits
        size_t alloc_chunks = (max_dim * sizeof(T) + 31) / alignment;
        auto aligned_buf_bytes = alignment * alloc_chunks;

        aligned_buffer = new T [aligned_buf_bytes / sizeof (T)];
    }

    void build_anchor_grid(T *data) {  // store anchor points. steplength: anchor_stride on each dimension
        std::array<size_t, N> strides;
        std::array<size_t, N> begins{0};
        std::fill(strides.begin(), strides.end(), anchor_stride);
        foreach
            <T, N>(data, 0, begins, original_dimensions, strides, original_dim_offsets,
                   [&](T *d) { quant_inds[quant_index++] = quantizer.force_save_unpred(*d); });
    }

    void recover_anchor_grid(T *data) {  // recover anchor points. steplength: anchor_stride on each dimension
        std::array<size_t, N> strides;
        std::array<size_t, N> begins{0};
        std::fill(strides.begin(), strides.end(), anchor_stride);
        foreach
            <T, N>(data, 0, begins, original_dimensions, strides, original_dim_offsets, [&](T *d) {
                *d = quantizer.recover_unpred();
                quant_index++;
            });
    }

    ALWAYS_INLINE void gather_base(const T* src, size_t len, T* dst) const
    {
        /*
    #ifdef __AVX2__
      const double* src_end = src + len;
      double* dst_evens = dst;
      double* dst_odds = dst + len - len / 2;

      // Process 8 elements at a time
      for (; src + 8 <= src_end; src += 8) {
        __m256d v0 = _mm256_loadu_pd(src);      // 0, 1, 2, 3
        __m256d v1 = _mm256_loadu_pd(src + 4);  // 4, 5, 6, 7

        __m256d evens = _mm256_unpacklo_pd(v0, v1);  // 0, 4, 2, 6
        __m256d odds = _mm256_unpackhi_pd(v0, v1);   // 1, 5, 3, 7

        __m256d result1 = _mm256_permute4x64_pd(evens, 0b11011000);  // 0, 2, 4, 6
        __m256d result2 = _mm256_permute4x64_pd(odds, 0b11011000);   // 1, 3, 5, 7

        _mm256_store_pd(dst_evens, result1);
        _mm256_storeu_pd(dst_odds, result2);

        dst_evens += 4;
        dst_odds += 4;
      }

      for (; src < src_end - 1; src += 2) {
        *(dst_evens++) = *src;
        *(dst_odds++) = *(src + 1);
      }

      if (src < src_end)
        *dst_evens = *src;
    #else*/
      size_t low_count = len - len / 2, high_count = len / 2;
      for (size_t i = 0; i < low_count; i++) {
        *dst = *(src + i * 2);
        ++dst;
      }
      for (size_t i = 0; i < high_count; i++) {
        *dst = *(src + i * 2 + 1);
        ++dst;
      }
    //#endif
    }

    ALWAYS_INLINE void scatter_base(const T* begin, size_t len, T* dst) const
    {
        /*
    #ifdef __AVX2__
      const double* even_end = begin + len - len / 2;
      const double* odd_beg = even_end;
      const double* dst_end = dst + len;

      // Process 8 elements at a time
      for (; begin + 4 < even_end; begin += 4) {
        __m256d v0 = _mm256_loadu_pd(begin);    // 0, 1, 2, 3
        __m256d v1 = _mm256_loadu_pd(odd_beg);  // 4, 5, 6, 7

        __m256d evens = _mm256_unpacklo_pd(v0, v1);  // 0, 4, 2, 6
        __m256d odds = _mm256_unpackhi_pd(v0, v1);   // 1, 5, 3, 7

        __m256d result1 = _mm256_permute2f128_pd(evens, odds, 0x20);  // 0, 4, 1, 5
        __m256d result2 = _mm256_permute2f128_pd(evens, odds, 0x31);  // 2, 6, 3, 7

        _mm256_store_pd(dst, result1);
        _mm256_store_pd(dst + 4, result2);

        dst += 8;
        odd_beg += 4;
      }

      for (; dst < dst_end - 1; dst += 2) {
        *dst = *(begin++);
        *(dst + 1) = *(odd_beg++);
      }

      if (dst < dst_end)
        *dst = *begin;
    #else*/
      size_t low_count = len - len / 2, high_count = len / 2;
      for (size_t i = 0; i < low_count; i++) {
        *(dst + i * 2) = *begin;
        ++begin;
      }
      for (size_t i = 0; i < high_count; i++) {
        *(dst + i * 2 + 1) = *begin;
        ++begin;
      }
    //#endif
    }

     //Gather the elements on each dimension that have coordinations divisible by 2*stride, using buffer
    ALWAYS_INLINE void gather_1D(T * data, size_t stride){
        //dim0. 

       
        if(stride == 1){
            gather_base(data, original_dimensions[N-1], aligned_buffer);
            std::copy(aligned_buffer, aligned_buffer + original_dimensions[N-1], data);
        }
        else{
            size_t buffer_idx = 0;
            for(size_t i = 0;i < original_dimensions[N - 1];i += stride)
                buffer[buffer_idx++] = * (data + i); 
            auto col_len = (original_dimensions[N - 1] - 1) / stride + 1;
            gather_base(buffer, col_len, aligned_buffer);
            buffer_idx = 0;
            for(size_t i = 0;i < original_dimensions[N - 1];i+= stride) 
                * (data + i) = aligned_buffer[buffer_idx++];
        }


    }


    //Scatter the elements on each dimension that have coordinations divisible by 2*stride, using buffer
    ALWAYS_INLINE void scatter_1D(T * data, size_t stride){
        ///dim1. 
       if(stride == 1){
                scatter_base(data, original_dimensions[N-1], aligned_buffer);
                std::copy(aligned_buffer, aligned_buffer + original_dimensions[N-1], data);
        }
        else{
            size_t buffer_idx = 0;
            for(size_t i = 0;i < original_dimensions[N - 1];i += stride)
                buffer[buffer_idx++] = * (data + i); 
            auto col_len = (original_dimensions[N - 1] - 1) / stride + 1;
            scatter_base(buffer, col_len, aligned_buffer);
            buffer_idx = 0;
            for(size_t i = 0;i < original_dimensions[N - 1];i+= stride) 
                * (data + i) = aligned_buffer[buffer_idx++];
        }


    }



    //Gather the elements on each dimension that have coordinations divisible by 2*stride, using buffer
    ALWAYS_INLINE void gather_2D(T * data, size_t stride){
        //dim1. 
        for(size_t i = 0;i < original_dimensions[N - 2];i += stride){
            auto pos = data + i * original_dim_offsets[N - 2] ;
            //just 1 column
            size_t buffer_idx = 0;
            if(stride == 1){
                gather_base(pos, original_dimensions[N-1], aligned_buffer);
                std::copy(aligned_buffer, aligned_buffer + original_dimensions[N-1], pos);
            }
            else{
                for(size_t j = 0;j < original_dimensions[N - 1];j += stride)
                    buffer[buffer_idx++] = * (pos + j); 
                auto col_len = (original_dimensions[N - 1] - 1) / stride + 1;
                gather_base(buffer, col_len, aligned_buffer);
                buffer_idx = 0;
                for(size_t j = 0;j < original_dimensions[N - 1];j+= stride) 
                    * (pos + j) = aligned_buffer[buffer_idx++];
            }


        }


        //dim0. Currently, only 1 column per iter. 
        for(size_t j = 0;j < original_dimensions[N-1];j += stride){
            auto pos = data  +  j;
            size_t buffer_idx = 0;
            //just 1 column
            for(size_t i = 0;i < original_dimensions[N-2];i += stride)
                buffer[buffer_idx++] = * (pos + i * original_dim_offsets[N-2]); 
            auto col_len = (original_dimensions[N-2] - 1) / stride + 1;
            gather_base(buffer, col_len, aligned_buffer);
            buffer_idx = 0;
            for(size_t i = 0;i < original_dimensions[N-2];i += stride)
                * (pos + i * original_dim_offsets[N-2]) = aligned_buffer[buffer_idx++];
        }

    }


    //Scatter the elements on each dimension that have coordinations divisible by 2*stride, using buffer
    ALWAYS_INLINE void scatter_2D(T * data, size_t stride){
        ///dim1. 
        for(size_t i = 0;i < original_dimensions[N - 2];i += stride){
            auto pos = data + i * original_dim_offsets[N - 2] ;
            //just 1 column
            size_t buffer_idx = 0;
            if(stride == 1){
                scatter_base(pos, original_dimensions[N-1], aligned_buffer);
                std::copy(aligned_buffer, aligned_buffer + original_dimensions[N-1], pos);
            }
            else{
                for(size_t j = 0;j < original_dimensions[N - 1];j += stride)
                    buffer[buffer_idx++] = * (pos + j); 
                auto col_len = (original_dimensions[N - 1] - 1) / stride + 1;
                scatter_base(buffer, col_len, aligned_buffer);
                buffer_idx = 0;
                for(size_t j = 0;j < original_dimensions[N - 1];j+= stride) 
                    * (pos + j) = aligned_buffer[buffer_idx++];
            }


        }


        //dim0. Currently, only 1 column per iter. 
        for(size_t j = 0;j < original_dimensions[N-1];j += stride){
            auto pos = data +  j;
            size_t buffer_idx = 0;
            //just 1 column
            for(size_t i = 0;i < original_dimensions[N-2];i += stride)
                buffer[buffer_idx++] = * (pos + i * original_dim_offsets[N-2]); 
            auto col_len = (original_dimensions[N-2] - 1) / stride + 1;
            scatter_base(buffer, col_len, aligned_buffer);
            buffer_idx = 0;
            for(size_t i = 0;i < original_dimensions[N-2];i += stride)
                * (pos + i * original_dim_offsets[N-2]) = aligned_buffer[buffer_idx++];
        }


    }




    //Gather the elements on each dimension that have coordinations divisible by 2*stride, using buffer
    ALWAYS_INLINE void gather_3D(T * data, size_t stride){
        //dim21
        for(size_t i = 0;i < original_dimensions[0];i += stride){
            auto pos = data + i * original_dim_offsets[0];
            gather_2D(pos,stride);
        }

        //dim0. Grouping columns
        for(size_t j = 0;j < original_dimensions[1];j += stride){
            for(size_t k = 0;k < original_dimensions[2];k += column_num){
                auto pos = data + j * original_dim_offsets[1] + k;
                const auto col_count = std::min((column_num - 1) / stride + 1, (original_dimensions[2] - k - 1) / stride + 1);
                
                auto col_len = (original_dimensions[0] - 1) / stride + 1;
                for(size_t i = 0; i < original_dimensions[0]; i += stride){
                    size_t buffer_idx = 0;
                    for(size_t kk = 0;kk < col_count * stride;kk += stride){
                        buffer[(buffer_idx++) * col_len + i] = * (pos + i * original_dim_offsets[0] + kk); 
                    }
                }
                for(size_t buffer_idx = 0;buffer_idx < col_count ;buffer_idx++){
                    auto buffer_pos = buffer + buffer_idx * col_len;
                    gather_base(buffer_pos, col_len, aligned_buffer);
                    std::copy(aligned_buffer, aligned_buffer + col_len, buffer_pos);
                }
                for(size_t i = 0; i < original_dimensions[0]; i += stride){
                     size_t buffer_idx = 0;
                    for(size_t kk = 0;kk < col_count * stride;kk += stride){
                        * (pos + i * original_dim_offsets[0] + kk) = buffer[(buffer_idx++) * col_len + i];
                    }
                }

            }
        }

    }


    //Scatter the elements on each dimension that have coordinations divisible by 2*stride, using buffer
    ALWAYS_INLINE void scatter_3D (T * data, size_t stride){
        //we don't need to reverse the dimension order of the gather one, so just keep the same
        //dim21
        for(size_t i = 0;i < original_dimensions[0];i += stride){
            auto pos = data + i * original_dim_offsets[0];
            scatter_2D(pos,stride);
        }

        //dim0. Grouping columns
        for(size_t j = 0;j < original_dimensions[1];j += stride){
            for(size_t k = 0;k < original_dimensions[2];k += column_num){
                auto pos = data + j * original_dim_offsets[1] + k;
                const auto col_count = std::min((column_num - 1) / stride + 1, (original_dimensions[2] - k - 1) / stride + 1);
                
                auto col_len = (original_dimensions[0] - 1) / stride + 1;
                for(size_t i = 0; i < original_dimensions[0]; i += stride){
                    size_t buffer_idx = 0;
                    for(size_t kk = 0;kk < col_count * stride;kk += stride){
                        buffer[(buffer_idx++) * col_len + i] = * (pos + i * original_dim_offsets[0] + kk); 
                    }
                }
                for(size_t buffer_idx = 0;buffer_idx < col_count ;buffer_idx++){
                    auto buffer_pos = buffer + buffer_idx * col_len;
                    scatter_base(buffer_pos, col_len, aligned_buffer);
                    std::copy(aligned_buffer, aligned_buffer + col_len, buffer_pos);
                }
                for(size_t i = 0; i < original_dimensions[0]; i += stride){
                    size_t buffer_idx = 0;
                    for(size_t kk = 0;kk < col_count * stride;kk += stride){
                        * (pos + i * original_dim_offsets[0] + kk) = buffer[(buffer_idx++) * col_len + i];
                    }
                }

            }
        }

    }

    void gather(T * data, size_t stride){
        if constexpr (N == 3) {
            gather_3D (data,stride);
        }
        else if constexpr (N == 2) {
            gather_2D (data,stride);
        }
        else if constexpr (N == 1) {
            gather_1D(data,stride);
        }
    }

    void scatter(T * data, size_t stride){
        if constexpr (N == 3) {
            scatter_3D (data,stride);
        }
        else if constexpr (N == 2) {
            scatter_2D (data,stride);
        }
        else if constexpr (N == 1) {
            scatter_1D (data,stride);
        }
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

        size_t stride3x = 3 * stride;
        size_t stride5x = 5 * stride;
        if (interp_func == "linear" || n < 5) {
            // if (pb == PB_predict_overwrite) {
            for (size_t i = 1; i + 1 < n; i += 2) {
                T *d = data + begin + i * stride;
                quantize_func(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
            }
            if (n % 2 == 0) {
                T *d = data + begin + (n - 1) * stride;
                if (n < 4) {
                    quantize_func(d - data, *d, *(d - stride));
                } else {
                    quantize_func(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                }
            }
            // }
        } else {
            T *d;
            size_t i;
            for (i = 3; i + 3 < n; i += 2) {
                d = data + begin + i * stride;
                quantize_func(d - data, *d,
                              interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
            }
            d = data + begin + stride;
            quantize_func(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

            d = data + begin + i * stride;
            quantize_func(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
            if (n % 2 == 0) {
                d = data + begin + (n - 1) * stride;
                quantize_func(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
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
        for (size_t i = 0; i < N; i++) {
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
            foreach
                <T, N>(data, offset, begins, ends, strides, dim_offsets,
                       [&](T *d) { quantize_func(d - data, *d, interp_linear(*(d - stride), *(d + stride))); });
            if (n % 2 == 0) {
                begins[direction] = n - 1;
                ends[direction] = n;
                foreach
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d) {
                        if (n < 3)
                            quantize_func(d - data, *d, *(d - stride));
                        else
                            quantize_func(d - data, *d, interp_linear1(*(d - stride2x), *(d - stride)));
                    });
            }
        } else {
            size_t stride3x = 3 * stride;
            size_t i_start = 3;
            begins[direction] = i_start;
            ends[direction] = (n >= 3) ? (n - 3) : 0;
            strides[direction] = 2;
            foreach
                <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d) {
                    quantize_func(d - data, *d,
                                  interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
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
                begins[direction] = boundary;
                ends[direction] = boundary + 1;
                foreach
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d) {
                        if (boundary >= 3) {
                            if (boundary + 3 < n)
                                quantize_func(
                                    d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            else if (boundary + 1 < n)
                                quantize_func(d - data, *d,
                                              interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                            else
                                quantize_func(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                        } else {
                            if (boundary + 3 < n)
                                quantize_func(d - data, *d,
                                              interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                            else if (boundary + 1 < n)
                                quantize_func(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                            else
                                quantize_func(d - data, *d, *(d - stride));
                        }
                    });
            }
        }
        return predict_error;
    }
    /*
    template <class QuantizeFunc>
    double interpolation_1d_xpack_3d(T *data, const std::array<size_t, N> &begin_idx,
                                              const std::array<size_t, N> &end_idx, const size_t &direction,
                                              std::array<size_t, N> &strides, const size_t &math_stride,
                                              const std::string &interp_func, QuantizeFunc &&quantize_func) {
        assert(N == 3 and direction == 0);
        for (size_t i = 0; i < N; i++) {
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
        for (size_t i = 0; i < N; i++) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            dim_offsets[i] = original_dim_offsets[i];
            offset += original_dim_offsets[i] * begin_idx[i];
        }
        dim_offsets[direction] = stride;
        size_t stride2x = 2 * stride;


        size_t pack_cluster_num = direction == 0 ? std::max( 1, 64 / sizeof(T) / strides[2] ): 1;

        if (interp_func == "linear") {
            begins[direction] = 1;
            ends[direction] = n - 1;
            strides[direction] = 2;
            foreach
                <T, N>(data, offset, begins, ends, strides, dim_offsets,
                       [&](T *d) { quantize_func(d - data, *d, interp_linear(*(d - stride), *(d + stride))); });
            if (n % 2 == 0) {
                begins[direction] = n - 1;
                ends[direction] = n;
                foreach
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d) {
                        if (n < 3)
                            quantize_func(d - data, *d, *(d - stride));
                        else
                            quantize_func(d - data, *d, interp_linear1(*(d - stride2x), *(d - stride)));
                    });
            }
        } else {
            size_t stride3x = 3 * stride;
            size_t i_start = 3;
            begins[direction] = i_start;
            ends[direction] = (n >= 3) ? (n - 3) : 0;
            strides[direction] = 2;
            foreach
                <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d) {
                    quantize_func(d - data, *d,
                                  interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
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
                begins[direction] = boundary;
                ends[direction] = boundary + 1;
                foreach
                    <T, N>(data, offset, begins, ends, strides, dim_offsets, [&](T *d) {
                        if (boundary >= 3) {
                            if (boundary + 3 < n)
                                quantize_func(
                                    d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            else if (boundary + 1 < n)
                                quantize_func(d - data, *d,
                                              interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                            else
                                quantize_func(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                        } else {
                            if (boundary + 3 < n)
                                quantize_func(d - data, *d,
                                              interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                            else if (boundary + 1 < n)
                                quantize_func(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                            else
                                quantize_func(d - data, *d, *(d - stride));
                        }
                    });
            }
        }
        return predict_error;
    }*/

    template <class QuantizeFunc>
    double interpolation_gathered_base(T * buffer,const size_t len, const std::string &interp_func,
                            QuantizeFunc &&quantize_func) {


        double predict_error = 0;
       
        size_t even_len = len - len / 2;
        //size_t odd_len = len / 2;
        if (interp_func == "linear" || len < 5) {
            // if (pb == PB_predict_overwrite) {
            auto d = buffer + even_len, pred_d = buffer;
            for (size_t i = 0; i  < even_len - 1; i ++) {
                quantize_func(d - buffer, *d, interp_linear(*pred_d, *(pred_d + 1)));
                d += 1;
                pred_d += 1;
            }
            if (len % 2 == 0) {
                if (len < 4) {
                    quantize_func(d - buffer, *d, *pred_d);
                } else {
                    quantize_func(d - buffer, *d, interp_linear1(*pred_d, *(pred_d - 1)));
                }
            }
        }
        else {
            auto d = buffer + (even_len+1) , pred_d = buffer + 1;
            for (size_t i = 1; i + 2 < even_len; i ++) {
                
                quantize_func(d - buffer, *d,
                              interp_cubic(*(pred_d- 1), *pred_d, *(pred_d + 1), *(pred_d + 2)));
            }
            d = buffer + even_len;
            pred_d = buffer;
            quantize_func(d - buffer, *d, interp_quad_1(*pred_d, *(pred_d + 1), *(pred_d + 2)));

            d = buffer + (even_len * 2 - 1) ;
            pred_d = buffer + (even_len - 1);
            quantize_func(d - buffer, *d, interp_quad_2(*(pred_d - 2), *(pred_d - 1), *pred_d));
            if (len % 2 == 0) {
                d += 1;
                quantize_func(d - buffer, *d, interp_quad_3(*(pred_d - 2), *(pred_d - 1), *pred_d));
            }
        }
        return predict_error;
       
        
    }


    template <class QuantizeFunc>
    double interpolation_gathered_1D(T * data,const size_t stride, const std::string &interp_func,
                            QuantizeFunc &&quantize_func) {

        double predict_error = 0;
        auto stride2x = stride * 2;
        size_t len = (original_dimensions[N -1] - 1) / stride + 1;
        if (stride2x < column_num || len < 16){//no buffer needed
           
            size_t even_len = len - len / 2;
            //size_t odd_len = len / 2;
            if (interp_func == "linear" || len < 5) {
                // if (pb == PB_predict_overwrite) {
                auto d = data + even_len * stride, pred_d = data;
                for (size_t i = 0; i  < even_len - 1; i ++) {
                    quantize_func(d - data, *d, interp_linear(*pred_d, *(pred_d + stride)));
                    d += stride;
                    pred_d += stride;
                }
                if (len % 2 == 0) {
                    if (len < 4) {
                        quantize_func(d - data, *d, *pred_d);
                    } else {
                        quantize_func(d - data, *d, interp_linear1(*pred_d, *(pred_d - stride)));
                    }
                }
            }
            else {
                auto d = data + (even_len+1) * stride, pred_d = data + stride;
                for (size_t i = 1; i + 2 < even_len; i ++) {
                    
                    quantize_func(d - data, *d,
                                  interp_cubic(*(pred_d- stride), *pred_d, *(pred_d + stride), *(pred_d + stride2x)));
                }
                d = data + even_len * stride;
                pred_d = data;
                quantize_func(d - data, *d, interp_quad_1(*pred_d, *(pred_d + stride), *(pred_d + stride2x)));

                d = data + (even_len * 2 - 1) * stride;
                pred_d = data + (even_len - 1) * stride;
                quantize_func(d - data, *d, interp_quad_2(*(pred_d - stride2x), *(pred_d - stride), *pred_d));
                if (len % 2 == 0) {
                    d += stride;
                    quantize_func(d - data, *d, interp_quad_3(*(pred_d - stride2x), *(pred_d - stride), *pred_d));
                }
            }

        }
        else{//pack
            size_t buffer_idx=0;
            for(size_t i = 0;i<original_dimensions[N-1];i+=stride){
                aligned_buffer[buffer_idx++] = *(data + i);
            }
            interpolation_gathered_base(aligned_buffer, len, interp_func, quantize_func);
            buffer_idx=0;
            for(size_t i = 0;i<original_dimensions[N-1];i+=stride){
                *(data + i) = aligned_buffer[buffer_idx++] ;
            }

        }
        return predict_error;
    }


    template <class QuantizeFunc>
    double interpolation_gathered_2D(T * data,const size_t stride, const std::string &interp_func,
                            QuantizeFunc &&quantize_func, const int direction) {

        double predict_error = 0;
        size_t len_x = (original_dimensions[N - 2] - 1) / stride + 1;
        size_t len_y = (original_dimensions[N - 1] - 1) / stride + 1;
        size_t even_len_x = len_x - len_x / 2;
        //size_t odd_len_x = len_x / 2;
        size_t even_len_y = len_y - len_y / 2;
        //size_t odd_len_y = len_y / 2;

        if(direction == 0){//slow (x) first
            for(size_t j = 0;j < even_len_y; j ++){
                auto pos = data +  j * stride;
                //just 1 column
                for(size_t i = 0;i < len_x;i ++)
                    aligned_buffer[i] = * (pos + i * stride * original_dim_offsets[N-2]); 
                interpolation_gathered_base(aligned_buffer, len_x, interp_func, quantize_func);
                 for(size_t i = 0;i < len_x;i ++)
                     *(pos + i * stride * original_dim_offsets[N-2]) = aligned_buffer[i]; 
            }

            for(size_t i = 0;i < len_x;i ++){
                auto pos = data + i * stride * original_dim_offsets[N-2];
                interpolation_gathered_1D(pos, stride, interp_func, quantize_func);
            }
        }
        else{//fast (y) first
            for(size_t i = 0;i < even_len_x;i ++){
                auto pos = data + i * stride * original_dim_offsets[N-2];
                interpolation_gathered_1D(pos, stride, interp_func, quantize_func);
            }

            for(size_t j = 0;j < len_y; j ++){
                auto pos = data +  j * stride;
                //just 1 column
                for(size_t i = 0;i < len_x;i ++)
                    aligned_buffer[i] = * (pos + i * stride * original_dim_offsets[N-2]); 
                interpolation_gathered_base(aligned_buffer, len_x, interp_func, quantize_func);
                 for(size_t i = 0;i < len_x;i ++)
                     *(pos + i * stride * original_dim_offsets[N-2]) = aligned_buffer[i]; 
            }

        }
        return predict_error;
    }


    template <class QuantizeFunc>
    double interpolation_gathered_3D(T * data,const size_t stride, const std::string &interp_func,
                            QuantizeFunc &&quantize_func, const int direction) {

        double predict_error = 0;

        auto offset_x = original_dim_offsets[N - 3];
        auto offset_y = original_dim_offsets[N - 2];
        size_t len_x = (original_dimensions[N - 3] - 1) / stride + 1;
        size_t len_y = (original_dimensions[N - 2] - 1) / stride + 1;
        size_t len_z = (original_dimensions[N - 1] - 1) / stride + 1;
        size_t even_len_x = len_x - len_x / 2;
        size_t even_len_y = len_y - len_y / 2;
        size_t even_len_z = len_z - len_z / 2;

        if(direction == 0){//slow (x) first


            

            //dim0. Grouping columns
            for(size_t j = 0;j < even_len_y;j ++){
                for(size_t k = 0;k < even_len_z;k += column_num){
                    auto pos = data + j * stride * offset_y + k * stride;
                    const auto col_count = std::min((column_num - 1) / stride + 1, (even_len_z - k) );
                    
                    for(size_t i = 0; i < original_dimensions[N - 3]; i += stride){
                        size_t buffer_idx = 0;
                        for(size_t kk = 0;kk < col_count * stride;kk += stride){
                            buffer[(buffer_idx++) * len_x  + i] = * (pos + i * offset_x + kk); 
                        }
                    }
                    for(size_t kk = 0;kk < col_count ;kk ++){
                        auto buffer_pos = buffer + kk * len_x;
                        interpolation_gathered_base(buffer_pos,len_x,interp_func,quantize_func);
                    }
                    for(size_t i = 0; i < original_dimensions[N-3]; i += stride){
                        size_t buffer_idx = 0;
                        for(size_t kk = 0;kk < col_count * stride;kk += stride){
                            * (pos + i * offset_x + kk) = buffer[(buffer_idx++) * len_x + i];
                        }
                    }

                }
            }

            for(size_t i = 0;i < original_dimensions[N-3];i += stride){
                auto pos = data + i * offset_x;
                interpolation_gathered_2D(pos,stride,interp_func,quantize_func, direction);
            }

        }
        else{//fast (z) first
            //dim0. Grouping columns

            for(size_t i = 0;i < even_len_x;i ++){
                auto pos = data + i * offset_x * stride;
                interpolation_gathered_2D(pos,stride,interp_func,quantize_func,direction);
            }

            for(size_t j = 0;j < len_y;j ++){
                for(size_t k = 0;k < len_z;k += column_num){
                    auto pos = data + j * stride * offset_y + k * stride;
                    const auto col_count = std::min((column_num - 1) / stride + 1, (len_z - k) );
                    
                    for(size_t i = 0; i < original_dimensions[N - 3]; i += stride){
                        size_t buffer_idx = 0;
                        for(size_t kk = 0;kk < col_count * stride;kk += stride){
                            buffer[(buffer_idx++) * len_x  + i] = * (pos + i * offset_x + kk); 
                        }
                    }
                    for(size_t kk = 0;kk < col_count ;kk ++){
                        auto buffer_pos = buffer + kk * len_x;
                        interpolation_gathered_base(buffer_pos,len_x,interp_func,quantize_func);
                    }
                    for(size_t i = 0; i < original_dimensions[N-3]; i += stride){
                        size_t buffer_idx = 0;
                        for(size_t kk = 0;kk < col_count * stride;kk += stride){
                            * (pos + i * offset_x + kk) = buffer[(buffer_idx++) * len_x + i];
                        }
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
            for (uint i = 1; i < N; i++) {
                begin_idx[dims[i]] = (begin[dims[i]] ? begin[dims[i]] + stride2x : 0);
                strides[dims[i]] = stride2x;
            }

            predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[0], strides, stride,
                                                                interp_func, quantize_func);
            for (uint i = 1; i < N; i++) {
                begin_idx[dims[i]] = begin[dims[i]];
                begin_idx[dims[i - 1]] = (begin[dims[i - 1]] ? begin[dims[i - 1]] + stride : 0);
                strides[dims[i - 1]] = stride;
                predict_error += interpolation_1d_fastest_dim_first(data, begin_idx, end_idx, dims[i], strides, stride,
                                                                    interp_func, quantize_func);
            }
            return predict_error;
        } else {
            throw std::runtime_error("Unsupported dimension in InterpolationDecomposition");
        }
    }

    //no begin and end idx. Just interp on the whole data array.
    template <class QuantizeFunc>
    double interpolation_gathered(T *data,
                         const std::string &interp_func, QuantizeFunc &&quantize_func, const int direction,
                         size_t stride = 1) {
        assert (N <= 3);
        if constexpr (N == 1) {  
            return interpolation_gathered_1D(data, stride, interp_func, quantize_func);
        } else if constexpr (N == 2) {  
            return interpolation_gathered_2D(data, stride, interp_func, quantize_func, direction);
        } else if constexpr (N == 3) {  // new API (for faster speed)
            return interpolation_gathered_3D(data, stride, interp_func, quantize_func, direction);
        } else {
            throw std::runtime_error("Unsupported dimension in Gathered Interpolation");
        }
    }



    int interp_level = -1;
    int interp_id;
    uint blocksize;
    std::vector<std::string> interpolators = {"linear", "cubic"};
    int *quant_inds;
    size_t quant_index = 0;
    double max_error;
    Quantizer quantizer;
    size_t num_elements;
    const size_t column_num = 64 / sizeof(T);
    size_t max_dim = 1;
    T * buffer, *aligned_buffer;

    std::array<size_t, N> original_dimensions;
    std::array<size_t, N> original_dim_offsets;
    std::vector<std::array<int, N>> dim_sequences;
    int direction_sequence_id;
    size_t anchor_stride = 0;
    double eb_alpha = -1;
    double eb_beta = -1;
    double eb_ratio = 0.5;  // To be deprecated
};

template <class T, uint N, class Quantizer>
InterpolationDecomposition<T, N, Quantizer> make_decomposition_interpolation(const Config &conf, Quantizer quantizer) {
    return InterpolationDecomposition<T, N, Quantizer>(conf, quantizer);
}

}  // namespace SZ3

#endif
