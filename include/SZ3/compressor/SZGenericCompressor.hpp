#ifndef SZ3_COMPRESSOR_TYPE_ONE_HPP
#define SZ3_COMPRESSOR_TYPE_ONE_HPP

#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {
/**
 * SZGenericCompressor glues together decomposition, encoder, and lossless modules to form the compressor.
 * It only takes Decomposition, not Predictor.
 * @tparam T original data type
 * @tparam N original data dimension
 * @tparam Decomposition decomposition module
 * @tparam Encoder encoder module
 * @tparam Lossless lossless module
 */
template <class T, uint N, class Decomposition, class Encoder, class Lossless>
class SZGenericCompressor : public concepts::CompressorInterface<T> {
   public:
    SZGenericCompressor(Decomposition decomposition, Encoder encoder, Lossless lossless)
        : decomposition(decomposition), encoder(encoder), lossless(lossless) {
        static_assert(std::is_base_of<concepts::DecompositionInterface<T, int, N>, Decomposition>::value,
                      "must implement the frontend interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) override {
        Timer timer(true);
        std::vector<int> quant_inds = decomposition.compress(conf, data);
        timer.stop("cmp interp");
        if (decomposition.get_out_range().first != 0) {
            throw std::runtime_error("The output range of the decomposition must start from 0 for this compressor");
        }
        //timer.start();
       // encoder.preprocess_encode(quant_inds, decomposition.get_out_range().second);
       // timer.stop("prehuff");
      

        auto * cmpData_pos = cmpData;

        decomposition.save(cmpData_pos);
        auto cmpSize = cmpData_pos - cmpData;
        //encoder.save(buffer_pos);

        //store the size of quant_inds is necessary as it is not always equal to conf.num
         timer.start();

        auto quant_inds_size =  quant_inds.size();
        write<size_t>(quant_inds_size, cmpData_pos);
        cmpSize += sizeof(size_t);

        #ifdef _OPENMP
        auto default_nthreads = omp_get_max_threads();
        //std::cout<<default_nthreads<<" "<<quant_inds_size<<std::endl;
        auto best_num_threads = std::min(default_nthreads, (int)(quant_inds_size / (1u<<16)));
        //std::cout<<best_num_threads<<std::endl;
        if (best_num_threads > 1) {
            omp_set_num_threads(best_num_threads);
            //uchar * offset_block_pos = buffer_pos + sizeof(int);
            //size_t bins_per_thread;
            auto quant_inds_data = quant_inds.data();
            std::vector<size_t>block_byte_offsets;
            size_t offset_chunk_size;
            size_t total_hz_size = 0;
            size_t nthreads;
            #pragma omp parallel
            {   

                #pragma omp single
                {
                    nthreads = omp_get_num_threads();
                    //std::cout<<nthreads<<std::endl;
                    block_byte_offsets.resize(nthreads);
                    write<int>(nthreads, cmpData_pos);
                    offset_chunk_size = nthreads * sizeof (size_t);
                }
                Encoder cur_encoder;
                auto tid = omp_get_thread_num();
                size_t start_idx = ((size_t)tid * quant_inds_size) / (size_t)nthreads, cur_len = ((size_t)(tid+1) * quant_inds_size) / (size_t)nthreads - start_idx;
                //#pragma omp critical
                //std::cout<<tid<<" "<<start_idx<<" "<<cur_len<<std::endl;
                size_t cur_bufferSize = std::max<size_t>(1000, 1.2 * sizeof(T) * cur_len);
                auto cur_buffer = static_cast<uchar *>(malloc(cur_bufferSize)); 
                auto cur_buffer_pos = cur_buffer;
                cur_encoder.preprocess_encode(quant_inds_data + start_idx, cur_len, decomposition.get_out_range().second);
                cur_encoder.save(cur_buffer_pos);
                //#pragma omp critical
                //std::cout<<tid<<" "<<cur_buffer_pos-cur_buffer <<std::endl;
                cur_encoder.encode(quant_inds_data + start_idx, cur_len, cur_buffer_pos);
                //std::cout<<tid<<" "<<encode_length<<std::endl;
                cur_encoder.postprocess_encode();

                auto cur_huffSize = cur_buffer_pos - cur_buffer;

                auto cur_cmp_buffer = static_cast<uchar *>(malloc(cur_bufferSize)); 
                auto cur_cmp_buffer_pos = cur_cmp_buffer;
                auto cur_cmpSize = lossless.compress(cur_buffer, cur_huffSize, cur_cmp_buffer_pos, cur_bufferSize);
                free(cur_buffer);
                //+ sizeof(size_t); //the original outsize doesn't contain the size header. Actually, since we already have the offset chunk, write the size in each block is a waste.
                                               //However, remove it will need to modify the huffman encoding api, which may bring compatability issue. So keep it now. 
                //#pragma omp critical
               // std::cout<<tid<<" "<<cur_outSize<<std::endl;

                block_byte_offsets[tid] = cur_cmpSize;
                // #pragma omp critical
                //std::cout<<"tid: "<<tid<<" outsize: "<<cur_outSize<<std::endl;
                #pragma omp barrier
                #pragma omp single
                {
                    size_t prefix_sum = 0;
                    for (size_t i = 0; i < nthreads; i++){
                        //std::cout<<"prefix, tid: "<<i<<", outsize: "<<block_byte_offsets[i]<<std::endl;
                        auto next_prefix_sum = prefix_sum + block_byte_offsets[i];
                        block_byte_offsets[i] = prefix_sum;
                        prefix_sum = next_prefix_sum;
                        //#pragma omp critical
                        //std::cout<<" offset: "<<block_byte_offsets[i]<<std::endl;

                    }
                }

                if(tid == nthreads - 1)
                    total_hz_size = block_byte_offsets[tid] + cur_cmpSize;
                auto temp_cmp_pos = cmpData_pos + tid * sizeof(size_t);
                write<size_t>(block_byte_offsets[tid],temp_cmp_pos);
                temp_cmp_pos = cmpData_pos + offset_chunk_size + block_byte_offsets[tid];
                write<uchar>(cur_cmp_buffer, cur_cmpSize, temp_cmp_pos);
                free(cur_cmp_buffer);

               


            }
            omp_set_num_threads(default_nthreads);
            //buffer_pos += offset_chunk_size + total_huffman_size;
            cmpSize += sizeof(int) + offset_chunk_size + total_hz_size;
        }
            
        else{
            
            write<int>(1, cmpData_pos); //1 thread
            write<size_t>(0, cmpData_pos); //offset = 0;
            size_t bufferSize = std::max<size_t>(
                1000, 1.2 * (decomposition.size_est() + encoder.size_est_without_init() + sizeof(T) * quant_inds.size()));
            std::cout<<bufferSize<<std::endl;
            auto buffer = static_cast<uchar *>(malloc(bufferSize));
            uchar *buffer_pos = buffer;
            encoder.preprocess_encode(quant_inds, decomposition.get_out_range().second);
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
             std::cout<<buffer_pos-buffer<<std::endl;
            encoder.postprocess_encode();
            cmpCap -= cmpData_pos - cmpData;
            std::cout<<cmpCap<<std::endl;
           // std::cout<<buffer_pos - buffer<<std::endl;
            auto zstdSize = lossless.compress(buffer, buffer_pos - buffer, cmpData_pos, cmpCap);
            std::cout<<zstdSize<<std::endl;
            cmpSize = cmpData_pos - cmpData + zstdSize;
            free(buffer);
        }


        #else
            write<int>(1, cmpData_pos); //1 thread
            write<size_t>(0, cmpData_pos); //offset = 0;
             size_t bufferSize = std::max<size_t>(
                1000, 1.2 * (decomposition.size_est() + encoder.size_est_without_init() + sizeof(T) * quant_inds.size()));

            auto buffer = static_cast<uchar *>(malloc(bufferSize));
            encoder.preprocess_encode(quant_inds, decomposition.get_out_range().second);
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
            cmpCap -= cmpData_pos - cmpData;
            auto zstdSize = lossless.compress(buffer, buffer_pos - buffer, cmpData_pos, cmpCap);
            cmpSize = cmpData_pos - cmpData + zstdSize;
            free(buffer);

        #endif
        timer.stop("huffzstd");
       //  timer.start();
      
        // timer.stop("zstd");

        return cmpSize;
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
       // uchar *buffer = nullptr;
        //size_t bufferSize = 0;
        //Timer timer(true);
        uchar const *cmpDataPos = cmpData;
         decomposition.load(cmpDataPos, cmpSize);
        //timer.stop("decmp zstd");

       
        

        size_t quant_inds_size = 0;
        read(quant_inds_size, cmpDataPos);
        //std::cout<<quant_inds_size<<std::endl;
        int compression_thread_num = 0;
        read(compression_thread_num, cmpDataPos);
        //std::cout<<compression_thread_num<<std::endl;
        std::vector<int> quant_inds; // todo: it should better match the encoder output type,
        if(compression_thread_num <=1){
            size_t offset;
            read(offset, cmpDataPos);
           // std::cout<<offset<<std::endl;
             uchar *buffer = nullptr;
            size_t bufferSize = 0;
            cmpSize -= cmpDataPos - cmpData;
            lossless.decompress(cmpDataPos, cmpSize, buffer, bufferSize);
             //std::cout<<bufferSize<<std::endl;
               uchar const *bufferPos = buffer;
            encoder.load(bufferPos, bufferSize);
            quant_inds = encoder.decode(bufferPos, quant_inds_size);
             encoder.postprocess_decode();
              free(buffer);
        }
        else{
            quant_inds.resize(quant_inds_size);
            std::vector<size_t>block_byte_offsets(compression_thread_num);

            read(block_byte_offsets.data(),compression_thread_num,cmpDataPos);
            cmpSize -= cmpDataPos - cmpData;
            #ifdef _OPENMP
            #pragma omp parallel for 
            #endif
            for(int tid=0; tid < compression_thread_num;tid++){
                
                size_t start_idx = ((size_t)tid * quant_inds_size) / (size_t)compression_thread_num, cur_len = ((size_t)(tid+1) * quant_inds_size) / (size_t)compression_thread_num - start_idx;
                size_t block_byte_offset = block_byte_offsets[tid];
                //std::cout<<tid<<" "<<block_byte_offset<<std::endl;

                size_t block_byte_length = (tid == compression_thread_num - 1) ? cmpSize-block_byte_offset: block_byte_offsets[tid + 1] - block_byte_offset;
                auto cur_cmpDataPos = cmpDataPos + block_byte_offset ;

                uchar *cur_buffer = nullptr;
                size_t cur_bufferSize = 0;
                lossless.decompress(cur_cmpDataPos, block_byte_length, cur_buffer, cur_bufferSize);
               // #pragma omp critical
                //std::cout<<"tid: "<<tid<<", prefix: "<<block_byte_offset<<std::endl;
                uchar const *cur_buffer_pos = cur_buffer;
                Encoder cur_encoder;
                cur_encoder.load(cur_buffer_pos,cur_bufferSize);

               //#pragma omp critical
               // std::cout<<"tid: "<<tid<<", loaded."<<std::endl;
               // #pragma omp critical
              // std::cout<<"tid: "<<tid<<", moved: "<<temp_buffer_pos - temp_start <<std::endl;
             
                auto cur_quant_inds = cur_encoder.decode(cur_buffer_pos, cur_len);
               // #pragma omp critical
               // std::cout<<"tid: "<<tid<<", ended at: "<<temp_buffer_pos - bufferPos<<" ,"<<cur_quant_inds.size()<<" bins extracted."<<std::endl;
                cur_encoder.postprocess_decode();
                free(cur_buffer);
               //  #pragma omp critical
               // std::cout<<"tid: "<<tid<<" postprocessed."<<std::endl;
                std::copy(cur_quant_inds.begin(), cur_quant_inds.end(), quant_inds.begin() + start_idx);
                cur_quant_inds.clear();
                


            }
             



        }

        //auto quant_inds = encoder.decode(bufferPos, quant_inds_size);
       

       

        decomposition.decompress(conf, quant_inds, decData);
        return decData;
    }

   private:
    Decomposition decomposition;
    Encoder encoder;
    Lossless lossless;
};

template <class T, uint N, class Decomposition, class Encoder, class Lossless>
std::shared_ptr<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>> make_compressor_sz_generic(
    Decomposition decomposition, Encoder encoder, Lossless lossless) {
    return std::make_shared<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>>(decomposition, encoder,
                                                                                         lossless);
}

}  // namespace SZ3
#endif