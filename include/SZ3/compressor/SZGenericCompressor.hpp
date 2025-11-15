#ifndef SZ3_COMPRESSOR_TYPE_ONE_HPP
#define SZ3_COMPRESSOR_TYPE_ONE_HPP

#include <cstring>

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
        //Timer timer(true);  
        std::vector<int> quant_inds = decomposition.compress(conf, data);
        std::vector<int> suffix_quant_inds;
        if(quant_inds_size > conf.num){
            suffix_quant_inds.insert(suffix_quant_inds.end(),quant_inds.begin()+conf.num,quant_inds.end);
            quant_inds.resize(conf.num;)
        }
        //timer.stop("interpquant");
        //timer.start();
        if (decomposition.get_out_range().first != 0) {
            throw std::runtime_error("The output range of the decomposition must start from 0 for this compressor");
        }
        encoder.preprocess_encode(quant_inds, decomposition.get_out_range().second);
        //timer.stop("prepro_huff");
        //timer.start();
        size_t bufferSize = std::max<size_t>(
            1000, 1.2 * (decomposition.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size()+ sizeof(T) * suffix_inds.size()));

        auto buffer = static_cast<uchar *>(malloc(bufferSize));
        uchar *buffer_pos = buffer;
        if(suffix_quant_inds > conf.num)
            *buffer_pos = 1;
        else
            *buffer_pos = 0;
        buffer_pos++;

        decomposition.save(buffer_pos);
        encoder.save(buffer_pos);

        //store the size of quant_inds is necessary as it is not always equal to conf.num
        write<size_t>(quant_inds.size(), buffer_pos);
        //timer.stop("memloc");
        //timer.start();
        encoder.encode(quant_inds, buffer_pos);
        encoder.postprocess_encode();
        

        size_t huff_size = buffer_pos - buffer;
        if(suffix_quant_inds > conf.num){
            auto old_pos = buffer_pos;
            encoder.preprocess_encode(suffix_quant_inds, decomposition.get_out_range().second);
            write<size_t>(quant_inds_suffix.size(), buffer_pos);
            
            encoder.encode(quant_inds_suffix, buffer_pos);
            auto suffix_huff_size = buffer_pos - old_pos;
            if (suffix_huff_size >= huff_size / 10){
                buffer[0] = 0;
                buffer_pos = old_pos;
            }
            encoder.postprocess_encode();

        }
        //timer.stop("huff");
       // timer.start();
        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
      //  timer.stop("zstd");
      //  timer.start();
        free(buffer);
        
        //std::cout<<"compress ended."<<std::endl; 

        return cmpSize;
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        //std::cout<<"dec star."<<std::endl; 
        uchar *buffer = nullptr;
        size_t bufferSize = 0;
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);

        uchar const *bufferPos = buffer;

        decomposition.load(bufferPos, bufferSize);
        encoder.load(bufferPos, bufferSize);

        bool have_suffix = *buffer_pos;
        buffer_pos++;

        size_t quant_inds_size = 0;
        read(quant_inds_size, bufferPos);
        auto quant_inds = encoder.decode(bufferPos, quant_inds_size);
        encoder.postprocess_decode();

        if(have_suffix){
            size_t suffix_quant_inds_size = 0;
            read(suffix_quant_inds_size, bufferPos);
            auto suffix_quant_inds = encoder.decode(bufferPos, quant_inds_size);
            encoder.postprocess_decode();
            quant_inds.resize(quant_inds.size()+suffix_quant_inds.size());
            quant_inds.insert(quant_inds.end(),suffix_quant_inds.begin(),suffix_quant_inds,end())
            suffix_quant_inds.clear();

        }


        free(buffer);

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
