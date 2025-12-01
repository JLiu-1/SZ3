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
        std::vector<int> quant_inds = decomposition.compress(conf, data);

        if (decomposition.get_out_range().first != 0) {
            throw std::runtime_error("The output range of the decomposition must start from 0 for this compressor");
        }
       
        size_t bufferSize = std::max<size_t>(
            1000, 1.2 * (decomposition.size_est() + encoder.size_est_without_init() + sizeof(T) * quant_inds.size()));

        auto buffer = static_cast<uchar *>(malloc(bufferSize));
        uchar *buffer_pos = buffer;

        decomposition.save(buffer_pos);
        

        //store the size of quant_inds is necessary as it is not always equal to conf.num
        size_t quant_inds_size = quant_inds.size();


        write<size_t>(quant_inds_size, buffer_pos);

        size_t quant_inds_size_one_eighth = quant_inds_size / 8;

        //auto buffer_pos_start = buffer_pos;
        //buffer_pos += sizeof(size_t);
        auto quant_inds_pos = quant_inds.data(); 


        encoder.preprocess_encode(quant_inds_pos,quant_inds_size_one_eighth, decomposition.get_out_range().second);
        encoder.save(buffer_pos);
        encoder.encode(quant_inds_pos,quant_inds_size_one_eighth, buffer_pos);
        encoder.postprocess_encode();
        encoder.preprocess_encode(quant_inds_pos+quant_inds_size_one_eighth, quant_inds_size-quant_inds_size_one_eighth, decomposition.get_out_range().second);
        encoder.save(buffer_pos);
        encoder.encode(quant_inds_pos+quant_inds_size_one_eighth, quant_inds_size-quant_inds_size_one_eighth, buffer_pos);
        encoder.postprocess_encode();
        
        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
        free(buffer);

        return cmpSize;
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        uchar *buffer = nullptr;
        size_t bufferSize = 0;
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);

        uchar const *bufferPos = buffer;

        decomposition.load(bufferPos, bufferSize);
      

        size_t quant_inds_size = 0;
        read(quant_inds_size, bufferPos);
        size_t quant_inds_size_one_eighth = quant_inds_size / 8;

        encoder.load(bufferPos, bufferSize);

        auto quant_inds = encoder.decode(bufferPos, quant_inds_size_one_eighth);
        encoder.postprocess_decode();

        encoder.load(bufferPos, bufferSize);

        auto quant_inds_2 = encoder.decode(bufferPos, quant_inds_size - quant_inds_size_one_eighth);
        encoder.postprocess_decode();


        free(buffer);
        quant_inds.reserve(quant_inds_size);
        quant_inds.insert(quant_inds.end(),std::make_move_iterator(quant_inds_2.begin()),std::make_move_iterator(quant_inds_2.end()));



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