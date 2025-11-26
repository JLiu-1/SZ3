#ifndef SZ3_SZALGO_INTERP_HPP
#define SZ3_SZALGO_INTERP_HPP

#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/decomposition/InterpolationDecomposition.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Sample.hpp"
#include "SZ3/utils/Statistic.hpp"


//#include "SZ3/sperr/SPERR3D_OMP_C.h"
#include "SZ3/sperr/SPECK2D_FLT.h"
#include "SZ3/sperr/SPECK3D_FLT.h"
//#include "SZ3/sperr/SPERR3D_OMP_D.h"
#include "SZ3/sperr/sperr_helper.h"


namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_Interp(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_INTERP);
    calAbsErrorBound(conf, data);
    /*
    if (conf.interpAnchorStride < 0) {  // set default anchor stride
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interpAnchorStride = anchor_strides[N - 1];
    }*/

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_Interp(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_INTERP);
    auto cmpDataPos = cmpData;

    if(N==3){
        auto ori_dims = conf.dims;
        std::array<size_t,3> downsampled_dims = {(conf.dims[0]+1)/2,(conf.dims[1]+1)/2,(conf.dims[2]+1)/2};
        size_t downsampled_num = downsampled_dims [0] * downsampled_dims [1] * downsampled_dims [2];
        size_t offset_y = ori_dims[2], offset_x = ori_dims[1] * offset_y;
        size_t downsampled_offset_y = downsampled_dims[2], downsampled_offset_x = downsampled_dims[1] * downsampled_offset_y;
  
        auto decompressor = std::make_unique<sperr::SPECK3D_FLT>();

          
        //const auto chunks = sperr::dims_type{1024,1024,1024};//ori 256^3, to tell the truth this is not large enough for scale but I just keep it, maybe set it large later.
        const auto sperr_dims = sperr::dims_type{downsampled_dims[2],downsampled_dims[1],downsampled_dims[0]};
        decompressor->set_dims(sperr_dims);
        size_t SPERR_cmpSize;
        read(SPERR_cmpSize,cmpDataPos,cmpSize);



        decompressor->use_bitstream(cmpDataPos, SPERR_cmpSize);
        cmpDataPos += SPERR_cmpSize;
        cmpSize -= SPERR_cmpSize;
        //compressor->set_tolerance(conf.absErrorBound);



        /*
        if (std::is_same<T, double>::value)
            rtn = compressor.copy_data<double>(reinterpret_cast<const double*>(data), conf.num,
                                    {conf.dims[2], conf.dims[1], conf.dims[0]}, {chunks[0], chunks[1], chunks[2]});
        else
            rtn = compressor.copy_data<float>(reinterpret_cast<const float*>(data), conf.num,
                                    {conf.dims[2], conf.dims[1], conf.dims[0]}, {chunks[0], chunks[1], chunks[2]});
        compressor.set_tolerance(conf.absErrorBound);*/
        /*
        if (std::is_same<T, double>::value)
            compressor->compress(reinterpret_cast<const double*>(data), conf.num);
        else{
            compressor->compress(reinterpret_cast<const float*>(data), conf.num);
        }*/
        decompressor->decompress();
        
        

        const auto decData_downsampled = decompressor->release_decoded_data();

        for(size_t i = 0; i < downsampled_dims[0]; i++){
            for(size_t j = 0; j < downsampled_dims[1]; j++){
                for(size_t k = 0; k < downsampled_dims[2]; k++){
                    auto downsampled_idx = i * downsampled_offset_x + j * downsampled_offset_y + k;
                    auto ori_idx = 2 * i * offset_x + 2 * j + offset_y + 2 * k;
                    decData[ori_idx] = decData_downsampled[downsampled_idx]; 

                }

            }
        }


        decompressor.reset();
        
        //return outData;
    }



   // std::cout<<"decomp started"<<std::endl; 
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

template <class T, uint N>
double interp_compress_test(
    const std::vector<std::vector<T>> sampled_blocks, const Config conf, int block_size, uchar *cmpData,
    size_t cmpCap) {  // test interp cmp on a set of sampled data blocks and return the compression ratio
    auto sz =
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));

    std::vector<int> total_quant_bins;
    for (size_t k = 0; k < sampled_blocks.size(); k++) {
        auto cur_block = sampled_blocks[k];
        auto quant_bins = sz.compress(conf, cur_block.data());
        total_quant_bins.insert(total_quant_bins.end(), quant_bins.begin(),
                                quant_bins.end());  // merge the quant bins. Lossless them together
    }

    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();

    encoder.preprocess_encode(total_quant_bins, sz.get_out_range().second);
    size_t bufferSize =
        std::max<size_t>(1000, 1.2 * (sz.size_est() + encoder.size_est() + sizeof(T) * total_quant_bins.size()));

    auto buffer = static_cast<uchar *>(malloc(bufferSize));
    uchar *buffer_pos = buffer;
    sz.save(buffer_pos);
    encoder.save(buffer_pos);
    // store the size of quant_inds is necessary as it is not always equal to conf.num
    write<size_t>(total_quant_bins.size(), buffer_pos);
    encoder.encode(total_quant_bins, buffer_pos);
    encoder.postprocess_encode();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    free(buffer);
    auto compression_ratio = conf.num * sampled_blocks.size() * sizeof(T) * 1.0 / cmpSize;
    return compression_ratio;
}

template <class T, uint N>
double lorenzo_compress_test(
    const std::vector<std::vector<T>> sampled_blocks, const Config &conf, uchar *cmpData,
    size_t cmpCap) {  // test lorenzo cmp on a set of sampled data blocks and return the compression ratio
    std::vector<int> total_quant_bins;
    // if ((N == 3 && !conf.regression2) || (N == 1 && !conf.regression && !conf.regression2)) {
    std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors;
    predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
    predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
    auto sz = make_decomposition_blockwise<T, N>(conf, ComposedPredictor<T, N>(predictors),
                                                 LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));
    // auto sz = make_decomposition_lorenzo_regression<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound,
    // conf.quantbinCnt / 2));
    for (size_t k = 0; k < sampled_blocks.size(); k++) {
        auto cur_block = sampled_blocks[k];
        auto quant_bins = sz.compress(conf, cur_block.data());
        total_quant_bins.insert(total_quant_bins.end(), quant_bins.begin(),
                                quant_bins.end());  // merge the quant bins. Lossless them together
    }
    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();
    encoder.preprocess_encode(total_quant_bins, conf.quantbinCnt);
    size_t bufferSize = std::max<size_t>(1000, 1.2 * (encoder.size_est() + sizeof(T) * total_quant_bins.size()));

    auto buffer = static_cast<uchar *>(malloc(bufferSize));
    uchar *buffer_pos = buffer;
    sz.save(buffer_pos);
    encoder.save(buffer_pos);

    // store the size of quant_inds is necessary as it is not always equal to conf.num
    write<size_t>(total_quant_bins.size(), buffer_pos);
    encoder.encode(total_quant_bins, buffer_pos);
    encoder.postprocess_encode();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    free(buffer);
    auto compression_ratio = conf.num * sampled_blocks.size() * sizeof(T) * 1.0 / cmpSize;
    return compression_ratio;
    // }
    // else{
    //     return 0.0;
    // }
}

template <class T, uint N>
size_t SZ_compress_Interp_lorenzo(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == ALGO_INTERP_LORENZO);

    //        Timer timer(true);
    calAbsErrorBound(conf, data);

    // sperr
    //uchar * SPERR_cmpData = nullptr;
    size_t SPERR_cmpSize = 0; 
    auto cmpDataPos = cmpData;
    if(N==3){
        auto ori_dims = conf.dims;
        std::array<size_t,3> downsampled_dims = {(conf.dims[0]+1)/2,(conf.dims[1]+1)/2,(conf.dims[2]+1)/2};
        size_t downsampled_num = downsampled_dims [0] * downsampled_dims [1] * downsampled_dims [2];
        size_t offset_y = ori_dims[2], offset_x = ori_dims[1] * offset_y;
        size_t downsampled_offset_y = downsampled_dims[2], downsampled_offset_x = downsampled_dims[1] * downsampled_offset_y;
        std::vector<double> downsampled_data(downsampled_num);
        for(size_t i = 0; i < downsampled_dims[0]; i++){
            for(size_t j = 0; j < downsampled_dims[1]; j++){
                for(size_t k = 0; k < downsampled_dims[2]; k++){
                    auto downsampled_idx = i * downsampled_offset_x + j * downsampled_offset_y + k;
                    auto ori_idx = 2 * i * offset_x + 2 * j + offset_y + 2 * k;
                    downsampled_data [downsampled_idx] = data [ori_idx]; 

                }

            }
        }

       

        double q_coeff = 1.5;
        auto compressor = std::make_unique<sperr::SPECK3D_FLT>();
        //compressor->set_num_threads(1);
        compressor->set_eb_coeff(q_coeff);
        compressor->take_data(std::move(downsampled_data));
        //auto rtn = sperr::RTNType::Good;
          
        //const auto chunks = sperr::dims_type{1024,1024,1024};//ori 256^3, to tell the truth this is not large enough for scale but I just keep it, maybe set it large later.
        const auto sperr_dims = sperr::dims_type{downsampled_dims[2],downsampled_dims[1],downsampled_dims[0]};
        compressor->set_dims(sperr_dims);
        compressor->set_tolerance(conf.absErrorBound);



        /*
        if (std::is_same<T, double>::value)
            rtn = compressor.copy_data<double>(reinterpret_cast<const double*>(data), conf.num,
                                    {conf.dims[2], conf.dims[1], conf.dims[0]}, {chunks[0], chunks[1], chunks[2]});
        else
            rtn = compressor.copy_data<float>(reinterpret_cast<const float*>(data), conf.num,
                                    {conf.dims[2], conf.dims[1], conf.dims[0]}, {chunks[0], chunks[1], chunks[2]});
        compressor.set_tolerance(conf.absErrorBound);*/
        /*
        if (std::is_same<T, double>::value)
            compressor->compress(reinterpret_cast<const double*>(data), conf.num);
        else{
            compressor->compress(reinterpret_cast<const float*>(data), conf.num);
        }*/
        compressor->compress();
        
        sperr::vec8_type stream(128);
        compressor->append_encoded_bitstream(stream);
        
            
        //SPERR_cmpData = new uchar[stream.size()];
        SPERR_cmpSize = stream.size();
        //std::cout<<outSize<<std::endl;

        //memcpy(outData,stream.data(),stream.size());//maybe not efficient
        

        write(SPERR_cmpSize,cmpDataPos);
        write(stream.data(),stream.size(),cmpDataPos);
        cmpCap -= sizeof(size_t) + SPERR_cmpSize;
        stream.clear();
        stream.shrink_to_fit();

        auto decData = compressor->release_decoded_data();

        for(size_t i = 0; i < downsampled_dims[0]; i++){
            for(size_t j = 0; j < downsampled_dims[1]; j++){
                for(size_t k = 0; k < downsampled_dims[2]; k++){
                    auto downsampled_idx = i * downsampled_offset_x + j * downsampled_offset_y + k;
                    auto ori_idx = 2 * i * offset_x + 2 * j + offset_y + 2 * k;
                    data[ori_idx] = decData[downsampled_idx]; 

                }

            }
        }


        compressor.reset();
        writefile("aftersperr.sz3.out",data,conf.num);
        //return outData;
    }






    /*
    if (conf.interpAnchorStride < 0) {  // set default anchor stride
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interpAnchorStride = anchor_strides[N - 1];
    }*/

    std::array<double, 4> sample_Rates = {0.005, 0.005, 0.005,
                                          0.005};  // default data sample rate. todo: add a config var to control
    auto sampleRate = sample_Rates[N - 1];
    std::array<size_t, 4> sampleBlock_Sizes = {4096, 128, 32,
                                               16};  // default sampled data block rate. Should better be no smaller
                                                     // than the anchor stride. todo: add a config var to control
    size_t sampleBlockSize = sampleBlock_Sizes[N - 1];
    size_t shortest_edge = conf.dims[0];
    for (size_t i = 0; i < N; i++) {
        shortest_edge = conf.dims[i] < shortest_edge ? conf.dims[i] : shortest_edge;
    }
    // Automatically adjust sampleblocksize.
    while (sampleBlockSize >= shortest_edge) sampleBlockSize /= 2;
    while (sampleBlockSize >= 16 && (pow(sampleBlockSize + 1, N) / conf.num) > 1.5 * sampleRate) sampleBlockSize /= 2;
    if (sampleBlockSize < 8) sampleBlockSize = 8;

    bool to_tune = pow(sampleBlockSize + 1, N) <= 0.05 * conf.num;  // to further revise
    for (auto &dim : conf.dims) {
        if (dim < sampleBlockSize) {
            to_tune = false;
            break;
        }
    }
    bool useInterp = true;
     size_t cmpSize = 0;


    

    if (!to_tune) {  // if the sampled data would be too many (currently it is 5% of the input), skip the tuning
        //conf.cmprAlgo = ALGO_INTERP;
        useInterp = true;
        //return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    }
    else{
        std::vector<std::vector<T>> sampled_blocks;
        size_t per_block_ele_num = pow(sampleBlockSize + 1, N);
        size_t sampling_num;
        std::vector<std::vector<size_t>> starts;
        auto profStride = sampleBlockSize / 4;  // larger is faster, smaller is better
        profiling_block<T, N>(data, conf.dims, starts, sampleBlockSize, conf.absErrorBound,
                              profStride);  // filter out the non-constant data blocks
        size_t num_filtered_blocks = starts.size();
        bool profiling = num_filtered_blocks * per_block_ele_num >= 0.5 * sampleRate * conf.num;  // temp. to refine
        // bool profiling = false;
        sampleBlocks<T, N>(data, conf.dims, sampleBlockSize, sampled_blocks, sampleRate, profiling,
                           starts);  // sample out same data blocks
        sampling_num = sampled_blocks.size() * per_block_ele_num;

        if (sampling_num == 0 || sampling_num >= conf.num * 0.2) {
            //conf.cmprAlgo = ALGO_INTERP;
            useInterp = true;
        }
        else{
            double best_lorenzo_ratio = 0, best_interp_ratio = 0, ratio;
            size_t bufferCap = conf.num * sizeof(T);
            auto buffer = static_cast<uchar *>(malloc(bufferCap));
            Config lorenzo_config = conf;

            {
                // tune interp
                conf.interpDirection = 0;
                conf.interpAlpha = 1.25;
                conf.interpBeta = 2.0;
                auto testConfig = conf;
                std::vector<size_t> dims(N, sampleBlockSize + 1);
                testConfig.setDims(dims.begin(), dims.end());
                for (auto &interp_op : {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC}) {
                    testConfig.interpAlgo = interp_op;
                    ratio = interp_compress_test<T, N>(sampled_blocks, testConfig, sampleBlockSize, buffer, bufferCap);
                    if (ratio > best_interp_ratio) {
                        best_interp_ratio = ratio;
                        conf.interpAlgo = interp_op;
                    }
                }

                testConfig.interpAlgo = conf.interpAlgo;
                testConfig.interpDirection = factorial(N) - 1;
                ratio = interp_compress_test<T, N>(sampled_blocks, testConfig, sampleBlockSize, buffer, bufferCap);
                if (ratio > best_interp_ratio * 1.02) {
                    best_interp_ratio = ratio;
                    conf.interpDirection = testConfig.interpDirection;
                }
                testConfig.interpDirection = conf.interpDirection;
                // test more alpha-beta pairs for best compression ratio,
                /*
                auto alphalist = std::vector<double>{1.0, 1.5, 2.0};
                auto betalist = std::vector<double>{1.0, 2.5, 3.0};
                for (size_t i = 0; i < alphalist.size(); i++) {
                    auto alpha = alphalist[i];
                    auto beta = betalist[i];
                    testConfig.interpAlpha = alpha;
                    testConfig.interpBeta = beta;
                    ratio = interp_compress_test<T, N>(sampled_blocks, testConfig, sampleBlockSize, buffer, bufferCap);
                    if (ratio > best_interp_ratio * 1.02) {
                        best_interp_ratio = ratio;
                        conf.interpAlpha = alpha;
                        conf.interpBeta = beta;
                    }
                }*/
            }
            {
                // only test lorenzo for 1D
                if (N == 1 && best_interp_ratio < 50) {
                    std::vector<size_t> sample_dims(N, sampleBlockSize + 1);
                    lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
                    lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
                    lorenzo_config.lorenzo = true;
                    lorenzo_config.lorenzo2 = true;
                    lorenzo_config.regression = false;
                    lorenzo_config.regression2 = false;
                    lorenzo_config.openmp = false;
                    lorenzo_config.blockSize = 5;
                    //        lorenzo_config.quantbinCnt = 65536 * 2;
                    best_lorenzo_ratio = lorenzo_compress_test<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
                    //            delete[]cmprData;
                    //    printf("Lorenzo ratio = %.2f\n", ratio);
                }
            }

            useInterp = !(best_lorenzo_ratio >= best_interp_ratio * 1.1 && best_lorenzo_ratio < 50 &&
                               best_interp_ratio < 50);  // 1.1 is a fix coefficient. subject to revise


            if(!useInterp){
                if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
                    auto quant_num = lorenzo_config.quantbinCnt;
                    lorenzo_config.quantbinCnt = 16384;
                    ratio = lorenzo_compress_test<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
                    if (ratio > best_lorenzo_ratio * 1.02) {
                        best_lorenzo_ratio = ratio;
                    } else {
                        lorenzo_config.quantbinCnt = quant_num;
                    }
                }
                lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
                conf = lorenzo_config;
            }
            free(buffer);
       
        }
 
    }
    if (useInterp) {
        conf.cmprAlgo = ALGO_INTERP;
        
        cmpSize = SZ_compress_Interp<T, N>(conf, data, cmpDataPos, cmpCap);
        if(N==3)
            cmpSize += sizeof(size_t) + SPERR_cmpSize;

        writefile("aftercmp.sz3.out",data,conf.num);


    } else {
        // no need to tune lorenzo for 3D anymore
        // if (N == 3) {
        //     float pred_freq, mean_freq;
        //     T mean_guess;
        //     lorenzo_config.quantbinCnt = optimize_quant_invl_3d<T>(
        //         data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound, pred_freq, mean_freq,
        //         mean_guess);
        //     lorenzo_config.pred_dim = 2;
        //     ratio  =
        //         lorenzo_compress_test<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
        //     if (ratio > best_lorenzo_ratio * 1.02) {
        //         best_lorenzo_ratio = ratio;
        //     } else {
        //         lorenzo_config.pred_dim = 3;
        //     }
        // }

        
        //            double tuning_time = timer.stop();
        cmpSize = SZ_compress_LorenzoReg<T, N>(conf, data, cmpData, cmpCap);
    }

    
    return cmpSize;
}
}  // namespace SZ3
#endif
