#ifndef SZ3_SZALGO_INTERP_HPP
#define SZ3_SZALGO_INTERP_HPP

#ifdef _OPENMP
#include "SZ3/decomposition/InterpolationDecomposition_Omp.hpp"
#include "SZ3/quantizer/LinearQuantizer_Omp.hpp"
#endif

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

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_Interp(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_INTERP);
    calAbsErrorBound(conf, data);
    if (conf.interpAnchorStride < 0) {  // set default anchor stride
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interpAnchorStride = anchor_strides[N - 1];
    }

    #ifdef _OPENMP

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation_omp<T, N>(conf, LinearQuantizerOMP<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());

    #else
        auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    #endif

    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_Interp(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_INTERP);
    auto cmpDataPos = cmpData;
   // std::cout<<"decomp started"<<std::endl;
    #ifdef _OPENMP

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation_omp<T, N>(conf, LinearQuantizerOMP<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());

    #else
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    #endif
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

template <class T, uint N>
double interp_compress_test(
    const std::vector<std::vector<T>> &sampled_blocks, const Config conf, int block_size, uchar *cmpData,
    size_t cmpCap) {  // test interp cmp on a set of sampled data blocks and return the compression ratio

    Timer timer(true);
    /*
    #ifdef _OPENMP
     auto sz =
        make_decomposition_interpolation_omp<T, N>(conf, LinearQuantizerOMP<T>(conf.absErrorBound, conf.quantbinCnt / 2));
     #else
         auto sz =
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));
     #endif
   
    */

    std::vector<int> quant_inds;
    std::vector<std::vector<int> > quant_inds_vec(sampled_blocks.size());

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t k = 0; k < sampled_blocks.size(); ++k) {
        auto sz =
            make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));
        auto cur_block = sampled_blocks[k];
        quant_inds_vec[k] = sz.compress(conf, cur_block.data());
        
    }
    for (size_t k = 0; k < sampled_blocks.size(); ++k)
        quant_inds.insert(quant_inds.end(), std::make_move_iterator(quant_inds_vec[k].begin()),
                                 std::make_move_iterator(quant_inds_vec[k].end()));  // merge the quant bins. Lossless them together
    timer.stop("att interp");
    timer.start();
    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();
    encoder.preprocess_encode(quant_inds, conf.quantbinCnt);
    timer.stop("att prehuff");
    timer.start();
    size_t bufferSize = std::max<size_t>(
        1000, 1.2 * (encoder.size_est() + sizeof(T) * quant_inds.size()));

    auto buffer = static_cast<uchar *>(malloc(bufferSize));
    uchar *buffer_pos = buffer;
    //store the size of quant_inds is necessary as it is not always equal to conf.num


    auto quant_inds_size =  quant_inds.size();
    write<size_t>(quant_inds_size, buffer_pos);

    #ifdef _OPENMP
    auto default_nthreads = omp_get_max_threads();
    //std::cout<<default_nthreads<<" "<<quant_inds_size<<std::endl;
    auto best_num_threads = std::min(default_nthreads, (int)(quant_inds_size / (1024)));
    //std::cout<<best_num_threads<<std::endl;
    if (best_num_threads > 1) {
        omp_set_num_threads(best_num_threads);
        //uchar * offset_block_pos = buffer_pos + sizeof(int);
        //size_t bins_per_thread;
        auto quant_inds_data = quant_inds.data();
        std::vector<size_t>block_byte_offsets;
        size_t offset_chunk_size;
        size_t total_huffman_size = 0;
        size_t nthreads;
        #pragma omp parallel
        {   

            #pragma omp single
            {
                nthreads = omp_get_num_threads();
                //std::cout<<nthreads<<std::endl;
                block_byte_offsets.resize(nthreads);
                write<int>(nthreads, buffer_pos);
                offset_chunk_size = nthreads * sizeof (size_t);
            }

            auto tid = omp_get_thread_num();
            size_t start_idx = ((size_t)tid * quant_inds_size) / (size_t)nthreads, cur_len = ((size_t)(tid+1) * quant_inds_size) / (size_t)nthreads - start_idx;
            //#pragma omp critical
            //std::cout<<tid<<" "<<start_idx<<" "<<cur_len<<std::endl;
            size_t cur_bufferSize = std::max<size_t>(1000, 1.2 * sizeof(T) * cur_len);
            auto cur_buffer = static_cast<uchar *>(malloc(cur_bufferSize)); 
            auto cur_buffer_pos = cur_buffer;
            auto cur_outSize = encoder.encode(quant_inds_data + start_idx, cur_len, cur_buffer_pos);

            cur_outSize += sizeof(size_t); //the original outsize doesn't contain the size header. Actually, since we already have the offset chunk, write the size in each block is a waste.
                                           //However, remove it will need to modify the huffman encoding api, which may bring compatability issue. So keep it now. 

            block_byte_offsets[tid] = cur_outSize;
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
                   // std::cout<<" offset: "<<block_byte_offsets[i]<<std::endl;

                }
            }

            if(tid == nthreads - 1)
                total_huffman_size = block_byte_offsets[tid] + cur_outSize;
            auto temp_buffer_pos = buffer_pos + tid * sizeof(size_t);
            write<size_t>(block_byte_offsets[tid],temp_buffer_pos);
            temp_buffer_pos = buffer_pos + offset_chunk_size + block_byte_offsets[tid];
            write<uchar>(cur_buffer, cur_outSize, temp_buffer_pos);

            free(cur_buffer);


        }
        omp_set_num_threads(default_nthreads);
        buffer_pos += offset_chunk_size + total_huffman_size;
    }
        
    else{
        write<int>(1, buffer_pos); //1 thread
        write<size_t>(0, buffer_pos); //offset = 0;
        encoder.encode(quant_inds, buffer_pos);
    }


    #else
        write<int>(1, buffer_pos); //1 thread
        write<size_t>(0, buffer_pos); //offset = 0;
        encoder.encode(quant_inds, buffer_pos);

    #endif
    timer.stop("att huff");
    timer.start();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    timer.stop("att zstd");
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

           Timer timer(true);

    calAbsErrorBound(conf, data);
    timer.stop("abseb compute");
    timer.start();
    if (conf.interpAnchorStride < 0) {  // set default anchor stride
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interpAnchorStride = anchor_strides[N - 1];
    }

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

    if (!to_tune) {  // if the sampled data would be too many (currently it is 5% of the input), skip the tuning
        conf.cmprAlgo = ALGO_INTERP;
        return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    }
    timer.stop("preparation");
    timer.start();
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
        conf.cmprAlgo = ALGO_INTERP;
        return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    }
    timer.stop("sampling");
    timer.start();
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
        for (auto &interp_op : {INTERP_ALGO_LINEAR,INTERP_ALGO_CUBIC}) {
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
        const int ablist_size = 3;
       auto alphalist = std::array<double,ablist_size>{1.0, 1.5, 2.0};
        auto betalist = std::array<double,ablist_size>{1.0, 2.5, 3.0};
        std::array<double,ablist_size> ratios;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < ablist_size; i++) {
            auto tempConfig = testConfig;
            tempConfig.interpAlpha = alphalist[i];
            tempConfig.interpBeta = betalist[i];
            ratios[i] = interp_compress_test<T, N>(sampled_blocks, tempConfig, sampleBlockSize, buffer, bufferCap);
            
        }
        for (size_t i = 0; i < ablist_size; i++) {
            auto ratio = ratios[i];
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpAlpha = alphalist[i];
                conf.interpBeta = betalist[i];
            }
        }
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


    bool useInterp = !(best_lorenzo_ratio >= best_interp_ratio * 1.1 && best_lorenzo_ratio < 50 &&
                       best_interp_ratio < 50);  // 1.1 is a fix coefficient. subject to revise
    size_t cmpSize = 0;
    timer.stop("interp tuning");
    timer.start();
    if (useInterp) {
        conf.cmprAlgo = ALGO_INTERP;
        cmpSize = SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
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
        //            double tuning_time = timer.stop();
        cmpSize = SZ_compress_LorenzoReg<T, N>(conf, data, cmpData, cmpCap);
    }

    free(buffer);
    return cmpSize;
}
}  // namespace SZ3
#endif
