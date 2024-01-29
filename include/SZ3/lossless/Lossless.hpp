//
// Created by Kai Zhao on 6/12/20.
//

#ifndef SZ_LOSSLESS_HPP
#define SZ_LOSSLESS_HPP


namespace SZ3::concepts {

    /**
     * Lossless compressors is used in addition to lossy compression to further reduce the data size
     * Usually this module calls into existing lossless compress APIs, instead of re-implementing the lossless algorithms.
     */
    class LosslessInterface {
    public:

        /**
         * compress data with lossless compressors
         * @param src  data to be compressed
         * @param srcLen length (in bytes) of the data to be compressed
         * @param dst compressed data
         * @param dstCap capacity to store compressed size (in bytes), will be overwritten by the actual size in return
         */
        virtual void compress(uchar *src, size_t srcLen, uchar *dst, size_t &dstCap) = 0;

        /**
         * reverse of compress(), decompress the data with lossless compressors
         * @param src data to be decompressed
         * @param srcLen length (in bytes) of the data to be decompressed (as input) or the data decompressed (as output).
         * @param dst decompressed data
         * @param dstCap capacity to store decompressed data (in bytes), will be overwritten by the actual size in return
         */
        virtual void decompress(const uchar *src, const size_t srcLen, uchar *dst, size_t &dstCap) = 0;
    };
}


#endif //SZ_LOSSLESS_HPP
