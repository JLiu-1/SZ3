#ifndef SZ3_HUFFMAN_ENCODER_HPP
#define SZ3_HUFFMAN_ENCODER_HPP

#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#if INTPTR_MAX == INT64_MAX  // 64bit system
#include "SZ3/utils/ska_hash/unordered_map.hpp"
#endif  // INTPTR_MAX == INT64_MAX
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace SZ3 {

template <class T>
class HuffmanEncoder : public concepts::EncoderInterface<T> {
   public:
    typedef struct node_t {
        struct node_t *left, *right;
        size_t freq;
        char t;  // in_node:0; otherwise:1
        T c;
    } *node;

    typedef struct HuffmanTree {
        unsigned int stateNum;
        unsigned int allNodes;
        struct node_t *pool;
        node *qqq, *qq;  // the root node of the HuffmanTree is qq[1]
        int n_nodes;     // n_nodes is for compression
        int qend;
        uint64_t *code;
        unsigned char *cout;
        int n_inode;  // n_inode is for decompression
        int maxBitCount;
    } HuffmanTree;

    HuffmanEncoder() {
        int x = 1;
        char *y = reinterpret_cast<char *>(&x);
        if (*y == 1)
            sysEndianType = 0;
        else  //=0
            sysEndianType = 1;
    }

    ~HuffmanEncoder() override { SZ_FreeHuffman(); }

    // build huffman tree
    HuffmanTree *createHuffmanTree(int stateNum) {
        HuffmanTree *huffmanTree = static_cast<HuffmanTree *>(malloc(sizeof(HuffmanTree)));
        memset(huffmanTree, 0, sizeof(HuffmanTree));
        huffmanTree->stateNum = stateNum;
        huffmanTree->allNodes = 2 * stateNum;

        huffmanTree->pool = static_cast<struct node_t *>(malloc(huffmanTree->allNodes * 2 * sizeof(struct node_t)));
        huffmanTree->qqq = static_cast<node *>(malloc(huffmanTree->allNodes * 2 * sizeof(node)));
        huffmanTree->code = static_cast<uint64_t *>(malloc(huffmanTree->stateNum * sizeof(uint64_t)));
        huffmanTree->cout = static_cast<unsigned char *>(malloc(huffmanTree->stateNum * sizeof(unsigned char)));

memset(huffmanTree->code, 0, huffmanTree->stateNum * sizeof(uint64_t));
memset(huffmanTree->cout, 0, huffmanTree->stateNum * sizeof(unsigned char));

        memset(huffmanTree->pool, 0, huffmanTree->allNodes * 2 * sizeof(struct node_t));
        memset(huffmanTree->qqq, 0, huffmanTree->allNodes * 2 * sizeof(node));
        memset(huffmanTree->code, 0, huffmanTree->stateNum * sizeof(uint64_t *));
        memset(huffmanTree->cout, 0, huffmanTree->stateNum * sizeof(unsigned char));
        huffmanTree->qq = huffmanTree->qqq - 1;
        huffmanTree->n_nodes = 0;
        huffmanTree->n_inode = 0;
        huffmanTree->qend = 1;

        return huffmanTree;
    }

    /**
     * build huffman tree using bins
     * @param bins
     * @param stateNum
     */
    void preprocess_encode(const std::vector<T> &bins, int stateNum) override {
        preprocess_encode(bins.data(), bins.size(), stateNum);
    }

    /**
     * build huffman tree using bins
     * @param bins
     * @param num_bin
     * @param stateNum
     */
    void preprocess_encode(const T *bins, size_t num_bin, int stateNum) {
        nodeCount = 0;
        if (num_bin == 0) {
            throw std::invalid_argument("Huffman bins should not be empty");
        }
        init(bins, num_bin);
        for (unsigned int i = 0; i < huffmanTree->stateNum; ++i)
            if (huffmanTree->cout[i] != 0) ++nodeCount;
        nodeCount = nodeCount * 2 - 1;
    }

    // save the huffman Tree in the compressed data
    void save(uchar *&c) override {
        // auto cc = c;
        write(offset, c);
        int32ToBytes_bigEndian(c, nodeCount);
        c += sizeof(int);
        int32ToBytes_bigEndian(c, huffmanTree->stateNum / 2);
        c += sizeof(int);
        uint totalSize = 0;  // = convert_HuffTree_to_bytes_anyStates(nodeCount, c);
        // std::cout << "nodeCount = " << nodeCount << std::endl;
        if (nodeCount <= 256)
            totalSize = convert_HuffTree_to_bytes_anyStates<unsigned char>(nodeCount, c);
        else if (nodeCount <= 65536)
            totalSize = convert_HuffTree_to_bytes_anyStates<unsigned short>(nodeCount, c);
        else
            totalSize = convert_HuffTree_to_bytes_anyStates<unsigned int>(nodeCount, c);
        c += totalSize;
        //            return c - cc;
    }

    size_t size_est() override {
        size_t b = (nodeCount <= 256) ? sizeof(unsigned char)
                                      : ((nodeCount <= 65536) ? sizeof(unsigned short) : sizeof(unsigned int));
        return 1 + 2 * nodeCount * b + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T) + sizeof(int) +
               sizeof(int) + sizeof(T);
    }

    // perform encoding
    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        return encode(bins.data(), bins.size(), bytes);
    }

    // perform encoding
    size_t encode(const T *bins, size_t num_bin, uchar *&bytes) {
        uchar *out_begin = bytes + sizeof(size_t);
        uchar *p = out_begin;

        uint64_t bitbuf = 0;   // LSB-first bit buffer
        unsigned nbits = 0;    // 当前 buffer 中已有多少 bit（低位开始）

        uint64_t *code = huffmanTree->code;
        unsigned char *len = huffmanTree->cout;

        for (size_t i = 0; i < num_bin; ++i) {
            int state = bins[i] - offset;
            assert(state >= 0 && static_cast<unsigned>(state) < huffmanTree->stateNum);

            uint64_t c = code[state];        // 低 len[state] bits 有效
            unsigned l = len[state];

            // 把当前码字的 bits 追加到 bitbuf 的高位（但表示成 LSB-first）
            // bitbuf = [低位是先来的 bit，逐渐往高位堆]
            bitbuf |= (c << nbits);
            nbits  += l;

            // 每凑够 8 bit，输出一个字节
            while (nbits >= 8) {
                unsigned char byte = static_cast<unsigned char>(bitbuf & 0xFFu);
                // 反转 bit 顺序，让存入字节是 MSB-first
                *p++ = reverse8(byte);
                bitbuf >>= 8;
                nbits  -= 8;
            }
        }

        // 剩余不到 8 bit 的部分，也输出一个字节
        if (nbits > 0) {
            unsigned char byte = static_cast<unsigned char>(bitbuf & 0xFFu);
            *p++ = reverse8(byte);
        }

        size_t outSize = static_cast<size_t>(p - out_begin);
        write(outSize, bytes);      // 保留原来的长度写入方式
        bytes += outSize;           // 移动到 bitstream 末尾
        return outSize;
    }

    void postprocess_encode() override { SZ_FreeHuffman(); }

    void preprocess_decode() override {}

    // perform decoding
    /*std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
        node t = treeRoot;
        std::vector<T> out(targetLength);
        size_t i = 0, byteIndex = 0, count = 0;
        int r;
        node n = treeRoot;
        size_t encodedLength = 0;
        read(encodedLength, bytes);
        if (n->t)  // root->t==1 means that all state values are the same (constant)
        {
            for (count = 0; count < targetLength; count++) out[count] = n->c + offset;
            return out;
        }

        for (i = 0; count < targetLength; i++) {
            byteIndex = i >> 3;  // i/8
            r = i % 8;
            if (((bytes[byteIndex] >> (7 - r)) & 0x01) == 0)
                n = n->left;
            else
                n = n->right;

            if (n->t) {
                out[count] = n->c + offset;
                n = t;
                count++;
            }
        }
        bytes += encodedLength;
        return out;
    }*/

    std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
        node root = treeRoot;
        std::vector<T> out(targetLength);
        size_t count = 0;

        // 先读出 bitstream 的字节长度
        size_t encodedLength = 0;
        read(encodedLength, bytes);

        node n = root;
        if (n->t) {  // root->t==1 表示只有一个符号（常数）
            T val = n->c + offset;
            for (size_t i = 0; i < targetLength; ++i) {
                out[i] = val;
            }
            // 注意：这里仍需把 bytes 向前挪 encodedLength 个字节，
            //       虽然其实不会用到这些 bit
            bytes += encodedLength;
            return out;
        }

        const uchar *bitptr = bytes;     // 指向 bitstream 起始位置
        size_t bytePos = 0;
        size_t byteLimit = encodedLength;

        if (byteLimit == 0) {
            // 理论上不应该出现：有 targetLength>0 但 encodedLength=0
            return out;
        }

        uchar curByte = bitptr[0];
        int bitsLeft = 8;

        // 遍历 bit，直到解出 targetLength 个符号或者没有更多 bit
        while (count < targetLength && bytePos < byteLimit) {
            // 取当前字节的最高 bit（MSB-first）
            int bit = (curByte & 0x80u) != 0;
            curByte <<= 1;
            --bitsLeft;

            // 走 Huffman 树
            n = bit ? n->right : n->left;

            if (n->t) {
                // 叶子：输出符号并回到根
                out[count++] = n->c + offset;
                n = root;
            }

            // 当前字节 bit 用完了，读取下一个字节
            if (bitsLeft == 0) {
                ++bytePos;
                if (bytePos >= byteLimit) {
                    break;  // 没有更多 bit 了
                }
                curByte = bitptr[bytePos];
                bitsLeft = 8;
            }
        }

        // 消费掉 encodedLength 个字节
        bytes += encodedLength;
        return out;
    }

    // empty function
    void postprocess_decode() override { SZ_FreeHuffman(); }

    // load Huffman tree
    void load(const uchar *&c, size_t &remaining_length) override {
        read(offset, c, remaining_length);
        nodeCount = bytesToInt32_bigEndian(c);
        int stateNum = bytesToInt32_bigEndian(c + sizeof(int)) * 2;
        size_t encodeStartIndex;
        if (nodeCount <= 256)
            encodeStartIndex = 1 + 3 * nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        else if (nodeCount <= 65536)
            encodeStartIndex =
                1 + 2 * nodeCount * sizeof(unsigned short) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        else
            encodeStartIndex =
                1 + 2 * nodeCount * sizeof(unsigned int) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);

        huffmanTree = createHuffmanTree(stateNum);
        treeRoot = reconstruct_HuffTree_from_bytes_anyStates(c + sizeof(int) + sizeof(int), nodeCount);
        c += sizeof(int) + sizeof(int) + encodeStartIndex;
        loaded = true;
    }

    bool isLoaded() const { return loaded; }

   private:
    HuffmanTree *huffmanTree = nullptr;
    node treeRoot;
    unsigned int nodeCount = 0;
    uchar sysEndianType;  // 0: little endian, 1: big endian
    bool loaded = false;
    T offset;

    static inline unsigned char reverse8(unsigned char x) {
        x = (unsigned char)((x & 0xF0u) >> 4 | (x & 0x0Fu) << 4);
        x = (unsigned char)((x & 0xCCu) >> 2 | (x & 0x33u) << 2);
        x = (unsigned char)((x & 0xAAu) >> 1 | (x & 0x55u) << 1);
        return x;
    }


    node reconstruct_HuffTree_from_bytes_anyStates(const unsigned char *bytes, uint nodeCount) {
        if (nodeCount <= 256) {
            unsigned char *L = static_cast<unsigned char *>(malloc(nodeCount * sizeof(unsigned char)));
            memset(L, 0, nodeCount * sizeof(unsigned char));
            unsigned char *R = static_cast<unsigned char *>(malloc(nodeCount * sizeof(unsigned char)));
            memset(R, 0, nodeCount * sizeof(unsigned char));
            T *C = static_cast<T *>(malloc(nodeCount * sizeof(T)));
            memset(C, 0, nodeCount * sizeof(T));
            unsigned char *t = static_cast<unsigned char *>(malloc(nodeCount * sizeof(unsigned char)));
            memset(t, 0, nodeCount * sizeof(unsigned char));
            // TODO: Endian type
            // unsigned char cmpSysEndianType = bytes[0];
            // if(cmpSysEndianType!=(unsigned char)sysEndianType)
            // {
            // 	unsigned char* p = (unsigned char*)(bytes+1+2*nodeCount*sizeof(unsigned char));
            // 	size_t i = 0, size = nodeCount*sizeof(unsigned int);
            // 	while(1)
            // 	{
            // 		symTransform_4bytes(p);
            // 		i+=sizeof(unsigned int);
            // 		if(i<size)
            // 			p+=sizeof(unsigned int);
            // 		else
            // 			break;
            // 	}
            // }
            memcpy(L, bytes + 1, nodeCount * sizeof(unsigned char));
            memcpy(R, bytes + 1 + nodeCount * sizeof(unsigned char), nodeCount * sizeof(unsigned char));
            memcpy(C, bytes + 1 + 2 * nodeCount * sizeof(unsigned char), nodeCount * sizeof(T));
            memcpy(t, bytes + 1 + 2 * nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T),
                   nodeCount * sizeof(unsigned char));
            node root = this->new_node2(C[0], t[0]);
            this->unpad_tree<uchar>(L, R, C, t, 0, root);
            free(L);
            free(R);
            free(C);
            free(t);
            return root;
        } else if (nodeCount <= 65536) {
            unsigned short *L = static_cast<unsigned short *>(malloc(nodeCount * sizeof(unsigned short)));
            memset(L, 0, nodeCount * sizeof(unsigned short));
            unsigned short *R = static_cast<unsigned short *>(malloc(nodeCount * sizeof(unsigned short)));
            memset(R, 0, nodeCount * sizeof(unsigned short));
            T *C = static_cast<T *>(malloc(nodeCount * sizeof(T)));
            memset(C, 0, nodeCount * sizeof(T));
            unsigned char *t = static_cast<unsigned char *>(malloc(nodeCount * sizeof(unsigned char)));
            memset(t, 0, nodeCount * sizeof(unsigned char));

            // TODO: Endian type
            // unsigned char cmpSysEndianType = bytes[0];
            // if(cmpSysEndianType!=(unsigned char)sysEndianType)
            // {
            // 	unsigned char* p = (unsigned char*)(bytes+1);
            // 	size_t i = 0, size = 3*nodeCount*sizeof(unsigned int);
            // 	while(1)
            // 	{
            // 		symTransform_4bytes(p);
            // 		i+=sizeof(unsigned int);
            // 		if(i<size)
            // 			p+=sizeof(unsigned int);
            // 		else
            // 			break;
            // 	}
            // }

            memcpy(L, bytes + 1, nodeCount * sizeof(unsigned short));
            memcpy(R, bytes + 1 + nodeCount * sizeof(unsigned short), nodeCount * sizeof(unsigned short));
            memcpy(C, bytes + 1 + 2 * nodeCount * sizeof(unsigned short), nodeCount * sizeof(T));

            memcpy(t, bytes + 1 + 2 * nodeCount * sizeof(unsigned short) + nodeCount * sizeof(T),
                   nodeCount * sizeof(unsigned char));

            node root = this->new_node2(0, 0);
            this->unpad_tree<unsigned short>(L, R, C, t, 0, root);
            free(L);
            free(R);
            free(C);
            free(t);
            return root;
        } else  // nodeCount>65536
        {
            unsigned int *L = static_cast<unsigned int *>(malloc(nodeCount * sizeof(unsigned int)));
            memset(L, 0, nodeCount * sizeof(unsigned int));
            unsigned int *R = static_cast<unsigned int *>(malloc(nodeCount * sizeof(unsigned int)));
            memset(R, 0, nodeCount * sizeof(unsigned int));
            T *C = static_cast<T *>(malloc(nodeCount * sizeof(T)));
            memset(C, 0, nodeCount * sizeof(T));
            unsigned char *t = static_cast<unsigned char *>(malloc(nodeCount * sizeof(unsigned char)));
            memset(t, 0, nodeCount * sizeof(unsigned char));
            // TODO: Endian type
            // unsigned char cmpSysEndianType = bytes[0];
            // if(cmpSysEndianType!=(unsigned char)sysEndianType)
            // {
            // 	unsigned char* p = (unsigned char*)(bytes+1);
            // 	size_t i = 0, size = 3*nodeCount*sizeof(unsigned int);
            // 	while(1)
            // 	{
            // 		symTransform_4bytes(p);
            // 		i+=sizeof(unsigned int);
            // 		if(i<size)
            // 			p+=sizeof(unsigned int);
            // 		else
            // 			break;
            // 	}
            // }

            memcpy(L, bytes + 1, nodeCount * sizeof(unsigned int));
            memcpy(R, bytes + 1 + nodeCount * sizeof(unsigned int), nodeCount * sizeof(unsigned int));
            memcpy(C, bytes + 1 + 2 * nodeCount * sizeof(unsigned int), nodeCount * sizeof(T));

            memcpy(t, bytes + 1 + 2 * nodeCount * sizeof(unsigned int) + nodeCount * sizeof(T),
                   nodeCount * sizeof(unsigned char));

            node root = this->new_node2(0, 0);
            this->unpad_tree<unsigned int>(L, R, C, t, 0, root);
            free(L);
            free(R);
            free(C);
            free(t);
            return root;
        }
    }

    node new_node(size_t freq, T c, node a, node b) {
        node n = huffmanTree->pool + huffmanTree->n_nodes++;
        if (freq) {
            n->c = c;
            n->freq = freq;
            n->t = 1;
            // printf("new_node: c = %d, freq = %zu, t = %d \n", n->c, n->freq, n->t);
        } else {
            n->left = a;
            n->right = b;
            n->freq = a->freq + b->freq;
            n->t = 0;
            // printf("new_node: c = %d, freq = %zu, t = %d, left = %d, right = %d \n", n->c, n->freq, n->t, n->left->c,
            // n->right->c);
            // n->c = 0;
        }
        return n;
    }

    node new_node2(T c, unsigned char t) {
        huffmanTree->pool[huffmanTree->n_nodes].c = c;
        huffmanTree->pool[huffmanTree->n_nodes].t = t;
        return huffmanTree->pool + huffmanTree->n_nodes++;
    }

    /* priority queue */
    void qinsert(node n) {
        int j, i = huffmanTree->qend++;
        while ((j = (i >> 1)))  // j=i/2
        {
            if (huffmanTree->qq[j]->freq <= n->freq) break;
            huffmanTree->qq[i] = huffmanTree->qq[j], i = j;
        }
        huffmanTree->qq[i] = n;
    }

    node qremove() {
        int i = 1, l;
        node n = huffmanTree->qq[i = 1];
        node p;
        if (huffmanTree->qend < 2) return nullptr;
        huffmanTree->qend--;
        huffmanTree->qq[i] = huffmanTree->qq[huffmanTree->qend];

        while ((l = (i << 1)) < huffmanTree->qend) {  // l=(i*2)
            if (l + 1 < huffmanTree->qend && huffmanTree->qq[l + 1]->freq < huffmanTree->qq[l]->freq) l++;
            if (huffmanTree->qq[i]->freq > huffmanTree->qq[l]->freq) {
                p = huffmanTree->qq[i];
                huffmanTree->qq[i] = huffmanTree->qq[l];
                huffmanTree->qq[l] = p;
                i = l;
            } else {
                break;
            }
        }
        return n;
    }

    /* walk the tree and put 0s and 1s */
    /**
     * @out1 should be set to 0.
     * @out2 should be 0 as well.
     * @index: the index of the byte
     * */
    void build_code(node n, uint64_t code_val, int len) {
        if (n->t) {
            assert(len <= 64);
            huffmanTree->code[n->c] = code_val;                  // 低 len bit 有效
            huffmanTree->cout[n->c] = static_cast<unsigned char>(len);
            if (len > huffmanTree->maxBitCount) {
                huffmanTree->maxBitCount = len;
            }
            return;
        }

        build_code(n->left, code_val, len + 1);

        uint64_t code_right = code_val | (uint64_t(1) << len);
        build_code(n->right, code_right, len + 1);
    }
    /**
     * Compute the frequency of the data and build the Huffman tree
     * @param HuffmanTree* huffmanTree (output)
     * @param int *s (input)
     * @param size_t length (input)
     * */
    void init(const T *s, size_t length) {

   



        T max = s[0];
        offset = 0;  // offset is min
         const size_t ui16_range = 1u << 16;
         std::vector<size_t> frequencyList(ui16_range, 0);
         auto frenqencies = frequencyList.data();



        #ifdef _OPENMP

            auto default_nthreads = omp_get_max_threads();
            if (default_nthreads > 1 && length >= 1u << 18) {
                auto best_num_threads = std::min(default_nthreads, (int)(length / ui16_range));
                omp_set_num_threads(best_num_threads);
                #pragma omp parallel
                {

                    //int tid = omp_get_thread_num();
                    
                    // 每个线程一个局部 freq
                    std::vector<size_t> local_freq(ui16_range, 0);

                    #pragma omp for
                    for (long long i = 0; i < (long long)length; ++i) {
                        auto v = s[i];
                        // 假设 v 已经在 [0, ui16_range)
                        local_freq[(unsigned)v]++;
                    }

                    // 归并到全局
                   
                        for (size_t k = 0; k < ui16_range; ++k) {
                            #pragma omp atomic
                            frenqencies[k] += local_freq[k];
                        }

                }
                omp_set_num_threads(default_nthreads);
            }
            else{
                for (size_t i = 0; i < length; ++i) {
                    ++frenqencies[s[i]];
                }
            }

        #else
            for (size_t i = 0; i < length; ++i) {
                ++frenqencies[s[i]];
            }
        #endif


/*
#if (SZ3_USE_SKA_HASH) && (INTPTR_MAX == INT64_MAX)  // use ska for 64bit system
        ska::unordered_map<T, size_t> frequency;
#else   // most likely 32bit system
        std::unordered_map<T, size_t> frequency;
#endif  // INTPTR_MAX == INT64_MAX
*/

/*
        for (const auto &kv : frequency) {
            auto k = kv.first;
            if (k > max) {
                max = k;
            }
            if (k < offset) {
                offset = k;
            }
        }
*/  
        //Timer timer(true);

        for (int i = 0; i < ui16_range; ++i) {
            if (frenqencies[i] != 0) {
                offset = i;
                break;
            }

        }
        for (int i = ui16_range - 1; i >= 0 ; --i) {
            if (frenqencies[i] != 0) {
                max = i;
                break;
            }

        }
       // std::cout<<offset<<" "<<max<<std::endl;


        int stateNum = max - offset + 2;
        //timer.stop("count");
        huffmanTree = createHuffmanTree(stateNum);
        // to produce the same huffman three on linux & win, we need to iterate through ordered_map in a fixed order
        
        for (int i = offset; i <= max; ++i) {
            if (frenqencies[i] != 0) {
                qinsert(new_node(frenqencies[i], i - offset, nullptr, nullptr));
            }
        }
        // for (const auto &f : frequency) {
        //     qinsert(new_node(f.second, f.first - offset, nullptr, nullptr));
        // }

        while (huffmanTree->qend > 2) {
            auto left = qremove();
            auto right = qremove();
            qinsert(new_node(0, 0, left, right));
        }


        build_code(huffmanTree->qq[1], 0ULL, 0);
        treeRoot = huffmanTree->qq[1];
    }

    template <class T1>
    void pad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root) {
        C[i] = root->c;
        t[i] = root->t;
        node lroot = root->left;
        if (lroot != nullptr) {
            huffmanTree->n_inode++;
            L[i] = huffmanTree->n_inode;
            pad_tree(L, R, C, t, huffmanTree->n_inode, lroot);
        }
        node rroot = root->right;
        if (rroot != nullptr) {
            huffmanTree->n_inode++;
            R[i] = huffmanTree->n_inode;
            pad_tree(L, R, C, t, huffmanTree->n_inode, rroot);
        }
    }

    template <class T1>
    void unpad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root) {
        // root->c = C[i];
        if (root->t == 0) {
            T1 l, r;
            l = L[i];
            if (l != 0) {
                node lroot = new_node2(C[l], t[l]);
                root->left = lroot;
                unpad_tree(L, R, C, t, l, lroot);
            }
            r = R[i];
            if (r != 0) {
                node rroot = new_node2(C[r], t[r]);
                root->right = rroot;
                unpad_tree(L, R, C, t, r, rroot);
            }
        }
    }

    template <class T1>
    unsigned int convert_HuffTree_to_bytes_anyStates(unsigned int nodeCount, unsigned char *out) {
        T1 *L = static_cast<T1 *>(malloc(nodeCount * sizeof(T1)));
        memset(L, 0, nodeCount * sizeof(T1));
        T1 *R = static_cast<T1 *>(malloc(nodeCount * sizeof(T1)));
        memset(R, 0, nodeCount * sizeof(T1));
        T *C = static_cast<T *>(malloc(nodeCount * sizeof(T)));
        memset(C, 0, nodeCount * sizeof(T));
        unsigned char *t = static_cast<unsigned char *>(malloc(nodeCount * sizeof(unsigned char)));
        memset(t, 0, nodeCount * sizeof(unsigned char));

        pad_tree(L, R, C, t, 0, huffmanTree->qq[1]);

        unsigned int totalSize =
            1 + 2 * nodeCount * sizeof(T1) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        //*out = (unsigned char*)malloc(totalSize);
        out[0] = sysEndianType;
        memcpy(out + 1, L, nodeCount * sizeof(T1));
        memcpy(out + 1 + nodeCount * sizeof(T1), R, nodeCount * sizeof(T1));
        memcpy(out + 1 + 2 * nodeCount * sizeof(T1), C, nodeCount * sizeof(T));
        memcpy(out + 1 + 2 * nodeCount * sizeof(T1) + nodeCount * sizeof(T), t, nodeCount * sizeof(unsigned char));

        free(L);
        free(R);
        free(C);
        free(t);
        return totalSize;
    }

    void SZ_FreeHuffman() {
    if (huffmanTree != nullptr) {
        free(huffmanTree->pool);
        huffmanTree->pool = nullptr;

        free(huffmanTree->qqq);
        huffmanTree->qqq = nullptr;

        if (huffmanTree->code != nullptr) {
            free(huffmanTree->code);
            huffmanTree->code = nullptr;
        }
        if (huffmanTree->cout != nullptr) {
            free(huffmanTree->cout);
            huffmanTree->cout = nullptr;
        }

        free(huffmanTree);
        huffmanTree = nullptr;
    }
}
};
}  // namespace SZ3

#endif
