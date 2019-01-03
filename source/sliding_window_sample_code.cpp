/*****************************************************************************
 *  Orbbec
 *  Copyright (C) 2019 by ORBBEC Technology., Inc.
 *
 *  This file is part of Orbbec DepthEngine.
 *
 *  This file belongs to ORBBEC Technology., Inc.
 *  It is considered a trade secret, and is not to be divulged or used by
 * parties who have NOT received written authorization from the owner.
 *
 *  Description:
 *
 *  Author: liuxianzhuo@orbbec.com
 *
 ****************************************************************************/


typedef struct _SlidingWindowData
{
    int rbuf_size;                 // Ring buffer size
    ring_buf_t *horiz_cost_rbuf;   // horizontal cost ring buffer

    int line_sad_align_offset;     // 行 SAD 对齐偏移（两边各自的扩展量）
    int line_sad_step;             // 行 SAD 的步长，为了方便求亚像素做了扩展，步长为 disp_range_ + 2*line_sad_align_offset
    int line_sad_nelem;            // 行 SAD 总的大小
    uint16_t *curr_line_sad;       // 当前行的sad，维度为[width][disp_range + 2*line_sad_align_offset]
} SlidingWindowData;


/**
 * 这个函数只摘抄了部分
 * @param helper
 * @return
 */
int DispSlidingWindow::AllocHorizSADRingBuf(SlidingWindowData *helper) const
{
    // 分配水平方向SAD循环缓冲区
    helper->rbuf_size = block_size_;
    helper->horiz_cost_rbuf = ring_buf_new(helper->rbuf_size);

    HorizWinRangeSAD *h_win_sad = NULL;
    for (int i = 0; i < helper->rbuf_size; ++i)
    {
        h_win_sad = (HorizWinRangeSAD *) malloc(sizeof(HorizWinRangeSAD));
        h_win_sad->disp_range = disp_range_;
        h_win_sad->width = width_;
        h_win_sad->sad_buffer = (uint16_t *) malloc(
            sizeof(h_win_sad->sad_buffer[0]) * h_win_sad->disp_range * h_win_sad->width);
        memset(h_win_sad->sad_buffer, 0, sizeof(h_win_sad->sad_buffer[0]) * h_win_sad->disp_range * h_win_sad->width);
        ring_buf_enque(helper->horiz_cost_rbuf, (void *) h_win_sad);
    }

    // 分配 SAD buffer
    helper->line_sad_align_offset = 1;
    helper->line_sad_step = disp_range_ + 2 * helper->line_sad_align_offset;
    helper->line_sad_nelem = helper->line_sad_step * width_;
    helper->curr_line_sad = (uint16_t *) malloc(sizeof(helper->curr_line_sad[0]) * helper->line_sad_nelem);

    return 0;
}


int FreeHorizSADRingBuf(SlidingWindowData *helper)
{
    assert(helper);

    HorizWinRangeSAD *h_win_sad = NULL;
    while (!ring_buf_empty(helper->horiz_cost_rbuf))
    {
        // 逐个释放ring_buffer内的元素
        h_win_sad = (HorizWinRangeSAD *) ring_buf_deque(helper->horiz_cost_rbuf);
        free(h_win_sad->sad_buffer);
        free(h_win_sad);
    }
    ring_buf_free(helper->horiz_cost_rbuf);


    if (helper->curr_line_sad)
    {
        free(helper->curr_line_sad);
        helper->curr_line_sad = NULL;
    }

    return 0;
}

/**
 * 更新横向窗口 SAD
 * @param prev_sad          前一列的SAD
 * @param curr_sad          当前列的SAD
 * @param disp_range        视差范围
 * @param tgt_sub           刚被移出窗口的target pixel
 * @param ref_sub           刚被移出窗口的reference pixels
 * @param tgt_add           刚被移入窗口的target pixel
 * @param ref_add           刚被移入窗口的reference pixels
 */
void UpdateHorizWinRangeSAD(const uint16_t *prev_sad, uint16_t *curr_sad, int disp_range,
                            uint8_t tgt_sub, const uint8_t *ref_sub,
                            uint8_t tgt_add, const uint8_t *ref_add)
{
    int d;

#if defined(__AVX2__)
    __m256i tgt_sub_vec = _mm256_set1_epi8(tgt_sub);    // broadcast
    __m256i tgt_add_vec = _mm256_set1_epi8(tgt_add);
    __m256i x1, y1, sub_vec1, add_vec1;
    __m256i x2, y2, sub_vec2, add_vec2;

    __m256i src1, src2, src3, src4;
    __m256i tmp1_low, tmp1_high, dst1_low, dst1_high;
    __m256i tmp2_low, tmp2_high, dst2_low, dst2_high;
    __m256i zero = _mm256_setzero_si256();
    __m256i x1_low_u16, x1_high_u16, y1_low_u16, y1_high_u16;
    __m256i x2_low_u16, x2_high_u16, y2_low_u16, y2_high_u16;

    for (d = 0; d < disp_range / 64; ++d)
    {
        //////////////////////////////////////////////////////////////////////////////
        // round 1
        x1 = _mm256_loadu_si256((__m256i *) (ref_sub + d * 64));
        y1 = _mm256_loadu_si256((__m256i *) (ref_add + d * 64));
        sub_vec1 = _mm256_abs_epi8(_mm256_sub_epi8(tgt_sub_vec, x1));
        add_vec1 = _mm256_abs_epi8(_mm256_sub_epi8(tgt_add_vec, y1));

        // 交换中间两个64位，方面调用_mm256_unpacklo_epi8/_mm256_unpackhi_epi8
        sub_vec1 = _mm256_permute4x64_epi64(sub_vec1, 0xD8);
        add_vec1 = _mm256_permute4x64_epi64(add_vec1, 0xD8);

        src1 = _mm256_loadu_si256((__m256i *) (prev_sad + d * 64));
        x1_low_u16 = _mm256_unpacklo_epi8(sub_vec1, zero);
        y1_low_u16 = _mm256_unpacklo_epi8(add_vec1, zero);
        tmp1_low = _mm256_sub_epi16(src1, x1_low_u16);
        dst1_low = _mm256_add_epi16(tmp1_low, y1_low_u16);
        _mm256_storeu_si256((__m256i *) (curr_sad + d * 64), dst1_low);

        src2 = _mm256_loadu_si256((__m256i *) (prev_sad + d * 64 + 16));
        x1_high_u16 = _mm256_unpackhi_epi8(sub_vec1, zero);
        y1_high_u16 = _mm256_unpackhi_epi8(add_vec1, zero);
        tmp1_high = _mm256_sub_epi16(src2, x1_high_u16);
        dst1_high = _mm256_add_epi16(tmp1_high, y1_high_u16);
        _mm256_storeu_si256((__m256i *) (curr_sad + d * 64 + 16), dst1_high);

        //////////////////////////////////////////////////////////////////////////////
        // round 2
        x2 = _mm256_loadu_si256((__m256i *) (ref_sub + d * 64 + 32));
        y2 = _mm256_loadu_si256((__m256i *) (ref_add + d * 64 + 32));
        sub_vec2 = _mm256_abs_epi8(_mm256_sub_epi8(tgt_sub_vec, x2));
        add_vec2 = _mm256_abs_epi8(_mm256_sub_epi8(tgt_add_vec, y2));

        // 交换中间两个64位，方面调用_mm256_unpacklo_epi8/_mm256_unpackhi_epi8
        sub_vec2 = _mm256_permute4x64_epi64(sub_vec2, 0xD8);
        add_vec2 = _mm256_permute4x64_epi64(add_vec2, 0xD8);

        src3 = _mm256_loadu_si256((__m256i *) (prev_sad + d * 64 + 32));
        x2_low_u16 = _mm256_unpacklo_epi8(sub_vec2, zero);
        y2_low_u16 = _mm256_unpacklo_epi8(add_vec2, zero);
        tmp2_low = _mm256_sub_epi16(src3, x2_low_u16);
        dst2_low = _mm256_add_epi16(tmp2_low, y2_low_u16);
        _mm256_storeu_si256((__m256i *) (curr_sad + d * 64 + 32), dst2_low);

        src4 = _mm256_loadu_si256((__m256i *) (prev_sad + d * 64 + 48));
        x2_high_u16 = _mm256_unpackhi_epi8(sub_vec2, zero);
        y2_high_u16 = _mm256_unpackhi_epi8(add_vec2, zero);
        tmp2_high = _mm256_sub_epi16(src4, x2_high_u16);
        dst2_high = _mm256_add_epi16(tmp2_high, y2_high_u16);
        _mm256_storeu_si256((__m256i *) (curr_sad + d * 64 + 48), dst2_high);
    }
#else
    // TODO: 改成向量化的方式
    for (d = 0; d < disp_range; ++d)
    {
        curr_sad[d] = prev_sad[d] - (uint16_t) abs(tgt_sub - ref_sub[d]) + (uint16_t) abs(tgt_add - ref_add[d]);
    }
#endif
}

/**
 * 向量减
 * @param minuend       被减数
 * @param subtrahend    减数
 * @param len           向量长度
 */
void VectorSubUInt16(uint16_t *minuend, const uint16_t *subtrahend, int16_t len)
{
    int i;
#if defined(__AVX2__)
    __m256i x0, x1, x2, x3;
    __m256i y0, y1, y2, y3;
    __m256i z0, z1, z2, z3;

    for (i = 0; i < len / 64; ++i)
    {
        x0 = _mm256_loadu_si256((__m256i *) (minuend + i * 64));
        y0 = _mm256_loadu_si256((__m256i *) (subtrahend + i * 64));
        z0 = _mm256_sub_epi16(x0, y0);
        _mm256_storeu_si256((__m256i *) (minuend + i * 64), z0);

        x1 = _mm256_loadu_si256((__m256i *) (minuend + i * 64 + 16));
        y1 = _mm256_loadu_si256((__m256i *) (subtrahend + i * 64 + 16));
        z1 = _mm256_sub_epi16(x1, y1);
        _mm256_storeu_si256((__m256i *) (minuend + i * 64 + 16), z1);

        x2 = _mm256_loadu_si256((__m256i *) (minuend + i * 64 + 32));
        y2 = _mm256_loadu_si256((__m256i *) (subtrahend + i * 64 + 32));
        z2 = _mm256_sub_epi16(x2, y2);
        _mm256_storeu_si256((__m256i *) (minuend + i * 64 + 32), z2);

        x3 = _mm256_loadu_si256((__m256i *) (minuend + i * 64 + 48));
        y3 = _mm256_loadu_si256((__m256i *) (subtrahend + i * 64 + 48));
        z3 = _mm256_sub_epi16(x3, y3);
        _mm256_storeu_si256((__m256i *) (minuend + i * 64 + 48), z3);
    }
#elif defined(__SSE__) || _M_IX86_FP >= 1
    __m128i x, y, z;

    for (i = 0; i < len/8; ++i)
    {
        x = _mm_loadu_si128((__m128i*)(minuend + i * 8));
        y = _mm_loadu_si128((__m128i*)(subtrahend + i * 8));
        z = _mm_sub_epi16(x, y);
        _mm_storeu_si128((__m128i*)(minuend + i * 8), z);
    }
#else
    // TODO: 改成向量化的方式
    for (i = 0; i < len; ++i)
        minuend[i] -= subtrahend[i];
#endif
}


/**
 * 向量加，并求最小值
 * @param[in]  src       源向量
 * @param[in]  tgt       目标向量
 * @param[in]  num       向量长度
 * @param[out] min_idx   最小值的索引（从0开始）
 * @param[out] min_val   最小值
 */
void VectorAddFindMinUInt16(const uint16_t *src, uint16_t *tgt, int16_t num, int16_t *min_idx, uint16_t *min_val)
{
    int16_t i;

#if defined(__AVX2__)

    __m256i x0, x1, x2, x3;
    __m256i y0, y1, y2, y3;
    __m256i z0, z1, z2, z3;
    __m256i vmin0, vmin1, vmin2, vmin3;
    __m256i vcmp0, vcmp1, vcmp2, vcmp3;
    uint32_t mask[4];
#ifdef _MSC_VER
    unsigned long idx[4];	// _BitScanForward(unsigned long * Index, unsigned long Mask)
#else
    uint8_t idx[4];
#endif
    uint16_t minval[4];

    int16_t min_idx_ = 0;
    uint16_t min_val_ = UINT16_MAX;
    // 16 × 16 bits = 256 bits，展开4级
    for (i = 0; i < num / 64; ++i)
    {
        x0 = _mm256_loadu_si256((__m256i *) (tgt + i * 64));
        y0 = _mm256_loadu_si256((__m256i *) (src + i * 64));
        z0 = _mm256_add_epi16(x0, y0);
        _mm256_storeu_si256((__m256i *) (tgt + i * 64), z0);

        x1 = _mm256_loadu_si256((__m256i *) (tgt + i * 64 + 16));
        y1 = _mm256_loadu_si256((__m256i *) (src + i * 64 + 16));
        z1 = _mm256_add_epi16(x1, y1);
        _mm256_storeu_si256((__m256i *) (tgt + i * 64 + 16), z1);

        x2 = _mm256_loadu_si256((__m256i *) (tgt + i * 64 + 32));
        y2 = _mm256_loadu_si256((__m256i *) (src + i * 64 + 32));
        z2 = _mm256_add_epi16(x2, y2);
        _mm256_storeu_si256((__m256i *) (tgt + i * 64 + 32), z2);

        x3 = _mm256_loadu_si256((__m256i *) (tgt + i * 64 + 48));
        y3 = _mm256_loadu_si256((__m256i *) (src + i * 64 + 48));
        z3 = _mm256_add_epi16(x3, y3);
        _mm256_storeu_si256((__m256i *) (tgt + i * 64 + 48), z3);

        // 求最小值
        vmin0 = _mm256_min_epu16(z0, _mm256_alignr_epi8(z0, z0, 2));
        vmin1 = _mm256_min_epu16(z1, _mm256_alignr_epi8(z1, z1, 2));
        vmin2 = _mm256_min_epu16(z2, _mm256_alignr_epi8(z2, z2, 2));
        vmin3 = _mm256_min_epu16(z3, _mm256_alignr_epi8(z3, z3, 2));

        vmin0 = _mm256_min_epu16(vmin0, _mm256_alignr_epi8(vmin0, vmin0, 4));
        vmin1 = _mm256_min_epu16(vmin1, _mm256_alignr_epi8(vmin1, vmin1, 4));
        vmin2 = _mm256_min_epu16(vmin2, _mm256_alignr_epi8(vmin2, vmin2, 4));
        vmin3 = _mm256_min_epu16(vmin3, _mm256_alignr_epi8(vmin3, vmin3, 4));

        vmin0 = _mm256_min_epu16(vmin0, _mm256_alignr_epi8(vmin0, vmin0, 8));
        vmin1 = _mm256_min_epu16(vmin1, _mm256_alignr_epi8(vmin1, vmin1, 8));
        vmin2 = _mm256_min_epu16(vmin2, _mm256_alignr_epi8(vmin2, vmin2, 8));
        vmin3 = _mm256_min_epu16(vmin3, _mm256_alignr_epi8(vmin3, vmin3, 8));

        vmin0 = _mm256_min_epu16(vmin0, _mm256_permute2x128_si256(vmin0, vmin0, 0x01));
        vmin1 = _mm256_min_epu16(vmin1, _mm256_permute2x128_si256(vmin1, vmin1, 0x01));
        vmin2 = _mm256_min_epu16(vmin2, _mm256_permute2x128_si256(vmin2, vmin2, 0x01));
        vmin3 = _mm256_min_epu16(vmin3, _mm256_permute2x128_si256(vmin3, vmin3, 0x01));

        vcmp0 = _mm256_cmpeq_epi16(z0, vmin0);  // 比较是否相等
        vcmp1 = _mm256_cmpeq_epi16(z1, vmin1);
        vcmp2 = _mm256_cmpeq_epi16(z2, vmin2);
        vcmp3 = _mm256_cmpeq_epi16(z3, vmin3);

        mask[0] = (uint32_t) _mm256_movemask_epi8(vcmp0); // 取每个8bit的最高有效位，组成一个mask
        mask[1] = (uint32_t) _mm256_movemask_epi8(vcmp1);
        mask[2] = (uint32_t) _mm256_movemask_epi8(vcmp2);
        mask[3] = (uint32_t) _mm256_movemask_epi8(vcmp3);

#if defined(__GNUC__)
        idx[0] = __builtin_ctz(mask[0]) >> 1; // 每两个字节
        idx[1] = __builtin_ctz(mask[1]) >> 1;
        idx[2] = __builtin_ctz(mask[2]) >> 1;
        idx[3] = __builtin_ctz(mask[3]) >> 1;
#elif defined(_MSC_VER)
        // 一定会有一个非0，所以这里不检查_BitScanForward()返回值
        _BitScanForward(idx, mask[0]);
        idx[0] >>= 1;
        _BitScanForward(idx + 1, mask[1]);
        idx[1] >>= 1;
        _BitScanForward(idx + 2, mask[2]);
        idx[2] >>= 1;
        _BitScanForward(idx + 3, mask[3]);
        idx[3] >>= 1;
#endif

        minval[0] = (uint16_t) _mm256_extract_epi16(vmin0, 0);  // 最后全是相同的值
        if (minval[0] < min_val_)
        {
            min_val_ = minval[0];
            min_idx_ = (int16_t) (idx[0] + i * 64);
        }
        minval[1] = (uint16_t) _mm256_extract_epi16(vmin1, 0);
        if (minval[1] < min_val_)
        {
            min_val_ = minval[1];
            min_idx_ = (int16_t) (idx[1] + 16 + i * 64);
        }
        minval[2] = (uint16_t) _mm256_extract_epi16(vmin2, 0);
        if (minval[2] < min_val_)
        {
            min_val_ = minval[2];
            min_idx_ = (int16_t) (idx[2] + 32 + i * 64);
        }
        minval[3] = (uint16_t) _mm256_extract_epi16(vmin3, 0);
        if (minval[3] < min_val_)
        {
            min_val_ = minval[3];
            min_idx_ = (int16_t) (idx[3] + 48 + i * 64);
        }
    }

    *min_idx = min_idx_;
    *min_val = min_val_;

#else
    tgt[0] += src[0];
    *min_idx = 0;
    *min_val = tgt[0];
    for (i = 1; i < num; ++i)
    {
        tgt[i] += src[i];
        if (tgt[i] < *min_val)
        {
            *min_idx = i;
            *min_val = tgt[i];
        }
    }
}


/**
 * @brief 更新最终的 disp_range_个SAD，它由 @block_size_ 行、同位置的disp_range_个横向的 @block_size_ 窗口的SAD累加得到
 *      并求最小cost，和最小cost对应的视差
 * @param[in]  hwin_sad      水平方向的sad
 * @param[in]  ext_sad       扩展SAD，左右两边各扩展1
 * @param[in]  disp_range    视差范围
 * @param[out] disp          最小视差
 * @param[out] min_cost      最小代价
 *
 */
void
UpdateVeritalAccumExtSAD(const uint16_t *hwin_sad, uint16_t *ext_sad, int16_t disp_range,
                         int16_t *disp, uint16_t *min_cost)
{
    uint16_t *regular_sad = ext_sad + 1;
    VectorAddFindMinUInt16(hwin_sad, regular_sad, disp_range, disp, min_cost);

    // TODO: 可以进一步做亚像素插值
}


/**
 * 滑窗匹配
 * @param[in]  ext_tgt   (扩展后的)目标图像
 * @param[in]  ext_ref   (扩展后的)参考图像
 * @param[out] disp      (原始分辨率的)视差图
 * @param[out] cost      (原始分辨率的)最小匹配代价图
 * @return
 *
 * @note    这里 @ext_tgt 和 @ext_ref 都是考虑到边界后扩展的，所以函数内部并不处理边界
 *      边界的扩展宽度，由类的成员变量决定
 *
 */
int DispSlidingWindow::SlidingWindowMatch(DEImageU8 ext_tgt, DEImageU8 ext_ref, DEImageS16 disp, DEImageU16 cost)
{
    SlidingWindowData helper;

    AllocHorizSADRingBuf(&helper);

    auto curr_line_sad = helper.curr_line_sad;
    const auto line_sad_nelem = helper.line_sad_nelem;
    const auto line_sad_step = helper.line_sad_step;
    const auto line_sad_align_offset = helper.line_sad_align_offset;

    auto horiz_cost_rbuf = helper.horiz_cost_rbuf;

    const uint8_t *tgt_ptr = ext_tgt.data;
    const uint8_t *ref_ptr = ext_ref.data;
    int16_t *disp_ptr = disp.data;
    uint16_t *cost_ptr = cost.data;

    const uint8_t *blk_tgt_ptr = NULL;
    const uint8_t *blk_ref_ptr = NULL;

    int extend_offset = std::max(0, min_disp_ - half_blk_sz_);

    // 清零SAD向量
    memset(curr_line_sad, 0, sizeof(curr_line_sad[0]) * line_sad_nelem);

    ///////////////////////////////////////////////////////////////////////////
    // 行启动阶段：初始化第一有效行 前面的 block_size - 1 行
    HorizWinRangeSAD *h_win_range_sads = NULL;

    int16_t disp_ignore;
    uint16_t cost_ignore;
    int col;    // 列号
    for (int r = 0; r < block_size_ - 1; ++r)
    {
        col = 0;
        h_win_range_sads = (HorizWinRangeSAD *) ring_buf_deque(horiz_cost_rbuf);
        assert(h_win_range_sads);

        blk_ref_ptr = ref_ptr + extend_offset;
        blk_tgt_ptr = tgt_ptr + extend_offset + max_disp_;

        // 启动行内SAD的计算：最开始的 @sliding_blk_sz_ 个长为 @disp_range_ 的AbsDiff向量相加
        memset(h_win_range_sads->sad_buffer, 0, sizeof(h_win_range_sads->sad_buffer[0]) * disp_range_);
        for (int i = 0; i < block_size_; ++i)
        {
            const uint8_t *ref_disp_range = blk_ref_ptr;
            // TODO: 修改成向量化操作
            for (int d = 0; d < disp_range_; ++d)
            {
                h_win_range_sads->sad_buffer[d] += abs(*blk_tgt_ptr - *ref_disp_range++);
            }

            blk_tgt_ptr++;
            blk_ref_ptr++;
        }
        UpdateVeritalAccumExtSAD(h_win_range_sads->sad_buffer + col * disp_range_,
                                 curr_line_sad + col * line_sad_step, disp_range_,
                                 &disp_ignore, &cost_ignore);

        for (col = 1; col < width_; ++col)
        {
            // 移到下一列时，更新h_win_range_sads
            UpdateHorizWinRangeSAD(h_win_range_sads->sad_buffer + (col - 1) * disp_range_,
                                   h_win_range_sads->sad_buffer + col * disp_range_,
                                   disp_range_,
                                   *(blk_tgt_ptr - block_size_), blk_ref_ptr - block_size_,   // substract
                                   *blk_tgt_ptr, blk_ref_ptr);                                // add
            // 累积到总的SAD中
            UpdateVeritalAccumExtSAD(h_win_range_sads->sad_buffer + col * disp_range_,
                                 curr_line_sad + col * line_sad_step, disp_range_,
                                 &disp_ignore, &cost_ignore);

            blk_ref_ptr++;
            blk_tgt_ptr++;
        }

        ring_buf_enque(horiz_cost_rbuf, h_win_range_sads);

        ref_ptr += ext_ref.width;
        tgt_ptr += ext_tgt.width;
    }

    ///////////////////////////////////////////////////////////////////////////
    // 真正开始计算输出
    for (int r = 0; r < height_; ++r)
    {
        col = 0;
        h_win_range_sads = (HorizWinRangeSAD *) ring_buf_deque(horiz_cost_rbuf);
        assert(h_win_range_sads);

        blk_ref_ptr = ref_ptr + extend_offset;
        blk_tgt_ptr = tgt_ptr + extend_offset + max_disp_;

        ///////////////////////////////////////////////////////////////////////
        // 减去被移出匹配窗口的行，这是与启动阶段的差别
        // curr_line_sad = prev_line_sad - (last) h_win_range_sads
        VectorSubUInt16(curr_line_sad + 0 * line_sad_step + line_sad_align_offset,
                        h_win_range_sads->sad_buffer + 0, disp_range_);
        ///////////////////////////////////////////////////////////////////////

        // 启动行内SAD的计算：最开始的 @sliding_blk_sz_ 个长为 @disp_range_ 的AbsDiff向量相加
        memset(h_win_range_sads->sad_buffer, 0, sizeof(h_win_range_sads->sad_buffer[0]) * disp_range_);
        for (int i = 0; i < block_size_; ++i)
        {
            const uint8_t *ref_disp_range = blk_ref_ptr;
            // TODO: 修改成向量化操作
            for (int d = 0; d < disp_range_; ++d)
            {
                h_win_range_sads->sad_buffer[d] += abs(*blk_tgt_ptr - *ref_disp_range++);
            }

            blk_tgt_ptr++;
            blk_ref_ptr++;
        }

        UpdateVeritalAccumExtSAD(h_win_range_sads->sad_buffer + col * disp_range_,
                                 curr_line_sad + col * line_sad_step, disp_range_,
                                 disp_ptr + col, cost_ptr + col);

        for (col = 1; col < width_; ++col)
        {
            ///////////////////////////////////////////////////////////////////
            // 减去被移出匹配窗口的行，这是与启动阶段的差别
            // curr_line_sad = prev_line_sad - (last) h_win_range_sads
            VectorSubUInt16(curr_line_sad + col * line_sad_step + line_sad_align_offset,
                            h_win_range_sads->sad_buffer + col * disp_range_, disp_range_);
            ///////////////////////////////////////////////////////////////////

            // 移到下一列时，更新h_win_range_sads
            UpdateHorizWinRangeSAD(h_win_range_sads->sad_buffer + (col - 1) * disp_range_,
                                   h_win_range_sads->sad_buffer + col * disp_range_,
                                   disp_range_,
                                   *(blk_tgt_ptr - block_size_), blk_ref_ptr - block_size_,   // substract
                                   *blk_tgt_ptr, blk_ref_ptr);                                // add

            // 累积到总的SAD中，并求最小cost，以及对应的视差
            UpdateVeritalAccumExtSAD(h_win_range_sads->sad_buffer + col * disp_range_,
                                 curr_line_sad + col * line_sad_step, disp_range_,
                                 disp_ptr + col, cost_ptr + col);


            blk_ref_ptr++;
            blk_tgt_ptr++;
        }

        ring_buf_enque(horiz_cost_rbuf, h_win_range_sads);

        // 输入输出图像的指针下移一行
        ref_ptr += ext_ref.width;
        tgt_ptr += ext_tgt.width;
        disp_ptr += disp.width;
        cost_ptr += cost.width;
    }

    FreeHorizSADRingBuf(&helper);

    return 0;
}
