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

#include "sliding_window_sample_code.hpp"


DispSlidingWindow::DispSlidingWindow()
{

}
DispSlidingWindow::~DispSlidingWindow()
{

}
void DispSlidingWindow::setParams(int block_size, int min_disp,
                    int disp_range,int width, int height, int max_disp)
{
    block_size_ = block_size;
    min_disp_ = min_disp;
    half_blk_sz_ = block_size/2;
    disp_range_ = disp_range;
    width_ = width;
    height_ = height;
    max_disp_ = max_disp;
}
/**
 * 这个函数只摘抄了部分
 * @param helper
 * @return
 */
int DispSlidingWindow::AllocHorizSADRingBuf(SlidingWindowData *helper) const
{
    // 分配水平方向SAD循环缓冲区
    helper->rbuf_size = DispSlidingWindow::block_size_;
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
    // TODO: 改成向量化的方式
    for (d = 0; d < disp_range; ++d)
    {
        curr_sad[d] = prev_sad[d] - (uint16_t) abs(tgt_sub - ref_sub[d]) + (uint16_t) abs(tgt_add - ref_add[d]);
    }

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

    // TODO: 改成向量化的方式
    for (i = 0; i < len; ++i)
        minuend[i] -= subtrahend[i];

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

    //int extend_offset = std::min(0, min_disp_ - half_blk_sz_);
    int extend_offset = 129;

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

        blk_ref_ptr = ref_ptr;
        blk_tgt_ptr = tgt_ptr + 127;

        // 启动行内SAD的计算：最开始的 @sliding_blk_sz_ 个长为 @disp_range_ 的AbsDiff向量相加
        memset(h_win_range_sads->sad_buffer, 0, sizeof(h_win_range_sads->sad_buffer[0]) * disp_range_);
        //计算第一个元素的256视差的代价
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
        //更新第一个元素的列
        UpdateVeritalAccumExtSAD(h_win_range_sads->sad_buffer + col * disp_range_,
                                 curr_line_sad + col * line_sad_step, disp_range_,
                                 &disp_ignore, &cost_ignore);
        //计算第一行其余元素的256个视差的代价
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
        h_win_range_sads = (HorizWinRangeSAD *) ring_buf_deque(horiz_cost_rbuf); //取出第一行
        assert(h_win_range_sads);

        //blk_ref_ptr = ref_ptr + extend_offset;
        //blk_tgt_ptr = tgt_ptr + extend_offset + max_disp_;
        blk_ref_ptr = ref_ptr;
        blk_tgt_ptr = tgt_ptr + 127;

        ///////////////////////////////////////////////////////////////////////
        // 减去被移出匹配窗口的行，这是与启动阶段的差别
        // curr_line_sad = prev_line_sad - (last) h_win_range_sads
        VectorSubUInt16(curr_line_sad + 0 * line_sad_step + line_sad_align_offset,
                        h_win_range_sads->sad_buffer + 0, disp_range_);  //减去第一行的第一个元素
        ///////////////////////////////////////////////////////////////////////

        // 启动行内SAD的计算：最开始的 @sliding_blk_sz_ 个长为 @disp_range_ 的AbsDiff向量相加
        memset(h_win_range_sads->sad_buffer, 0, sizeof(h_win_range_sads->sad_buffer[0]) * disp_range_);
        //计算新一行的第一个元素的横向SAD
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
        //计算新一行的第一个元素的总SAD
        UpdateVeritalAccumExtSAD(h_win_range_sads->sad_buffer + col * disp_range_,
                                 curr_line_sad + col * line_sad_step, disp_range_,
                                 disp_ptr + col, cost_ptr + col);  //加上新计算出来的行sad,即为总的SAD
        //计算新一行的其余元素的SAD
        for (col = 1; col < width_; ++col)
        {
            ///////////////////////////////////////////////////////////////////
            // 减去被移出匹配窗口的行，这是与启动阶段的差别
            // curr_line_sad = prev_line_sad - (last) h_win_range_sads
            VectorSubUInt16(curr_line_sad + col * line_sad_step + line_sad_align_offset,
                            h_win_range_sads->sad_buffer + col * disp_range_, disp_range_);  //需要减去原来行
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
