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
#include <iostream>
#include <cstring>
#include "ring_buffer.h"
#include "depth_engine_img.h"
#include <assert.h>
#include <emmintrin.h>

typedef struct _SlidingWindowData
{
    int rbuf_size;                 // Ring buffer size
    ring_buf_t *horiz_cost_rbuf;   // horizontal cost ring buffer

    int line_sad_align_offset;     // 行 SAD 对齐偏移（两边各自的扩展量）
    int line_sad_step;             // 行 SAD 的步长，为了方便求亚像素做了扩展，步长为 disp_range_ + 2*line_sad_align_offset
    int line_sad_nelem;            // 行 SAD 总的大小
    uint16_t *curr_line_sad;       // 当前行的sad，维度为[width][disp_range + 2*line_sad_align_offset]
} SlidingWindowData;


typedef struct _HorizWinRangeSAD
{
    int disp_range;
    int width;
    uint16_t *sad_buffer;


}HorizWinRangeSAD;


class DispSlidingWindow
{
public:
    DispSlidingWindow();
    ~DispSlidingWindow();

    void setParams(int block_size, int min_disp, int disp_range,
                   int width, int height, int max_disp);

    int AllocHorizSADRingBuf(SlidingWindowData *helper) const;
    int SlidingWindowMatch(DEImageU8 ext_tgt, DEImageU8 ext_ref, DEImageS16 disp, DEImageU16 cost);


private:

    int block_size_;
    int min_disp_;
    int half_blk_sz_;
    int disp_range_;
    int width_;
    int height_;
    int max_disp_;

};


int FreeHorizSADRingBuf(SlidingWindowData *helper);

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
                            uint8_t tgt_add, const uint8_t *ref_add);

/**
 * 向量减
 * @param minuend       被减数
 * @param subtrahend    减数
 * @param len           向量长度
 */
void VectorSubUInt16(uint16_t *minuend, const uint16_t *subtrahend, int16_t len);

/**
 * 向量加，并求最小值
 * @param[in]  src       源向量
 * @param[in]  tgt       目标向量
 * @param[in]  num       向量长度
 * @param[out] min_idx   最小值的索引（从0开始）
 * @param[out] min_val   最小值
 */
void VectorAddFindMinUInt16(const uint16_t *src, uint16_t *tgt, int16_t num, int16_t *min_idx, uint16_t *min_val);


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
                         int16_t *disp, uint16_t *min_cost);