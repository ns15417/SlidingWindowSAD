/*****************************************************************************
 *  Orbbec
 *  Copyright (C) 2018 by ORBBEC Technology., Inc.
 *
 *  This file is part of Orbbec DepthEngine.
 *
 *  This file belongs to ORBBEC Technology., Inc.
 *  It is considered a trade secret, and is not to be divulged or used by
 * parties who have NOT received written authorization from the owner.
 *
 *  Description:
 *
 ****************************************************************************/
#ifndef __OB_DEPTH_ENGINE_MX6500_RING_BUFFER_H__
#define __OB_DEPTH_ENGINE_MX6500_RING_BUFFER_H__

#include <stddef.h>

/**
 * @brief Ring buffer of void * (so we can hold any thing)
 */
typedef struct ring_buf
{
    void **buffer;
    int head;
    int tail;
    int size;
} ring_buf_t;

#ifdef __cplusplus
extern "C" {
#endif


/**
 * 创建循环缓冲区
 * @param capacity  容量
 * @return
 */
ring_buf_t *ring_buf_new(int capacity);

/**
 * 释放循环缓冲区
 * @param rbuf
 */
void ring_buf_free(ring_buf_t *rbuf);

/**
 * 获取循环缓冲区的当前容量
 * @param rbuf
 */
int ring_buf_num_item(const ring_buf_t *rbuf);

/**
 * 是否已满
 * @param rbuf
 * @return
 */
int ring_buf_full(const ring_buf_t *rbuf);

/**
 * 是否为空
 * @param rbuf
 * @return
 */
int ring_buf_empty(const ring_buf_t *rbuf);

/**
 * 从头部添加一个数据
 * @param rbuf
 * @param data
 */
void ring_buf_enque(ring_buf_t *rbuf, void *data);

/**
 * 从尾部取出一个数据
 * @param rbuf
 * @param data
 */
void *ring_buf_deque(ring_buf_t *rbuf);

#ifdef __cplusplus
}
#endif


#endif //__OB_DEPTH_ENGINE_MX6500_CIRCLE_BUFFER_H__
