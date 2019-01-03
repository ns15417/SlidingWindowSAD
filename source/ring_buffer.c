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
 *  Author: liuxianzhuo@orbbec.com
 *
 ****************************************************************************/
#include <malloc.h>
#include "ring_buffer.h"


ring_buf_t *ring_buf_new(int capacity)
{
    ring_buf_t *rb = malloc(sizeof(struct ring_buf));

    if (!rb)
        return rb;

    rb->size = capacity + 1;
    rb->buffer = malloc(rb->size * sizeof(void *));  // 要多申请一个，用于区分满和空
    if (!rb->buffer)
    {
        free(rb);
        return NULL;
    }

    rb->head = 0;
    rb->tail = 0;

    return rb;
}

void ring_buf_free(ring_buf_t *rbuf)
{
    if (!rbuf)
        return;

    if (rbuf->buffer)
        free(rbuf->buffer);

    free(rbuf);
}

int ring_buf_num_item(const ring_buf_t *rbuf)
{
    return (rbuf->head - rbuf->tail) % rbuf->size;
}

int ring_buf_full(const ring_buf_t *rbuf)
{
    return (rbuf->size - 1) == ((rbuf->head - rbuf->tail) % rbuf->size);
}

int ring_buf_empty(const ring_buf_t *rbuf)
{
    return rbuf->head == rbuf->tail;
}

void ring_buf_enque(ring_buf_t *rbuf, void *data)
{
    if (ring_buf_full(rbuf))
    {
        rbuf->tail = (rbuf->tail + 1) % rbuf->size;
    }

    rbuf->buffer[rbuf->head] = data;
    rbuf->head = (rbuf->head + 1) % rbuf->size;
//    printf("rbuf->head = %d\n", rbuf->head);
}

void *ring_buf_deque(ring_buf_t *rbuf)
{
    if (ring_buf_empty(rbuf))
    {
        return NULL;
    }

    void *result = rbuf->buffer[rbuf->tail];
    rbuf->tail = (rbuf->tail + 1) % rbuf->size;
//    printf("rbuf->tail = %d\n", rbuf->tail);

    return result;
}
//
//void * ring_buf_head(ring_buf_t *rbuf)
//{
//    if (ring_buf_empty(rbuf))
//    {
//        return NULL;
//    }
//
//    return rbuf->buffer[rbuf->head];
//}
//
//void * ring_buf_tail(ring_buf_t *rbuf)
//{
//    if (ring_buf_empty(rbuf))
//    {
//        return NULL;
//    }
//
//    return rbuf->buffer[rbuf->tail];
//}
