/*****************************************************************************
 *  Orbbec DepthEngine MX6500
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
#ifndef __OB_DEPTH_ENGINE_MX6500_DEPTH_ENGINE_IMG_H__
#define __OB_DEPTH_ENGINE_MX6500_DEPTH_ENGINE_IMG_H__

#include "stdint.h"

typedef struct DEImageU8_s
{
    uint8_t *data;      /// 图像数据，行优先存储
    int16_t width;      /// 图像宽度
    int16_t height;     /// 图像高度
} DEImageU8;

typedef struct DEImageS8_s
{
    int8_t *data;      /// 图像数据，行优先存储
    int16_t width;      /// 图像宽度
    int16_t height;     /// 图像高度
} DEImageS8;

typedef struct DEImageU16_s
{
    uint16_t *data;     /// 图像数据，行优先存储
    int16_t width;      /// 图像宽度
    int16_t height;     /// 图像高度
} DEImageU16;

typedef struct DEImageS16_s
{
    int16_t *data;      /// 图像数据，行优先存储
    int16_t width;      /// 图像宽度
    int16_t height;     /// 图像高度
} DEImageS16;


typedef struct DEPoint_s
{
    int16_t x;
    int16_t y;
} DEPoint;

#endif //__OB_DEPTH_ENGINE_MX6500_DEPTH_ENGINE_IMG_H__
