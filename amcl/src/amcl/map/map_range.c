/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
/**************************************************************************
 * Desc: Range routines
 * Author: Andrew Howard
 * Date: 18 Jan 2003
 * CVS: $Id: map_range.c 1347 2003-05-05 06:24:33Z inspectorg $
**************************************************************************/

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "amcl/map/map.h"

// Extract a single range reading from the map.  Unknown cells and/or
// out-of-bound cells are treated as occupied, which makes it easy to
// use Stage bitmap files.
// map是我们的地图对象，ox、oy、oa分别指示了机器人的位置坐标和方向角，max_range则是激光传感器的量程。
double map_calc_range(map_t *map, double ox, double oy, double oa, double max_range)
{
  // Bresenham raytracing
  int x0,x1,y0,y1;
  int x,y;
  int xstep, ystep;
  char steep;
  int tmp;
  int deltax, deltay, error, deltaerr;
  // 首先通过宏MAP_GXWX和MAP_GYWY获取扫描波束的起点和终点所对应的删格索引
  // 通过 MAP_GXWX 和 MAP_GYWY 宏来将实际坐标转换为地图中的栅格坐标，以获取起点 x0、y0
  x0 = MAP_GXWX(map,ox);
  y0 = MAP_GYWY(map,oy);

  // 通过 MAP_GXWX 和 MAP_GYWY 宏来将实际坐标转换为地图中的栅格坐标，以获取终点 x1、y1。
  x1 = MAP_GXWX(map,ox + max_range * cos(oa));
  y1 = MAP_GYWY(map,oy + max_range * sin(oa));

  // 当扫描波束的y轴差异大于x轴差异时，将该方向标记为陡峭的(steep = 1)

  // 如果 abs(y1 - y0) 大于 abs(x1 - x0)，表示在 y 方向上的差异比在 x 方向上的差异更大，即斜率大于 1。这种情况下，将 steep 设置为 1，表示将进行坐标轴交换。否则，如果斜率不大于 1，将 steep 设置为 0，表示不需要坐标轴交换。

  // 这个操作的目的是为了保证 Bresenham 射线跟踪算法能够正确地在各种方向上进行射线追踪，无论是斜率大于 1 还是斜率小于等于 1。在交换坐标轴后，算法将能够正确地处理不同方向上的射线追踪，从而准确地预测机器人可能的碰撞距离。
  if(abs(y1-y0) > abs(x1-x0))
    steep = 1;
  else
    steep = 0;
  // 如果是陡峭的，则交换x和y的坐标。
  if(steep)
  {
    tmp = x0;
    x0 = y0;
    y0 = tmp;

    tmp = x1;
    x1 = y1;
    y1 = tmp;
  }
  // 下面分别计算x和y轴的索引扫描偏差量，并为error和deltaerr赋予初值。
  deltax = abs(x1-x0);
  deltay = abs(y1-y0);
  error = 0;
  deltaerr = deltay;

  x = x0;
  y = y0;
  // 然后根据x和y轴的发展方向，设定搜索方向。
  if(x0 < x1)
    xstep = 1;
  else
    xstep = -1;
  if(y0 < y1)
    ystep = 1;
  else
    ystep = -1;
  // 我们先检查机器人的位置坐标，确保它在地图的上，并且没有处于未知区域或者占用删格中。
  if(steep)
  {
    if(!MAP_VALID(map,y,x) || map->cells[MAP_INDEX(map,y,x)].occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
  }
  else
  {
    if(!MAP_VALID(map,x,y) || map->cells[MAP_INDEX(map,x,y)].occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
  }

  while(x != (x1 + xstep * 1))
  {
    x += xstep;
    error += deltaerr;
    if(2*error >= deltax)
    {
      y += ystep;
      error -= deltax;
    }

    if(steep)
    {
      if(!MAP_VALID(map,y,x) || map->cells[MAP_INDEX(map,y,x)].occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
    }
    else
    {
      if(!MAP_VALID(map,x,y) || map->cells[MAP_INDEX(map,x,y)].occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
    }
  }
  // 如果迭代完成，没有发现障碍物，直接返回预测的最大距离 max_range。
  return max_range;
}
