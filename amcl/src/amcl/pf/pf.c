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
 * Desc: Simple particle filter for localization.
 * Author: Andrew Howard
 * Date: 10 Dec 2002
 * CVS: $Id: pf.c 6345 2008-04-17 01:36:39Z gerkey $
 *************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "amcl/pf/pf.h"
#include "amcl/pf/pf_pdf.h"
#include "amcl/pf/pf_kdtree.h"
#include "portable_utils.hpp"


// Compute the required number of samples, given that there are k bins
// with samples in them.
static int pf_resample_limit(pf_t *pf, int k);



// Create a new filter
pf_t *pf_alloc(int min_samples, int max_samples,
               double alpha_slow, double alpha_fast,
               pf_init_model_fn_t random_pose_fn, void *random_pose_data)
{
  int i, j;
  pf_t *pf;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  // 调用了C语言的库函数srand48完成随机数种子的初始化
  srand48(time(NULL));

  // 通过calloc申请滤波器对象的内存。
  pf = calloc(1, sizeof(pf_t));

  pf->random_pose_fn = random_pose_fn;
  pf->random_pose_data = random_pose_data;

  pf->min_samples = min_samples;
  pf->max_samples = max_samples;

  // Control parameters for the population size calculation.  [err] is
  // the max error between the true distribution and the estimated
  // distribution.  [z] is the upper standard normal quantile for (1 -
  // p), where p is the probability that the error on the estimated
  // distrubition will be less than [err].
  pf->pop_err = 0.01;
  pf->pop_z = 3;
  pf->dist_threshold = 0.5; 

  // Number of leaf nodes is never higher than the max number of samples
  pf->limit_cache = calloc(max_samples, sizeof(int));
  
  pf->current_set = 0;
  for (j = 0; j < 2; j++)
  {
    set = pf->sets + j;
      
    set->sample_count = max_samples;
    set->samples = calloc(max_samples, sizeof(pf_sample_t));

    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;
      sample->pose.v[0] = 0.0;
      sample->pose.v[1] = 0.0;
      sample->pose.v[2] = 0.0;
      sample->weight = 1.0 / max_samples;
    }

    // HACK: is 3 times max_samples enough?
    // 构建了三倍于样本集合尺寸的直方图对象kdtree
    set->kdtree = pf_kdtree_alloc(3 * max_samples);

    set->cluster_count = 0;
    set->cluster_max_count = max_samples;
    set->clusters = calloc(set->cluster_max_count, sizeof(pf_cluster_t));

    set->mean = pf_vector_zero();
    set->cov = pf_matrix_zero();
  }

  pf->w_slow = 0.0;
  pf->w_fast = 0.0;

  pf->alpha_slow = alpha_slow;
  pf->alpha_fast = alpha_fast;

  //set converged to 0
  pf_init_converged(pf);

  return pf;
}

// Free an existing filter
void pf_free(pf_t *pf)
{
  int i;

  free(pf->limit_cache);

  for (i = 0; i < 2; i++)
  {
    free(pf->sets[i].clusters);
    pf_kdtree_free(pf->sets[i].kdtree);
    free(pf->sets[i].samples);
  }
  free(pf);
  
  return;
}

// Initialize the filter using a gaussian
// pf是将要初始化的滤波器对象
// mean和cov分别是机器人初始位姿和协方差描述
void pf_init(pf_t *pf, pf_vector_t mean, pf_matrix_t cov)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  pf_pdf_gaussian_t *pdf;
  
  set = pf->sets + pf->current_set;
  
  // Create the kd tree for adaptive sampling
  // 先清空直方图，并根据刚刚构建的正态分布采样，填充粒子集合更新统计直方图。
  pf_kdtree_clear(set->kdtree);

  set->sample_count = pf->max_samples;

  pdf = pf_pdf_gaussian_alloc(mean, cov);
  //通过一个for循环完成两个粒子集合的初始化工作。并在7到13行中的为每个粒子集申请内存，并赋予初值。在第14行中构建了三倍于样本集合尺寸的直方图对象kdtree。 最后完成粒子簇的内存申请和均值方差的初始化。
  // Compute the new sample poses
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;
    sample->weight = 1.0 / pf->max_samples;
    sample->pose = pf_pdf_gaussian_sample(pdf);

    // Add sample to histogram
    pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
  }

  pf->w_slow = pf->w_fast = 0.0;

  pf_pdf_gaussian_free(pdf);
    
  // Re-compute cluster statistics
  pf_cluster_stats(pf, set); 

  //通过函数pf_init_converged设定滤波器未收敛，并返回滤波器对象。
  //set converged to 0
  pf_init_converged(pf);

  return;
}


// Initialize the filter using some model
void pf_init_model(pf_t *pf, pf_init_model_fn_t init_fn, void *init_data)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;

  set = pf->sets + pf->current_set;

  // Create the kd tree for adaptive sampling
  pf_kdtree_clear(set->kdtree);

  set->sample_count = pf->max_samples;

  // Compute the new sample poses
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;
    sample->weight = 1.0 / pf->max_samples;
    sample->pose = (*init_fn) (init_data);

    // Add sample to histogram
    pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
  }

  pf->w_slow = pf->w_fast = 0.0;

  // Re-compute cluster statistics
  pf_cluster_stats(pf, set);
  
  //set converged to 0
  pf_init_converged(pf);

  return;
}

void pf_init_converged(pf_t *pf){
  pf_sample_set_t *set;
  set = pf->sets + pf->current_set;
  set->converged = 0; 
  pf->converged = 0; 
}
// 检测粒子滤波算法是否收敛
int pf_update_converged(pf_t *pf)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  double total;

  set = pf->sets + pf->current_set;
  double mean_x = 0, mean_y = 0;
  // 计算样本的位置均值
  for (i = 0; i < set->sample_count; i++){
    sample = set->samples + i;

    mean_x += sample->pose.v[0];
    mean_y += sample->pose.v[1];
  }
  mean_x /= set->sample_count;
  mean_y /= set->sample_count;
  // 检查每个样本的位置是否接近均值，如果不是则认为没有收敛
  for (i = 0; i < set->sample_count; i++){
    sample = set->samples + i;
    if(fabs(sample->pose.v[0] - mean_x) > pf->dist_threshold || 
       fabs(sample->pose.v[1] - mean_y) > pf->dist_threshold){
      set->converged = 0; 
      pf->converged = 0; 
      return 0;
    }
  }
  // 设置集合和算法为已收敛状态
  set->converged = 1; 
  pf->converged = 1; 
  return 1; 
}

// Update the filter with some new action
void pf_update_action(pf_t *pf, pf_action_model_fn_t action_fn, void *action_data)
{
  pf_sample_set_t *set;

  set = pf->sets + pf->current_set;

  (*action_fn) (action_data, set);
  
  return;
}


#include <float.h>
//Augmented_MCL算法。概算法的核心思想就是， 通过对比短期和长期样本的均值估计来确定注入随机粒子的数量，并在重采样过程中注入新的随机粒子。

// 完成的是短期和长期样本均值的估计， 也就是字段w_slow和w_fast的计算。该函数在测量更新的过程中被调用。
// 它有三个参数，其中pf是滤波器对象，sensor_fn则是一个函数指针实现了传感器的模型，用于计算各个粒子的概
// 率权重。参数sensor_data则是用于更新的传感器数据。
// Update the filter with some new sensor observation
void pf_update_sensor(pf_t *pf, pf_sensor_model_fn_t sensor_fn, void *sensor_data)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  double total;
  // 首先获取当前激活的样本集合，并通过函数指针sensor_fn计算各个样本的概率权重。
  set = pf->sets + pf->current_set;

  // Compute the sample weights
  // 计算样本权重总和
  total = (*sensor_fn) (sensor_data, set);

  // set->n_effective 的值越接近 1，表示粒子集合越均匀，越适合进行下一步的状态估计操作。如果 n_effective 的值较小，说明粒子权重分布不均衡，可能需要进行重采样来重新调整权重，以保持粒子集合的有效性。
  set->n_effective = 0;
  
  if (total > 0.0)
  {
    // 权重归一化并计算平均权重
    // Normalize weights
    // 计算粒子的平均权重
    double w_avg=0.0;
    // 完成样本权重的归一化。
    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;
      w_avg += sample->weight;
      sample->weight /= total;
      // 每个粒子的权重都会被平方并加到 n_effective 上。这个操作实际上在强调权重较大的粒子的影响，因为权重较大的粒子会对 n_effective 有更大的贡献。
      // 当权重分布不均衡时，一些粒子的权重很小，它们的平方项会减小 n_effective 的值，从而提示算法可能需要进行重采样来调整权重分布，增加较重要粒子的数量。
      set->n_effective += sample->weight*sample->weight;
    }
    // 按照右边算法伪代码中的第10和11行的公式，计算统计量w_slow,w_fast。
    // Update running averages of likelihood of samples (Prob Rob p258)
    w_avg /= set->sample_count;

    // 当粒子质量下降时，平均权重随之下降，w_slow、w_fast也会随之下降，但是显然w_fast下降的速度要快于w_slow
    // 计算w_slow(长期似然平均估计)
    // w_slow 利用一个比较小的权重更新系数 alpha_slow 来逐步调整自身，以适应权重的变化趋势。如果权重变化趋势比较缓慢
    if(pf->w_slow == 0.0)
      pf->w_slow = w_avg;
    else
      pf->w_slow += pf->alpha_slow * (w_avg - pf->w_slow);

    // 计算w_fast(短期似然平均估计)
    // w_fast 也会根据权重的变化进行调整，但是使用了一个比较大的权重更新系数 alpha_fast，因此它对权重变化的敏感性较高
    if(pf->w_fast == 0.0)
      pf->w_fast = w_avg;
    else
      pf->w_fast += pf->alpha_fast * (w_avg - pf->w_fast);
    //printf("w_avg: %e slow: %e fast: %e\n", 
           //w_avg, pf->w_slow, pf->w_fast);
  }
  else
  {
    // 处理总权重为零的情况
    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;
      sample->weight = 1.0 / set->sample_count;
    }
  }

  set->n_effective = 1.0/set->n_effective;
  return;
}

// copy set a to set b

void copy_set(pf_sample_set_t* set_a, pf_sample_set_t* set_b)
{
  int i;
  double total;
  pf_sample_t *sample_a, *sample_b;

  // Clean set b's kdtree
  // 清空 set_b 的 kd 树
  pf_kdtree_clear(set_b->kdtree);

  // Copy samples from set a to create set b
  // 复制样本从 set_a 到 set_b
  total = 0;
  set_b->sample_count = 0;

  for(i = 0; i < set_a->sample_count; i++)
  {
    sample_b = set_b->samples + set_b->sample_count++;

    sample_a = set_a->samples + i;

    assert(sample_a->weight > 0);

    // Copy sample a to sample b
    // 复制样本 a 到样本 b
    sample_b->pose = sample_a->pose;
    sample_b->weight = sample_a->weight;

    total += sample_b->weight;

    // Add sample to histogram
    // 将样本添加到 kd 树中，用于自适应采样
    pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);
  }

  // Normalize weights
  for (i = 0; i < set_b->sample_count; i++)
  {
    sample_b = set_b->samples + i;
    sample_b->weight /= total;
  }

  set_b->converged = set_a->converged;
}
// 具体实现了重采样操作，相比于传统的粒子滤波器的重采样操作，它多了一个根据统计量w_slow,w_fast插入随机样本的过程。
// 它只有一个参数pf用于指示进行重采样更新的滤波器对象。在函数的一开始用两个指针set_a和set_b分别获取了当前激活的粒子集合和待切换的集合。
// Resample the distribution
void pf_update_resample(pf_t *pf)
{
  int i;
  double total;
  pf_sample_set_t *set_a, *set_b;
  pf_sample_t *sample_a, *sample_b;

  //double r,c,U;
  //int m;
  //double count_inv;
  double* c;

  double w_diff;

  set_a = pf->sets + pf->current_set;
  set_b = pf->sets + (pf->current_set + 1) % 2;

  if (pf->selective_resampling != 0)
  {
    // 如果有效样本数量（n_effective）超过了总样本数量的一半,就不执行重采样，直接复制。有效样本数量通常是根据粒子权重来计算的，表示有多少个粒子在重采样过程中会被保留下来。
    // 0.5这个阈值实际上是一个经验性的选择，通常用于决定是否进行重采样。这个阈值并没有一个严格的理论依据，而是在实际应用中根据经验确定的。
    if (set_a->n_effective > 0.5*(set_a->sample_count))
    {
      // copy set a to b
      copy_set(set_a,set_b);

      // Re-compute cluster statistics
      pf_cluster_stats(pf, set_b);

      // Use the newly created sample set
      // 这种切换样本集合的操作通常会在滤波的不同迭代步骤中进行，以便在每一步中使用新的样本进行状态预测和更新。
      pf->current_set = (pf->current_set + 1) % 2;
      return;
    }
  }

  // Build up cumulative probability table for resampling.
  // TODO: Replace this with a more efficient procedure
  // (e.g., http://www.network-theory.co.uk/docs/gslref/GeneralDiscreteDistributions.html)
  // 对粒子的权重进行积分，获得分布函数，后面用于生成样本。
  c = (double*)malloc(sizeof(double)*(set_a->sample_count+1));
  c[0] = 0.0;
  for(i=0;i<set_a->sample_count;i++)
    // 权重累积
    c[i+1] = c[i]+set_a->samples[i].weight;
  // 我们重置更新粒子集合set_b，它将保存重采样后的粒子。
  // Create the kd tree for adaptive sampling
  pf_kdtree_clear(set_b->kdtree);
  
  // Draw samples from set a to create set b.
  total = 0;
  set_b->sample_count = 0;
  // 接着计算算法Augmented_MCL中第13行的判据，该判据是用于判定是否需要插入新粒子。
  w_diff = 1.0 - pf->w_fast / pf->w_slow;
  if(w_diff < 0.0)
    w_diff = 0.0;
  //printf("w_diff: %9.6f\n", w_diff);

  // Can't (easily) combine low-variance sampler with KLD adaptive
  // sampling, so we'll take the more traditional route.
  /*
  // Low-variance resampler, taken from Probabilistic Robotics, p110
  count_inv = 1.0/set_a->sample_count;
  r = drand48() * count_inv;
  c = set_a->samples[0].weight;
  i = 0;
  m = 0;
  */
  // 在循环中先生成一个随机数，并根据算法的第13行的判据插入新粒子。
  while(set_b->sample_count < pf->max_samples)
  {
    sample_b = set_b->samples + set_b->sample_count++;

    if(drand48() < w_diff)
      // 插入随机粒子
      sample_b->pose = (pf->random_pose_fn)(pf->random_pose_data);
    // 如果判据不生效，那么我们就根据原粒子集合所描述的样本分布进行采样。
    else
    {
      // Can't (easily) combine low-variance sampler with KLD adaptive
      // sampling, so we'll take the more traditional route.
      /*
      // Low-variance resampler, taken from Probabilistic Robotics, p110
      U = r + m * count_inv;
      while(U>c)
      {
        i++;
        // Handle wrap-around by resetting counters and picking a new random
        // number
        if(i >= set_a->sample_count)
        {
          r = drand48() * count_inv;
          c = set_a->samples[0].weight;
          i = 0;
          m = 0;
          U = r + m * count_inv;
          continue;
        }
        c += set_a->samples[i].weight;
      }
      m++;
      */

      // Naive discrete event sampler
      double r;
      r = drand48();
      for(i=0;i<set_a->sample_count;i++)
      {
        // 轮盘赌（独立随机采样） https://blog.csdn.net/HERO_CJN/article/details/93458674
        // 权重大的粒子被选取的概率更大
        if((c[i] <= r) && (r < c[i+1]))
          break;
      }
      assert(i<set_a->sample_count);

      sample_a = set_a->samples + i;

      assert(sample_a->weight > 0);

      // Add sample to list
      sample_b->pose = sample_a->pose;
    }
    // 为重采样后的粒子赋予相同的权重，并更新统计直方图。当满足样本限制的时候，就退出循环，结束重采样。
    sample_b->weight = 1.0;
    total += sample_b->weight;

    // Add sample to histogram
    pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);

    // See if we have enough samples yet
    if (set_b->sample_count > pf_resample_limit(pf, set_b->kdtree->leaf_count))
      break;
  }
  // 在函数的最后，重置统计量w_slow,w_fast，归一化样本权重，更新粒子簇，交换当前粒子集，检查粒子集是否收敛。
  // Reset averages, to avoid spiraling off into complete randomness.
  if(w_diff > 0.0)
    pf->w_slow = pf->w_fast = 0.0;

  //fprintf(stderr, "\n\n");
  // // 归一化样本权重
  // Normalize weights
  for (i = 0; i < set_b->sample_count; i++)
  {
    sample_b = set_b->samples + i;
    sample_b->weight /= total;
  }
  
  // Re-compute cluster statistics
  pf_cluster_stats(pf, set_b);

  // Use the newly created sample set
  pf->current_set = (pf->current_set + 1) % 2;

  // 检测粒子滤波算法是否收敛 -->通过检查每个粒子是不是接近均值
  pf_update_converged(pf);

  free(c);
  return;
}


// Compute the required number of samples, given that there are k bins
// with samples in them.  This is taken directly from Fox et al.
int pf_resample_limit(pf_t *pf, int k)
{
  double a, b, c, x;
  int n;

  // Return max_samples in case k is outside expected range, this shouldn't
  // happen, but is added to prevent any runtime errors
  if (k < 1 || k > pf->max_samples)
      return pf->max_samples;

  // Return value if cache is valid, which means value is non-zero positive
  if (pf->limit_cache[k-1] > 0)
    return pf->limit_cache[k-1];

  if (k == 1)
  {
    pf->limit_cache[k-1] = pf->max_samples;
    return pf->max_samples;
  }

  a = 1;
  b = 2 / (9 * ((double) k - 1));
  c = sqrt(2 / (9 * ((double) k - 1))) * pf->pop_z;
  x = a - b + c;

  n = (int) ceil((k - 1) / (2 * pf->pop_err) * x * x * x);

  if (n < pf->min_samples)
  {
    pf->limit_cache[k-1] = pf->min_samples;
    return pf->min_samples;
  }
  if (n > pf->max_samples)
  {
    pf->limit_cache[k-1] = pf->max_samples;
    return pf->max_samples;
  }
  
  pf->limit_cache[k-1] = n;
  return n;
}


// Re-compute the cluster statistics for a sample set
// 根据样本的位置和权重计算聚类的统计信息，包括均值、权重和协方差等，以及整体滤波器的统计信息。
// 这些信息在后续的粒子滤波算法中可能会被用来进行权重更新、重采样或预测步骤，以更好地估计系统状态的概率分布。
void pf_cluster_stats(pf_t *pf, pf_sample_set_t *set)
{
  int i, j, k, cidx;
  pf_sample_t *sample;
  pf_cluster_t *cluster;
  
  // Workspace
  double m[4], c[2][2];
  size_t count;
  double weight;

  // Cluster the samples
  pf_kdtree_cluster(set->kdtree);
  
  // Initialize cluster stats
  set->cluster_count = 0;

  for (i = 0; i < set->cluster_max_count; i++)
  {
    cluster = set->clusters + i;
    cluster->count = 0;
    cluster->weight = 0;
    cluster->mean = pf_vector_zero();
    cluster->cov = pf_matrix_zero();

    for (j = 0; j < 4; j++)
      cluster->m[j] = 0.0;
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        cluster->c[j][k] = 0.0;
  }

  // Initialize overall filter stats
  count = 0;
  weight = 0.0;
  set->mean = pf_vector_zero();
  set->cov = pf_matrix_zero();
  for (j = 0; j < 4; j++)
    m[j] = 0.0;
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
      c[j][k] = 0.0;
  
  // Compute cluster stats
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;

    //printf("%d %f %f %f\n", i, sample->pose.v[0], sample->pose.v[1], sample->pose.v[2]);

    // Get the cluster label for this sample
    cidx = pf_kdtree_get_cluster(set->kdtree, sample->pose);
    assert(cidx >= 0);
    if (cidx >= set->cluster_max_count)
      continue;
    if (cidx + 1 > set->cluster_count)
      set->cluster_count = cidx + 1;
    
    cluster = set->clusters + cidx;

    cluster->count += 1;
    cluster->weight += sample->weight;

    count += 1;
    weight += sample->weight;

    // Compute mean
    cluster->m[0] += sample->weight * sample->pose.v[0];
    cluster->m[1] += sample->weight * sample->pose.v[1];
    cluster->m[2] += sample->weight * cos(sample->pose.v[2]);
    cluster->m[3] += sample->weight * sin(sample->pose.v[2]);

    m[0] += sample->weight * sample->pose.v[0];
    m[1] += sample->weight * sample->pose.v[1];
    m[2] += sample->weight * cos(sample->pose.v[2]);
    m[3] += sample->weight * sin(sample->pose.v[2]);

    // Compute covariance in linear components
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
      {
        cluster->c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
        c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
      }
  }

  // Normalize
  for (i = 0; i < set->cluster_count; i++)
  {
    cluster = set->clusters + i;
        
    cluster->mean.v[0] = cluster->m[0] / cluster->weight;
    cluster->mean.v[1] = cluster->m[1] / cluster->weight;
    cluster->mean.v[2] = atan2(cluster->m[3], cluster->m[2]);

    cluster->cov = pf_matrix_zero();

    // Covariance in linear components
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        cluster->cov.m[j][k] = cluster->c[j][k] / cluster->weight -
          cluster->mean.v[j] * cluster->mean.v[k];

    // Covariance in angular components; I think this is the correct
    // formula for circular statistics.
    cluster->cov.m[2][2] = -2 * log(sqrt(cluster->m[2] * cluster->m[2] +
                                         cluster->m[3] * cluster->m[3]) / cluster->weight);

    //printf("cluster %d %d %f (%f %f %f)\n", i, cluster->count, cluster->weight,
           //cluster->mean.v[0], cluster->mean.v[1], cluster->mean.v[2]);
    //pf_matrix_fprintf(cluster->cov, stdout, "%e");
  }

  assert(fabs(weight) >= DBL_EPSILON);
  if (fabs(weight) < DBL_EPSILON)
  {
    printf("ERROR : divide-by-zero exception : weight is zero\n");
    return;
  }
  // Compute overall filter stats
  set->mean.v[0] = m[0] / weight;
  set->mean.v[1] = m[1] / weight;
  set->mean.v[2] = atan2(m[3], m[2]);

  // Covariance in linear components
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
      set->cov.m[j][k] = c[j][k] / weight - set->mean.v[j] * set->mean.v[k];

  // Covariance in angular components; I think this is the correct
  // formula for circular statistics.
  set->cov.m[2][2] = -2 * log(sqrt(m[2] * m[2] + m[3] * m[3]));

  return;
}

void pf_set_selective_resampling(pf_t *pf, int selective_resampling)
{
  pf->selective_resampling = selective_resampling;
}

// Compute the CEP statistics (mean and variance).
void pf_get_cep_stats(pf_t *pf, pf_vector_t *mean, double *var)
{
  int i;
  double mn, mx, my, mrr;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  
  set = pf->sets + pf->current_set;

  mn = 0.0;
  mx = 0.0;
  my = 0.0;
  mrr = 0.0;
  
  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;

    mn += sample->weight;
    mx += sample->weight * sample->pose.v[0];
    my += sample->weight * sample->pose.v[1];
    mrr += sample->weight * sample->pose.v[0] * sample->pose.v[0];
    mrr += sample->weight * sample->pose.v[1] * sample->pose.v[1];
  }

  assert(fabs(mn) >= DBL_EPSILON);
  if (fabs(mn) < DBL_EPSILON)
  {
    printf("ERROR : divide-by-zero exception : mn is zero\n");
    return;
  }

  mean->v[0] = mx / mn;
  mean->v[1] = my / mn;
  mean->v[2] = 0.0;

  *var = mrr / mn - (mx * mx / (mn * mn) + my * my / (mn * mn));

  return;
}


// Get the statistics for a particular cluster.
int pf_get_cluster_stats(pf_t *pf, int clabel, double *weight,
                         pf_vector_t *mean, pf_matrix_t *cov)
{
  pf_sample_set_t *set;
  pf_cluster_t *cluster;

  set = pf->sets + pf->current_set;

  if (clabel >= set->cluster_count)
    return 0;
  cluster = set->clusters + clabel;

  *weight = cluster->weight;
  *mean = cluster->mean;
  *cov = cluster->cov;

  return 1;
}


