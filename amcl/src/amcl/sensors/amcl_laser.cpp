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
///////////////////////////////////////////////////////////////////////////
//
// Desc: AMCL laser routines
// Author: Andrew Howard
// Date: 6 Feb 2003
// CVS: $Id: amcl_laser.cc 7057 2008-10-02 00:44:06Z gbiggs $
//
///////////////////////////////////////////////////////////////////////////

#include <sys/types.h> // required by Darwin
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#ifdef HAVE_UNISTD_H

#include <unistd.h>

#endif

#include "amcl/sensors/amcl_laser.h"

using namespace amcl;

////////////////////////////////////////////////////////////////////////////////
// Default constructor
AMCLLaser::AMCLLaser(size_t max_beams, map_t *map) : AMCLSensor(),
                                                     max_samples(0), max_obs(0),
                                                     temp_obs(NULL) {
    this->time = 0.0;

    this->max_beams = max_beams;
    this->map = map;

    return;
}

AMCLLaser::~AMCLLaser() {
    if (temp_obs) {
        for (int k = 0; k < max_samples; k++) {
            delete[] temp_obs[k];
        }
        delete[]temp_obs;
    }
}

void
AMCLLaser::SetModelBeam(double z_hit,
                        double z_short,
                        double z_max,
                        double z_rand,
                        double sigma_hit,
                        double lambda_short,
                        double chi_outlier) {
    this->model_type = LASER_MODEL_BEAM;
    this->z_hit = z_hit;
    this->z_short = z_short;
    this->z_max = z_max;
    this->z_rand = z_rand;
    this->sigma_hit = sigma_hit;
    this->lambda_short = lambda_short;
    this->chi_outlier = chi_outlier;
}

void
AMCLLaser::SetModelLikelihoodField(double z_hit,
                                   double z_rand,
                                   double sigma_hit,
                                   double max_occ_dist) {
    this->model_type = LASER_MODEL_LIKELIHOOD_FIELD;
    this->z_hit = z_hit;
    this->z_rand = z_rand;
    this->sigma_hit = sigma_hit;

    map_update_cspace(this->map, max_occ_dist);
}

void
AMCLLaser::SetModelLikelihoodFieldProb(double z_hit,
                                       double z_rand,
                                       double sigma_hit,
                                       double max_occ_dist,
                                       bool do_beamskip,
                                       double beam_skip_distance,
                                       double beam_skip_threshold,
                                       double beam_skip_error_threshold) {
    this->model_type = LASER_MODEL_LIKELIHOOD_FIELD_PROB;
    this->z_hit = z_hit;
    this->z_rand = z_rand;
    this->sigma_hit = sigma_hit;
    this->do_beamskip = do_beamskip;
    this->beam_skip_distance = beam_skip_distance;
    this->beam_skip_threshold = beam_skip_threshold;
    this->beam_skip_error_threshold = beam_skip_error_threshold;
    map_update_cspace(this->map, max_occ_dist);
}


////////////////////////////////////////////////////////////////////////////////
// Apply the laser sensor model
bool AMCLLaser::UpdateSensor(pf_t *pf, AMCLSensorData *data) {
    if (this->max_beams < 2)
        return false;

    // Apply the laser sensor model
    if (this->model_type == LASER_MODEL_BEAM)
        pf_update_sensor(pf, (pf_sensor_model_fn_t) BeamModel, data);
    else if (this->model_type == LASER_MODEL_LIKELIHOOD_FIELD)
        // LikelihoodFieldModel 基于激光束的几何模型，将激光测量与地图栅格的占用概率进行比较，从而计算权重。
        // 它提供了更丰富的信息，能够考虑激光束与障碍物的相对位置和占用状态，从而在一定程度上提高了准确性。
        pf_update_sensor(pf, (pf_sensor_model_fn_t) LikelihoodFieldModel, data);
    else if (this->model_type == LASER_MODEL_LIKELIHOOD_FIELD_PROB)
        // 主要是针对动态障碍物概率去除,选通过判断每个粒子相同方向的激光点,在地图内与障碍物的最近距离,是否小于一定距离阈值,
        // 并且超过一定数量阈值, 就不进行计算权重,如果在符合范围内就权重叠加.
        pf_update_sensor(pf, (pf_sensor_model_fn_t) LikelihoodFieldModelProb, data);
    else
        pf_update_sensor(pf, (pf_sensor_model_fn_t) BeamModel, data);

    return true;
}


////////////////////////////////////////////////////////////////////////////////
// Determine the probability for the given pose
// 激光测量(波束)模型
double AMCLLaser::BeamModel(AMCLLaserData *data, pf_sample_set_t *set) {
    AMCLLaser *self;
    int i, j, step;
    double z, pz;
    double p;
    double map_range;
    double obs_range, obs_bearing;
    double total_weight;
    pf_sample_t *sample;
    pf_vector_t pose;

    self = (AMCLLaser *) data->sensor;

    total_weight = 0.0;

    // Compute the sample weights
    // 计算样本的权重
    for (j = 0; j < set->sample_count; j++) {
        sample = set->samples + j;
        pose = sample->pose;

        // Take account of the laser pose relative to the robot
        // 考虑激光传感器在机器人上的姿态偏移
        pose = pf_vector_coord_add(self->laser_pose, pose);

        p = 1.0;

        step = (data->range_count - 1) / (self->max_beams - 1);
        for (i = 0; i < data->range_count; i += step) {
            // 获取激光测量数据
            // 代表激光束测量得到的距离值
            obs_range = data->ranges[i][0];
            // 代表激光束的方位角
            obs_bearing = data->ranges[i][1];

            // Compute the range according to the map
            // 根据地图计算测量的距离
            // pose.v[0] 和 pose.v[1]：是机器人当前位置的 x 和 y 坐标。
            map_range = map_calc_range(self->map, pose.v[0], pose.v[1],
                                       pose.v[2] + obs_bearing, data->range_max);
            pz = 0.0;

            // Part 1: good, but noisy, hit
            // 第一部分：测量和地图预测相符，但有噪声
            // 计算了实际测量距离 obs_range 与地图预测距离 map_range 之间的差值，用变量 z 来表示。
            z = obs_range - map_range;

            // 从权重的角度来看，z 的值越接近零，粒子的权重越大，而 z 的值越大，粒子的权重越小。这种权重调整是根据激光测量数据与地图匹配程度的一种表现
            pz += self->z_hit * exp(-(z * z) / (2 * self->sigma_hit * self->sigma_hit));

            // Part 2: short reading from unexpected obstacle (e.g., a person)
            // 第二部分：预测障碍物但测量到的距离较短（例如，人的影响）
            // self->lambda_short 是一个负值，用于在短测量模型中控制概率的衰减速率。这有助于模拟短测量情况下的不确定性，以便在权重更新中适当地减小粒子的权重
            if (z < 0)
                pz += self->z_short * self->lambda_short * exp(-self->lambda_short * obs_range);

            // Part 3: Failure to detect obstacle, reported as max-range
            // 第三部分：未能探测到障碍物，测量报告为最大距离
            if (obs_range == data->range_max)
                // 将最大距离处的观测数据的权重增加了 --> 因为这边不确定性多，需要更多例粒子
                pz += self->z_max * 1.0;

            // Part 4: Random measurements
            // 第四部分：随机测量
            if (obs_range < data->range_max)
                pz += self->z_rand * 1.0 / data->range_max;

            // TODO: outlier rejection for short readings
            // TODO：对于短测量，可能需要排除异常值

            assert(pz <= 1.0);
            assert(pz >= 0.0);
            //      p *= pz;
            // here we have an ad-hoc weighting scheme for combining beam probs
            // works well, though...
            // 这里使用了一种自定义的权重方案来组合激光束的概率
            // 尽管这是一种自定义的方法，但效果很好...
            p += pz * pz * pz;
        }

        // 在概率定位中，通过将不同激光测量模型的部分结合，可以更好地考虑实际测量与地图预测之间的差异，以及各种不同情况下的权重影响。
        // 这有助于更准确地更新粒子的权重，以反映机器人在环境中的位置和姿态。在这里，使用了一个自定义的权重组合方法，将各部分的权重
        // 增量相乘，并用于更新粒子的权重。
        sample->weight *= p;
        total_weight += sample->weight;
    }

    return (total_weight);
}

double AMCLLaser::LikelihoodFieldModel(AMCLLaserData *data, pf_sample_set_t *set) {
    AMCLLaser *self;
    int i, j, step;
    // z表示从激光束的测量位置（endpoint）到地图上最近障碍物的距离，也就是击中模型所使用的障碍物距离。这个距离是根据激光测量数据和机器人的位姿计算出来的。
    double z, pz;
    double p;
    double obs_range, obs_bearing;
    double total_weight;
    pf_sample_t *sample;
    pf_vector_t pose;
    // 向量存储了激光束射线与地图上的物体或障碍物相交的位置，也就是激光束测量的终点。
    pf_vector_t hit;

    self = (AMCLLaser *) data->sensor;
    // ... 初始化和变量声明 ...

    total_weight = 0.0;

    // Compute the sample weights
    // 计算样本的权重
    for (j = 0; j < set->sample_count; j++) {
        sample = set->samples + j;
        // 遍历每个粒子，这是粒子对应的位姿，是经运动模型更新后的先验位姿
        pose = sample->pose;

        // Take account of the laser pose relative to the robot
        // 考虑激光在机器人上的相对姿态
        // self->laser_pose 进行pose选择和变换 ,使得激光点变换到世界坐标
        pose = pf_vector_coord_add(self->laser_pose, pose);

        p = 1.0;

        // Pre-compute a couple of things
        // 测量噪声的方差
        double z_hit_denom = 2 * self->sigma_hit * self->sigma_hit;
        // 随机测量相关，1/雷达最大距离
        double z_rand_mult = 1.0 / data->range_max;

        step = (data->range_count - 1) / (self->max_beams - 1);

        // Step size must be at least 1
        // 步长至少为1
        if (step < 1)
            step = 1;

        for (i = 0; i < data->range_count; i += step) {
            // ... 获取激光测量数据 ...
            obs_range = data->ranges[i][0];
            obs_bearing = data->ranges[i][1];

            // This model ignores max range readings
            // 此模型忽略了最大测量范围的读数
            if (obs_range >= data->range_max)
                continue;

            // Check for NaN
            // 检查NaN值
            if (obs_range != obs_range)
                continue;

            pz = 0.0;

            // Compute the endpoint of the beam
            // pose.v[0] 和 pose.v[1] 表示机器人的当前位置（x 和 y 坐标），pose.v[2] 表示机器人的当前方向（姿态角）。
            // 通过将距离 obs_range 与方向 obs_bearing 转换为笛卡尔坐标系的增量，可以计算出激光束射线的击中点。
            hit.v[0] = pose.v[0] + obs_range * cos(pose.v[2] + obs_bearing);
            hit.v[1] = pose.v[1] + obs_range * sin(pose.v[2] + obs_bearing);

            // Convert to map grid coords.
            // 转换为地图栅格坐标
            int mi, mj;
            // 用于将世界坐标系中的X坐标转换为地图网格坐标系中的X坐标
            mi = MAP_GXWX(self->map, hit.v[0]);
            // 用于将世界坐标系中的Y坐标转换为地图网格坐标系中的Y坐标
            mj = MAP_GYWY(self->map, hit.v[1]);

            // Part 1: Get distance from the hit to closest obstacle.
            // Off-map penalized as max distance
            if (!MAP_VALID(self->map, mi, mj))
                // 如果激光束射线射出的位置在地图之外，那么 z 被设置为地图的最大障碍物距离。这是一个默认值，表示激光射线射出地图的范围之外。
                z = self->map->max_occ_dist;
            else
                // 如果激光束射线射出的位置在地图内，那么 z 被设置为地图网格中对应位置的障碍物距离。这个距离可能是根据地图中的障碍物分布计算得到的，用于描述地图上各个点到最近障碍物的距离。
                z = self->map->cells[MAP_INDEX(self->map, mi, mj)].occ_dist;
            // Gaussian model
            // 高斯模型
            // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
            // 这里使用了高斯分布来建模，z_hit 参数影响了高斯分布的峰值和形状。较大的z_hit 值会使高斯分布更加尖锐，对于匹配较好的激光测量，权重增加。较小的z_hit 值会导致高斯分布更加平坦，对于距离不太匹配的激光测量，权重减少。
            // self->z_hit：这是一个参数，用于控制激光测量的击中模型的权重。它表示激光测量击中物体的概率。
            // z：这是机器人到地图上最近障碍物的距离，是之前计算得到的值。
            // z_hit_denom：这是一个常数，用于计算高斯分布的分母。它与激光测量噪声有关。

            // 整个计算表达的意思是：当激光射出点到地图上最近障碍物的距离越接近，击中模型的概率就会增加。这符合我们对于激光测量的预期：测量越接近障碍物，就越有可能是准确的测量，从而机器人的位姿估计也会更准确。
            pz += self->z_hit * exp(-(z * z) / z_hit_denom);
            // Part 2: random measurements
            // 第二部分：随机测量
            // 整个计算表达的意思是，通过将 self->z_rand 与 z_rand_mult 相乘，你得到了随机测量模型的概率。这个概率部分描述了激光测量可能是由于传感器的随机噪声引起，而不是实际的击中物体。
            // 这个部分通常会为每个激光束测量都添加一个相同的小的概率，以考虑随机误差对位姿估计的影响。
            pz += self->z_rand * z_rand_mult;

            // TODO: outlier rejection for short readings

            assert(pz <= 1.0);
            assert(pz >= 0.0);
            //      p *= pz;
            // here we have an ad-hoc weighting scheme for combining beam probs
            // works well, though...
            // 这里使用了一种自定义的权重方案来组合激光束的概率
            // 尽管这是一种自定义的方法，但效果很好...
            p += pz * pz * pz;
        }

        sample->weight *= p;
        // 计算每个粒子的权重后，通过累加所有粒子的权重，得到了total_weight。这个值可以用来进行粒子权重的归一化，以便在后续的重采样步骤中更好地选择具有更高权重的粒子，从而提高定位的准确性。
        total_weight += sample->weight;
    }

    return (total_weight);
}
// 主要是针对动态障碍物概率去除,选通过判断每个粒子相同方向的激光点,在地图内与障碍物的最近距离,是否小于一定距离阈值,
// 并且超过一定数量阈值, 就不进行计算权重,如果在符合范围内就权重叠加.
double AMCLLaser::LikelihoodFieldModelProb(AMCLLaserData *data, pf_sample_set_t *set) {
    AMCLLaser *self;
    int i, j, step;//点云序号,粒子序号,一个区域有多少个点
    double z, pz;//障碍物距离,概率值
    double log_p;//概率对数值
    double obs_range, obs_bearing;//数据距离,和角度
    double total_weight;//总权重
    pf_sample_t *sample;//单粒子
    pf_vector_t pose;//机器人位姿->雷达世界坐标
    pf_vector_t hit;//世界坐标
    // 最后都是在创建对象AMCLSensor
    self = (AMCLLaser *) data->sensor;
    // 初始化权重
    total_weight = 0.0;
    // 算要分成多少个区域
    // data->range_count: 这是激光测量数据中的点云数量，即一帧中激光传感器测量得到的激光点的数量。
    // self->max_beams: 这是算法中设定的一个值，表示每个区域中最多包含的激光束数量。
    // 整个表达式的目的是将激光测量数据分成一定数量的区域
    step = ceil((data->range_count) / static_cast<double>(self->max_beams));

    // Step size must be at least 1
    if (step < 1)
        step = 1;

    // Pre-compute a couple of things
    // 预算噪声协方差和随机概率
    double z_hit_denom = 2 * self->sigma_hit * self->sigma_hit;
    double z_rand_mult = 1.0 / data->range_max;
    // 高斯分布模型,概率分布,取最大距离的概率
    double max_dist_prob = exp(-(self->map->max_occ_dist * self->map->max_occ_dist) / z_hit_denom);

    //Beam skipping - ignores beams for which a majority of particles do not agree with the map
    //prevents correct particles from getting down weighted because of unexpected obstacles
    //such as humans
    //光束跳过-忽略大多数粒子与贴图不一致的光束
    //防止正确的粒子由于意外的障碍物而降低权重
    //比如人类
    bool do_beamskip = self->do_beamskip;//如果需要进行去除动态障碍运算
    double beam_skip_distance = self->beam_skip_distance;//障碍物激光点栅格坐标与地图障碍最近距离距离阈值
    double beam_skip_threshold = self->beam_skip_threshold;//光束跳过阈值

    //we only do beam skipping if the filter has converged
    //如果粒子群收敛就不进行设置跳过光束,就不进行去除动态障碍物
    if (do_beamskip && !set->converged) {
        do_beamskip = false;
    }

    //we need a count the no of particles for which the beam agreed with the map
    // 我们需要计算出光束与地图一致的粒子数
    int *obs_count = new int[self->max_beams]();//区域内粒子数量

    //we also need a mask of which observations to integrate (to decide which beams to integrate to all particles)
    // 我们还需要一个obs_mask来合并观测结果决定哪些光束要整合到所有粒子上
    bool *obs_mask = new bool[self->max_beams]();

    //一桢数据中采样 点云序号
    int beam_ind = 0;

    //realloc indicates if we need to reallocate the temp data structure needed to do beamskipping
    //realloc(重新分配)表示是否需要重新分配执行beamskipping所需的temp数据结构
    bool realloc = false;

    if (do_beamskip) {
        if (self->max_obs < self->max_beams) {
            realloc = true;
        }
        //如果粒子群粒子粒子数量大于粒子群设定最大值,重新分配
        if (self->max_samples < set->sample_count) {
            realloc = true;
        }
        // 重新分配
        if (realloc) {
            self->reallocTempData(set->sample_count, self->max_beams);
            fprintf(stderr, "Reallocing temp weights %d - %d\n", self->max_samples, self->max_obs);
        }
    }

    // Compute the sample weights
    // 遍历所有激光雷达数据点
    for (j = 0; j < set->sample_count; j++) {
        //粒子序号
        sample = set->samples + j;
        //获取机器人位姿
        pose = sample->pose;

        // Take account of the laser pose relative to the robot
        //结合机器人相对雷达的坐标,求出雷达在世界的坐标
        pose = pf_vector_coord_add(self->laser_pose, pose);

        log_p = 0;

        beam_ind = 0;
        // beam_ind 检测点序号,每个区域为step个粒子
        for (i = 0; i < data->range_count; i += step, beam_ind++) {
            // 当前点云的测量距离
            obs_range = data->ranges[i][0];
            // 当前点云的测量角度
            obs_bearing = data->ranges[i][1];

            // This model ignores max range readings
            // 点云距离合理判断
            if (obs_range >= data->range_max) {
                continue;
            }

            // Check for NaN
            // NaN 值的一个特点。NaN 与任何值都不相等，包括自身。
            if (obs_range != obs_range) {
                continue;
            }

            pz = 0.0;

            // Compute the endpoint of the beam
            // 激光点云结合雷达坐标求出点云在世界的坐标

            // 激光点云在世界坐标系中的x坐标
            hit.v[0] = pose.v[0] + obs_range * cos(pose.v[2] + obs_bearing);
            // 激光点云在世界坐标系中的y坐标
            hit.v[1] = pose.v[1] + obs_range * sin(pose.v[2] + obs_bearing);

            // Convert to map grid coords.
            // 世界转栅格坐标
            int mi, mj;
            mi = MAP_GXWX(self->map, hit.v[0]);
            mj = MAP_GYWY(self->map, hit.v[1]);

            // Part 1: Get distance from the hit to closest obstacle.
            // Off-map penalized as max distance

            if (!MAP_VALID(self->map, mi, mj)) {
                //如果超出地图范围,根据最大距离概率,计算概率值
                pz += self->z_hit * max_dist_prob;
            } else {
                z = self->map->cells[MAP_INDEX(self->map, mi, mj)].occ_dist;
                //否则求出该激光点离障碍物最短的距离.可以理解为当前束激光测的距离和地图匹配的距离做差值
                //但是它实际是在激光点在附近每个方向都进行搜索最近的障碍的距离,取最小值

                //激光点与地图障碍距离小于阈值
                if (z < beam_skip_distance) {
                    //一帧内第 beam_ind 个点云,与地图一致的粒子数量obs_count加一
                    obs_count[beam_ind] += 1;
                }
                //障碍权重乘以高斯概率模型
                pz += self->z_hit * exp(-(z * z) / z_hit_denom);
            }

            // Gaussian model
            // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)

            // Part 2: random measurements
            //随机权重乘以随机概率模型
            pz += self->z_rand * z_rand_mult;

            assert(pz <= 1.0);
            assert(pz >= 0.0);

            // TODO: outlier rejection for short readings
            // 如果不跳过,就算他的概率对数,并求和各个粒子的
            // 不做跳变进行独立处理，不做跳变时权重正常相加
            if (!do_beamskip) {
                log_p += log(pz);
            } else {
                //第j个粒子的第beam_ind点云 的距离值
                self->temp_obs[j][beam_ind] = pz;
            }
        }
        if (!do_beamskip) {
            //更新粒子权重
            sample->weight *= exp(log_p);
            //权重求和
            total_weight += sample->weight;
        }
    }
    // 做跳变时，如果跳变的个数占总个数一定权重时(即大于光束跳过阈值)，设定障碍物遮挡
    if (do_beamskip) {
        int skipped_beam_count = 0;
        for (beam_ind = 0; beam_ind < self->max_beams; beam_ind++) {
            if ((obs_count[beam_ind] / static_cast<double>(set->sample_count)) > beam_skip_threshold) {
                obs_mask[beam_ind] = true;
            } else {
                // 否则不进行遮挡
                obs_mask[beam_ind] = false;
                skipped_beam_count++;
            }
        }

        //we check if there is at least a critical number of beams that agreed with the map
        //otherwise it probably indicates that the filter converged to a wrong solution
        //if that's the case we integrate all the beams and hope the filter might converge to
        //the right solution
        bool error = false;
        //跳过的数量 大于 无效波束比率的阈值-则定位收敛出错
        //如果skipped_beam_count/beam_ind >= beam_skip_error_threshold 则收敛到一个错误的位姿
        //beam_ind为一桢数据中采样个数
        if (skipped_beam_count >= (beam_ind * self->beam_skip_error_threshold)) {
            fprintf(stderr,
                    "Over %f%% of the observations were not in the map - pf may have converged to wrong pose - integrating all observations\n",
                    (100 * self->beam_skip_error_threshold));
            error = true;
        }

        for (j = 0; j < set->sample_count; j++) {
            //由set指向第一个粒子的地址开始
            sample = set->samples + j;
            pose = sample->pose;

            log_p = 0;

            for (beam_ind = 0; beam_ind < self->max_beams; beam_ind++) {
                //收敛错误或者第j个粒子的第beam_ind点云不是是遮挡光束则被用于计算概率和
                if (error || obs_mask[beam_ind]) {
                    log_p += log(self->temp_obs[j][beam_ind]);
                }
            }

            sample->weight *= exp(log_p);

            total_weight += sample->weight;//权重总和
        }
    }

    delete[] obs_count;
    delete[] obs_mask;
    return (total_weight);
}

void AMCLLaser::reallocTempData(int new_max_samples, int new_max_obs) {
    if (temp_obs) {
        for (int k = 0; k < max_samples; k++) {
            delete[] temp_obs[k];
        }
        delete[]temp_obs;
    }
    max_obs = new_max_obs;
    max_samples = fmax(max_samples, new_max_samples);

    temp_obs = new double *[max_samples]();
    for (int k = 0; k < max_samples; k++) {
        temp_obs[k] = new double[max_obs]();
    }
}
