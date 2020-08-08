/*
 * cuda_export.hpp
 *
 *  Created on: May 28, 2020
 *      Author: dmarce1
 */

#ifndef NTIGER_CUDA_EXPORT_HPP_
#define NTIGER_CUDA_EXPORT_HPP_




#ifdef __CUDA_ARCH__
#define CUDA_EXPORT __device__ __host__
#define CONSTEXPR
#else
#define CUDA_EXPORT
#define CONSTEXPR constexpr
#endif



#endif /* NTIGER_CUDA_EXPORT_HPP_ */
