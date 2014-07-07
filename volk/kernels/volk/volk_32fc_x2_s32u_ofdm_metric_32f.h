#include <inttypes.h>
#include <stdio.h>
#include "volk_32fc_32f_dot_prod_32fc.h"
#include "volk_32fc_conjugate_32fc.h"
#include "volk_32fc_x2_multiply_32fc.h"
#include "volk_32fc_magnitude_squared_32f.h"
#include "volk_32f_x2_add_32f.h"
#include "volk_32f_x2_dot_prod_32f.h"
#include "volk_32f_x2_divide_32f.h"
#include <string.h>
#include <volk/volk_malloc.h>

#define SSE_ALIGNMENT 16
#define AVX_ALIGNMENT 32

#ifndef INCLUDED_volk_32fc_x2_s32u_ofdm_metric_32f_a_H
#define INCLUDED_volk_32fc_x2_s32u_ofdm_metric_32f_a_H

#ifdef LV_HAVE_GENERIC
/*!
  \brief Computes timing utility metric for ofdm frame detection
  \param cVector The output vector of utility metric
  \param bVector The input vector that is to be delayed
  \param aVector The input vector that is to be conjugated
  \param num_points Number of points for which metric is to be computed
*/

static inline void volk_32fc_x2_s32u_ofdm_metric_32f_a_generic(float* cVector, lv_32fc_t* bVector, lv_32fc_t* aVector, uint32_t fftLen, unsigned int num_points){

	lv_32fc_t* aPtr = aVector;
	lv_32fc_t* bPtr = bVector;
	float* cPtr = cVector;

	//lv_32fc_t* aConj = (lv_32fc_t *)malloc(num_points*sizeof(lv_32fc_t));
	lv_32fc_t* aConj = (lv_32fc_t *)volk_malloc(num_points*sizeof(lv_32fc_t), SSE_ALIGNMENT);
	lv_32fc_t* aConjPtr = aConj;
	//lv_32fc_t* bDelay = (lv_32fc_t *)malloc(num_points*sizeof(lv_32fc_t));
	lv_32fc_t* bDelay = (lv_32fc_t *)volk_malloc(num_points*sizeof(lv_32fc_t), SSE_ALIGNMENT);
	lv_32fc_t* bDelayPtr = bDelay;
	//lv_32fc_t* abProd = (lv_32fc_t *)malloc(num_points*sizeof(lv_32fc_t));
	lv_32fc_t* abProd = (lv_32fc_t *)volk_malloc(num_points*sizeof(lv_32fc_t), SSE_ALIGNMENT);
	lv_32fc_t* abProdPtr = abProd;
	//float* abSum = (float *)malloc(num_points*sizeof(float));
	float* abSum = (float *)volk_malloc(num_points*sizeof(float), SSE_ALIGNMENT);
	float* abSumPtr = abSum;
	//float* normalizer = (float *)malloc(num_points*sizeof(float));
	float* normalizer = (float *)volk_malloc(num_points*sizeof(float), SSE_ALIGNMENT);
	float* normPtr = normalizer;

	uint32_t L = fftLen/2;
	uint32_t j = 0;
	
	//float* taps = (float *)malloc(L*sizeof(float));
	float* taps = (float *)volk_malloc(L*sizeof(float), SSE_ALIGNMENT);
	for(j = 0; j < L; j++)	taps[j] = 1.0;

	unsigned int outputPoints = 0;

	volk_32fc_conjugate_32fc_generic(aConjPtr, aPtr, num_points);	
	
	for(j = 0; j < L; j++)	*bDelayPtr++ = lv_cmake(0.0, 0.0);
	memcpy(bDelayPtr, bPtr, (num_points-L)*sizeof(lv_32fc_t));
	bDelayPtr = bDelay;

	aConjPtr = aConj;
	bDelayPtr = bDelay;

	volk_32fc_x2_multiply_32fc_generic(abProd, bDelayPtr, aConjPtr, num_points);

	//float* aConjSqr = (float *)malloc(num_points*sizeof(float));
    //float* bDelaySqr = (float *)malloc(num_points*sizeof(float));
	float* aConjSqr = (float *)volk_malloc(num_points*sizeof(float), SSE_ALIGNMENT);
    float* bDelaySqr = (float *)volk_malloc(num_points*sizeof(float), SSE_ALIGNMENT);
	
	volk_32fc_magnitude_squared_32f_a_generic(aConjSqr, aConjPtr, num_points);
	volk_32fc_magnitude_squared_32f_a_generic(bDelaySqr, bDelayPtr, num_points);
	volk_32f_x2_add_32f_generic(abSumPtr, aConjSqr, bDelaySqr, num_points);	

	// Reuse aConj, find a better way
	for(; outputPoints < L; outputPoints++){
		volk_32fc_32f_dot_prod_32fc_generic(aConjPtr++, abProdPtr, taps, outputPoints+1);
		volk_32f_x2_dot_prod_32f_generic(normPtr++, abSumPtr, taps, outputPoints+1);
	}

		abProdPtr++;
		abSumPtr++;
	
	for(; outputPoints < num_points; outputPoints++){
		volk_32fc_32f_dot_prod_32fc_generic(aConjPtr++, abProdPtr++, taps, L);
		volk_32f_x2_dot_prod_32f_generic(normPtr++, abSumPtr++, taps, L);
	}


	normPtr = normalizer;	
	aConjPtr = aConj;
	
	volk_32fc_magnitude_squared_32f_a_generic(cPtr, aConjPtr, num_points);	
	volk_32f_x2_divide_32f_generic(cPtr, cPtr, normPtr, num_points);	

/*	
	free(taps);
	free(abSum);
	free(abProd);
	free(aConj);
	free(bDelay);
	free(normalizer);
	free(aConjSqr);
	free(bDelaySqr);
*/
	volk_free(taps);
	volk_free(abSum);
	volk_free(abProd);
	volk_free(aConj);
	volk_free(bDelay);
	volk_free(normalizer);
	volk_free(aConjSqr);
	volk_free(bDelaySqr);

}

#endif /* LV_HAVE_GENERIC */

#endif /* INCLUDED_volk_32fc_x2_s32u_ofdm_metric_32f_a_H */
