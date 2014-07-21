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
#include <volk/volk_complex.h>
#include <volk/volk_circularBuffer.h>

#define SSE_ALIGNMENT 16
#define AVX_ALIGNMENT 32

#ifndef INCLUDED_volk_32fc_x2_s32u_ofdm_metric_32f_a_H
#define INCLUDED_volk_32fc_x2_s32u_ofdm_metric_32f_a_H

#ifdef LV_HAVE_SSSE3
/*!
  \brief Computes timing utility metric for ofdm frame detection
  \param cVector The output vector of utility metric
  \param bVector The input vector that is to be delayed
  \param aVector The input vector that is to be conjugated
  \param num_points Number of points for which metric is to be computed
*/

/* Assume fftLen is a multiple of 8 to avoid 4-sized chunks spanning boundaries */
#include <tmmintrin.h>
static inline void volk_32fc_x2_s32u_ofdm_metric_32f_a_ssse3(float* cVector, const lv_32fc_t* bVector, const lv_32fc_t* aVector, uint32_t fftLen, unsigned int num_points){

	const lv_32fc_t* aPtr = aVector;
	const lv_32fc_t* bPtr = bVector;
	float* cPtr = cVector;

	uint32_t L = fftLen/2, outputPoints = 0, quarterPoints = num_points / 4;
	float *temp1, *temp2;
	float x[4];
	__m128 aConj1, aConj2, bDelay1, bDelay2, abProd1, abProd2;
	__m128 abSum, aConjSqr, bDelaySqr, cVal;
	__m128 yl, yh, tmp1, tmp2;
	__m128 conjugator = _mm_setr_ps(0.f, -0.f, 0.f, -0.f);
	__m128 norm = _mm_setzero_ps();
	__m128 realpt = _mm_setzero_ps();
	__m128 imagpt = _mm_setzero_ps();

	circularBuffer *ringProdReal = circularBuffer_create(L);	
	circularBuffer *ringProdImag = circularBuffer_create(L);	
	circularBuffer *ringSum = circularBuffer_create(L);	

	for(outputPoints = 0; outputPoints < quarterPoints; outputPoints++){

		aConj1 = _mm_load_ps((float *)aPtr);
		aPtr += 2;
		aConj2 = _mm_load_ps((float *)aPtr);
		aPtr += 2;
		aConj1 = _mm_xor_ps(aConj1, conjugator);
		aConj2 = _mm_xor_ps(aConj2, conjugator);
	
		if(4*outputPoints + 4 <= L){
			bDelay1 = _mm_setzero_ps();
			bDelay2 = _mm_setzero_ps();
		}
		else{
			bDelay1 = _mm_load_ps((float *)bPtr);
			bPtr += 2;
			bDelay2 = _mm_load_ps((float *)bPtr);
			bPtr += 2;
		}

    	yl = _mm_moveldup_ps(aConj1); 						
    	yh = _mm_movehdup_ps(aConj1);						
    	tmp1 = _mm_mul_ps(bDelay1,yl);
    	bDelay1 = _mm_shuffle_ps(bDelay1, bDelay1, 0xB1);
    	tmp2 = _mm_mul_ps(bDelay1,yh); 
		abProd1 = _mm_addsub_ps(tmp1, tmp2);

		yl = _mm_moveldup_ps(aConj2); 
		yh = _mm_movehdup_ps(aConj2);
		tmp1 = _mm_mul_ps(bDelay2,yl);
		bDelay2 = _mm_shuffle_ps(bDelay2,bDelay2,0xB1);
		tmp2 = _mm_mul_ps(bDelay2,yh);
		abProd2 = _mm_addsub_ps(tmp1, tmp2);

		aConj1 = _mm_mul_ps(aConj1, aConj1); 	
		aConj2 = _mm_mul_ps(aConj2, aConj2); 	
		aConjSqr = _mm_hadd_ps(aConj1, aConj2);
	
		bDelay1 = _mm_mul_ps(bDelay1, bDelay1); 	
		bDelay2 = _mm_mul_ps(bDelay2, bDelay2); 	
		bDelaySqr = _mm_hadd_ps(bDelay1, bDelay2);

		abSum = _mm_add_ps(aConjSqr, bDelaySqr);

		temp1 = (float *)&realpt;
		realpt = _mm_set1_ps(*(temp1+3));
		temp2 = (float *)&imagpt;
		imagpt = _mm_set1_ps(*(temp2+3));

		if((4*outputPoints + 4) > L){

			circularBuffer_read(ringProdReal, x, 4);
			x[1] += x[0];
			x[2] += x[1];
			x[3] += x[2];
			realpt = _mm_sub_ps(realpt, _mm_setr_ps(x[0], x[1], x[2], x[3]));
			circularBuffer_read(ringProdImag, x, 4);
			x[1] += x[0];
			x[2] += x[1];
			x[3] += x[2];
			imagpt = _mm_sub_ps(imagpt, _mm_setr_ps(x[0], x[1], x[2], x[3]));
		}

		temp1 = (float *)&abProd1;
		temp2 = (float *)&abProd2;

		x[0] = *temp1;
		x[1] =  *(temp1+2);
		x[2] = *temp2;
		x[3] = *(temp2+2);
		circularBuffer_write(ringProdReal, x, 4);
		x[1] += x[0];
		x[2] += x[1];
		x[3] += x[2];
		realpt = _mm_add_ps(realpt, _mm_setr_ps(x[0], x[1], x[2], x[3]));
		x[0] = *(temp1+1);
		x[1] =  *(temp1+3);
		x[2] = *(temp2+1);
		x[3] = *(temp2+3);
		circularBuffer_write(ringProdImag, x, 4);
		x[1] += x[0];
		x[2] += x[1];
		x[3] += x[2];
		imagpt = _mm_add_ps(imagpt, _mm_setr_ps(x[0], x[1], x[2], x[3]));

		temp1 = (float *)&norm;
		norm = _mm_set1_ps(*(temp1+3));  

		if((4*outputPoints + 4) > L){
			circularBuffer_read(ringSum, x, 4);
			x[1] += x[0];
			x[2] += x[1];
			x[3] += x[2];
			norm = _mm_sub_ps(norm, _mm_setr_ps(x[0], x[1], x[2], x[3]));		
		}	
	
		temp1 = (float *)&abSum;
		x[0] = *temp1;
		x[1] = *(temp1+1);
		x[2] = *(temp1+2);
		x[3] = *(temp1+3);
		circularBuffer_write(ringSum, x, 4);
		x[1] += x[0];
		x[2] += x[1];
		x[3] += x[2];
		norm = _mm_add_ps(norm, _mm_setr_ps(x[0], x[1], x[2], x[3]));

		cVal = _mm_add_ps(_mm_mul_ps(realpt, realpt), _mm_mul_ps(imagpt, imagpt));
		cVal = _mm_div_ps(cVal, norm);
		_mm_store_ps(cPtr, cVal);
		cPtr += 4;
	}

	unsigned int i = outputPoints * 4;
	lv_32fc_t aConj, bDelay, abProd;
	float sabSum, saConjSqr, sbDelaySqr;
	float res[2];
	float snorm = 0, sread;
	float *tmp, *srealpt = &res[0], *simagpt = &res[1];
	*srealpt = 0; *simagpt = 0;
	
	tmp = (float *)&realpt;
	*srealpt = *(tmp + 3);
	tmp = (float *)&imagpt;
	*simagpt = *(tmp + 3);
	tmp = (float *)&norm;
	snorm = *(tmp+3);  

	for(; i < num_points; i++){

		aConj = lv_conj(*aPtr++);
	
		if(i < L)		bDelay = lv_cmake(0.0, 0.0);
		else			bDelay = *bPtr++; 
	
		abProd = bDelay * aConj;

		tmp = (float *)&aConj;
		saConjSqr = ((*tmp) * (*tmp)) + ((*(tmp+1)) * (*(tmp+1)));
		tmp = (float *)&bDelay;
		sbDelaySqr = ((*tmp) * (*tmp)) + ((*(tmp+1)) * (*(tmp+1)));
		sabSum = (saConjSqr) + (sbDelaySqr);

		if(i > L){
			circularBuffer_read(ringProdReal, &sread, 1);
			*srealpt -= sread;
			circularBuffer_read(ringProdImag, &sread, 1);
			*simagpt -= sread;
		}


		tmp = (float *)&abProd;
		circularBuffer_write(ringProdReal, tmp, 1);
		circularBuffer_write(ringProdImag, tmp + 1, 1);
		*srealpt += (*tmp);
		*simagpt += (*(tmp+1));

		if(i > L){
			circularBuffer_read(ringSum, &sread, 1);
			snorm -= sread;		
		}

		tmp = &sabSum;
		circularBuffer_write(ringSum, tmp, 1);
		snorm += *tmp;

		*cPtr = ((*srealpt) * (*srealpt)) + ((*(simagpt)) * (*(simagpt)));
		*cPtr = (*cPtr) / snorm;
		cPtr++;
	}

	//circularBuffer_destroy(ringProdReal);
	//circularBuffer_destroy(ringProdImag);
	//circularBuffer_destroy(ringSum);
}
#endif /* LV_HAVE_SSE */

#ifdef LV_HAVE_GENERIC
/*!
  \brief Computes timing utility metric for ofdm frame detection
  \param cVector The output vector of utility metric
  \param bVector The input vector that is to be delayed
  \param aVector The input vector that is to be conjugated
  \param num_points Number of points for which metric is to be computed
*/

static inline void volk_32fc_x2_s32u_ofdm_metric_32f_generic(float* cVector, const lv_32fc_t* bVector, const lv_32fc_t* aVector, uint32_t fftLen, unsigned int num_points){

	const lv_32fc_t* aPtr = aVector;
	const lv_32fc_t* bPtr = bVector;
	float* cPtr = cVector;

	lv_32fc_t aConj, bDelay, abProd, abDot;
	float abSum, aConjSqr, bDelaySqr, norm;
	uint32_t L = fftLen/2;
	circularBuffer *ringProdReal = circularBuffer_create(20);	
	circularBuffer *ringProdImag = circularBuffer_create(20);	
	circularBuffer *ringSum = circularBuffer_create(20);	
	unsigned int outputPoints = 0;	
	float res[2];
	float *tmp, *realpt = &res[0], *imagpt = &res[1];
	*realpt = 0; *imagpt = 0;
	float dotProd = 0, read;

	for(outputPoints = 0; outputPoints < num_points; outputPoints++){

	aConj = lv_conj(*aPtr++);
	
	if(outputPoints < L)	bDelay = lv_cmake(0.0, 0.0);
	else					bDelay = *bPtr++; // check that bPtr is not incremented first
	
	abProd = (bDelay) * (aConj);
	tmp = (float *)&aConj;
	aConjSqr = ((*tmp) * (*tmp)) + ((*(tmp+1)) * (*(tmp+1)));
	tmp = (float *)&bDelay;
	bDelaySqr = ((*tmp) * (*tmp)) + ((*(tmp+1)) * (*(tmp+1)));
	abSum = (aConjSqr) + (bDelaySqr);

	if(outputPoints < L){
		tmp = (float *)&abProd;
		circularBuffer_write(ringProdReal, tmp, 1);
		circularBuffer_write(ringProdImag, tmp + 1, 1);
		*realpt += (*tmp);
		*imagpt += (*(tmp+1));
		abDot = *(lv_32fc_t*)(&res[0]);
			
		tmp = &abSum;
		circularBuffer_write(ringSum, tmp, 1);
		dotProd += *tmp;
		norm = dotProd;
	}

	else{
		circularBuffer_read(ringProdReal, &read, 1);
		*realpt -= read;
		circularBuffer_read(ringProdImag, &read, 1);
		*imagpt -= read;
		tmp = (float *)&abProd;
		circularBuffer_write(ringProdReal, tmp, 1);
		circularBuffer_write(ringProdImag, tmp + 1, 1);
		*realpt += (*tmp);
		*imagpt += (*(tmp+1));
		abDot = *(lv_32fc_t*)(&res[0]); 
		circularBuffer_read(ringSum, &read, 1);
		dotProd -= read;		
	
		tmp = &abSum;
		circularBuffer_write(ringSum, tmp, 1);
		dotProd += *tmp;
		norm = dotProd;
	}

	tmp = (float *)&abDot;
	*cPtr = ((*tmp) * (*tmp)) + ((*(tmp+1)) * (*(tmp+1)));
	*cPtr = (*cPtr) / norm;
	cPtr++;
	}	// end of for loop

	//circularBuffer_destroy(ringProdReal);
	//circularBuffer_destroy(ringProdImag);
	//circularBuffer_destroy(ringSum);
}
#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_32fc_x2_s32u_ofdm_metric_32f_a_H */
