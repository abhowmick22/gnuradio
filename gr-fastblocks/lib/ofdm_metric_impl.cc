/* -*- c++ -*- */
/* 
 * Copyright 2014 Free Software Foundation, Inc.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "ofdm_metric_impl.h"
#include <volk/volk.h>
#include <assert.h>

namespace gr {
  namespace fastblocks {

    ofdm_metric::sptr
    ofdm_metric::make(int fft_len, int cp_len, bool use_even_carriers)
    {
      return gnuradio::get_initial_sptr
        (new ofdm_metric_impl(fft_len, cp_len, use_even_carriers));
    }

    /*
     * The private constructor
     */
    ofdm_metric_impl::ofdm_metric_impl(int fft_len, int cp_len, bool use_even_carriers)
      : gr::block("ofdm_metric",
              gr::io_signature::make(2, 2, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(float))),
		d_fft_len(fft_len),
		d_cp_len(cp_len)
    {
	// TODO : Check its purpose
      const int alignment_multiple =
		volk_get_alignment() / sizeof(float);
      set_alignment(std::max(1, alignment_multiple));
	}

    /*
     * Our virtual destructor.
     */
    ofdm_metric_impl::~ofdm_metric_impl()
    {
    }


// TODO : Rely on the default forecast function

    void
    ofdm_metric_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items; 
        ninput_items_required[1] = noutput_items; 
    }

 

    int
    ofdm_metric_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
		assert(input_items.size() == 2 && "More than 2 input vectors!");
        const gr_complex *in1 = (const gr_complex *) input_items[0];
        const gr_complex *in2 = (const gr_complex *) input_items[1];
        float *out = (float *) output_items[0];
		unsigned int noi = noutput_items;
		
        // Do <+signal processing+>
		
		volk_32fc_x2_s32u_ofdm_metric_32f(out, in1, in2, d_fft_len, noi);
        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace fastblocks */
} /* namespace gr */

