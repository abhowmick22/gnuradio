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


#ifndef INCLUDED_FASTBLOCKS_OFDM_METRIC_H
#define INCLUDED_FASTBLOCKS_OFDM_METRIC_H

#include <fastblocks/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace fastblocks {

    /*!
     * \brief Calculation of timing utility metric for OFDM frame detection
     * \ingroup fastblocks
     *
	 * \details
	 * Enter more details about the block here
     */
    class FASTBLOCKS_API ofdm_metric : virtual public gr::block
	// TODO : Check if this should be hier_block2 instead of gr::block
    {
     public:
      typedef boost::shared_ptr<ofdm_metric> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of fastblocks::ofdm_metric.
       *
       * To avoid accidental use of raw pointers, fastblocks::ofdm_metric's
       * constructor is in a private implementation
       * class. fastblocks::ofdm_metric::make is the public interface for
       * creating new instances.
       */
      static sptr make(int fft_len, int cp_len, bool use_even_carriers=false);
    };

  } // namespace fastblocks
} // namespace gr

#endif /* INCLUDED_FASTBLOCKS_OFDM_METRIC_H */

