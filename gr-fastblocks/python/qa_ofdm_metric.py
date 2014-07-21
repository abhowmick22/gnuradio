#!/usr/bin/env python
# 
# Copyright 2014 Free Software Foundation, Inc.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

import numpy
import random

from gnuradio import gr, gr_unittest, blocks, analog, channels
from gnuradio import digital
from gnuradio.digital.utils import tagged_streams
from gnuradio.digital.ofdm_txrx import ofdm_tx
import fastblocks_swig as fastblocks


class qa_ofdm_metric (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_detect (self):
        """ Send two bursts, with zeros in between, and check
        they are both detected at the correct position and no
        false alarms occur """
        n_zeros = 15
        fft_len = 32
        cp_len = 4
        sig_len = (fft_len + cp_len) * 10
        sync_symbol = [(random.randint(0, 1)*2)-1 for x in range(fft_len/2)] * 2
        tx_signal = [0,] * n_zeros + \
                    sync_symbol[-cp_len:] + \
                    sync_symbol + \
                    [(random.randint(0, 1)*2)-1 for x in range(sig_len)]
        tx_signal = tx_signal * 2
        add = blocks.add_cc()
        #sync = digital.ofdm_sync_sc_cfb(fft_len, cp_len)
        metric = fastblocks.ofdm_metric(fft_len, cp_len)
        detector = blocks.plateau_detector_fb(cp_len)
        #sink_freq   = blocks.vector_sink_f()
        sink_detect = blocks.vector_sink_b()
        self.tb.connect(blocks.vector_source_c(tx_signal), (add, 0))
        self.tb.connect(analog.noise_source_c(analog.GR_GAUSSIAN, .01), (add, 1))
        self.tb.connect(add, (metric, 0))
        self.tb.connect(add, (metric, 1))
        self.tb.connect(metric, detector)
        self.tb.connect(detector, sink_detect)
        self.tb.run ()
        sig1_detect = sink_detect.data()[0:len(tx_signal)/2]
        sig2_detect = sink_detect.data()[len(tx_signal)/2:]
        self.assertTrue(abs(sig1_detect.index(1) - (n_zeros + fft_len + cp_len)) < cp_len)
        self.assertTrue(abs(sig2_detect.index(1) - (n_zeros + fft_len + cp_len)) < cp_len)
        self.assertEqual(numpy.sum(sig1_detect), 1)
        self.assertEqual(numpy.sum(sig2_detect), 1)
        #self.assertEqual(1,1)



if __name__ == '__main__':
    gr_unittest.run(qa_ofdm_metric, "qa_ofdm_metric.xml")
