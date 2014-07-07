/* -*- c++ -*- */

#define FASTBLOCKS_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "fastblocks_swig_doc.i"

%{
#include "fastblocks/ofdm_metric.h"
%}


%include "fastblocks/ofdm_metric.h"
GR_SWIG_BLOCK_MAGIC2(fastblocks, ofdm_metric);
