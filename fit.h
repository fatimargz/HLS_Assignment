#ifndef FIT_H
#define FIT_H

#include "hls_stream.h"
extern "C" {
    void fit(hls::stream<int>& x_stream, hls::stream<int>& y_stream, hls::stream<int>& y_unc_stream, hls::stream<int>& last_stream, int size);
}

#endif // FIT_H
