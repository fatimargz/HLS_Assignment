#include "fit.h"
#include "hls_math.h"
#include <iostream>
#include <ap_int.h>
#include <ap_fixed.h> // Using fixed-point arithmetic for resource optimization
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include <iostream>  // For debug prints
#include <bitset>

extern "C" {
	void fit(hls::stream<int>& x_stream, hls::stream<int>& y_stream, hls::stream<int>& y_unc_stream, hls::stream<int>& last_stream, int size){

        #pragma HLS INTERFACE axis port=x_stream
        #pragma HLS INTERFACE axis port=y_stream
        #pragma HLS INTERFACE axis port=y_unc_stream
        #pragma HLS INTERFACE axis port=last_stream
	    #pragma HLS INTERFACE ap_none port=size
        #pragma HLS PIPELINE II=1
		
		#define MAX 500
		ap_fixed<32, 16> sx = 0;
		ap_fixed<32, 16> sy = 0;
		ap_fixed<32, 16> s = 0;
		ap_fixed<32, 16> st = 0;
		ap_fixed<32, 16> sty = 0;
		ap_uint<32> index = 0;
		ap_fixed<32, 16> x_event[MAX];
		ap_fixed<32, 16> y_event[MAX];
		ap_fixed<32, 16> y_unc_event[MAX];
		ap_uint<32> event_size = 0;
		ap_uint<32> event_index =0;

		for (int i = 0; i < size; i++) {
			#pragma HLS PIPELINE II=1
			// Read from all input streams
			ap_fixed<32, 16> x = x_stream.read();
			ap_fixed<32, 16> y = y_stream.read();
			ap_fixed<32, 16> y_unc = y_unc_stream.read();
			ap_fixed<32, 16> last = last_stream.read();

			// calculations
			ap_fixed<32, 16> y_unc_sq = y_unc * y_unc;
			ap_fixed<32, 16> inv_y_unc_sq = ap_fixed<32,16>(1.0) / y_unc_sq;
			if (last == 0) {
				sx += x * inv_y_unc_sq;
				sy += y * inv_y_unc_sq;
				s += inv_y_unc_sq;
				x_event[index] = x;
				y_event[index] = y;
				y_unc_event[index] = y_unc;
				index++;
				event_size++;
			}else if (last == 1){
				sx += x * inv_y_unc_sq;
				sy += y * inv_y_unc_sq;
				s += inv_y_unc_sq;
				x_event[index] = x;
                y_event[index] = y;
				y_unc_event[index] = y_unc;
				event_size++;

				// a, b, a_err, b_err, chi2 calculations

				#pragma HLS ARRAY_PARTITION variable=x_event cyclic factor=2
				#pragma HLS ARRAY_PARTITION variable=y_event cyclic factor=2
				#pragma HLS ARRAY_PARTITION variable=y_unc_event cyclic factor=2

				
				for (int j = 0; j < event_size ; j++){
					ap_fixed<32, 16> inv_y_unc = (ap_fixed<32,16>(1.0)/y_unc_event[j]);
					ap_fixed<32,16> sx_s = sx/s;
					ap_fixed<32, 16> ti = (inv_y_unc * (x_event[j] - sx_s));
					st += ti * ti;
					ap_fixed<32, 16> tiy = ti*y_event[j];
					sty += tiy * inv_y_unc;
				}
				// dubugging statements for simulation. something is wrong with sty
				//std::cout << "sx: " << sx << " sy: " << sy << " s: " << s << std::endl;
				//std::cout << "st: " << st << " sty: " << sty << std::endl;

				ap_fixed<32, 16> inv_st = (ap_fixed<32,16>(1.0)/st);
				ap_fixed<32, 16> b = inv_st * sty;
				ap_fixed<32, 16> a = (sy - sx * b)/s;
				ap_fixed<32, 16> sig_a_sq = (ap_fixed<32,16>(1.0) / s) * (ap_fixed<32,16>(1.0) + ((sx * sx) / (s * st)));
				ap_fixed<32, 16> sig_a = hls::sqrt(sig_a_sq);
				ap_fixed<32, 16> sig_b = hls::sqrt(ap_fixed<32,16>(1.0)/st);

				ap_fixed<32, 16> sr = 0;

				//#pragma HLS UNROLL factor=2 - talked to hayden he said this is useless
				for (int k = 0; k < event_size; k++){
					ap_fixed<32, 16> residuals = (y_event[k] - a - b * x_event[k]) / y_unc_event[k];
					sr += residuals * residuals;
				}

				ap_fixed<32, 16> chi2 = (ap_fixed<32,16>(1.0)/(event_size-2)) * sr;

				// print results (for simulation only)
				#ifndef __SYNTHESIS__
				std::cout << "For event " << event_index << ":\n";
				std::cout << "  a = " << a << " +/- " << sig_a << "\n";
				std::cout << "  b = " << b << " +/- " << sig_b << "\n";
				std::cout << "  chi2/ndf = " << chi2 << "\n\n";
				#endif

				// reset everything
				sx = 0;
				sy = 0;
				s = 0;
				st = 0;
				sty = 0;
				event_size = 0;
				index = 0;
				event_index++;
			}
		}
	}
}
