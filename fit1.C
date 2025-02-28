#include "fit.h"
#include "hls_math.h"
#include <iostream>

extern "C" {
	void fit(hls::stream<int>& x_stream, hls::stream<int>& y_stream, hls::stream<int>& y_unc_stream, hls::stream<int>& last_stream, int size){
        	#pragma HLS INTERFACE axis port=x_stream
        	#pragma HLS INTERFACE axis port=y_stream
        	#pragma HLS INTERFACE axis port=y_unc_stream
        	#pragma HLS INTERFACE axis port=last_stream
        	#pragma HLS PIPELINE II=1

		float sx = 0;
		float sy = 0;
		float s = 0;
		float st = 0;
		float sty = 0;
		float sum_chi = 0; //sum needed to calculate chi2
		int event_size = 0;
		int index = 0;
		for (int i = 0; i < size; i++) {
			// Read from all input streams
			int x = x_stream.read();
			int y = y_stream.read();
			int y_unc = y_unc_stream.read();
			int last = last_stream.read();

			// calculations
			float y_unc_sq = y_unc * y_unc;
			float inv_y_unc_sq = 1.0f / y_unc_sq;
			if (last == 0) {
				sx += x * inv_y_unc_sq;
				sy += y * inv_y_unc_sq;
				s += inv_y_unc_sq;
				float ti = (1.0f / y_unc) * (x - sx/s);
				st += ti * ti;
				sty += ti * y / y_unc;
				event_size++;
				float residual = (y - (sy - sx * sty /st) / s - (sty / st) * x) /y_unc;
				sum_chi += residual * residual;
			}else if (last == 1){
				sx += x * inv_y_unc_sq;
				sy += y * inv_y_unc_sq;
				s += inv_y_unc_sq;
				float ti = (1.0f / y_unc) * (x - sx/s);
				st += ti * ti;
				sty += ti * y / y_unc;
				event_size++;

				float residual = (y - (sy - sx * sty / st) / s - (sty / st) * x) / y_unc;
				sum_chi += residual * residual;

				// calculate parameters
				float b = (1.0f / st) * sty;
				float a = (sy - sx * b)/s;
				float sig_a = hls::sqrt((1.0f / s) * (1.0f + ((sx * sx) / (s * st))));
				float sig_b = hls::sqrt(1.0f/st);
				float chi2 = (1.0f/(event_size-2)) * sum_chi;

				// print results (for simulation only)
				std::cout << "For event " << index << ":\n";
				std::cout << "  a = " << a << " +/- " << sig_a << "\n";
				std::cout << "  b = " << b << " +/- " << sig_b << "\n";
				std::cout << "  chi2/ndf = " << chi2 << "\n\n";

				// reset everything
				float sx = 0;
				float sy = 0;
				float s = 0;
				float st = 0;
				float sty = 0;
				float sum_chi = 0; //sum needed to calculate chi2
				event_size = 0;
				index++;
			}
		}
	}
}
