#include "main.h"
#include <ap_int.h>
#include <ap_fixed.h> // Using fixed-point arithmetic for resource optimization
#include <hls_stream.h> // input size is unknown at compile time, this allows for data to process
#include <iostream>  // For debug prints
#include "hls_math.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>
#include "fit.h"

using namespace std;
			
int main() {
	// Declare streams
    	hls::stream<int> x_stream, y_stream, y_unc_stream, last_stream;
    	
	for (int i = 0; i < tempN; i++) {
        	x_stream.write(tempxs[i]);
        	y_stream.write(tempys[i]);
        	y_unc_stream.write(tempsigmas[i]);
        	last_stream.write(templasts[i]);
    }
	
	fit(x_stream, y_stream, y_unc_stream, last_stream, tempN);
    return 0;
}

