#include <math.h>
#include <stdio.h>
#include "fft4g_h.h"
#include "spectrogram.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))

/* random number generator, 0 <= RND < 1 */
#define RND(p) ((*(p) = (*(p) * 7141 + 54773) % 259200) * (1.0 / 259200.0))

#ifndef NMAX
#define NMAX 8192
#define NMAXSQRT 64
#endif

#define audio_length 5632
#define frame_size 256

const int spect_time = audio_length / frame_size; // 22
const int spect_y = frame_size / 2; // 128

double frame[256];


void stft(double audio[audio_length], float spectrogram[spect_time * spect_y]) {
    // does simplified stft
    // no overlap
    // no windowing function

    int total_length = spect_time * spect_y;

    // allocating space
    int ip[NMAXSQRT + 2];
    double w[NMAX * 5 / 4];

    // n_fft == hop_length for this stft (no overlap)
    int i, j, start_audio, start_spect;
    double a[frame_size];

    int thrown_out = 0;
    float amin = 1e-10;

    double val, max_val;
    max_val = 0;
    for (i = 0; i < spect_time; i++) {
        start_audio = (i - thrown_out) * 256;
        for (j = 0; j < frame_size; j++) { // copy section into array a
            a[j] = audio[start_audio + j];
        }

        rdft(frame_size, 1, a, ip, w);      // do fft

        start_spect = i * spect_y;

        for (j = 0; j < spect_y; j++) {
        	if (j == 0) val = fabs(a[0]); // a[1] is some weird number, ignore
        	else val = sqrt(pow(a[2 * j], 2) + pow(a[2 * j + 1], 2));    // magnitude of real and imaginary vector at point

            if (val != val) {
				thrown_out ++;
				break;
		    }
            else if (val > 10000) {
            	val = 20;
            }

            if (val > max_val) {
                max_val = val;
            }

            spectrogram[start_spect + j] = (float) val;
        }
    }


    float log_max = log10(max_val);

    int i2 = (spect_time - thrown_out) * 256;
    for (i = i2; i < total_length; i++) {
    	spectrogram[i] = amin;
    }

    float new_val;
    float mean = 0;     // to calculate mean after

    int i_start;
    double temp_sum;
    double temp_mean;
    for (i = 0; i < spect_time; i++) {
        temp_sum = 0;
        i_start = i * spect_y;
		for (j = 1; j < spect_y; j++) {
			val = spectrogram[i_start + j];

			if (val < amin || val != val) val = amin;			// if zero or nan

			new_val = log10(val) - log_max;

			if (new_val < -80) new_val = -80;

			spectrogram[i_start + j] = new_val;
			temp_sum = temp_sum + new_val;
		}
		temp_mean = (float) (temp_sum / spect_y);
		mean = mean + temp_mean;
	}
    mean = mean / spect_time;

    // calculate std of spectrogram
    double sum = 0;
    for (i = 0; i < total_length; i++) {
        sum = sum + fabs(spectrogram[i] - mean) / total_length;
    }
    float std = (float) sqrt(sum);

    for (i = 0; i < total_length; i++) {
        val = (spectrogram[i] - mean) / std;
        spectrogram[i] = val;
    }
}
