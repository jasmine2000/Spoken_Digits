#include <math.h>
#include <stdio.h>
#include "fft4g.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))

/* random number generator, 0 <= RND < 1 */
#define RND(p) ((*(p) = (*(p) * 7141 + 54773) % 259200) * (1.0 / 259200.0))

#ifndef NMAX
#define NMAX 8192
#define NMAXSQRT 64
#endif

// #define audio_length 5888
#define audio_length 256
#define frame_size 256

const int spect_time = audio_length / frame_size; // 23
const int spect_y = frame_size / 2; // 128

double frame[frame_size];

void stft(double audio[audio_length], double spectogram[spect_time * spect_y]);

int main()
{
    printf("\n");
    double spectogram[spect_time * spect_y];
    stft(frame, spectogram);

    int i;
    for (i = 0; i < spect_y; i++) {
        printf("%f\n", spectogram[i]);
    }
}


void stft(double audio[audio_length], double spectogram[spect_time * spect_y]) {
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

    double val, max_val;
    max_val = 0;
    for (i = 0; i < spect_time; i++) {
        start_audio = i * 256;
        for (j = 0; j < frame_size; j++) { // copy section into array a
            a[j] = audio[start_audio + j];
        }

        rdft(frame_size, 1, a, ip, w);      // do fft

        start_spect = i * spect_y;
        val = fabs(a[0]);
        spectogram[start_spect] = val; // a[1] is some weird number, ignore

        if (val > max_val) {
            max_val = val;
        }
        for (j = 1; j < spect_y; j++) {
            val = sqrt(pow(a[2 * j], 2) + pow(a[2 * j + 1], 2));    // magnitude of real and imaginary vector at point
            // val = pow(a[2 * j], 2) + pow(a[2 * j + 1], 2);
            spectogram[start_spect + j] = val;

            if (val > max_val) {
                max_val = val;
            }
        }
    }
    double amin = 1e-10;
    double log_max = log10(max_val);

    double new_val;
    double sum;     // to calculate mean after
    for (i = 0; i < total_length; i++) {
        val = spectogram[i];
        if (val < amin) val = amin;

        new_val = log10(val) - log_max;
        sum = sum + new_val;

        spectogram[i] = new_val;
    }

    double mean = sum / total_length;
    
    // calculate std of spectogram
    sum = 0;    // reuse sum var
    for (i = 0; i < total_length; i++) {
        sum += pow((spectogram[i] - mean), 2);
    }
    double std = sqrt(sum / total_length);
    
    for (i = 0; i < total_length; i++) {
        val = spectogram[i];
        spectogram[i] = (val - mean) / std;
    }
}

