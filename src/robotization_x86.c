/*
Author: Kim Radmacher

Date: 03.09.2022

Description:
  FFT based robotization. Algorithm is based on DAFX - U. Zoelzer page 287 ff.
  Project is deployed on x86
*/

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

//#define DEBUG

#if defined(_MSC_VER)
#include <getopt.h>
#else
#include <unistd.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "wavreader.h"
#include "wavwriter.h"
#ifdef __cplusplus
}
#endif

#include "kissfft/kiss_fft.h"

#define BUFFER_SIZE (16384)

void cleanup(void);

// set the frequency of the oscillators
int16_t gInputBuffer[BUFFER_SIZE];
int gInputBufferPointer = 0;
int16_t gOutputBuffer[BUFFER_SIZE];
int32_t inL, inR, outL, outR;
int gOutputBufferWritePointer = 0;
int gOutputBufferReadPointer = 0;
int gSampleCount = 0;
float modPhase = 0.f; // good value is -180...0...180

// These variables used internally in the example:
static const int gFFTSize = 10;
// phase parameters
static const int Ha = 441;//(1<<gFFTSize)/4-0; /* analysis hopsize */
static const int Hs = 441;//(1<<gFFTSize)/4;   /* synthisis hopsize */
int gHopSize = Hs;
int gPeriod = gHopSize;
float gFFTScaleFactor = 0;
float *gWindowBuffer;

// FFT vars
kiss_fft_cpx* timeDomainIn;
kiss_fft_cpx* timeDomainOut;
kiss_fft_cpx* frequencyDomain;

kiss_fft_cfg fft_cfg;
kiss_fft_cfg ifft_cfg;

// phase processing vars
float *amplitude;
float *phi;

// Set the analog channels to read from
const int gAnalogIn = 0;
int gAudioFramesPerAnalogFrame = 0;
int gFFTInputBufferPointer = 0;
int gFFTOutputBufferPointer = 0;

bool setup(void)
{
    fft_cfg = kiss_fft_alloc(1<<gFFTSize, 0, (void*)0, 0);
    ifft_cfg = kiss_fft_alloc(1<<gFFTSize, 1, (void*)0, 0);

    gFFTScaleFactor = 1.0f / (float)(1<<gFFTSize);
    gOutputBufferWritePointer = Hs;
    gOutputBufferReadPointer = 0;

    timeDomainIn = (kiss_fft_cpx*) malloc ((1<<gFFTSize) * sizeof (kiss_fft_cpx));
    if(timeDomainIn == NULL)
    {
        return false;
    }
    timeDomainOut = (kiss_fft_cpx*) malloc ((1<<gFFTSize) * sizeof (kiss_fft_cpx));
    if(timeDomainOut == NULL)
    {
        return false;
    }
    frequencyDomain = (kiss_fft_cpx*) malloc ((1<<gFFTSize) * sizeof (kiss_fft_cpx));
    if(frequencyDomain == NULL)
    {
        return false;
    }
    
    // list of fft c-implementations
    // https://community.vcvrack.com/t/complete-list-of-c-c-fft-libraries/9153

    memset(timeDomainIn, 0, (1<<gFFTSize) * sizeof (kiss_fft_cpx));
    memset(timeDomainOut, 0, (1<<gFFTSize) * sizeof (kiss_fft_cpx));

    memset(gOutputBuffer, 0, BUFFER_SIZE * sizeof(int16_t));
    memset(gInputBuffer, 0, BUFFER_SIZE * sizeof(int16_t));
    
    // Allocate processing buffer and init vars
    amplitude = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(amplitude == NULL)
    {
        return false;
    }
    phi = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(phi == NULL)
    {
        return false;
    }
    // Allocate the window buffer based on the FFT size
    gWindowBuffer = (float *)malloc((1<<gFFTSize) * sizeof(float));

    if(gWindowBuffer == NULL)
    {
        return false;
    }

    // Calculate a Hann window
    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        gWindowBuffer[n] = 0.5f * (1.0f - cosf(2.0f * M_PI * (float)n / (float)((1<<gFFTSize))));
    }

    return true;

}

// This function handles the FFT based robotization
void process_robotization(int16_t *inBuffer, int inWritePointer, int16_t *outBuffer, int outWritePointer)
{
    // Copy buffer into FFT input
    int pointer = (inWritePointer - (1<<gFFTSize) + BUFFER_SIZE) % BUFFER_SIZE;

    // interleave real imag parts
    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        timeDomainIn[n].r = (int16_t)((float)inBuffer[pointer] * gWindowBuffer[n]);
        timeDomainIn[n].i = 0;

        pointer++;
        if(pointer >= BUFFER_SIZE)
        {
            pointer = 0;
        }
    }
    
    kiss_fft(fft_cfg, timeDomainIn , frequencyDomain);

    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        amplitude[n] = sqrtf((float)frequencyDomain[n].r * frequencyDomain[n].r + (float)frequencyDomain[n].i * frequencyDomain[n].i);
        phi[n] = (modPhase)*(float)n;
    }

    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        frequencyDomain[n].r = (1<<(gFFTSize-1))*(int16_t)(cosf(phi[n])*amplitude[n]);
        frequencyDomain[n].i = (1<<(gFFTSize-1))*(int16_t)(cosf(phi[n])*amplitude[n]);
    }

    //ne10_fft_c2c_1d_float32_neon (timeDomainOut, frequencyDomain, cfg, 1);
    kiss_fft(ifft_cfg, frequencyDomain, timeDomainIn);

    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        timeDomainIn[n].r = (int16_t)((float)timeDomainIn[n].r * gWindowBuffer[n]);
    }

#ifdef DEBUG
    float absMax = 0.f;
    uint32_t ueqZero = 0;
#endif
    // Overlap-and-add timeDomainOut into the output buffer
    pointer = outWritePointer;

    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        outBuffer[pointer] += (int16_t)((float)timeDomainIn[n].r);

#ifdef DEBUG
        if(timeDomainIn[n].i != 0)
        {
            absMax = (absMax < fabs(timeDomainIn[n].i)) ? fabs(timeDomainIn[n].i) : absMax;
            ueqZero++;
        }

        if(isnan(outBuffer[pointer]))
        {
            printf("outBuffer OLA\n");
        }
#endif
        pointer++;
        if(pointer >= BUFFER_SIZE)
        {
            pointer = 0;
        }
    }

#ifdef DEBUG
    if(ueqZero)
    {
        printf("WARNING: timeDomainIn[N].i not zero %d times (max = %f)\n", ueqZero, absMax);
    }
#endif
}

void usage(const char* name)
{
    fprintf(stderr, "%s in.wav out.wav <-200...200>\n", name);
}

int main(int argc, char *argv[])
{
    const char *infile, *outfile;
    FILE *out;
    void *wavIn;
    void *wavOut;
    int format, sample_rate, channels, bits_per_sample;
    uint32_t data_length;
    int input_size;
    uint8_t* input_buf;
    int16_t* convert_buf;

    if(!setup())
    {
        fprintf(stderr, "setup failed\n");
    }
    
    if (argc - optind < 2)
    {
        fprintf(stderr, "Error: not enough parameter provided\n");
        usage(argv[0]);
        return 1;
    }
    
    infile = argv[optind];
    outfile = argv[optind + 1];

    if (argc - optind > 2)
    {
        modPhase = atof(argv[optind + 2]);
        if(modPhase > 180) modPhase = 180;
        if(modPhase < -180) modPhase = -180;
    }

    wavIn = wav_read_open(infile);
    if (!wavIn)
    {
        fprintf(stderr, "Unable to open wav file %s\n", infile);
        return 1;
    }
    if (!wav_get_header(wavIn, &format, &channels, &sample_rate, &bits_per_sample, &data_length))
    {
        fprintf(stderr, "Bad wav file %s\n", infile);
        return 1;
    }
    if (format != 1)
    {
        fprintf(stderr, "Unsupported WAV format %d\n", format);
        return 1;
    }

    wavOut = wav_write_open(outfile, sample_rate, bits_per_sample, channels);

    if (!wavOut)
    {
        fprintf(stderr, "Unable to open wav file for writing %s\n", infile);
        return 1;
    }

    input_size = data_length;
    input_buf = (uint8_t*) malloc(input_size);
    convert_buf = (int16_t*) malloc(input_size);

    if (input_buf == NULL || convert_buf == NULL)
    {
        fprintf(stderr, "Unable to allocate memory for buffer\n");
        return 1;
    }

    int read = wav_read_data(wavIn, input_buf, input_size);

    printf("using modPhase = %f\n", modPhase);
    
    printf("data_length = %d\tread = %d\tinput_size = %d \n", data_length, read, input_size);
    printf("sample_rate = %d\tbits_per_sample = %d\tchannels = %d \n", sample_rate, bits_per_sample, channels);

    int numSamples = read/2;
    for(unsigned int n = 0; n < numSamples; n++)
    {
        const uint8_t* in = &input_buf[2*n];
        convert_buf[n] = in[0] | (in[1] << 8);
    }
    
    // iterate over the audio frames and create three oscillators, seperated in phase by PI/2
    for(unsigned int n = 0; n < numSamples; n+=channels)
    {
        // Read audio inputs
        if(channels == 1)
        {
            inL = (int32_t)convert_buf[n];
            inR = inL;
        }
        else if(channels == 2)
        {
            // interleaved left right channel
            inL = (int32_t)convert_buf[n];
            inR = (int32_t)convert_buf[n+1];
        }
        else
        {
            fprintf(stderr, "channel = %d\n", channels);
            return -1;
        }

#if 1 // apply robotization if defined. otherwise it's bypassed
        gInputBuffer[gInputBufferPointer] = (int16_t)((inR+inL)/2);

        outL = gOutputBuffer[gOutputBufferReadPointer];
        outR = outL;
        
        // Clear the output sample in the buffer so it is ready for the next overlap-add
        gOutputBuffer[gOutputBufferReadPointer] = 0;
        
        gOutputBufferReadPointer++;
        if(gOutputBufferReadPointer >= (BUFFER_SIZE))
        {
            gOutputBufferReadPointer = 0;
        }   
        gOutputBufferWritePointer++;
        if(gOutputBufferWritePointer >= (BUFFER_SIZE))
        {
            gOutputBufferWritePointer = 0;
        }
        gInputBufferPointer++;
        if(gInputBufferPointer >= (BUFFER_SIZE))
        {
            gInputBufferPointer = 0;
        }    
        gSampleCount++;
        if(gSampleCount >= Hs)
        {
            process_robotization(gInputBuffer, gInputBufferPointer, gOutputBuffer, gOutputBufferWritePointer);
            gSampleCount = 0;
        }    
        
#else
        outL = inL;
        outR = inR;
#endif        

        int16_t oL = (int16_t)outL;
        int16_t oR = (int16_t)outR;
        wav_write_data(wavOut, (unsigned char*)&oL, 2);
        if(channels > 1)
        {
            wav_write_data(wavOut, (unsigned char*)&oR, 2);
        }
    }    

    free(convert_buf);
    free(input_buf);
    
    cleanup();

    wav_write_close(wavOut);
    wav_read_close(wavIn);

    return 0;
}
// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup()
{

    kiss_fft_free(fft_cfg);

    free(timeDomainIn);
    free(timeDomainOut);
    free(frequencyDomain);

    free(gWindowBuffer);
    free(phi);
    free(amplitude);
}
