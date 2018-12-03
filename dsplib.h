/*
   dsplib.h - plik nagłówkowy biblioteki DSP "dsplib"
*/
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define TABLEN(x) (sizeof(x)/sizeof(x[0]))

/*
  Typ danych pomocny przy zwracaniu z funkcji
  tablicy utworzonej dynamicznie, oraz jej długości
*/
typedef struct {
    double * data;
    long length;
}Signal;


/*
  Typ danych pomocny przy zwracaniu z funkcji
  dymacznej tablicy liczb zespolonych
*/
typedef struct {
    double complex * data;
    long length;
}complexSignal;

typedef enum {
    HORIZONTAL_INT,
    HORIZONTAL_FP,
    VERTICAL_INT,
    VERTICAL_FP
} printStyle;

typedef enum {
    INT_TXT,
    INT_BIN,
    FP_TXT,
    FP_BIN
}writeStyle;

typedef enum{
    FALSE = 0,
    TRUE  = 1
}BOOLEAN;

typedef enum {  // Type of normalization in correlation function
    WITHOUT,    // Without any normalization
    BIASED,     // Biased estimate of the cross-correlation
    UNBIASED,   // Unbiased estimate of the cross-correlation
    COEFF       //  Normalizes the sequence so that the autocorrelations at zero lag equal 1
}corNormalizationType;

// GNUPLOT

// Useful functions
void tabcpy(const double * source, double * target, long len);
void printTab1D(const double tab[], long len, printStyle style);
void printComplexTab(const double complex tab[], long len);
BOOLEAN saveSignalInFile(const char * filePath, const double sig[], long len, writeStyle style); // TODO
double * readSignalFromFIle(const char * filePath);
Signal frequency_vector(long len, double fs);
Signal dft_abs(const double complex tab[], long len);

// Signal parameters
double getMax(const double tab[], long len);
double getMin(const double tab[], long len);
double mean(const double tab[], long len);
double energy(const double tab[], long len);
double power(const double tab[], long len);
double rms(const double tab[], long len);
double variance(const double tab[], long len);
double stddev(const double tab[], long len); // standard deviation

// DSP functions
Signal conv(const double x[], long xlen, const double y[], long ylen);
Signal xcorr(const double x[], long xlen,  const double y[], long ylen, corNormalizationType normType);
Signal random_noise(long length);                                // do poprawy
Signal normal_noise(long length, double mean, double variance);  // do poprawy
complexSignal dft(const double x[], long xlen);
