/*
   main.c - plik testuj¹cy funkcje biblioteki DSP "dsplib"
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "dsplib.h"

#define TEST_SPLOTU_1 0
#define TEST_SPLOTU_2 0
#define TEST_SPLOTU_3 0
#define TEST_ZAPISU_DO_PLIKU 0
#define TEST_ODCZYTU_Z_PLIKU 0
#define TEST_KORELACJI 0
#define TEST_PARAMETROW 0
#define TEST_SZUMOW 0
#define TEST_DFT 1


int main()
{
    #if TEST_SPLOTU_1

    double x[] = {1,2};
    double y[] = {1,2};

    Signal conv1 = conv(x, TABLEN(x), y, TABLEN(y));

    printf("X: ");      printTab1D(x, TABLEN(x), HORIZONTAL_INT); printf("\n");
    printf("Y: ");      printTab1D(y, TABLEN(y), HORIZONTAL_INT); printf("\n");
    printf("convXY: "); printTab1D(conv1.data,conv1.length,HORIZONTAL_INT); printf("\n\n");


    #endif // TEST_SPLOTU_1

    #if TEST_SPLOTU_2

    double xx[] = {1,2,3,4};
    double yy[] = {1,2,3};

    Signal conv2 = conv(xx, TABLEN(xx), yy, TABLEN(yy));

    printf("X: ");      printTab1D(xx, TABLEN(xx), HORIZONTAL_INT);         printf("\n");
    printf("Y: ");      printTab1D(yy, TABLEN(yy), HORIZONTAL_INT);         printf("\n");
    printf("convXY: "); printTab1D(conv2.data,conv2.length,HORIZONTAL_INT); printf("\n\n");

    #endif // TEST_SPLOTU_2

    #if TEST_SPLOTU_3

    double xxx[] = {1,2,3,4};
    double yyy[] = {1,2,3};

    Signal conv3 = conv(xxx, TABLEN(xxx), yyy, TABLEN(yyy));

    printf("X: ");      printTab1D(xxx, TABLEN(xxx), HORIZONTAL_FP);       printf("\n");
    printf("Y: ");      printTab1D(yyy, TABLEN(yyy), HORIZONTAL_FP);       printf("\n");
    printf("convXY: "); printTab1D(conv3.data,conv3.length,HORIZONTAL_FP); printf("\n\n");

    #endif // TEST_SPLOTU_3

    #if TEST_ZAPISU_DO_PLIKU

    double longTab[] = {1,2,3,4,5,6,7,8,9,10};
    double doubleTab[] = {1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,10.1};

    char * file1 = "dane_calkowite.txt";
    char * file2 = "dane_rzeczywiste.txt";
    char * file3 = "dane_calkowite_binarne.bin";

    saveSignalInFile(file1, longTab, 10, INT_TXT);
    saveSignalInFile(file2, doubleTab, 10, FP_TXT);
    saveSignalInFile(file3, longTab, 10, INT_BIN);

    #endif // TEST_ZAPISU_DO_PLIKU

    #if TEST_ODCZYTU_Z_PLIKU

    printTab1D(readSignalFromFIle("test_odczytu.txt"), 10, VERTICAL_FP);
    printf("\n");
    printTab1D(readSignalFromFIle("dane_rzeczywiste.txt"), 10, HORIZONTAL_FP);
    printf("\n");

    #endif // TEST_ODCZYTU_Z_PLIKU

    #if TEST_KORELACJI

    //double xc[] = {1,2,3,4,5,6};
    //double yc[] = {3,2,1};

    //double xc[] = {3,2,1};
    //double yc[] = {1,2,3,4,5,6};

    double xc[] = {1,2,3,4,5,6};
    double yc[] = {3,2,1,1,2,3};

    Signal corrXY;

    printf("X: ");  printTab1D(xc, TABLEN(xc), HORIZONTAL_INT);  printf("\n");
    printf("Y: ");  printTab1D(yc, TABLEN(yc), HORIZONTAL_INT);  printf("\n");

    corrXY = xcorr(xc, TABLEN(xc), yc, TABLEN(yc), WITHOUT);
    printf("corrXY  [WITHOUT]: "); printTab1D(corrXY.data,corrXY.length,HORIZONTAL_INT); printf("\n");
    free(corrXY.data);

    corrXY = xcorr(xc, TABLEN(xc), yc, TABLEN(yc), BIASED);
    printf("corrXY   [BIASED]: "); printTab1D(corrXY.data,corrXY.length,HORIZONTAL_FP); printf("\n");
    free(corrXY.data);

    corrXY = xcorr(xc, TABLEN(xc), yc, TABLEN(yc), UNBIASED);
    printf("corrXY [UNBIASED]: "); printTab1D(corrXY.data,corrXY.length,HORIZONTAL_FP); printf("\n\n");
    free(corrXY.data);


    #endif // TEST_KORELACJI

    #if TEST_PARAMETROW

    double forParams[] = {2, 2, 2, 2, 4, 4, 4, 4};

    printf("Sygnal: ");
    printTab1D(forParams, TABLEN(forParams), HORIZONTAL_INT);
    printf("\nWart. srednia: %f", mean(forParams, TABLEN(forParams)));
    printf("\nEnergia: %f", energy(forParams, TABLEN(forParams)));
    printf("\nMoc: %f", power(forParams, TABLEN(forParams)));
    printf("\nRMS: %f", rms(forParams, TABLEN(forParams)));
    printf("\n");


    #endif // TEST_PARAMETROW

    #if TEST_SZUMOW

    Signal noise = random_noise(10);
    Signal normal = normal_noise(9,0,1);

    printf("\n\n");
    printf("Szum o rozkladzie rownomiernym: \n");
    printTab1D(noise.data, noise.length, VERTICAL_FP);
    printf("\n\n");

    printf("\n\n");
    printf("Szum o rozkladzie normalnym: \n");
    printTab1D(normal.data, normal.length, VERTICAL_FP);
    printf("\n\n");


    free(noise.data);
    free(normal.data);

    #endif // TEST_SZUMOW

    #if TEST_DFT

    clock_t start, stop;
    double execution_time;
    long fs = 100;
    double start_time = 0;
    double end_time = 1;
    long signal_length =  (long) (end_time*fs - start_time*fs);

    double * time = (double *) calloc(signal_length, sizeof(* time));
    double * sig = (double *) calloc(signal_length, sizeof(* sig));
    double s1, s2 ,s3;
    long i;

    // Utworzenie wektora czasu i sygnału
    for(i = 0; i < signal_length; i++){

        time[i] = (double) i/fs;

        s1 = sin(2*M_PI*10*time[i]);
        s2 = 3*sin(2*M_PI*20*time[i] + M_PI/2);
        s3 = 2*sin(2*M_PI*40*time[i] + M_PI/3);

        sig[i] = s1 + s2 + s3;
    }

    start = clock();
    complexSignal sig_dft = dft(sig, signal_length);
    stop = clock();

    execution_time = (double) (stop - start) / CLOCKS_PER_SEC;

    printf("\n");
    printf("DFT\n\n");
    printf("Signal length: %ld\n", signal_length);
    printf("Czas obliczenia DFT: %f\n\n", execution_time);

    printf("DFT utworzonego sygnalu:\n");
    printComplexTab(sig_dft.data, sig_dft.length);
    printf("\nABS z DFT utworzonego sygnalu:\n");
    Signal abs_dft = dft_abs(sig_dft.data,sig_dft.length);
    printTab1D(abs_dft.data, abs_dft.length,VERTICAL_FP);

    free(time);
    free(sig);
    free(sig_dft.data);

    #endif // TEST_DFT

    return 0;
}
