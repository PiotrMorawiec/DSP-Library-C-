/*
   dsplib.c - plik z  cia³ami funkcji nale¹cych do biblioteki DSP "dsplib"
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include "dsplib.h"

// -------------------------------------------------------------------------------------------------
void tabcpy(const double * source, double * target, long len){

    long i;

    for(i = 0; i < len; i++){
        target[i] = source[i];
    }
}

// -------------------------------------------------------------------------------------------------
void printTab1D(const double tab[], long len, printStyle style){

    int i;

    switch(style){

    case VERTICAL_FP:

        for(i = 0; i < len; i++){
            printf("%.6f\n",tab[i]);
        }
        break;

    case VERTICAL_INT:

        for(i = 0; i < len; i++){
            printf("%ld\n",(long)tab[i]);
        }
        break;

    case HORIZONTAL_FP:

        printf("{");
        for(i = 0; i < len; i++){
            if(i == len-1)  printf("%.6f",tab[i]);
            else            printf("%.6f, ",tab[i]);
        }
        printf("}");
        break;

    case HORIZONTAL_INT:

        printf("{");
        for(i = 0; i < len; i++){
            if(i == len-1)  printf("%ld",(long)tab[i]);
            else            printf("%ld, ",(long)tab[i]);
        }
        printf("}");
        break;
    }
}

// -------------------------------------------------------------------------------------------------
void printComplexTab(const double complex tab[], long len){

    long i;
    for(i = 0; i < len; i++){
        double re = creal(tab[i]);
        double im = cimag(tab[i]);
        printf("%.6f + %.6fi\n", re, im);
    }
}


// -------------------------------------------------------------------------------------------------
BOOLEAN saveSignalInFile(const char * filePath, const double sig[], long len, writeStyle style){

    FILE * file;

    switch(style){

    case INT_TXT:
        file = fopen(filePath,"w");
        if(!file){
            perror("Cannot to open the file");
            return FALSE;
        }
        else{
            int i;
            for(i = 0; i < len; i++){
                fprintf(file,"%ld\n",(long)sig[i]);
            }

            fclose(file);
            return TRUE;
        }
        break;

    case INT_BIN:
        file = fopen(filePath,"wb");
        if(!file){
            perror("Cannot to open the file");
            return FALSE;
        }
        else{
            int i;
            for(i = 0; i < len; i++){
                fprintf(file,"%ld\n",(long)sig[i]);
            }

            fclose(file);
            return TRUE;
        }
        break;

    case FP_TXT:
        file = fopen(filePath,"w");
        if(!file){
            perror("Cannot to open the file");
            return FALSE;
        }
        else{
            int i;
            for(i = 0; i < len; i++){
                fprintf(file,"%f\n",sig[i]);
            }

            fclose(file);
            return TRUE;
        }
        break;

    case FP_BIN:
        file = fopen(filePath,"wb");
        if(!file){
            perror("Cannot to open the file");
            return FALSE;
        }
        else{
            int i;
            for(i = 0; i < len; i++){
                fprintf(file,"%f\n",sig[i]);
            }

            fclose(file);
            return TRUE;
        }
        break;
    }
}

// -------------------------------------------------------------------------------------------------
double * readSignalFromFIle(const char * filePath){

    FILE * file;
    file = fopen(filePath, "r");
    if(!file){
        fprintf(stderr,"Cannot to open the file \"%s\"",filePath);
        return NULL;
    }
    else{

        long long  dataLength;
        long long  fileLength;

        /* Sprawdzenie długości pliku */
        fseek(file,0,SEEK_END);
        fgetpos(file,&fileLength);

        /*
           Długość pliku podana jest w bajtach.
           Dlugosc zmiennej rzeczywistej w pliku: 4
           Dlugosc znaku \n : 1

           Calkowita ilosc linii w pliku = (dlugosc_w_bajtach + 1)/ (4 + 1)
        */

        dataLength = (fileLength + 1) / 5;

        /* Utworzenie tablicy o odpowiedniej długości */
        double * data = (double *) malloc(dataLength * sizeof(*data));

        /* Odczyt danych z pliku */
        fseek(file,0,SEEK_SET);
        long long i = 0;
        while(1){
            fscanf(file,"%lf\n",&data[i++]);
            if(feof(file)) break;
        }

        fclose(file);
        return data;
    }
}

// -------------------------------------------------------------------------------------------------
Signal frequency_vector(long len, double fs){

    double * freq = (double *) calloc(len, sizeof(* freq));
    long i;

    for(i = 0; i < len; i++){
        freq[i] = i*fs/len;
    }

    Signal freq_vector;
    freq_vector.data = freq;
    freq_vector.length = len;

    return freq_vector;
}

// -------------------------------------------------------------------------------------------------
Signal dft_abs(const double complex tab[], long len){

    double * absolute_values = (double *) calloc(len, sizeof(* absolute_values));
    long i;

    for(i = 0; i < len; i ++){
        absolute_values[i] = cabs(tab[i]);
    }

    Signal abs_dft;
    abs_dft.data = absolute_values;
    abs_dft.length = len;

    return abs_dft;
}

// -------------------------------------------------------------------------------------------------
double getMax(const double tab[], long len){

    long idx;
    double max = tab[0];

    for(idx = 0; idx < len; idx++){
        if(tab[idx] > max) max = tab[idx];
    }

    return max;
}

// -------------------------------------------------------------------------------------------------
double getMin(const double tab[], long len){

    long idx;
    double min = tab[0];

    for(idx = 0; idx < len; idx++){
        if(tab[idx] < min) min = tab[idx];
    }

    return min;
}

// -------------------------------------------------------------------------------------------------
double mean(const double tab[], long len){

    long i;
    double sum = 0;

    for(i = 0; i < len; i++){
        sum += tab[i];
    }

    return sum/len;
}

// -------------------------------------------------------------------------------------------------
double energy(const double tab[], long len){

    long i;
    double en = 0;

    for(i = 0; i < len; i++){
        en += tab[i] * tab[i];
    }

    return en;
}

// -------------------------------------------------------------------------------------------------
double power(const double tab[], long len){

    return energy(tab, len)/len;
}

// -------------------------------------------------------------------------------------------------
double rms(const double tab[], long len){

    return sqrt(power(tab, len));
}

// -------------------------------------------------------------------------------------------------
double variance(const double tab[], long len){

    double avr = mean(tab, len);
    double sum = 0;
    long i;

    for(i = 0; i < len; i++){
        sum += pow((tab[i] - avr),2);
    }

    return sum/len;
}

// -------------------------------------------------------------------------------------------------
double stddev(const double tab[], long len){

    return sqrt(variance(tab, len));
}

// -------------------------------------------------------------------------------------------------
Signal conv(const double x[], long xlen, const double y[], long ylen){

    long convLen = xlen + ylen - 1;

    long i,j;

    // Utworzenie i zerowanie tablicy convolution
    double * convolution = (double *) calloc(convLen, sizeof(*convolution));

    long count_from = 0;
    long count_to = 0;

    for(i = 0; i < convLen; i++){

        for(j = count_from; j <= count_to; j++){
            convolution[i] = convolution[i] + x[j]*y[i - j];
        }

        if(i + 1 < xlen) count_to += 1;
        if(i + 1 > ylen - 1) count_from += 1;
    }

    Signal sig;
    sig.data = convolution;
    sig.length = convLen;

    return sig;
}

// -------------------------------------------------------------------------------------------------
Signal xcorr(const double x[], long xlen,  const double y[], long ylen, corNormalizationType normType){

    long len;
    long width;
    double remain;
    long i,j;
    char equalLengths = 0;
    double * xx;
    double * yy;

    // Uzupelnianie zerami krótszego z wektorów
    if(ylen > xlen){
        len = ylen;
        remain = ylen - xlen;

        xx = (double *) malloc(len * sizeof(*xx));
        yy = (double *) malloc(len * sizeof(*yy));

        tabcpy(y,yy,ylen);
        for(i = 0; i < xlen; i++)   xx[i] = x[i];
        for(i = 0; i < remain; i++) xx[xlen + i] = 0;
    }
    else if(xlen > ylen){
        len = xlen;
        remain = xlen - ylen;
        xx = (double *) malloc(len * sizeof(*xx));
        yy = (double *) malloc(len * sizeof(*yy));

        tabcpy(x,xx,xlen);
        for(i = 0; i < ylen; i++)   yy[i] = y[i];
        for(i = 0; i < remain; i++) yy[ylen + i] = 0;
    }
    else if(xlen == ylen){
        len = xlen;
        xx = (double *) malloc(len * sizeof(*xx));
        yy = (double *) malloc(len * sizeof(*yy));

        tabcpy(y,yy,ylen);
        tabcpy(x,xx,xlen);

        equalLengths = 1;
    }

    long corrLen = 2*len - 1;

    // Utworzenie i zerowanie tablicy correlation
    double * correlation = (double *) calloc(corrLen,sizeof(*correlation));

    // Obliczenie korelacji
    long count_from = 0;
    long count_to = 0;

    double coeff;

    for(i = 0; i < corrLen; i++){

        for(j = count_from; j <= count_to; j++){
            correlation[i] = correlation[i] + yy[len - 1 - j]*xx[i - j];  // bylo xlen zamiast len
        }

        if(i + 1 < len) count_to += 1;
        if(i + 1 > len - 1) count_from += 1;

        if (i == len) coeff = correlation[i];
    }


    // Normalizacja
    if(equalLengths > 0){
        switch(normType){

        case BIASED:
            for(i = 0; i < corrLen; i++) correlation[i] /= len;
            break;

        case UNBIASED:
            width = len - 1;
            j = 0;
            for(i = -width; i <= width; i++){
                if(i < 0) correlation[j] /= fabs(len + i);
                else      correlation[j] /= fabs(len - i);
                j += 1;
            }
            break;

        case COEFF:
            for(i = 0; i < corrLen; i++) correlation[i] /= coeff;
            break;

        case WITHOUT:
            break;
        }
    }

    Signal sig;
    sig.data = correlation;
    sig.length = corrLen;

    free(xx);
    free(yy);

    return sig;
}

// -------------------------------------------------------------------------------------------------
Signal random_noise(long length){

    unsigned long seed = (long) time(NULL);
    unsigned long m = (unsigned long) pow(2,32);
    unsigned long a = 1103515245;
    unsigned long c = 12345;
    unsigned long i;

    unsigned long * random_numbers = (unsigned long *) calloc(length, sizeof(unsigned long));
    double * noise = (double *) calloc(length, sizeof(double));

    random_numbers[0] = (a*seed + c)%m;

    for(i = 1; i < length; i++){
        random_numbers[i] = (a*random_numbers[i - 1] + c)%m;
    }

    for(i = 0; i < length; i++){
        noise[i] = (double) random_numbers[i] / m;
    }

    Signal result;
    result.length = length;
    result.data = noise;

    free(random_numbers);
    free(noise);

    return result;
}

// -------------------------------------------------------------------------------------------------
Signal normal_noise(long length, double mean, double variance){

    long length2;
    long i;
    double X1, X2, V1, V2, R, Y1, Y2;

    if(length%2){
        length2 = length + 1;
    }
    else{
        length2 = length;
    }

    Signal noise = random_noise(length2);
    double * normal = (double *) calloc(length, sizeof(double));

    for(i = 1; i < length2; i += 2){

        X1 = noise.data[i - 1];
        X2 = noise.data[i];

        V1 = 2*X1 - 1;
        V2 = 2*X2 - 1;

        R = sqrt(V1*V1 + V2*V2);

        Y1 = sqrt(-2*log(X1))*(V1/R);
        Y2 = sqrt(-2*log(X1))*(V2/R);

        normal[i - 1] = mean + variance*Y1;
        if(i < length2 - 1) normal[i] = mean + variance*Y2;
    }

    Signal result;
    result.length = length;
    result.data = normal;

    free(noise.data);

    return result;
}
// -------------------------------------------------------------------------------------------------
complexSignal dft(const double x[], long xlen){

    double complex * dft_tab = (double complex *) calloc(xlen, sizeof(*dft_tab));
    long k;
    long n;
    long N = xlen;
    double ft;
    double Re = 0;
    double Im = 0;

    // Nie można wykorzystać exp(wykladnik) bo exp() przyjmuje
    // wartosci double a nie double coplex
    for(k = 0; k < N; k++){
        for(n = 0; n < N; n++){
            //Re += x[n]*cos(-2*M_PI*k*n/N);
            //Im += x[n]*sin(-2*M_PI*k*n/N);
            dft_tab[k] +=   x[n] * (cos(-2*M_PI*k*n/N) + I*sin(-2*M_PI*k*n/N));
        }
        //dft_tab[k] = Re + I*Im;
        //Re = 0;
        //Im = 0;
    }

    complexSignal digital_fourier_transform;
    digital_fourier_transform.data = dft_tab;
    digital_fourier_transform.length = N;

    return digital_fourier_transform;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
