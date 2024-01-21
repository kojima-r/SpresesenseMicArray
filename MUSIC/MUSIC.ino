#include <Audio.h>

#include <cstdio>
#include <cmath>
#include <cstring>

//#include <stlport.h>
//#include <Eigen30.h>

#include "FFT.h"

/*-----------------------------------------------------------------*/
/*
 * FFT parameters
 */
/* Select FFT length */

// #define FFT_LEN 32
// #define FFT_LEN 64
// #define FFT_LEN 128
// #define FFT_LEN 256
//#define FFT_LEN 512

#define FFT_LEN 1024

// #define FFT_LEN 2048
// #define FFT_LEN 4096

/* Number of channels*/
// #define MAX_CHANNEL_NUM 1
//#define MAX_CHANNEL_NUM 2
#define MAX_CHANNEL_NUM 4
#define SAMPLING_RATE 48000
#define RESOLUTION 20
#define MAX_ANGLE 360
#define MUSIC_FREQ_BIN 1

FFTClass<MAX_CHANNEL_NUM, FFT_LEN> FFT;
/* Parameters */
// 円形マイク
const float radius = 0.1; // m
const float offset = 45;  // degree
// 線形マイク:マイク間距離
const float dx = 0.05;     // m

const float c = 340;      // m/sec

const int delta_freq_MUSIC = 10;


/* Select mic channel number */
//const int mic_channel_num = 1;
//const int mic_channel_num = 2;
const int mic_channel_num = 4;



struct Complex {
    float real;
    float imag;

    Complex() : real(0), imag(0) {}
    Complex(float r, float i) : real(r), imag(i) {}

    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator*(const float& other) const {
    	return Complex(real*other, imag*other);
    }
    Complex operator*(const Complex& other) const {
        return Complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    Complex operator/(const float& scalar) const {
        return Complex(real / scalar, imag / scalar);
    }
};

//typename struct Complex_ Complex;

inline struct Complex conj(const struct Complex& a) {
    return Complex(a.real,-a.imag);
};
inline float norm2(const struct Complex& a) {
    return std::sqrt(a.real*a.real+a.imag*a.imag);
};

// Function to print a 2D array
void printArray(Complex* arr, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%2.2f + %2.2fi ", arr[i * cols + j].real, arr[i * cols + j].imag);
        }
        printf("\n");
    }
};



//ステアリングベクトル
//      角度[0~360)、周波数ビン、チャンネル、実数、虚数
Complex stearing_vec[MAX_ANGLE/RESOLUTION][FFT_LEN/2][MAX_CHANNEL_NUM];

/*
// 円形マイク
void calc_stearing_vec(){
  for(int theta=0; theta*RESOLUTION<MAX_ANGLE; theta++){
    // theta: target direction
    float radian_theta = theta*RESOLUTION*2*M_PI/360.0;
    for(int ch=0; ch<MAX_CHANNEL_NUM; ch++){
      // ch: mic direction
      float radian_ch = (offset+360/MAX_CHANNEL_NUM*ch)*2*M_PI/360.0;
      //円の中心を時刻０としたときの平面波の到達時間
      float tau = -radius*arm_cos_f32(radian_theta-radian_ch)/c;
      for(int k=0;k<FFT_LEN/2;k++){
        stearing_vec[theta][k][ch][0] = arm_cos_f32(2.0*M_PI*k/(float)FFT_LEN*SAMPLING_RATE*tau)/(float)MAX_CHANNEL_NUM;
        stearing_vec[theta][k][ch][1] = arm_sin_f32(2.0*M_PI*k/(float)FFT_LEN*SAMPLING_RATE*tau)/(float)MAX_CHANNEL_NUM;
      }
    }
  }
}*/

// 線形マイク
// https://setoti.hatenablog.com/entry/beamformer
void calc_stearing_vec(){
  for(int theta=0; theta*RESOLUTION<MAX_ANGLE; theta++){
    // theta: target direction
    float radian_theta = theta*RESOLUTION*2*M_PI/360.0;
    for(int ch=0; ch<MAX_CHANNEL_NUM; ch++){
      // ch: mic index
      // 円の中心を時刻０としたときの平面波の到達時間
      float tau=(ch-(MAX_CHANNEL_NUM-1)/2.0f)*radius*arm_sin_f32(radian_theta)/c;
      for(int k=0;k<FFT_LEN/2;k++){
        stearing_vec[theta][k][ch].real = arm_cos_f32(2.0*M_PI*k/(float)FFT_LEN*SAMPLING_RATE*tau)/(float)MAX_CHANNEL_NUM;
        stearing_vec[theta][k][ch].imag = arm_sin_f32(2.0*M_PI*k/(float)FFT_LEN*SAMPLING_RATE*tau)/(float)MAX_CHANNEL_NUM;
      }
    }
  }
}



struct Complex corrMat[MAX_ANGLE/RESOLUTION][MUSIC_FREQ_BIN*MAX_CHANNEL_NUM][MAX_CHANNEL_NUM];
struct Complex X[MAX_CHANNEL_NUM * MAX_CHANNEL_NUM];
struct Complex An[MAX_CHANNEL_NUM * MAX_CHANNEL_NUM];

AudioClass *theAudio;
arm_rfft_fast_instance_f32 iS;

void setup()
{
  Serial.begin(115200);
  while (!Serial);

  Serial.println("Init Audio Library");
  theAudio = AudioClass::getInstance();
  theAudio->begin();

  Serial.println("Init Audio Recorder");
  /* Select input device as AMIC */
  theAudio->setRecorderMode(AS_SETRECDR_STS_INPUTDEVICE_MIC, 210);

  /* Set PCM capture */
  uint8_t channel;
  switch (mic_channel_num) {
  case 1: channel = AS_CHANNEL_MONO;   break;
  case 2: channel = AS_CHANNEL_STEREO; break;
  case 4: channel = AS_CHANNEL_4CH;    break;
  }
  theAudio->initRecorder(AS_CODECTYPE_PCM, "/mnt/sd0/BIN", AS_SAMPLINGRATE_48000, channel);
    /* begin FFT */
  FFT.begin();

  arm_rfft_fast_init_f32(&iS, FFT_LEN);
  calc_stearing_vec();
  Serial.println("Rec start!");
  theAudio->startRecorder();
}

int result_size = 4;

void loop()
{
  static int pos = 0;
  static const int32_t buffer_sample = 768 * mic_channel_num;//768 * mic_channel_num;
  static const int32_t buffer_size = buffer_sample * sizeof(int16_t);
  static char  buffer[buffer_size];
  static float result[4][MAX_ANGLE/RESOLUTION];

  uint32_t read_size;
  static uint32_t succ_count=0;
  static uint32_t corr_count=0;
  
  /* Read frames to record in buffer */
  int err = theAudio->readFrames(buffer, buffer_size, &read_size);
  if (err != AUDIOLIB_ECODE_OK && err != AUDIOLIB_ECODE_INSUFFICIENT_BUFFER_AREA) {
    printf("Error err = %d\n", err);
    sleep(1);
    theAudio->stopRecorder();
    exit(1);
  }

  if ((read_size != 0) && (read_size >= buffer_size)) {
    succ_count+=1;
    if(succ_count%10==0){
      static float pTmp[MAX_CHANNEL_NUM][FFT_LEN];
      static float pDst[FFT_LEN];
      static float pPower[FFT_LEN/2];
      corr_count+=1;
      FFT.put((q15_t*)buffer,buffer_sample/mic_channel_num);
      while (!FFT.empty(0)) {
        for (int i = 0; i < MAX_CHANNEL_NUM; i++) {
          FFT.get_raw(&pTmp[i][0], i);
        }
        for(int theta=0;theta*RESOLUTION<MAX_ANGLE;theta++){
          /*
          for(int k=0;k<FFT_LEN;k++) pDst[k] = 0;
          for(int k=0;k<FFT_LEN/2;k++){
            for(int ch=0;ch<MAX_CHANNEL_NUM;ch++){
              pDst[2*k] += stearing_vec[theta][k][ch][0] * pTmp[ch][2*k];
              pDst[2*k] += -stearing_vec[theta][k][ch][1] * pTmp[ch][2*k+1];
              pDst[2*k+1] += stearing_vec[theta][k][ch][0] * pTmp[ch][2*k+1];
              pDst[2*k+1] += stearing_vec[theta][k][ch][1] * pTmp[ch][2*k]; 
            }
          }*/
          for(int k=0;k<MUSIC_FREQ_BIN;k++){//2*(MUSIC_FREQ_BIN*delta_freq_MUSIC)<FFT_LEN
            for(int c1=0;c1<MAX_CHANNEL_NUM;c1++){
              for(int c2=0;c2<MAX_CHANNEL_NUM;c2++){ // pTemp*conj(pTemp)
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].real +=  pTmp[c1][2*(k*delta_freq_MUSIC)]  *pTmp[c2][2*(k*delta_freq_MUSIC)];
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].real +=  pTmp[c1][2*(k*delta_freq_MUSIC)+1]*pTmp[c2][2*(k*delta_freq_MUSIC)+1];
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].imag += -pTmp[c1][2*(k*delta_freq_MUSIC)]  *pTmp[c2][2*(k*delta_freq_MUSIC)+1];
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].imag +=  pTmp[c1][2*(k*delta_freq_MUSIC)+1]*pTmp[c2][2*(k*delta_freq_MUSIC)];
              }
            }
          }
          /*
          result[pos][theta] = 0;
          arm_cmplx_mag_f32(pDst, pPower, FFT_LEN/2);
          for(int k=0;k<FFT_LEN/2;k++){
            result[pos][theta] += pPower[k];
          }
          printf("%d, %f\n", theta*RESOLUTION, result[pos][theta]);
          */
        }
      }
      //
      if(corr_count > 10){
        for(int theta=0;theta*RESOLUTION<MAX_ANGLE;theta++){
          result[pos][theta] = 0;
          for(int k=0;k<MUSIC_FREQ_BIN;k++){
            for(int c1=0;c1<MAX_CHANNEL_NUM;c1++){
              for(int c2=0;c2<MAX_CHANNEL_NUM;c2++){
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].real/=10.0f;
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].imag/=10.0f;
              }
            }
            //
            //corrMat[theta][k]
            get_eigen_val(corrMat[theta][k],MAX_CHANNEL_NUM);
            //
            float denom=0;
            for(int noise_ch=2;noise_ch<MAX_CHANNEL_NUM;noise_ch++){
              Complex denom_dot;
              for(int c1=0;c1<MAX_CHANNEL_NUM;c1++){
                denom_dot = denom_dot + conj(stearing_vec[theta][k*delta_freq_MUSIC][c1])*X[c1,noise_ch];
              }
              denom+=norm2(denom_dot);
            }
            ///
            float s=0;
            for(int c1=0;c1<MAX_CHANNEL_NUM;c1++){
              s+=norm2(stearing_vec[theta][k*delta_freq_MUSIC][c1]);
            }
            result[pos][theta] += s/denom;
          }
          printf("%d, %f\n", theta*RESOLUTION, result[pos][theta]);
          pos = (pos+1)%result_size;
        }
        for(int theta=0;theta*RESOLUTION<MAX_ANGLE;theta++){
          for(int k=0;k<MUSIC_FREQ_BIN;k++){
            for(int c1=0;c1<MAX_CHANNEL_NUM;c1++){
              for(int c2=0;c2<MAX_CHANNEL_NUM;c2++){
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].real=0;
                corrMat[theta][k][c1*MAX_CHANNEL_NUM+c2].imag=0;
              }
            }
          }
        }
        corr_count=0;
      }
    }
  } else {
  }
}

void errorLoop(int num)
{
  int i;

  while (1) {
    for (i = 0; i < num; i++) {
      ledOn(LED0);
      delay(300);
      ledOff(LED0);
      delay(300);
    }
    delay(1000);
  }
}


int get_eigen_val(Complex* A, int N) {
    
    // Initialize A and X matrices
    for (int i = 0; i < N; ++i) {
      X[i*N+i] = Complex(1.0, 0.0); // Identity matrix
    }
    
    for (int iter = 0; iter < 10; ++iter) {
        float max_A = 0.0;
        int p = 0, q = 0;
        // Find p, q
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (i != j) {
                    float a = A[i * N + j].real * A[i * N + j].real + A[i * N + j].imag * A[i * N + j].imag;
                    if (a > max_A) {
                        max_A = a;
                        p = i;
                        q = j;
                    }
                }
            }
        }

        //printf("max:%f (%d,%d)\n", max_A , p ,q);
      	if(max_A<1.0e-10){
          break;
	      }
        float theta = M_PI / 4.0 - 0.5 * std::atan2((A[p * N + p].real - A[q * N + q].real), 2.0 * std::sqrt(max_A));
        float phi = -std::atan2(A[p * N + q].imag, A[p * N + q].real);
        //printf("%f %f\n", theta, phi);

        float cos_theta = std::cos(theta);
        float sin_theta = std::sin(theta);
        Complex e_neg_phi(std::cos(phi),-std::sin(phi));
        Complex e_phi(std::cos(phi),std::sin(phi));

        // Update An matrix
        for (int i = 0; i < N; ++i) {
            for (int j = i; j < N; ++j) {
                if (i == p && j != q && j != p) {
                    An[p * N + j] = A[p * N + j] * cos_theta + A[q * N + j] * sin_theta * e_neg_phi;
                    An[j * N + p] = conj(An[p * N + j]);
                } else if (i != q && i != p && j == p) {
                    An[p * N + i] = A[p * N + i] * cos_theta + A[q * N + i] * sin_theta * e_neg_phi;
                    An[i * N + p] = conj(An[p * N + i]);

                } else if (i == q && j != q) {
                    An[q * N + j] = A[p * N + j] * sin_theta*e_phi - A[q * N + j] * cos_theta;
                    An[j * N + q] = conj(An[q * N + j]);
                } else if (i != q && i != p && j == q) {
                    An[q * N + i] = A[p * N + i] * sin_theta*e_phi - A[q * N + i] * cos_theta;
                    An[i * N + q] = conj(An[q * N + i]);

                } else if (i == p && j == p) {
                    An[p * N + p] = A[p * N + p] * cos_theta * cos_theta + A[q * N + q] * sin_theta * sin_theta +
                                    A[p * N + q] * sin_theta * cos_theta * e_phi + A[q * N + p] * sin_theta * cos_theta *e_neg_phi;
                } else if (i == q && j == q) {
                    An[q * N + q] = A[p * N + p] * sin_theta * sin_theta + A[q * N + q] * cos_theta * cos_theta -
                                    A[p * N + q] * sin_theta * cos_theta * e_phi - A[q * N + p] * sin_theta * cos_theta * e_neg_phi;
                } else if (i == p && j == q) {
                    An[p * N + q] = A[p * N + p] * sin_theta * cos_theta - A[q * N + q] * sin_theta * cos_theta -
                                    A[p * N + q] * cos_theta * cos_theta * e_phi + A[q * N + p] * sin_theta * sin_theta * e_neg_phi;
                    An[p * N + q] = e_neg_phi*An[p * N + q];
                    An[q * N + p] = conj(An[p * N + q]);
                } else {
                    An[i * N + j] = A[i * N + j];
                    An[j * N + i] = A[j * N + i];
                }
            }
        }

        // Update X matrix
        for (int i = 0; i < N; ++i) {
            Complex bip = X[i * N + p];
            Complex biq = X[i * N + q];
            X[i * N + p] = bip * cos_theta + biq * sin_theta*e_phi;
            X[i * N + q] = bip * sin_theta*e_neg_phi - biq * cos_theta;
        }

        // Display An matrix
        //printArray(An, N, N);
	      memcpy(A,An,sizeof(Complex)*N*N);
        // Display diagonal elements of An
        //for (int i = 0; i < N; ++i) {
        //    printf(">> %f %f\n", An[i * N + i].real,An[i * N + i].imag);
	      //}
    }
    /*
    for (int i = 0; i < N; ++i) {
      printf(">> %f %f\n", An[i * N + i].real,An[i * N + i].imag);
    }
    printArray(X,N,N);
    */
    return 0;
}


