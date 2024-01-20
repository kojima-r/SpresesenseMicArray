#include <Audio.h>


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

FFTClass<MAX_CHANNEL_NUM, FFT_LEN> FFT;
/* Parameters */
// 円形マイク
const float radius = 0.1; // m
const float offset = 45;  // degree
// 線形マイク:マイク間距離
const float dx = 0.05;     // m

const float c = 340;      // m/sec

//ステアリングベクトル
//      角度[0~360)、周波数ビン、チャンネル、実数、虚数
float stearing_vec[MAX_ANGLE/RESOLUTION][FFT_LEN/2][MAX_CHANNEL_NUM][2];


AudioClass *theAudio;
arm_rfft_fast_instance_f32 iS;

/* Select mic channel number */
//const int mic_channel_num = 1;
//const int mic_channel_num = 2;
const int mic_channel_num = 4;

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
        stearing_vec[theta][k][ch][0] = arm_cos_f32(2.0*M_PI*k/(float)FFT_LEN*SAMPLING_RATE*tau)/(float)MAX_CHANNEL_NUM;
        stearing_vec[theta][k][ch][1] = arm_sin_f32(2.0*M_PI*k/(float)FFT_LEN*SAMPLING_RATE*tau)/(float)MAX_CHANNEL_NUM;
      }
    }
  }
}


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
  static const int32_t buffer_sample = 768 * mic_channel_num;
  static const int32_t buffer_size = buffer_sample * sizeof(int16_t);
  static char  buffer[buffer_size];
  static float result[4][MAX_ANGLE/RESOLUTION];
  uint32_t read_size;
  static uint32_t succ_count=0;
  
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
      FFT.put((q15_t*)buffer,buffer_sample/mic_channel_num);
      while (!FFT.empty(0)) {
        for (int i = 0; i < MAX_CHANNEL_NUM; i++) {
          FFT.get_raw(&pTmp[i][0], i);
        }
        for(int theta=0;theta*RESOLUTION<MAX_ANGLE;theta++){
          for(int k=0;k<FFT_LEN;k++) pDst[k] = 0;
          for(int k=0;k<FFT_LEN/2;k++){
            for(int ch=0;ch<MAX_CHANNEL_NUM;ch++){
              pDst[2*k] += stearing_vec[theta][k][ch][0] * pTmp[ch][2*k];
              pDst[2*k] += -stearing_vec[theta][k][ch][1] * pTmp[ch][2*k+1];
              pDst[2*k+1] += stearing_vec[theta][k][ch][0] * pTmp[ch][2*k+1];
              pDst[2*k+1] += stearing_vec[theta][k][ch][1] * pTmp[ch][2*k]; 
            }
          }
          result[pos][theta] = 0;
          arm_cmplx_mag_f32(pDst, pPower, FFT_LEN/2);
          for(int k=0;k<FFT_LEN/2;k++){
            result[pos][theta] += pPower[k];
          }
          printf("%d, %f\n", theta*RESOLUTION, result[pos][theta]);
        }
        pos = (pos+1)%result_size;
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
