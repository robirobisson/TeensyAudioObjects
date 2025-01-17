#include <stdint.h>
#include<vector>
#include <SerialFlash.h>
#include "tubeAmp_obj.h"

std::vector<float> dataBlock;


// initialisation of tubeAmp object
void TubeAmp_Obj::init_tubeAmp_processing(float fs) {
  m_fs = fs;
  m_gain = 1.0;
  m_lowPoint = -1.0;
  m_highPoint = 1.0;
  m_tubeMem[0] = 0.0;
  m_tubeMem[1] = 0.0;
  m_RC = 0.0;
  m_feedback = 0.0;
  m_nOversamp = 1;
  
  m_firCoeff = new float[m_firOrder];
  float firCoeff[m_firOrder] = {
    #include "fir_coeff.h"
  };
  for (int i = 0; i < m_firOrder; ++i) {
    m_firCoeff[i] = firCoeff[i];
  }
  
  m_coeffsMat = new float*[m_nOversamp];
  for (int i = 0; i < m_nOversamp; ++i) {
    m_coeffsMat[i] = new float[m_firOrder/m_nOversamp];
  }

  // Matrix mit den Polyphasen-Koeffizienten
  for (int j=0; j<m_nOversamp; j++) {
    for (int i=0; i<m_firOrder/m_nOversamp; i++)
      m_coeffsMat[j][i]= m_firCoeff[m_nOversamp*i+j];
  }
  
  m_upMem = new float[m_firOrder/m_nOversamp];
  m_downMem = new float[m_firOrder];
}

// update function (necessary for teensy audio-object)
void TubeAmp_Obj::update(void) {
  dataBlock.resize(AUDIO_BLOCK_SAMPLES);
  
  audio_block_t *cur_block;
	cur_block = receiveWritable(0);
 
  if (!cur_block) return;

  for (int i=0; i<AUDIO_BLOCK_SAMPLES; i++)
  {
    dataBlock[i] = cur_block->data[i] / 32768.;          // conversion to float32 [-1, 1)
  }
    
  getData(dataBlock);  // processing of audio block

  for (int i=0; i<AUDIO_BLOCK_SAMPLES; i++)
  {
    dataBlock[i] = dataBlock[i] * 32768;          // conversion to int16

    if (dataBlock[i] > 32767.)                       // hard clipper
        dataBlock[i] = 32767.;

    if (dataBlock[i] < -32768.)                      // hard clipper
        dataBlock[i] = -32768.;

    cur_block->data[i] = (int16_t) dataBlock[i];         // truncate to int16
  }

	transmit(cur_block, 0);
  release(cur_block);
}

// get data function to process one full audio block
int TubeAmp_Obj::getData(std::vector<float>& data)
{
  float up_buffer[m_nOversamp];
  for (auto kk = 0U; kk < data.size(); ++kk)
  { 
    // processing
    float curIn = data[kk];
    upsampling(curIn, up_buffer);
    for (int ii=0; ii<m_nOversamp; ii++) {
      up_buffer[ii] = vacTube(up_buffer[ii], m_fs*m_nOversamp);
    }
    data[kk] = downsampling(up_buffer);
  }
  return 0;
}


// transfer function of vacuum tube
float TubeAmp_Obj::transFunc(float x, float a, float b, float g, float d, float o)
{
  if ( x < (a/g)) {
    float k1 = a*a;
    float k2 = 1 + 2*a;

    float out = ((k1 + g*x)/(k2 - g*x)) -o;
    return out;
  }
  else if ((b/g) < x){
    float k1 = b*b;
    float k2 = 1 - 2*b;

    float out = ((g*x - k1)/(g*x + k2)) -o;
    return out;
  }
  else {
  float out = g*x +d -o;
  return out;
  }
}

// tube amp simulation
float TubeAmp_Obj::vacTube(float input, double fs)
{
  float o = 0;  // shifting parameter
  float d = -0.4;  // shifting parameter
  
  float alpha = 1/(2*m_RC*fs);
  float C1 = ((1-alpha)/(1+alpha)) * m_tubeMem[1] + (alpha/(1+alpha)) * m_tubeMem[0];
  float C2 = alpha/(1+alpha);

  float v_a = ((input/m_feedback) - (m_lowPoint/(m_gain*m_feedback)));
  float f_a = C1 + C2 * (m_lowPoint +d -o);
  float v_b = ((input/m_feedback) - (m_highPoint/(m_gain*m_feedback)));
  float f_b = C1 + C2 * (m_highPoint +d -o);


// computing v_n
  float v_n = 0;
    if ( f_a > v_a )  // A-Section
  {
    float k1 = m_lowPoint * m_lowPoint;
    float k2 = 1 + 2 * m_lowPoint;

    float A = m_gain * m_feedback;
    float B = k2 - m_gain*input - C1*m_gain*m_feedback + C2*m_gain*m_feedback*(1+o);
    float C = C1*m_gain*input - C1*k2 - C2*k1 + C2*k2*o - C2*m_gain*input*(1+o);
        
    v_n = ((-1)*B + sqrt(B*B - 4*A*C))/(2 * A);
  }
    else if ( f_b < v_b )  // B-Section
  {
    float k1 = m_highPoint*m_highPoint;
    float k2 = 1 - 2*m_highPoint;

    float A = (-1)*m_gain*m_feedback;
    float B = k2 + m_gain*input + C1*m_gain*m_feedback + C2*m_gain*m_feedback*(1-o);
    float C = (-1)*C1*m_gain*input - C1*k2 + C2*k1 + C2*k2*o - C2*m_gain*input*(1-o);
        
    v_n = ((-1)*B - sqrt(B*B - 4*A*C))/(2 * A);
  }
    else  // linear section
  {
    v_n = (C1 + C2*m_gain*input + C2*d - C2*o) / (1 + C2*m_gain*m_feedback);
  }

  m_tubeMem[1] = v_n; // vn-1 for next step/sample
  float x = input - m_feedback*v_n;
  float out = transFunc(x, m_lowPoint, m_highPoint, m_gain, d, o);
  m_tubeMem[0] = out; // yn-1 for next step/sample

  return out;
}

void TubeAmp_Obj::upsampling(float input, float buffer[])
{
  float acc;
  int numCoeffs = m_firOrder/m_nOversamp;

  // copy input to shift register
  m_upMem[0] = input;

  // convolution
  for (int j=0; j<m_nOversamp; j++)
  {
      acc = 0;
      for (int i=0; i<numCoeffs; i++)
          acc += m_coeffsMat[j][i] * m_upMem[i];

      buffer[j] = acc * m_nOversamp;
  }

    // shift register
    for (int i=numCoeffs-1; i>0; i--)
        m_upMem[i] = m_upMem[i-1];
}

float TubeAmp_Obj::downsampling(float buffer[]) {
  float acc;

  // copy samples to shift register
  for (int i=0; i<m_nOversamp; i++)
      m_downMem[i] = buffer[i];

  // convolution
  acc = 0;
  for (int i=0; i<m_firOrder; i++)
      acc += m_firCoeff[i] * m_downMem[i];

  // shift register
  for (int i=m_firOrder-1; i>=m_nOversamp; i--)
      m_downMem[i] = m_downMem[i-m_nOversamp];
  return acc;
}

void TubeAmp_Obj::setGain(float gain)
{
  // input value boundaries
  if (gain < 1){
    gain = 1;
  }
  if (gain > 60){
    gain = 60;
  }
  m_gain = gain; 
}

void TubeAmp_Obj::reset()
{
  m_tubeMem[0] = 0.0;
  m_tubeMem[1] = 0.0;
  
  m_downMem = new float[m_firOrder];
  for (int i=0; i<m_firOrder; i++) {
    //m_upBuffer[i] = 0;
    m_downMem[i] = 0;
  }
  m_upMem = new float[m_firOrder/m_nOversamp];
  for (int i=0; i<m_firOrder/m_nOversamp; i++)
    m_upMem[i] = 0;
  m_coeffsMat = new float*[m_nOversamp];
  for (int i = 0; i < m_nOversamp; ++i) {
    m_coeffsMat[i] = new float[m_firOrder/m_nOversamp];
  }
  // write coefficients into matrix
  for (int j=0; j<m_nOversamp; j++) {
    for (int i=0; i<m_firOrder/m_nOversamp; i++)
      m_coeffsMat[j][i]= m_firCoeff[m_nOversamp*i+j];
  }
}
