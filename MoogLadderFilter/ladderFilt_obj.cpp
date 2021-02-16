#include <stdint.h>
#include<vector>

#include "ladderFilt_obj.h"

std::vector<float> dataBlock;


// initialisation of filter object
void LadderFilt_Obj::init_ladder_processing(float fs) {
  m_fs = fs;
  m_cutoff = 1.0;
  m_resonance = 0.707;

  m_alpha_0 = 0.0;
  m_alpha = 0.0;
  m_beta1 = 0.0;
  m_beta2 = 0.0;
  m_beta3 = 0.0;
  m_beta4 = 0.0;
  m_state1 = 0.0;
  m_state2 = 0.0;
  m_state3 = 0.0;
  m_state4 = 0.0;
}

// update function (necessary for teensy audioobject)
void LadderFilt_Obj::update(void) {
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
int LadderFilt_Obj::getData(std::vector<float>& data)
{
  for (auto kk = 0U; kk < data.size(); ++kk)
  { 
    // processing
    double curIn = data[kk];
    double u = m_alpha_0 * (curIn - m_resonance * (m_state1 * m_beta1 +
      m_state2 * m_beta2 + m_state3 * m_beta3 + m_state4 * m_beta4));

    double v1 = m_alpha * (u - m_state1);
    double y1 = v1 + m_state1;

    double v2 = m_alpha * (y1 - m_state2);
    double y2 = v2 + m_state2;

    double v3 = m_alpha * (y2 - m_state3);
    double y3 = v3 + m_state3;

    double v4 = m_alpha * (y3 - m_state4);
    double y4 = v4 + m_state4;

    data[kk] = y4;

    m_state1 = v1 + y1;
    m_state2 = v2 + y2;
    m_state3 = v3 + y3;
    m_state4 = v4 + y4;
  }
  return 0;
}

// reset function to clear internal states
void LadderFilt_Obj::reset()
{
  m_state1 = 0.0;
  m_state2 = 0.0;
  m_state3 = 0.0;
  m_state4 = 0.0;
}

// compute internal parameters if freq/res is changed
void LadderFilt_Obj::changeParameter()
{
  double g = tan(M_PI * m_cutoff / m_fs);
  m_alpha = g / (1.0 + g);
  m_alpha_0 = 1.0 / (1.0 + m_resonance * m_alpha * m_alpha * m_alpha * m_alpha);
  m_beta4 = 1.0 / (1.0 + g);
  m_beta3 = m_alpha / (1.0 + g);
  m_beta2 = m_alpha * m_alpha / (1.0 + g);
  m_beta1 = m_alpha * m_alpha * m_alpha / (1.0 + g);
}
