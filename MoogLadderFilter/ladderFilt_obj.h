#include <Audio.h>
#include <vector>

class LadderFilt_Obj : public AudioStream
{
  public:
    // 1 channel input
	  LadderFilt_Obj(void) : AudioStream(1, inputQueueArray) { }

    void init_ladder_processing(float fs);
	  virtual void update(void);

    int getData(std::vector<float>& data);
    void setResonance(double resonance) { m_resonance = resonance; changeParameter(); };
    void setFrequency(double cutoff) { m_cutoff = cutoff; changeParameter(); };
    void setSamplingrate(double samplerate) { m_fs = samplerate; changeParameter(); reset();};
    void reset();

  private:
    // 1 channel input
    audio_block_t *inputQueueArray[1];

    double m_fs;
    double m_cutoff;
    double m_resonance;

    double m_alpha_0;
    double m_alpha;
    double m_beta1, m_beta2, m_beta3, m_beta4;
    double m_state1, m_state2, m_state3, m_state4;

    void changeParameter();
};
