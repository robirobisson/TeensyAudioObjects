#include <Audio.h>
#include <vector>

class LadderFilt_Obj : public AudioStream
{
  public:
    // 1 channel input
	  LadderFilt_Obj(void) : AudioStream(1, inputQueueArray) { }

	  virtual void update(void);

    int getData(std::vector<float>& data);
    void setResonance(double resonance) { if(resonance <= 4.0) {m_resonance = resonance;} changeParameter(); };  // 4.0 is self-osc border
    void setFrequency(double cutoff) { m_cutoff = cutoff; changeParameter(); };
    void setSamplingrate(double samplerate) { m_fs = samplerate; changeParameter(); reset();};
    void reset();

  private:
    // 1 channel input
    audio_block_t *inputQueueArray[1];

    double m_fs{AUDIO_SAMPLE_RATE_EXACT};
    double m_cutoff{1.0};
    double m_resonance{0.707};

    double m_alpha_0{0.0};
    double m_alpha{0.0};
    double m_beta1{0.0};
    double m_beta2{0.0};
    double m_beta3{0.0};
    double m_beta4{0.0};
    double m_state1{0.0};
    double m_state2{0.0};
    double m_state3{0.0};
    double m_state4{0.0};

    void changeParameter();
};
