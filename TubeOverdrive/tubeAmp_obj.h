#include <Audio.h>
#include <array>
#include <vector>

class TubeAmp_Obj : public AudioStream
{
  public:
    // 1 channel input
	  TubeAmp_Obj(void)
     : AudioStream(1, inputQueueArray)
     , k_fs{AUDIO_SAMPLE_RATE_EXACT}
    {}

    void init_tubeAmp_processing();
	  virtual void update(void);

    int getData(std::vector<float>& data);
    void setGain(float gain);
    void setLowPoint(float lowPoint) { m_lowPoint = lowPoint; };
    void setHighPoint(float highPoint) { m_highPoint = highPoint; };
    void setFcrit(float f0) { m_RC = 1 / (2 * M_PI * f0); };
    void setFeedback(float fb) { m_feedback = fb; };
    void reset();

  private:
    // 1 channel input
    audio_block_t *inputQueueArray[1];

    static constexpr int k_firOrder = 64;
    static constexpr int k_nOversamp = 4;
    const double k_fs;
    float m_gain;
    float m_lowPoint;
    float m_highPoint;
    float m_tubeMem[2];
    float m_RC;
    float m_feedback;
    
    std::array<float, k_firOrder> m_firCoeffs;
    std::array<std::array<float, (k_firOrder / k_nOversamp)>, k_nOversamp> m_coeffsMat;
    std::array<float, (k_firOrder / k_nOversamp)> m_upMem;
    std::array<float, k_firOrder> m_downMem;
    std::array<float, k_nOversamp> m_upBuffer;
    
    float transFunc(float x, float a, float b, float g, float d, float o);
    float vacTube(float input, double fs);
    void upsampling(float input);
    float downsampling();
    void updateMatrix();
};
