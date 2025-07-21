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
     , m_gain{1.f}
     , m_lowPoint{-0.9f}
     , m_highPoint{0.9f}
     , m_RC{0.f}
     , m_feedback{0.f}
    {}

    void init();
	  virtual void update(void);

    int getData(std::vector<float>& data);
    void setGain(float gain) { gain = std::max(std::min(gain, 60.f),1.f); m_gain = gain; };
    void setLowPoint(float lowPoint) { lowPoint = std::max(std::min(lowPoint, -0.01f), -0.99f); m_lowPoint = lowPoint; };
    void setHighPoint(float highPoint) { highPoint = std::max(std::min(highPoint, 0.99f), 0.01f); m_highPoint = highPoint; };
    void setFcrit(float f0) { m_RC = 1.f / (2.f * M_PI * f0); };
    void setFeedback(float fb) { m_feedback = fb; };
    void reset();

  private:
    // 1 channel input
    audio_block_t *inputQueueArray[1];

    static constexpr int k_firOrder = 64;
    static constexpr int k_nOversamp = 4;
    static constexpr int k_coeffsPerStage = k_firOrder / k_nOversamp;
    const double k_fs;
    float m_gain;
    float m_lowPoint;
    float m_highPoint;
    float m_RC;
    float m_feedback;
    std::array<float, 2> m_tubeMem;
    
    std::array<float, k_firOrder> m_firCoeffs;
    std::array<float, k_coeffsPerStage> m_upMem;
    std::array<float, k_firOrder> m_downMem;
    std::array<float, k_nOversamp> m_upBuffer;
    
    float transFunc(float x, float a, float b, float g, float d, float o);
    float vacTube(float input, double fs);
    void upsampling(float input);
    float downsampling();
};
