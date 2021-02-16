#include <Audio.h>
#include <vector>

class TubeAmp_Obj : public AudioStream
{
  public:
    // 1 channel input
	  TubeAmp_Obj(void) : AudioStream(1, inputQueueArray) { }

    void init_tubeAmp_processing(float fs);
	  virtual void update(void);

    int getData(std::vector<float>& data);
    void setGain(float gain);
    void setLowPoint(float lowPoint) {m_lowPoint = lowPoint; };
    void setHighPoint(float highPoint) {m_highPoint = highPoint; };
    void setFcrit(float f0) {m_RC = 1/(2*M_PI*f0); };
    void setFeedback(float fb) {m_feedback = fb; };
    void setSamplingrate(double samplerate) {m_fs = samplerate; reset();};
    void setOversampling(int oversampFac) {m_nOversamp = oversampFac; reset();};
    void reset();

  private:
    // 1 channel input
    audio_block_t *inputQueueArray[1];

    double m_fs;
    float m_gain;
    float m_lowPoint;
    float m_highPoint;
    float m_tubeMem[2];
    float m_RC;
    float m_feedback;
    int m_nOversamp;
    int m_firOrder = 64;
    float* m_firCoeff;
    float** m_coeffsMat;
    float* m_upMem;
    float* m_downMem;
    
    float transFunc(float x, float a, float b, float g, float d, float o);
    float vacTube(float input, double fs);
    void upsampling(float input, float buffer[]);
    float downsampling(float buffer[]);
};
