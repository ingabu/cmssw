#ifndef HcalTimingParam_h
#define HcalTimingParam_h

#include <boost/cstdint.hpp>

class HcalTimingParam{
 public:
  HcalTimingParam():mId(0),m_nhits(0),m_phase(0.0),m_rms(0.0){}
  HcalTimingParam(unsigned long fId, unsigned int nhits, float phase,float rms):
      mId(fId), m_nhits(nhits),m_phase(phase),m_rms(rms){}
  uint32_t rawId () const {return mId;}
  float phase() const {return m_phase;}
  float rms() const {return m_rms;}
  unsigned int nhits() const {return m_nhits;}
 private:
  uint32_t mId;
  uint32_t m_nhits;
  float m_phase;
  float m_rms;
};
#endif
