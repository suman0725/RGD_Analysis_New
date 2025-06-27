#ifndef PID_FUNCTIONS_H
#define PID_FUNCTIONS_H

#include "ParticleData.h"
#include "reader.h"

ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO, hipo::bank& SCIN,
                            IndexMap& cherMap, IndexMap& caloMap, IndexMap& scinMap);

IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName);

bool isSimpleElectron(const ParticleData& pd);

int bestPidFromTiming(const ParticleData& p, float eventStartTime, int eventNum, size_t particleIdx);

void assignPids(std::vector<ParticleData>& particles, float eventStartTime, int eventNum);

float PIDQuality(const ParticleData& p, int pid);

float computeDeltaT(const ParticleData& p, int pid);

void loadCCDBParams();

float getSamplingFractionMean(int sector, float measuredEnergy);

float getSamplingFractionSigma(int sector, float measuredEnergy);

float getSamplingFractionNSigma(float samplingFraction, float mean, float sigma);

float getTheoryBeta(float p, float mass);

bool hasHit(const ParticleData& p, int detector, int layer = -1);

#endif