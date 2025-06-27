#ifndef CONTAMINATION_UTILS_H
#define CONTAMINATION_UTILS_H

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

double calculateContamination(TH1F* histPions, TH1F* histKaons, double min, double max,
                              TCanvas* canvas, const char* binName, TLegend* legend);

#endif