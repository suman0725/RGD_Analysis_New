#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include "reader.h"
#include "ParticleData.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: ./skim <target> <bending> (e.g., CxC IB or CxC OB)" << endl;
        return 1;
    }
    string target = argv[1];
    string bending = argv[2];
    if (bending != "IB" && bending != "OB") {
        cerr << "Error: Bending must be 'IB' or 'OB'" << endl;
        return 1;
    }
    float torusValue = (bending == "IB") ? -1.0f : 1.0f;

    string baseDir = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11";
    string targetDir = baseDir + "/" + target;

    if (!fs::exists(targetDir)) {
        cerr << "Error: Input directory not found: " << targetDir << endl;
        return 1;
    }

    string outputBasePath = "/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/Contamination/Data-driven/skim_data/skim_pass0v11";
    try {
        fs::create_directories(outputBasePath);
        fs::create_directories(outputBasePath + "/" + target);
        if (target == "CuSn") {
            fs::create_directories(outputBasePath + "/CuSn/Cu/root");
            fs::create_directories(outputBasePath + "/CuSn/Sn/root");
        } else {
            fs::create_directories(outputBasePath + "/" + target + "/root");
        }
    } catch (const fs::filesystem_error& e) {
        cerr << "Error: Failed to create output directories: " << e.what() << endl;
        return 1;
    }

    loadCCDBParams();

    TFile* outputFile = nullptr;
    TFile* outputFileCu = nullptr;
    TFile* outputFileSn = nullptr;

    if (target != "CuSn") {
        string outputPath = outputBasePath + "/" + target + "/root/pkptree_" + target + "_" + bending + ".root";
        outputFile = new TFile(outputPath.c_str(), "RECREATE");
        if (!outputFile || outputFile->IsZombie()) {
            cerr << "Error: Failed to create output file: " << outputPath << endl;
            return 1;
        }
    } else {
        string outputPathCu = outputBasePath + "/CuSn/Cu/root/pkptree_Cu_" + bending + ".root";
        string outputPathSn = outputBasePath + "/CuSn/Sn/root/pkptree_Sn_" + bending + ".root";
        outputFileCu = new TFile(outputPathCu.c_str(), "RECREATE");
        outputFileSn = new TFile(outputPathSn.c_str(), "RECREATE");
        if (!outputFileCu || outputFileCu->IsZombie() || !outputFileSn || outputFileSn->IsZombie()) {
            cerr << "Error: Failed to create output files: " << outputPathCu << " or " << outputPathSn << endl;
            return 1;
        }
    }

    TTree* treeEBPidPions = new TTree("EB_pid_pions", "Trigger particles with pid == 211");
    TTree* treeEBPidKaons = new TTree("EB_pid_kaons", "Trigger particles with pid == 321");
    TTree* treeEBPidProtons = new TTree("EB_pid_protons", "Trigger particles with pid == 2212");
    TTree* treeEBPidNegPions = new TTree("EB_pid_neg_pions", "Trigger particles with pid == -211");
    TTree* treeEBPidNegKaons = new TTree("EB_pid_neg_kaons", "Trigger particles with pid == -321");
    TTree* treePionHypothesis = new TTree("pion_hypothesis", "Pion hypothesis");
    TTree* treeKaonHypothesis = new TTree("kaon_hypothesis", "Kaon hypothesis");
    TTree* treeProtonHypothesis = new TTree("proton_hypothesis", "Proton hypothesis");
    TTree* treeNegPionHypothesis = new TTree("neg_pion_hypothesis", "Negative pion hypothesis");
    TTree* treeNegKaonHypothesis = new TTree("neg_kaon_hypothesis", "Negative kaon hypothesis");

    TTree* treeEBPidPionsCu = nullptr;
    TTree* treeEBPidKaonsCu = nullptr;
    TTree* treeEBPidProtonsCu = nullptr;
    TTree* treeEBPidNegPionsCu = nullptr;
    TTree* treeEBPidNegKaonsCu = nullptr;
    TTree* treePionHypothesisCu = nullptr;
    TTree* treeKaonHypothesisCu = nullptr;
    TTree* treeProtonHypothesisCu = nullptr;
    TTree* treeNegPionHypothesisCu = nullptr;
    TTree* treeNegKaonHypothesisCu = nullptr;
    TTree* treeEBPidPionsSn = nullptr;
    TTree* treeEBPidKaonsSn = nullptr;
    TTree* treeEBPidProtonsSn = nullptr;
    TTree* treeEBPidNegPionsSn = nullptr;
    TTree* treeEBPidNegKaonsSn = nullptr;
    TTree* treePionHypothesisSn = nullptr;
    TTree* treeKaonHypothesisSn = nullptr;
    TTree* treeProtonHypothesisSn = nullptr;
    TTree* treeNegPionHypothesisSn = nullptr;
    TTree* treeNegKaonHypothesisSn = nullptr;

    if (target == "CuSn") {
        treeEBPidPionsCu = new TTree("EB_pid_pions", "Trigger particles with pid == 211");
        treeEBPidKaonsCu = new TTree("EB_pid_kaons", "Trigger particles with pid == 321");
        treeEBPidProtonsCu = new TTree("EB_pid_protons", "Trigger particles with pid == 2212");
        treeEBPidNegPionsCu = new TTree("EB_pid_neg_pions", "Trigger particles with pid == -211");
        treeEBPidNegKaonsCu = new TTree("EB_pid_neg_kaons", "Trigger particles with pid == -321");
        treePionHypothesisCu = new TTree("pion_hypothesis", "Pion hypothesis");
        treeKaonHypothesisCu = new TTree("kaon_hypothesis", "Kaon hypothesis");
        treeProtonHypothesisCu = new TTree("proton_hypothesis", "Proton hypothesis");
        treeNegPionHypothesisCu = new TTree("neg_pion_hypothesis", "Negative pion hypothesis");
        treeNegKaonHypothesisCu = new TTree("neg_kaon_hypothesis", "Negative kaon hypothesis");
        treeEBPidPionsSn = new TTree("EB_pid_pions", "Trigger particles with pid == 211");
        treeEBPidKaonsSn = new TTree("EB_pid_kaons", "Trigger particles with pid == 321");
        treeEBPidProtonsSn = new TTree("EB_pid_protons", "Trigger particles with pid == 2212");
        treeEBPidNegPionsSn = new TTree("EB_pid_neg_pions", "Trigger particles with pid == -211");
        treeEBPidNegKaonsSn = new TTree("EB_pid_neg_kaons", "Trigger particles with pid == -321");
        treePionHypothesisSn = new TTree("pion_hypothesis", "Pion hypothesis");
        treeKaonHypothesisSn = new TTree("kaon_hypothesis", "Kaon hypothesis");
        treeProtonHypothesisSn = new TTree("proton_hypothesis", "Proton hypothesis");
        treeNegPionHypothesisSn = new TTree("neg_pion_hypothesis", "Negative pion hypothesis");
        treeNegKaonHypothesisSn = new TTree("neg_kaon_hypothesis", "Negative kaon hypothesis");
    }

    float chi2pid, momentum, beta;
    int pid;
    float chi2pid_kaon = 99999.0f;
    float chi2pid_proton = 99999.0f;
    float chi2pid_neg_pion = 99999.0f;
    float chi2pid_neg_kaon = 99999.0f;
    float beta_pion = 99999.0f;
    float beta_kaon = 99999.0f;
    float beta_proton = 99999.0f;
    float beta_neg_pion = 99999.0f;
    float beta_neg_kaon = 99999.0f;

    treeEBPidPions->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidPions->Branch("p", &momentum, "p/F");
    treeEBPidPions->Branch("beta", &beta, "beta/F");
    treeEBPidKaons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidKaons->Branch("p", &momentum, "p/F");
    treeEBPidKaons->Branch("beta", &beta, "beta/F");
    treeEBPidProtons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidProtons->Branch("p", &momentum, "p/F");
    treeEBPidProtons->Branch("beta", &beta, "beta/F");
    treeEBPidNegPions->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidNegPions->Branch("p", &momentum, "p/F");
    treeEBPidNegPions->Branch("beta", &beta, "beta/F");
    treeEBPidNegKaons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidNegKaons->Branch("p", &momentum, "p/F");
    treeEBPidNegKaons->Branch("beta", &beta, "beta/F");
    treePionHypothesis->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treePionHypothesis->Branch("p", &momentum, "p/F");
    treePionHypothesis->Branch("beta", &beta_pion, "beta/F");
    treePionHypothesis->Branch("pid", &pid, "pid/I");
    treeKaonHypothesis->Branch("chi2pid_kaon", &chi2pid_kaon, "chi2pid_kaon/F");
    treeKaonHypothesis->Branch("p", &momentum, "p/F");
    treeKaonHypothesis->Branch("beta", &beta_kaon, "beta/F");
    treeProtonHypothesis->Branch("chi2pid_proton", &chi2pid_proton, "chi2pid_proton/F");
    treeProtonHypothesis->Branch("p", &momentum, "p/F");
    treeProtonHypothesis->Branch("beta", &beta_proton, "beta/F");
    treeNegPionHypothesis->Branch("chi2pid_neg_pion", &chi2pid_neg_pion, "chi2pid_neg_pion/F");
    treeNegPionHypothesis->Branch("p", &momentum, "p/F");
    treeNegPionHypothesis->Branch("beta", &beta_neg_pion, "beta/F");
    treeNegKaonHypothesis->Branch("chi2pid_neg_kaon", &chi2pid_neg_kaon, "chi2pid_neg_kaon/F");
    treeNegKaonHypothesis->Branch("p", &momentum, "p/F");
    treeNegKaonHypothesis->Branch("beta", &beta_neg_kaon, "beta/F");

    if (target == "CuSn") {
        treeEBPidPionsCu->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidPionsCu->Branch("p", &momentum, "p/F");
        treeEBPidPionsCu->Branch("beta", &beta, "beta/F");
        treeEBPidKaonsCu->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidKaonsCu->Branch("p", &momentum, "p/F");
        treeEBPidKaonsCu->Branch("beta", &beta, "beta/F");
        treeEBPidProtonsCu->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidProtonsCu->Branch("p", &momentum, "p/F");
        treeEBPidProtonsCu->Branch("beta", &beta, "beta/F");
        treeEBPidNegPionsCu->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidNegPionsCu->Branch("p", &momentum, "p/F");
        treeEBPidNegPionsCu->Branch("beta", &beta, "beta/F");
        treeEBPidNegKaonsCu->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidNegKaonsCu->Branch("p", &momentum, "p/F");
        treeEBPidNegKaonsCu->Branch("beta", &beta, "beta/F");
        treePionHypothesisCu->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treePionHypothesisCu->Branch("p", &momentum, "p/F");
        treePionHypothesisCu->Branch("beta", &beta_pion, "beta/F");
        treePionHypothesisCu->Branch("pid", &pid, "pid/I");
        treeKaonHypothesisCu->Branch("chi2pid_kaon", &chi2pid_kaon, "chi2pid_kaon/F");
        treeKaonHypothesisCu->Branch("p", &momentum, "p/F");
        treeKaonHypothesisCu->Branch("beta", &beta_kaon, "beta/F");
        treeProtonHypothesisCu->Branch("chi2pid_proton", &chi2pid_proton, "chi2pid_proton/F");
        treeProtonHypothesisCu->Branch("p", &momentum, "p/F");
        treeProtonHypothesisCu->Branch("beta", &beta_proton, "beta/F");
        treeNegPionHypothesisCu->Branch("chi2pid_neg_pion", &chi2pid_neg_pion, "chi2pid_neg_pion/F");
        treeNegPionHypothesisCu->Branch("p", &momentum, "p/F");
        treeNegPionHypothesisCu->Branch("beta", &beta_neg_pion, "beta/F");
        treeNegKaonHypothesisCu->Branch("chi2pid_neg_kaon", &chi2pid_neg_kaon, "chi2pid_neg_kaon/F");
        treeNegKaonHypothesisCu->Branch("p", &momentum, "p/F");
        treeNegKaonHypothesisCu->Branch("beta", &beta_neg_kaon, "beta/F");

        treeEBPidPionsSn->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidPionsSn->Branch("p", &momentum, "p/F");
        treeEBPidPionsSn->Branch("beta", &beta, "beta/F");
        treeEBPidKaonsSn->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidKaonsSn->Branch("p", &momentum, "p/F");
        treeEBPidKaonsSn->Branch("beta", &beta, "beta/F");
        treeEBPidProtonsSn->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidProtonsSn->Branch("p", &momentum, "p/F");
        treeEBPidProtonsSn->Branch("beta", &beta, "beta/F");
        treeEBPidNegPionsSn->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidNegPionsSn->Branch("p", &momentum, "p/F");
        treeEBPidNegPionsSn->Branch("beta", &beta, "beta/F");
        treeEBPidNegKaonsSn->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treeEBPidNegKaonsSn->Branch("p", &momentum, "p/F");
        treeEBPidNegKaonsSn->Branch("beta", &beta, "beta/F");
        treePionHypothesisSn->Branch("chi2pid", &chi2pid, "chi2pid/F");
        treePionHypothesisSn->Branch("p", &momentum, "p/F");
        treePionHypothesisSn->Branch("beta", &beta_pion, "beta/F");
        treePionHypothesisSn->Branch("pid", &pid, "pid/I");
        treeKaonHypothesisSn->Branch("chi2pid_kaon", &chi2pid_kaon, "chi2pid_kaon/F");
        treeKaonHypothesisSn->Branch("p", &momentum, "p/F");
        treeKaonHypothesisSn->Branch("beta", &beta_kaon, "beta/F");
        treeProtonHypothesisSn->Branch("chi2pid_proton", &chi2pid_proton, "chi2pid_proton/F");
        treeProtonHypothesisSn->Branch("p", &momentum, "p/F");
        treeProtonHypothesisSn->Branch("beta", &beta_proton, "beta/F");
        treeNegPionHypothesisSn->Branch("chi2pid_neg_pion", &chi2pid_neg_pion, "chi2pid_neg_pion/F");
        treeNegPionHypothesisSn->Branch("p", &momentum, "p/F");
        treeNegPionHypothesisSn->Branch("beta", &beta_neg_pion, "beta/F");
        treeNegKaonHypothesisSn->Branch("chi2pid_neg_kaon", &chi2pid_neg_kaon, "chi2pid_neg_kaon/F");
        treeNegKaonHypothesisSn->Branch("p", &momentum, "p/F");
        treeNegKaonHypothesisSn->Branch("beta", &beta_neg_kaon, "beta/F");
    }

    vector<string> hipoFiles;
    for (const auto& entry : fs::directory_iterator(targetDir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
            hipoFiles.push_back(entry.path().string());
        }
    }
    if (hipoFiles.empty()) {
        cerr << "Error: No .hipo files found in " << targetDir << endl;
        return 1;
    }

    int event_count = 0;
    int maxEvents = 10000;
    int total_EBelectron_count = 0;
    int total_electron_count = 0;
    int total_EBpositron_count = 0;
    int trigger_electron_count = 0;
    int events_with_trigger = 0;
    int skipped_events = 0;

    for (const auto& file : hipoFiles) {
        cout << "Processing file: " << file << endl;
        hipo::reader reader;
        reader.open(file.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        while (reader.next() && event_count < maxEvents) {
            event_count++;
            hipo::event event;
            reader.read(event);

            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::bank EVENT(factory.getSchema("REC::Event"));
            hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
            hipo::bank CHER(factory.getSchema("REC::Cherenkov"));
            hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
            hipo::bank RUN(factory.getSchema("RUN::config"));
            event.getStructure(PART);
            event.getStructure(EVENT);
            event.getStructure(SCIN);
            event.getStructure(CHER);
            event.getStructure(CALO);
            event.getStructure(RUN);

            // Initialize index maps
            IndexMap cherMap = loadMapByIndex(CHER, "pindex");
            IndexMap caloMap = loadMapByIndex(CALO, "pindex");
            IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

            // Extract particles and count electrons/positrons
            vector<ParticleData> particles(PART.getRows());
            for (int i = 0; i < PART.getRows(); ++i) {
                particles[i] = getParticleData(i, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
            }
            for (const auto& p : particles) {
                if (p.pid == 11) {total_EBelectron_count++;}
                if (p.pid == -11) {total_EBpositron_count++;}
            }

            // Torus check
            float torus = RUN.getFloat("torus", 0);
            /* if (event_count <= 10) {
                cout << "Event " << event_count-1 << ": torus = " << torus << endl;
            } */
            if (torus != torusValue) {
                cout << "Warning: Skipped event " << event_count-1 << " with torus = " << torus 
                     << ", expected " << torusValue << " for " << bending << endl;
                skipped_events++;
                continue;
            }

            bool hadTrigger = false;
            string effectiveTarget = target;
            for (const auto& p : particles) {
                if (isSimpleElectron(p)) {
                    if (p.charge == -1) total_electron_count++;
                }

                if (isTriggerElectron(p)) {
                    float sf = p.energy_total / p.p;
                    float mean = getSamplingFractionMean(p.sector, p.energy_total);
                    float sigma = getSamplingFractionSigma(p.sector, p.energy_total);
                    float chi2pid_e = getSamplingFractionNSigma(sf, mean, sigma);
                    if (abs(chi2pid_e) >= 5 || abs(p.status) / 1000 != 2) continue;

                    pair<float, float> vzRange;
                    if (target == "CuSn") {
                        vzRange = getTriggerElectronVzRange("Cu", bending);
                        if (p.vz >= vzRange.first && p.vz < vzRange.second) {
                            effectiveTarget = "Cu";
                            trigger_electron_count++;
                            hadTrigger = true;
                            
                        }
                        vzRange = getTriggerElectronVzRange("Sn", bending);
                        if (p.vz >= vzRange.first && p.vz <= vzRange.second) {
                            effectiveTarget = "Sn";
                            trigger_electron_count++;
                            hadTrigger = true;
                            
                        }
                    } else {
                        vzRange = getTriggerElectronVzRange(target, bending);
                        if (p.vz >= vzRange.first && p.vz <= vzRange.second) {
                            trigger_electron_count++;
                            hadTrigger = true;
                            
                        }
                    }
                }
            }

            if (!hadTrigger) continue;
            events_with_trigger++;

            for (const auto& p : particles) {
                if (abs(p.status) / 2000 != 1) continue;
                if (isSimpleElectron(p)) continue;

                momentum = p.p;
                beta = p.beta;
                pid = p.pid;
                chi2pid = p.chi2pid;
                float beta_calculated = getCalculatedBeta(p, particles);
                // Print both beta values for comparison
    std::cout << "Particle ID: " << p.pid << std::endl;
    std::cout << "Event Builder Beta: " << beta << std::endl;
    std::cout << "Calculated Beta: " << beta_calculated << std::endl;
    
    // Optional: You can also calculate the difference between the two
    std::cout << "Beta Difference: " << (beta - beta_calculated) << std::endl;

                chi2pid = computeChi2pid(p, PION_MASS, beta_pion);
                chi2pid_kaon = computeChi2pid(p, KAON_MASS, beta_kaon);
                chi2pid_proton = computeChi2pid(p, PROTON_MASS, beta_proton);
                chi2pid_neg_pion = computeChi2pid(p, PION_MASS, beta_neg_pion);
                chi2pid_neg_kaon = computeChi2pid(p, KAON_MASS, beta_neg_kaon);

                if (effectiveTarget == "Cu") {
                    if (p.charge > 0) {
                        if (p.pid == 211) treeEBPidPionsCu->Fill();
                        if (p.pid == 321) treeEBPidKaonsCu->Fill();
                        if (p.pid == 2212) treeEBPidProtonsCu->Fill();
                    } else if (p.charge < 0) {
                        if (p.pid == -211) treeEBPidNegPionsCu->Fill();
                        if (p.pid == -321) treeEBPidNegKaonsCu->Fill();
                    }
                    if (chi2pid != 99999.0f) {
                        treePionHypothesisCu->Fill();
                        treeKaonHypothesisCu->Fill();
                        treeProtonHypothesisCu->Fill();
                        treeNegPionHypothesisCu->Fill();
                        treeNegKaonHypothesisCu->Fill();
                    }
                } else if (effectiveTarget == "Sn") {
                    if (p.charge > 0) {
                        if (p.pid == 211) treeEBPidPionsSn->Fill();
                        if (p.pid == 321) treeEBPidKaonsSn->Fill();
                        if (p.pid == 2212) treeEBPidProtonsSn->Fill();
                    } else if (p.charge < 0) {
                        if (p.pid == -211) treeEBPidNegPionsSn->Fill();
                        if (p.pid == -321) treeEBPidNegKaonsSn->Fill();
                    }
                    if (chi2pid != 99999.0f) {
                        treePionHypothesisSn->Fill();
                        treeKaonHypothesisSn->Fill();
                        treeProtonHypothesisSn->Fill();
                        treeNegPionHypothesisSn->Fill();
                        treeNegKaonHypothesisSn->Fill();
                    }
                } else {
                    if (p.charge > 0) {
                        if (p.pid == 211) treeEBPidPions->Fill();
                        if (p.pid == 321) treeEBPidKaons->Fill();
                        if (p.pid == 2212) treeEBPidProtons->Fill();
                    } else if (p.charge < 0) {
                        if (p.pid == -211) treeEBPidNegPions->Fill();
                        if (p.pid == -321) treeEBPidNegKaons->Fill();
                    }
                    if (chi2pid != 99999.0f) {
                        treePionHypothesis->Fill();
                        treeKaonHypothesis->Fill();
                        treeProtonHypothesis->Fill();
                        treeNegPionHypothesis->Fill();
                        treeNegKaonHypothesis->Fill();
                    }
                }
            }
        }
        if (event_count >= maxEvents) break;
    }

    cout << "Total events processed: " << event_count << endl;
    cout << "Events skipped due to torus mismatch: " << skipped_events << endl;
    cout << "Total EB electron count: " << total_EBelectron_count << endl;
    cout << "Total electron count: " << total_electron_count << endl;
    cout << "Total positron count: " << total_EBpositron_count << endl;
    cout << "Total trigger electron count: " << trigger_electron_count << endl;
    cout << "Events with trigger electrons: " << events_with_trigger << endl;

    if (outputFile) {
        outputFile->Write();
        outputFile->Close();
        delete outputFile;
    }
    if (outputFileCu) {
        outputFileCu->Write();
        outputFileCu->Close();
        delete outputFileCu;
    }
    if (outputFileSn) {
        outputFileSn->Write();
        outputFileSn->Close();
        delete outputFileSn;
    }

    return 0;
}