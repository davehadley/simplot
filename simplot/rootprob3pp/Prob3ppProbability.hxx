#ifndef ROOTPROB3PP_PROBABILITY_H
#define ROOTPROB3PP_PROBABILITY_H

#include "BargerPropagator.h"
#include <string>
#include "TMath.h"
#include <stdexcept>

namespace crootprob3pp
{

class Flavour
{
public:
    static const int NU_E = 1;
    static const int NU_MU = 2;
    static const int NU_TAU = 3;

    static int frompdg(int pdg);
};

class CP
{
public:
    static const int MATTER = 1;
    static const int ANTI_MATTER = -1;
    static int frompdg(int pdg);
};

class Probability
{
private:
    BargerPropagator* bargerprop;
    double theta12;
    double theta13;
    double theta23;
    double deltacp;
    double sdm;
    double ldm;
    double length;

    static const bool kSquared = 1;

public:

    static void Init(std::string fname);

    Probability();

    ~Probability();

    void setAll(double theta12_, double theta23_, double theta13_, double deltacp_, double sdm_, double ldm_);

    void setTheta12(double theta12_);

    void setTheta23(double theta23_);

    void setTheta13(double theta13_);

    void setDeltaCP(double deltacp_);

    void setSmallDeltaMassSquared(double sdm_);

    void setLargeDeltaMassSquared(double ldm_);

    void setBaseline(double l);

    void update();

    double prob(int initFlavour, int finalFlavour, double energy, int cp=1);

    double averageprob(int initFlavour, int finalFlavour, double energylow, double energyhigh, int cp=1, int nsamples=10);

    double getVacuumProbability(int initFlavour, int finalFlavour, double energy, int cp=1);

private:
    bool istouched;

    inline void touch()
    {
        this->istouched = 1;
    }

    inline double sinsq(double x) {
        x = TMath::Sin(x);
        return x*x;
    }

};

}

#endif

