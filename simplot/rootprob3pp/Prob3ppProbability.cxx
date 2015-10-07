#include "Prob3ppProbability.hxx"

int crootprob3pp::Flavour::frompdg(int pdg)
{
	if(abs(pdg) == 12)
	{
		return NU_E;
	}
	else if(abs(pdg) == 14)
	{
		return NU_MU;
	}
	else if(abs(pdg) == 16)
	{
		return NU_TAU;
	}
	throw std::exception();
}

int crootprob3pp::CP::frompdg(int pdg)
{
	if(pdg < 0)
	{
		return ANTI_MATTER;
	}
	return MATTER;
}


void crootprob3pp::Probability::Init(std::string fname)
{
	//only allow Init to be called once (otherwise there are seg faults)
	static bool hasbeencalled = 0;
	//osc values
	if (!hasbeencalled)
	{
		hasbeencalled  =1;
	}
}

crootprob3pp::Probability::Probability()
{
    bargerprop = new BargerPropagator();
    bargerprop->UseMassEigenstates(0);
    //bargerprop->SetOneMassScaleMode(0);
    //bargerprop->SetWarningSuppression(1);

	double sinSq_2theta23 = 1.0;
	double sinSq_2theta13 = 0.1;
	double sinSq_2theta12 = 0.8495;
	double deltaCP = TMath::Pi()/2.0;
	double delta_MSq32 = 2.4e-3;
	double delta_MSq12 = 7.6e-5;

	//compute angles
	this->theta12 = TMath::ASin(TMath::Sqrt(sinSq_2theta12))/2.0; //in radians
	this->theta13 = TMath::ASin(TMath::Sqrt(sinSq_2theta13))/2.0; //in radians
	this->theta23 = TMath::ASin(TMath::Sqrt(sinSq_2theta23))/2.0; //in radians
	this->deltacp = deltaCP; //in radians
	this->sdm = delta_MSq12; //in eV^2
	this->ldm = delta_MSq32; //in eV^2
	this->length = 295.0; // in km
	this->update();
}

crootprob3pp::Probability::~Probability()
{
	return;
}

void crootprob3pp::Probability::setAll(double theta12_, double theta23_, double theta13_, double deltacp_, double sdm_, double ldm_)
{
	this->theta12 = theta12_;
	this->theta23 = theta23_;
	this->theta13 = theta13_;
	this->deltacp = deltacp_;
	this->sdm = sdm_;
	this->ldm = ldm_;
	touch();
}

void crootprob3pp::Probability::setTheta12(double theta12_)
{
	this->theta12 = theta12_;
	touch();
}

void crootprob3pp::Probability::setTheta23(double theta23_)
{
	this->theta23 = theta23_;
	touch();
}

void crootprob3pp::Probability::setTheta13(double theta13_)
{
	this->theta13 = theta13_;
	touch();
}

void crootprob3pp::Probability::setDeltaCP(double deltacp_)
{
	this->deltacp = deltacp_;
	touch();
}

void crootprob3pp::Probability::setSmallDeltaMassSquared(double sdm_)
{
	this->sdm = sdm_;
	touch();
}

void crootprob3pp::Probability::setLargeDeltaMassSquared(double ldm_)
{
	this->ldm = ldm_;
	touch();
}

void crootprob3pp::Probability::setBaseline(double l)
{
	this->length = l;
}

void crootprob3pp::Probability::update()
{
//    double sinsq_theta12 = sinsq(theta12);
//    double sinsq_theta13 = sinsq(theta13);
//    double sinsq_theta23 = sinsq(theta23);
//    double dm32 = ldm - sdm;
//	this->bargerprop->SetMNS(sinsq_theta12, sinsq_theta13, sinsq_theta23, sdm, dm32, deltacp, 1.0, kSquared, 1);
	this->istouched = 0;
	return;
}

double crootprob3pp::Probability::prob(int initFlavour, int finalFlavour, double energy, int cp)
{
	return this->getVacuumProbability(initFlavour, finalFlavour, energy, cp);
}

double crootprob3pp::Probability::averageprob(int initFlavour, int finalFlavour, double energylow, double energyhigh, int cp, int nsamples)
{
	double psum = 0.0;
	int count = 0;
	for (int i = 0; i < nsamples; ++i)
	{
		double frac = double(i) / double(nsamples);
		double enu = ((1.0 - frac) * energylow) + (frac * energyhigh);
		psum += prob(initFlavour, finalFlavour, enu, cp);
		count += 1;
	}
	return psum / double(count);
}

double crootprob3pp::Probability::getVacuumProbability(int initFlavour, int finalFlavour, double energy, int cp)
{
	/***************************************************************************
	 * Function glbVacuumProbability                                           *
	 ***************************************************************************
	 * Returns the vacuum oscillation probability for one specific oscillation *
	 * channel                                                                 *
	 ***************************************************************************
	 * Parameters:                                                             *
	 *   initial_flavour/final_flavour: The desired oscillation channel        *
	 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
	 *   E: Neutrino energy in GeV                                             *
	 *   L: Baseline in km                                                     *
	 ***************************************************************************/
	if (this->istouched)
	{
		throw std::exception();
	}
    double sinsq_theta12 = sinsq(theta12);
    double sinsq_theta13 = sinsq(theta13);
    double sinsq_theta23 = sinsq(theta23);
    double dm32 = ldm - sdm;
	this->bargerprop->SetMNS(sinsq_theta12, sinsq_theta13, sinsq_theta23, sdm, dm32, this->deltacp, energy, kSquared, cp);
	this->bargerprop->propagateLinear(cp, this->length, 2.6);
	//double p = this->bargerprop->GetVacuumProb(cp*initFlavour, cp*finalFlavour, energy, this->length);
	double p = this->bargerprop->GetProb(cp*initFlavour, cp*finalFlavour);
	return p;
}


