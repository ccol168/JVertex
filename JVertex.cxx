#include "JVertex.h"

std::tuple <float,float,float> JVertex::GetEventPosition() {
	
	double xiqi = 0,yiqi = 0, ziqi = 0;
	double XCb,YCb,ZCb;

	double RecoPositionX, RecoPositionY, RecoPositionZ, Jt0, LogLikelihood, MinuitFlag;

	bool pmtCut = true;

	// Calculate charge barycenter

	double totalcharge = 0.;

	for (int j = 0; j < m_Charge.size(); j++) {

		double x_pmt = m_JUNO_PMTs.GetX(m_PMTID[j]);
		double y_pmt = m_JUNO_PMTs.GetY(m_PMTID[j]);
		double z_pmt = m_JUNO_PMTs.GetZ(m_PMTID[j]);

		double charge = m_Charge[j];
		xiqi += x_pmt*charge;
		yiqi += y_pmt*charge;
		ziqi += z_pmt*charge;

		totalcharge += charge;
		
	}

	//calculating the center of mass coordinates for the first guess 
	XCb = xiqi/totalcharge;
	YCb = yiqi/totalcharge;
	ZCb = ziqi/totalcharge;

	//Minimizing function
	InitializeMinuit(XCb,YCb,ZCb,RecoPositionX,RecoPositionY,RecoPositionZ, Jt0, totalcharge, LogLikelihood, MinuitFlag);

	//std::cout << "Reconstructed event position = ( " << RecoPositionX << " , " << RecoPositionY << " , " << RecoPositionZ << " )" << std::endl; 

	return std::make_tuple ((float) RecoPositionX, (float) RecoPositionY, (float) RecoPositionZ); 
}


void JVertex::FillPDFs() {
	TFile* Time_PDFs = TFile::Open(m_PdfPath.c_str());

	if (!Time_PDFs || Time_PDFs->IsZombie()) {
		delete Time_PDFs;
		throw std::runtime_error(
            "JVertex::FillPDFs - cannot open PDF file: " + m_PdfPath
        );
	}

	for (int i=0; i < 10; i++) {
		TString name = Form("funcPdf%d", i + 1);
		TF1* temp = (TF1*) Time_PDFs->Get(name);

		if (!temp) {
			Time_PDFs -> Close();
			delete Time_PDFs;
			throw std::runtime_error("Invalid name inside the PDF file " + m_PdfPath );
    	}
		
		m_PDFs_perhit[i] = std::unique_ptr<TF1>( static_cast<TF1*>( temp->Clone() ) );
	}


	Time_PDFs -> Close();
	delete Time_PDFs;

	return;
}

double JVertex::NLL (const Double_t *par) {

	double recoX = par[0];
	double recoY = par[1];
	double recoZ = par[2];
	double T0 = par[3];

	double logL = 0.0;

	for (int i = 0; i < m_Times.size(); i++) {

		

		double traveledDistanceReco = sqrt(pow(m_JUNO_PMTs.GetX(m_PMTID[i]) - recoX, 2) + 
									pow(m_JUNO_PMTs.GetY(m_PMTID[i]) - recoY, 2) + pow(m_JUNO_PMTs.GetZ(m_PMTID[i]) - recoZ, 2));
		double ToFReco = traveledDistanceReco * 1e-3 * m_nEff / 2.99792458e8 * 1e9; //ns

		double t_residual_data = m_Times[i] - ToFReco - T0;
		double y_mc = 0;
		
		if      (int(m_Charge[i]) < 2)       y_mc = m_PDFs_perhit[0] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 2)      y_mc = m_PDFs_perhit[1] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 3)      y_mc = m_PDFs_perhit[2] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 4)      y_mc = m_PDFs_perhit[3] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 5)      y_mc = m_PDFs_perhit[4] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 6)      y_mc = m_PDFs_perhit[5] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 7)      y_mc = m_PDFs_perhit[6] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 8)      y_mc = m_PDFs_perhit[7] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) == 9)      y_mc = m_PDFs_perhit[8] -> Eval(t_residual_data);
		else if (int(m_Charge[i]) > 10)       y_mc = m_PDFs_perhit[9] -> Eval(t_residual_data);

		constexpr double EPS = 1e-9;
		double y = std::max(y_mc, 0.0);
		y = std::max(y, EPS);
		logL += std::log(y);

		
	}
	
	return -logL;
}


void JVertex::InitializeMinuit(double StartingPointX, double StartingPointY, double StartingPointZ, 
			double& RecoPositionX, double& RecoPositionY, double& RecoPositionZ, double& Jt0, double totalPE, double& Jlkl, double& Jflag) {

  	const double xyz_range  = 20000.0;
  	const double t0_range  = 500.0;
  	const double errdef     = 0.5;
  	const double step_xyz_s = 1000.0;       //SIMPLEX steps
  	const double step_t0_s  = 10.0;        // ns
  	const double step_xyz_m = 100.0;       //MIGRAD steps
  	const double step_t0_m  = 1.0;        // ns

  	const int    N_ATTEMPTS = 6;
  	const double tol_list[N_ATTEMPTS]    = {0.01, 0.005, 0.0025, 0.0000001, 0.5, 1.0};
  	const int    maxFcn_list[N_ATTEMPTS] = {200000, 400000, 600000, 10000000, 20000, 200000};

  	ROOT::Math::Functor f(this, &JVertex::NLL, 4);

  	// --- First step: SIMPLEX (rough search for minimum region) ---
  	ROOT::Math::Minimizer* minSimplex = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  	minSimplex->SetFunction(f);
  	minSimplex->SetErrorDef(errdef);
  	minSimplex->SetStrategy(1);
  	minSimplex->SetPrecision(1.E-12);
  	minSimplex->SetTolerance(0.1);
  	minSimplex->SetMaxIterations(20000);

  	minSimplex->SetVariable(0, "x",  StartingPointX, step_xyz_s);
  	minSimplex->SetVariable(1, "y",  StartingPointY, step_xyz_s);
  	minSimplex->SetVariable(2, "z",  StartingPointZ, step_xyz_s);
  	minSimplex->SetVariable(3, "t0", 0.0,            step_t0_s);

  	minSimplex->SetVariableLimits(0, -xyz_range,  xyz_range);
  	minSimplex->SetVariableLimits(1, -xyz_range,  xyz_range);
  	minSimplex->SetVariableLimits(2, -xyz_range,  xyz_range);
  	minSimplex->SetVariableLimits(3, -t0_range,  t0_range);

  	minSimplex->Minimize();
  	const double* x_simplex = minSimplex->X();  // seed for MIGRAD

  	// --- Second step: MIGRAD (fine scan around the minimum) ---
  	ROOT::Math::Minimizer* minMigrad = nullptr;
  	bool validity_check = false;
  	double bestF = std::numeric_limits<double>::infinity();
  	double bestP[4] = {x_simplex[0], x_simplex[1], x_simplex[2], x_simplex[3]};

  	for (int iteration = 0; iteration < N_ATTEMPTS && !validity_check; ++iteration) {
  		if (minMigrad) {
      		delete minMigrad;
      		minMigrad = nullptr;
    	}

    	minMigrad = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    	minMigrad->SetFunction(f);
    	minMigrad->SetErrorDef(errdef);
    	minMigrad->SetStrategy(2);
    	minMigrad->SetPrecision(1.E-12);
    	minMigrad->SetTolerance(tol_list[iteration]);
    	minMigrad->SetMaxFunctionCalls(maxFcn_list[iteration]);

    	minMigrad->SetVariable(0, "x",  x_simplex[0], step_xyz_m);
    	minMigrad->SetVariable(1, "y",  x_simplex[1], step_xyz_m);
    	minMigrad->SetVariable(2, "z",  x_simplex[2], step_xyz_m);
    	minMigrad->SetVariable(3, "t0", x_simplex[3], step_t0_m);

    	minMigrad->SetVariableLimits(0, -xyz_range,  xyz_range);
    	minMigrad->SetVariableLimits(1, -xyz_range,  xyz_range);
    	minMigrad->SetVariableLimits(2, -xyz_range,  xyz_range);
    	minMigrad->SetVariableLimits(3, -t0_range,  t0_range);

    	minMigrad->Minimize();

    	int status = minMigrad->Status();
    	double fmin = minMigrad->MinValue();

    	if (fmin < bestF) {
      		//LogInfo << "The fmin has varied " << fmin << std::endl;
      		bestF = fmin;
      		const double* ptmp = minMigrad->X();
      		bestP[0] = ptmp[0]; bestP[1] = ptmp[1]; bestP[2] = ptmp[2]; bestP[3] = ptmp[3];
    	}

    	validity_check = (status == 0 || status == 1);
    	//std::cout << "Validity check = " << validity_check << " with tolerance " << tol_list[iteration] << std::endl;
    	if (validity_check) break;
    
  	}

    	minMigrad->SetVariableValue(0, bestP[0]);
    	minMigrad->SetVariableValue(1, bestP[1]);
    	minMigrad->SetVariableValue(2, bestP[2]);
    	minMigrad->SetVariableValue(3, bestP[3]);

  	minMigrad->Hesse();

  	const double* best = minMigrad->X();
  	RecoPositionX = best[0];
  	RecoPositionY = best[1];
  	RecoPositionZ = best[2];
  	Jt0           = best[3];
  	Jlkl          = minMigrad->MinValue();
  	Jflag         = minMigrad->Status();

  	delete minSimplex;
  	delete minMigrad;
}

