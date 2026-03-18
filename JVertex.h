#include <cstdlib>
#include <TF1.h>
#include <TFile.h>
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "JUNO_PMTs.h"

#ifndef JVERTEX_CXX
#define JVERTEX_CXX

class JVertex {

	private :
		std::vector <float> m_Times;
		std::vector <int> m_PMTID;
		std::vector <float> m_Charge;
		float m_nEff;
		std::string m_PdfPath;
		JUNO_PMTs m_JUNO_PMTs;

		std::array<std::unique_ptr<TF1>, 10> m_PDFs_perhit{};

	public :

		JVertex(float nEff, std::string PdfPath, std::string PMTFile) :
			m_nEff{nEff},
			m_PdfPath{PdfPath},
			m_JUNO_PMTs {PMTFile} {
				JVertex::FillPDFs();
			};

		JVertex(float nEff, std::string PdfPath, std::string PMTFile, std::vector <float> Times,
				 std::vector<int> PMTID, std::vector <float> Charge) :
			m_Times{Times},
			m_PMTID{PMTID},
			m_nEff{nEff},
			m_PdfPath{PdfPath},
			m_Charge{Charge} ,
			m_JUNO_PMTs{PMTFile} {
				JVertex::FillPDFs();
			};

		~JVertex() = default;
		JVertex(JVertex const& other) = delete;
		JVertex(JVertex && other) noexcept = default;
		JVertex& operator=(JVertex const& other) = delete;
		JVertex& operator=(JVertex&& other) noexcept = default;

		void FillPDFs();
		double NLL(const Double_t *par);
		void InitializeMinuit (double StartingPointX, double StartingPointY, double StartingPointZ, double& RecoPositionX,
							 double& RecoPositionY, double& RecoPositionZ, double& Jt0, double totalPE, double& Jlkl, double& Jflag);

		void ChangeEvent (std::vector <float> Times, std::vector <int> PMTID, std::vector <float> Charge) {
			m_Times = Times;
			m_PMTID = PMTID;
			m_Charge = Charge;
			return;
		};

		std::tuple <float,float,float> GetEventPosition();
		
};

#endif

