#include <cstdlib>
#include <TF1.h>
#include <TFile.h>

class JVertex {

	private :
		std::vector <float> m_Times;
		std::vector <int> m_PMTID;
		float m_nEff;
		std::string m_PdfPath;

		std::array<std::unique_ptr<TF1>, 10> m_PDFs_perhit{};

	public :

		JVertex(std::vector <float> Times, std::vector<int> PMTID, float nEff, std::string PdfPath) :
			m_Times{Times},
			m_PMTID{PMTID},
			m_nEff{nEff},
			m_PdfPath{PdfPath} {return;};

		~JVertex() = default;
		JVertex(JVertex const& other) = delete;
		JVertex(JVertex && other) noexcept = default;
		JVertex& operator=(JVertex const& other) = delete;
		JVertex& operator=(JVertex&& other) noexcept = default;

		
		void FillPDFs();

		void ChangeEvent (std::vector <float> Times, std::vector <int> PMTID) {
			m_Times = Times;
			m_PMTID = PMTID;
			return;
		};

		std::tuple <float,float,float> GetEventPosition(std::vector <float> Times, std::vector <int> PMTID);
};

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
