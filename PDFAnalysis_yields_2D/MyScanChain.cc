
//
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"
#include "../CORE/CMS2.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TTreeCache.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <cmath>
#include <cassert>

// LHAPDF
#include "/code/osgcode/imacneill/lhapdf-5.8.9/include/LHAPDF/LHAPDF.h"

MyScanChain::MyScanChain(std::string genPdfName, unsigned int genPdfSubset)
{                   
    // gen set parameters
    genPdfName_     = genPdfName;
    genPdfSubset_   = genPdfSubset;

    // set up LHAPDF
    LHAPDF::setPDFPath("/code/osgcode/imacneill/lhapdf-5.8.9/PDFSets/");
    LHAPDF::initPDFSetM(genset_, genPdfName_);
    LHAPDF::initPDFM(genset_, genPdfSubset_);

}   

//
// Main function
//

int MyScanChain::ScanChain(std::string sampleName, TChain *chain, std::string pdfName)
{

    TObjArray *listOfFiles = chain->GetListOfFiles();
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[MyScanChain::ScanChain] " << sampleName << " is not defined" << std::endl;
        return 1;
    }
    else {
        std::cout << "[MyScanChain::ScanChain] " << sampleName << std::endl;
    }

    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    if (rootdir == 0){
        std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
        if (rootdir){
            std::cout<<"OK: Got Rint:"<<std::endl;
        } else {
            std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
        }
    } 
    rootdir->cd(); 

    //
    // setup pdf stuff
    //

    // cteq6ll is only available in LHpdf format
    if (pdfName != "cteq6ll")   LHAPDF::initPDFSetM(set_, pdfName + ".LHgrid");
    else                        LHAPDF::initPDFSetM(set_, pdfName + ".LHpdf");

    unsigned int nsets = 1 + LHAPDF::numberPDFM(set_);
    if (pdfName == "NNPDF20_as_0116_100" || pdfName == "NNPDF20_as_0122_100") nsets = 5;
    if (pdfName == "NNPDF20_as_0117_100" || pdfName == "NNPDF20_as_0121_100") nsets = 27;
    if (pdfName == "NNPDF20_as_0118_100" || pdfName == "NNPDF20_as_0120_100") nsets = 72;
    if (pdfName == "cteq6mE") nsets = 1;
    if (pdfName == "cteq6ll") nsets = 1;

    // are we calculating the central value of 
    // the observable? 
    // - e.g. no pdf re-weighting needed
    bool doingGenSet = false;
    if (pdfName == genPdfName_) doingGenSet = true;

    //
    // setup histograms
    //

//     Int_t nbins = 40;
//     Float_t min = 0.0;
//     Float_t max = 200.0;
    std::vector <TH2F*> histArr;
    for (unsigned int i = 0; i < nsets; ++i) {
        histArr.push_back(new TH2F(Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i),
				   Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i), 41,-12.5,1012.5, 41,-12.5,1012.5));
    }

    //
    // branch variables
    //

    float   x1          = 0.0;      // momentum fraction for parton1
    float   x2          = 0.0;
    int     id1         = 0;        // pdgid of parton1
    int     id2         = 0;
    float   Q           = 0.0;      // event momentum scale
    float   scale1fb    = 0.0;      // event weight
    float   var1         = 0.0;      // variable to study pdf uncertainty wrt
    float   var2         = 0.0;      // variable to study pdf uncertainty wrt

    //
    // loop over pdf subsets
    //

    for (unsigned int subset = 0; subset < nsets; ++subset)
    {

        std::cout << "doing set, subset: " << set_ << ", " << subset << std::endl;
        LHAPDF::initPDFM(set_, subset);

        //
        // loop over content of sample
        //

        TIter fileIter(listOfFiles);
        while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

            // get the tree
            TFile *f = TFile::Open(currentFile->GetTitle()); 
            assert(f);
            TTree *tree = (TTree*)f->Get("t");
            assert(tree);
            TTreeCache::SetLearnEntries(10);
            tree->SetCacheSize(128*1024*1024);
            cms2.Init(tree);

            //
            // loop over events in file
            //

            ULong64_t nEvents = tree->GetEntries();
            for(ULong64_t event = 0; event < nEvents; ++event) {

                cms2.GetEntry(event);
                x1          = cms2.pdfx1();
                x2          = cms2.pdfx2();
                id1         = cms2.pdfid1();
                id2         = cms2.pdfid2();
                Q           = cms2.pdfQ();
                scale1fb    = cms2.mini_weight();
                var1        = cms2.mg();
                var2        = cms2.ml();

                // check if event passes kinematic selection
                if (!Cuts()) continue;

                // calculate the pdf weight
                double pdf_weight = 1.0;
		//                double experimental_weight = scale1fb * cms2.mini_sltrigeff() * cms2.mini_isrweight();
                double experimental_weight = scale1fb * cms2.mini_dltrigeff() * cms2.mini_isrweight();

                // if looper has been invoked to calculate the same pdf the sample was
                // generated with, then we know the pdf_weight will always be 1.0
                // so no need to actually do any work
                // e.g. this will represent the "central value" of the observable

                if (!doingGenSet) {
                    // generated pdf values
                    double fx1Q0gen = LHAPDF::xfxM(genset_, x1, Q, id1) / x1;
                    double fx2Q0gen = LHAPDF::xfxM(genset_, x2, Q, id2) / x2;
                    // subset pdf values
                    double fx1Qi = LHAPDF::xfxM(set_, x1, Q, id1) / x1;
                    double fx2Qi = LHAPDF::xfxM(set_, x2, Q, id2) / x2;
                    // calculate weight and fill histogram
                    pdf_weight = ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen));
                }

                // inclusive uncertainty, one fixed bin
                // but could equally easily fill with a physical variable...
                histArr[subset]->Fill(var2,var1, pdf_weight);
		//                histArr[subset]->Fill(var2,var1, pdf_weight * experimental_weight);

            } // end loop on events

            delete tree;
            f->Close();
            delete f;

        } // end loop on files in chain

    } // end loop on subsets

    //
    // make sure we're back in the right root dir
    //

    rootdir = gROOT->GetDirectory("root:");
    if (rootdir) rootdir->cd();
    else{
        std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
        rootdir = gROOT->GetDirectory("Rint:");
        if (rootdir){
            std::cout<<"OK, got Rint: "<<std::endl;
            rootdir->cd();
        } else {
            std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
        }
    }

    return 0;

}

// test event selection
bool MyScanChain::Cuts()
{
    return true;
}

