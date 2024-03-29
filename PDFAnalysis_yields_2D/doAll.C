
void doAll() {

    //
    // the looper
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
    gSystem->Load("/code/osgcode/imacneill/lhapdf-5.8.9/lib/.libs/libLHAPDF.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");

    //
    // create looper
    //

    // create a looper for a sample that was generated by using
    // CTEQ6LL.  Check if this is the case for your sample!
    // generally, LO samples are generated with CTEQ6LL and 
    // NLO samples are generated with CT10.  Note that in the 
    // CT10 case the central subset is at index 5, not 0.
    MyScanChain *looper_cteq6ll = new MyScanChain("cteq6ll.LHpdf", 0);

    //
    // run all pdf sets
    //

    std::vector<std::string> pdfSets;

    // CT10
    pdfSets.push_back("CT10");
    pdfSets.push_back("CT10as");
    // MSTW
    pdfSets.push_back("MSTW2008nlo68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz+68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz+68clhalf");
    pdfSets.push_back("MSTW2008nlo68cl_asmz-68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz-68clhalf");
    // NNPDF
    pdfSets.push_back("NNPDF20_as_0116_100");
    pdfSets.push_back("NNPDF20_as_0117_100");
    pdfSets.push_back("NNPDF20_as_0118_100");
    pdfSets.push_back("NNPDF20_100");
    pdfSets.push_back("NNPDF20_as_0120_100");
    pdfSets.push_back("NNPDF20_as_0121_100");
    pdfSets.push_back("NNPDF20_as_0122_100");

    // data sample
    TChain *chain_ttbar = new TChain("Events");
    chain_ttbar->Add("/nfs-6/userdata/vimartin/mariasigtest/scan/skim/superskim_1l4b.root");
 
    // do gensets 
    // the variation of the genset with respect to itself
    // is by definition 1.0, but we'll need this number later
    // stored just like all the others
    looper_cteq6ll->ScanChain("ttbar", chain_ttbar, "cteq6ll");

    // do other sets
    for (unsigned int i = 0; i < pdfSets.size(); ++i) {
        std::cout << "===== Doing =====> " << pdfSets[i] << std::endl;
        looper_cteq6ll->ScanChain("ttbar", chain_ttbar, pdfSets[i]);
    }

    //
    // write histograms
    // 

    const std::string outFile = "results_1l4b.root";
    saveHist(outFile.c_str());
    deleteHistos();

    //
    // tidy up
    //

    delete looper_cteq6ll;
    delete chain_ttbar;
}

