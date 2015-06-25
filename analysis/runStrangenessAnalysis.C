void runStrangenessAnalysis()
{
//    gStyle->SetOptStat(kFALSE);

    gROOT->ProcessLine(".L /Users/macbook/alice/simpleAnalysis/commonTools/QAforWindows.cxx+");
    gROOT->ProcessLine(".L /Users/macbook/alice/simpleAnalysis/analysers/AnalyserBase.cxx+");
//    gROOT->ProcessLine(".L ../../analysers/diHadronMethod/DiHadronAnalyser.cxx+");
    gROOT->ProcessLine(".L /Users/macbook/alice/simpleAnalysis/analysers/tileCorrelations/TileCorrelations.cxx+");

    gROOT->ProcessLine(".L StrangenessAnalysis.cxx+");

    //    int nEvents = 200000;
    StrangenessAnalysis an;
//    an.PrepareAnalysisAndRun();
    an.RunAnalysis(); // nEvents, i*0.05+0.0001, 0 );//, 0.99 );

    //    gROOT->ProcessLine( ".q" );

}


