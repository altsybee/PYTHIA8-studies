//#include "/Users/macbook/alice/simpleAnalysis/analysers/diHadronMethod/DiHadronAnalyser.h"
//#include "/Users/macbook/alice/simpleAnalysis/commonTools/MiniEvent.cxx"
#include "/Users/macbook/alice/simpleAnalysis/commonTools/Tools.cxx"
#include "/Users/macbook/alice/simpleAnalysis/analysers/tileCorrelations/TileCorrelations.h"
//#include "../NuclearStructure/StringDecayer/DecayInTwo.h"
//#include "../NuclearStructure/StringDecayer/StringFragmentation.h"
#include "/Users/macbook/alice/simpleAnalysis/commonTools/SimpleTrack.cxx"

#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include <TArrayF.h>
#include <TObjString.h>


#include "StrangenessAnalysis.h"
#include <iostream>

using namespace std;



//TH2D *fHist2DEtaPhi;
//TH2D *fHist2DEtaPhiWeightMinv;





//const int nAnalysers = 4;

//const double mRho = 0.77;// fRand->Gaus(0.77,0.05);
//const double mPion = 0.14;


const double mRho = 0.775;// fRand->Gaus(0.77,0.05);
const double mRhoWidth = 0.16;
const double mPion = 0.14;




StrangenessAnalysis::StrangenessAnalysis()
{
    //    strFr = new StringFragmentation;
}

void StrangenessAnalysis::extractPidMap(TFile *file, map <int,string> &pidMap)
{
    TArrayF *arrPid = (TArrayF*)file->Get( "arrPid" );
    TH1D *histArrPidNames  = (TH1D*)file->Get( "histArrPidNames" );
    for (int i = 0; i < arrPid->GetSize(); i++)
    {
        pidMap.insert ( pair<int,string>(  arrPid->At(i),histArrPidNames->GetXaxis()->GetBinLabel( i+1 )  ) );
    }

    int iPid = 0;
    for (auto it = pidMap.begin(); it != pidMap.end(); ++it)
    {
        cout << iPid << ": " << (*it).first << " : " << (*it).second << endl;
        iPid++;
    }
}


void StrangenessAnalysis::setVariablesForTrees()
{
    // ##### particle tree
    // particle tree vars
    treeParticles->SetBranchAddress("eventId"     ,&eventId    );
    treeParticles->SetBranchAddress("particleId"  ,&particleId );
    treeParticles->SetBranchAddress("pt"          ,&pt         );
    treeParticles->SetBranchAddress("eta"         ,&eta        );
    treeParticles->SetBranchAddress("phi"         ,&phi        );
    treeParticles->SetBranchAddress("charge"      ,&charge     );
    treeParticles->SetBranchAddress("y"           ,&y     );

    treeParticles->SetBranchAddress("id"          ,&id         );
    treeParticles->SetBranchAddress("status"      ,&status     );
    treeParticles->SetBranchAddress("mother1"     ,&mother1    );
    treeParticles->SetBranchAddress("mother2"     ,&mother2    );
    treeParticles->SetBranchAddress("daughter1"   ,&daughter1  );
    treeParticles->SetBranchAddress("daughter2"   ,&daughter2  );
    treeParticles->SetBranchAddress("px"          ,&px         );
    treeParticles->SetBranchAddress("py"          ,&py         );
    treeParticles->SetBranchAddress("pz"          ,&pz         );
    treeParticles->SetBranchAddress("e"           ,&e          );
    treeParticles->SetBranchAddress("m"           ,&m          );
    treeParticles->SetBranchAddress("hasVertex"   ,&hasVertex  );
    treeParticles->SetBranchAddress("xProd"       ,&xProd      );
    treeParticles->SetBranchAddress("yProd"       ,&yProd      );
    treeParticles->SetBranchAddress("zProd"       ,&zProd      );
    treeParticles->SetBranchAddress("tProd"       ,&tProd      );
    treeParticles->SetBranchAddress("tau"         ,&tau        );

    treeParticles->SetBranchAddress("theta"       ,&theta      );

    treeParticles->SetBranchAddress("nDaughters"  ,&nDaughters);
    treeParticles->SetBranchAddress("nMothers"    ,&nMothers  );
    treeParticles->SetBranchAddress("nSisters"    ,&nSisters  );

    treeParticles->SetBranchAddress("isFinal"   ,&isFinal        );

    treeParticles->SetBranchAddress("isCharged"   ,&isCharged    );
    treeParticles->SetBranchAddress("isNeutral"   ,&isNeutral    );
    treeParticles->SetBranchAddress("tau0"        ,&tau0         );
    treeParticles->SetBranchAddress("mayDecay"    ,&mayDecay     );
    treeParticles->SetBranchAddress("canDecay"    ,&canDecay     );
    treeParticles->SetBranchAddress("isResonance" ,&isResonance  );
    treeParticles->SetBranchAddress("isVisible"   ,&isVisible    );
    treeParticles->SetBranchAddress("isLepton"    ,&isLepton     );
    treeParticles->SetBranchAddress("isQuark"     ,&isQuark      );
    treeParticles->SetBranchAddress("isGluon"     ,&isGluon      );
    treeParticles->SetBranchAddress("isDiquark"   ,&isDiquark    );
    treeParticles->SetBranchAddress("isParton"    ,&isParton     );
    treeParticles->SetBranchAddress("isHadron"    ,&isHadron     );


    // event info tree vars
    treeEvent->SetBranchAddress("nCharged",                &nChargedFromTree );
    //    treeEvent->SetBranchAddress("nPionsPrimary",           &nPionsPrimary );
    //    treeEvent->SetBranchAddress("nPionsFromResonances",    &nPionsFromResonances );
    //    treeEvent->SetBranchAddress("nPionsFromRho0",          &nPionsFromRho0 );
    //    treeEvent->SetBranchAddress("nPionsFromRhoPlusMinus",  &nPionsFromRhoPlusMinus );
    //    treeEvent->SetBranchAddress("pionsFromRho0_bothDetected",&pionsFromRho0_bothDetected );

}


void StrangenessAnalysis::tuneAnalysers() //int nEvents )
{
    // ##### set analysers
    const int nPhiRotations = 1;//16;
    //    double extraPhiShift[nPtBins];

    const int nPtBins = 1;//16;
    double ptEdges[nPtBins][2] = {
        //        { 0.2, 1.0 },
        //        { 1.0, 10.0 },
        //        { 0.2, 1.0 },
        { 0.2, 4.0 },
        //        { 0.2, 0.8 },
        //        { 0.5, 1.5 },
        //        { 0.8, 4.0 },
        //TMP!!!! to be done!
        //        { 0.2, 1.0 },
        //        { 0.5, 2.0 },
        //        { 0.8, 4.0 },
        //        { 1.0, 5.0 },
        //        { 0.8, 4.0 },

    };

    nAnalysers = nPhiRotations*nPtBins;
    analysersArray = new TileCorrelations*[nAnalysers];

    for ( int phiShiftId = 0; phiShiftId < nPhiRotations; phiShiftId++ )
    {
        //        extraPhiShift[anId] = TMath::Pi()/8 * anId;
        double extraPhiShift = TMath::Pi()/8 * phiShiftId;
        cout << "prepare extra phi shift: " << extraPhiShift << endl;

        //        double tileSizeCoeffs[nAnalysers];
        //        for ( int anId = 0; anId < nAnalysers; anId++ )
        //        {
        //            //        tileSizeCoeffs[anId] = 1 * sqrt( anId + 1 );
        //            tileSizeCoeffs[anId] = ( anId + 1 );
        //        }


        for ( int ptBinId = 0; ptBinId < nPtBins; ptBinId++ )
        {
            TileCorrelations *analyser = new TileCorrelations;
            analysersArray[phiShiftId*nPtBins + ptBinId] = analyser;
            analyser->SetShortDef( Form( "list%d", ptBinId ));
            analyser->SetOutputFileId( phiShiftId );
            analyser->SetPtRange( ptEdges[ptBinId][0], ptEdges[ptBinId][1] );
            //        analyser->SetPtRange( ptEdges[0][0], ptEdges[0][1] );

            //        analyser->SetUseEllipseTiles(true);
            analyser->SetUseEllipseTiles(false);

            analyser->SetEtaRangeSymmetric( 0.8 );

            double tileEtaWidth = 0.4;//1.6;//0.8;//0.4 * tileSizeCoeffs[anId];
            double tilePhiWidth = TMath::Pi()/4; // * tileSizeCoeffs[anId];
            //        double tileEtaWidth = 0.1 * tileSizeCoeffs[anId];
            //        double tilePhiWidth = 2*TMath::Pi() * tileSizeCoeffs[anId];
            if ( tilePhiWidth > 2*TMath::Pi() )
                tilePhiWidth = 2*TMath::Pi();
            analyser->SetEtaTileWidth( tileEtaWidth );
            analyser->SetPhiTileWidth( tilePhiWidth );
            int nTilesEta = 1;//7;//0.8 / tileEtaWidth * 4 - 1;
            //        int nTilesPhi = 2*TMath::Pi() / tilePhiWidth * 2;
            int nTilesPhi = 1;//2*TMath::Pi() / tilePhiWidth * 2;
            analyser->SetNumberOfEtaTiles( nTilesEta );
            analyser->SetNumberOfPhiTiles( nTilesPhi );
            analyser->SetExtraPhiShift( extraPhiShift );
            cout << "##### analyzer " << ptBinId << ": " << endl;
            cout << "tileEtaWidth=" << tileEtaWidth << endl;
            cout << "tilePhiWidth=" << tilePhiWidth << endl;
            cout << "nTilesEta=" << nTilesEta << endl;
            cout << "nTilesPhi=" << nTilesPhi << endl;
        }
    }
}

void StrangenessAnalysis::RunAnalysis() //int nEvents )
{
    tuneAnalysers();

    //    fHist2DEtaPhi = new TH2D("fHist2DEtaPhi",";#eta;#phi",20,-2,2,20,-2,2);
    //    fHist2DEtaPhiWeightMinv = new TH2D("fHist2DEtaPhiWeightMinv",";#eta;#phi",20,-2,2,25,-2,TMath::TwoPi());

    TFile *file = new TFile( "/opt/mygit/PYTHIA_resonances/output.root" );

    treeEvent = (TTree*)file->Get( "eventMultTree" );


    // ##### trees
    treeParticles = (TTree*)file->Get( "particleTree" );
    treeEvent = (TTree*)file->Get( "eventMultTree" );
    setVariablesForTrees();

    //    //special tree with pid and names
    //    //    std::vector<int> *vPid = 0;
    //    //    std::vector<string> *vNames = 0;
    //    //    treePidNames = (TTree*)file->Get( "pidNames" );
    //    //    treePidNames->SetBranchAddress("vPid",&vPid);
    //    //    treePidNames->SetBranchAddress("vNames",&vNames);

    //    //    for (int i = 0; i < vPid->size(); i++)
    //    //    {
    //    //        cout << i << ": " << vPid->at(i) << " : " << vNames->at(i) << endl;
    //    //    }
    //    arrPidNames->ls();
    //    for (int i = 0; i < arrPid->GetSize(); i++)
    //    {
    //        cout << i << ": " << arrPid->At(i) << " : " << ((TObjString*)arrPidNames->At(i))->GetString() << endl;
    //    }

    map <int,string> pidMap;
    extractPidMap(file, pidMap);

    //init analysers
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
        analysersArray[iAn]->InitDataMembers();

    TH1D *histEta = new TH1D("histEta", ";#eta;entries", 200, -10, 10);
    TH1D *histPhi = new TH1D("histPhi", ";#phi;entries", 200, -4, 4 );
    TH1D *histPt = new TH1D("histPt", ";p_{T} (GeV/c);entries", 400, 0, 20);


    TH1F *histMult_charged = new TH1F("histMult_charged",";N_{ch}", 400, -0.5, 399.5);
    TH1F *histMult_Lambda0 = new TH1F("histMult_Lambda0",";N_{#Lambda}", 50, -0.5, 49.5);
    TH1F *histMult_K0 = new TH1F("histMult_K0",";N_{K_{0}}", 50, -0.5, 49.5);

    SimpleTrack *simpleTrack = new SimpleTrack;


    int prevEvId = -1;
    int nAnalysedEvents = 0;
    int nentries = treeParticles->GetEntries();
    //    cout << "nentries=" << nentries << endl;

    const int nEventsToAnalyse = 1500;
    int nCharged_inCuts = 0;
    int nLambda0_inCuts = 0;
    int nK0_inCuts = 0;
    //loop over particles in the tree
    for ( int iEntry = 0; iEntry < nentries; iEntry++ )
    {
        treeParticles->GetEntry(iEntry);

        //        cout <<
        //                "iEntry=" << iEntry <<
        //                ", myEventId =" << myEventId <<
        //                ", prevEvId =" << prevEvId <<
        //                ", nCharged =" << nCharged <<
        //             endl;

        if ( prevEvId != eventId ) //i.e. this is the next event
        {
            if ( eventId >= nEventsToAnalyse ) //we analysed enough events
                break;

            // read treeEvent to get info about this new event!
            treeEvent->GetEntry(eventId);

            // !!! cut on event using event info from treeEvent
            //            if ( nChargedFromTree < 4 || nChargedFromTree > 12 ) //if this event is not interesting - go to the next particle! (i.e. wait for next event...)
            //                continue;

            //            histMult_charged->Fill(nChargedFromTree);
            histMult_charged->Fill(nCharged_inCuts);
            histMult_Lambda0->Fill(nLambda0_inCuts);
            histMult_K0->Fill(nK0_inCuts);

            //finish prev event if was opened and start new
            if ( prevEvId >= 0 ) //i.e. we have previous event to be closed
                for ( int iAn = 0; iAn < nAnalysers; iAn++ )
                    analysersArray[iAn]->FinishEvent();

            prevEvId = eventId;
            nAnalysedEvents++;

            // start first event in analysers...
            for ( int iAn = 0; iAn < nAnalysers; iAn++ )
                analysersArray[iAn]->StartEvent();


            //print event counter
            if ( eventId % 100 == 0 )
                printf("Processing %d event...\n", eventId );

            //refresh
            nCharged_inCuts = 0;
            nLambda0_inCuts = 0;
            nK0_inCuts = 0;
        }


        if ( fabs(eta) < 2.5 ) // eta cut
        {
            if (id == 3122)
                ++nLambda0_inCuts;

            if (id == 311)
                ++nK0_inCuts;

            if ( isFinal && isCharged )
            {
                ++nCharged_inCuts;

                simpleTrack->id       = particleId;
                simpleTrack->eta      = eta;
                simpleTrack->phi      = phi;//    fixPhi(track->phi);
                simpleTrack->pt       = pt;
                simpleTrack->charge   = charge;
                simpleTrack->pid      = id;

                for ( int iAn = 0; iAn < nAnalysers; iAn++ )
                    analysersArray[iAn]->AddTrack( simpleTrack );

            }
        }
        histEta->Fill( eta );
        histPhi->Fill( phi );
        histPt->Fill( pt );

    }
    //finish last event
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
        analysersArray[iAn]->FinishEvent();

    //last actions in analysers
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
        analysersArray[iAn]->Terminate();

    printf("nAnalysedEvents=%d\n", nAnalysedEvents );

    TCanvas *canvChecks = new TCanvas("canvChecks","check plots",10,10,600,600 );
    canvChecks->Divide(2,2);
    canvChecks->SetTopMargin(0.05);
    canvChecks->SetRightMargin(0.15);
    canvChecks->SetBottomMargin(0.1);
    canvChecks->SetLeftMargin(0.15);
    //    pad->SetLogz(1);

    //    histPtInvM->DrawCopy();
    //    histPtMixedInvM->SetLineColor( kRed );
    //    histPtMixedInvM->DrawCopy("same");

    //    fHist2DEtaPhi->Draw("colz");
    //    histEta
    canvChecks->cd(1);
    histMult_charged->DrawCopy();//"surf1");

    canvChecks->cd(2);
    histMult_Lambda0->DrawCopy();

    canvChecks->cd(3);
    histMult_K0->DrawCopy();

    TFile *f = new TFile( "analysisResults.root","RECREATE");
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
    {
        //        if ( analysersArray[iAn]->GetOutputFileId() == fileId )
        //        {
        TList *list = analysersArray[iAn]->GetOutputList();
        //        list->ls();
        list->Write( list->GetName(), TObject::kSingleKey );
        cout << " list " << analysersArray[iAn]->GetShortDef() << " is written" << endl;
        //        }
    }
    //    fQAhistos.fHistNch->SetName("hist70-80");
    //    fQAhistos.fHistNch->Write();
    f->Close();

    TCanvas *canvBcorr = new TCanvas("canvBcorr","b plots",200,80,600,600 );
    canvBcorr->SetTheta(50);
    canvBcorr->SetPhi(45);
    canvBcorr->SetTopMargin(0.05);
    canvBcorr->SetRightMargin(0.15);
    canvBcorr->SetBottomMargin(0.1);
    canvBcorr->SetLeftMargin(0.15);
    TH2D *histB = (TH2D *)analysersArray[0]->GetOutputList()->At(2);

    const int phiShiftNbins = 4;
    shiftPhiIn2Dplot(histB,phiShiftNbins);
    histB = flipPlotInEtaAndFixPhiAxis(histB,phiShiftNbins);
    tuneHistogram2D(histB);
    histB->DrawCopy("surf1");


    //    file->Close();

}

