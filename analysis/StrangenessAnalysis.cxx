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
#include "TProfile.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include <TArrayF.h>
#include <TObjString.h>


#include "StrangenessAnalysis.h"
#include <iostream>

using namespace std;


map <int,string> mapPID;
map <string,double> mapTau0_byName;

//TH2D *fHist2DEtaPhi;
//TH2D *fHist2DEtaPhiWeightMinv;





//const int nAnalysers = 4;

//const double mRho = 0.77;// fRand->Gaus(0.77,0.05);
//const double mPion = 0.14;

/*
 * Find out the decay chain of the final state particle after the hadronization
 * This function should be only used for final state particles, since it searches for the decay chain recursively
 * If used for non-final particles, it will find a part of a decay chain for some final state particle
 *
 *  Pythia8::Event &event = PYTHIA8 event information for the current event
 *  Pythia8::Particle &currentParticle = PYTHIA8 particle information for the considered particle
 *  TString decayChain = The decay chain of the particle until this point
 */

TString StrangenessAnalysis::getDecayChain(Particle &currentParticle, TString decayChain )
{
    // The mother of the current particle is born before hadronization. Stop the decay chain here.
    if(partArr[currentParticle.mother1].status > -80)
    {
        return decayChain;
    }

    // The particles that are formed after hadronization can have only one mother. If we are not looking
    // at the particles before the hadronization, mother1 is the only mother we need to care about
    TString newAddition = Form( "|%d",partArr[currentParticle.mother1].id );
    decayChain += newAddition;
    return getDecayChain(partArr[currentParticle.mother1], decayChain);
}



void StrangenessAnalysis::getDecayChainById(Particle &part)//, TString decayChain )
{

    //    int curPartPID = part.id;   // << ", name: " << mapPID[part.id]
    //    int iMother = part.mother1;
    //    int motherPID = partArr[part.mother1].id;

    bool flagPrint = 1;

    //    if (1) cout << mapPID[part.id] << "  <---   ";

    if (flagPrint) cout << "particle: " << part.id << ", name: " << mapPID[part.id]
                        << ", tau = " << part.tau << ", tau0 = " << part.tau0;
    if ( part.nDaughters > 1 )
    {
        //        for (int i = 0; i < part.nDaughters; i++)
        //        {
        if (flagPrint) cout << ", nDaughters=" << part.nDaughters;
        int d1 = part.daughter1;
        int d2 = part.daughter2;
        if (flagPrint) cout << "     ---  d1: " << partArr[d1].id << ", name: " << mapPID[partArr[d1].id];
        if (flagPrint) cout << "     ---  d2: " << partArr[d2].id << ", name: " << mapPID[partArr[d2].id];
        //        }
    }

    if (flagPrint) cout << endl;
    //         << "  ---  mother: " << partArr[part.mother1].id << ", name: " << mapPID[partArr[part.mother1].id]
    //            //<< endl;
    // //                cout << "decay chain: " <<getDecayChain(part) << endl;
    //    //cout
    //            << "  ---  grandmother: " << partArr[partArr[part.mother1].mother1].id << ", name: " << mapPID[partArr[partArr[part.mother1].mother1].id] << endl;


    if (part.isParton )// || part.isDiquark)
    {
        if (flagPrint) cout << "####################" << endl;
        return;
    }

    getDecayChainById(partArr[part.mother1]);
}

bool StrangenessAnalysis::hasPrimaryMotherByALICEdefinition(Particle &part)
{
    if (partArr[part.mother1].tau0 > 10) // 1cm
        return false;

    if (part.isParton )// || part.isDiquark)
    {
        //        if (flagPring) cout << "####################" << endl;
        return true;
    }

    return hasPrimaryMotherByALICEdefinition(partArr[part.mother1]);
}

int idInArray(int *arr, int size, int lookFor)
{
    for ( int i = 0; i < size; i++ )
        if ( arr[i] == lookFor )
            return i;
    return -1;
}

StrangenessAnalysis::StrangenessAnalysis()
{
    //    strFr = new StringFragmentation;
}

void StrangenessAnalysis::extractPidMap(TFile *file)
{
    TArrayF *arrPid = (TArrayF*)file->Get( "arrPid" );
    TH1D *histArrPidNames  = (TH1D*)file->Get( "histArrPidNames" );
    for (int i = 0; i < arrPid->GetSize(); i++)
    {
        mapPID.insert ( pair<int,string>(  arrPid->At(i),histArrPidNames->GetXaxis()->GetBinLabel( i+1 )  ) );
    }

    int iPid = 0;
    for (auto it = mapPID.begin(); it != mapPID.end(); ++it)
    {
        cout << iPid << ": " << (*it).first << " : " << (*it).second << endl;
        iPid++;
    }
}


void StrangenessAnalysis::setVariablesForTrees()
{
    // ##### particle tree
    // particle tree vars
    treeParticles->SetBranchAddress("eventId"     ,&part.eventId    );
    treeParticles->SetBranchAddress("particleId"  ,&part.particleId );
    treeParticles->SetBranchAddress("pt"          ,&part.pt         );
    treeParticles->SetBranchAddress("eta"         ,&part.eta        );
    treeParticles->SetBranchAddress("phi"         ,&part.phi        );
    treeParticles->SetBranchAddress("charge"      ,&part.charge     );
    treeParticles->SetBranchAddress("y"           ,&part.y     );

    treeParticles->SetBranchAddress("id"          ,&part.id         );
    treeParticles->SetBranchAddress("status"      ,&part.status     );
    treeParticles->SetBranchAddress("mother1"     ,&part.mother1    );
    treeParticles->SetBranchAddress("mother2"     ,&part.mother2    );
    treeParticles->SetBranchAddress("daughter1"   ,&part.daughter1  );
    treeParticles->SetBranchAddress("daughter2"   ,&part.daughter2  );
    //    treeParticles->SetBranchAddress("px"          ,&part.px         );
    //    treeParticles->SetBranchAddress("py"          ,&part.py         );
    //    treeParticles->SetBranchAddress("pz"          ,&part.pz         );
    //    treeParticles->SetBranchAddress("e"           ,&part.e          );
    //    treeParticles->SetBranchAddress("m"           ,&part.m          );
    //    treeParticles->SetBranchAddress("hasVertex"   ,&part.hasVertex  );
    //    treeParticles->SetBranchAddress("xProd"       ,&part.xProd      );
    //    treeParticles->SetBranchAddress("yProd"       ,&part.yProd      );
    //    treeParticles->SetBranchAddress("zProd"       ,&part.zProd      );
    //    treeParticles->SetBranchAddress("tProd"       ,&part.tProd      );
    treeParticles->SetBranchAddress("tau"         ,&part.tau        );

    treeParticles->SetBranchAddress("theta"       ,&part.theta      );

    treeParticles->SetBranchAddress("nDaughters"  ,&part.nDaughters);
    treeParticles->SetBranchAddress("nMothers"    ,&part.nMothers  );
    treeParticles->SetBranchAddress("nSisters"    ,&part.nSisters  );

    treeParticles->SetBranchAddress("isFinal"   ,&part.isFinal        );

    treeParticles->SetBranchAddress("isCharged"   ,&part.isCharged    );
    //    treeParticles->SetBranchAddress("isNeutral"   ,&part.isNeutral    );
    treeParticles->SetBranchAddress("tau0"        ,&part.tau0         );
    //    treeParticles->SetBranchAddress("mayDecay"    ,&part.mayDecay     );
    //    treeParticles->SetBranchAddress("canDecay"    ,&part.canDecay     );
    //    treeParticles->SetBranchAddress("isResonance" ,&part.isResonance  );
    //    treeParticles->SetBranchAddress("isVisible"   ,&part.isVisible    );
    //    treeParticles->SetBranchAddress("isLepton"    ,&part.isLepton     );
    //    treeParticles->SetBranchAddress("isQuark"     ,&part.isQuark      );
    //    treeParticles->SetBranchAddress("isGluon"     ,&part.isGluon      );
    //    treeParticles->SetBranchAddress("isDiquark"   ,&part.isDiquark    );
    treeParticles->SetBranchAddress("isParton"    ,&part.isParton     );
    treeParticles->SetBranchAddress("isHadron"    ,&part.isHadron     );


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

void StrangenessAnalysis::RunAnalysis(int nEventsToAnalyse )
{
    tuneAnalysers();

    //    fHist2DEtaPhi = new TH2D("fHist2DEtaPhi",";#eta;#phi",20,-2,2,20,-2,2);
    //    fHist2DEtaPhiWeightMinv = new TH2D("fHist2DEtaPhiWeightMinv",";#eta;#phi",20,-2,2,25,-2,TMath::TwoPi());

    //    TString strTuneName("100k_ct1cm_Monash2013_pure");
    //    TString strTuneName("10kTEST_ct1cm_Monash2013_newCR_mode0");

    //        TString strTuneName("100k_ct1cm_Monash2013_newCR_mode0");
    //        TString strTuneName("100k_ct1cm_Monash2013_newCR_mode2");

//    TString strTuneName("1000k_ct1cm_Monash2013_pure");
//    TString strTuneName("1000k_ct1cm_Monash2013_newCR_mode2");
            TString strTuneName("1000k_ct1cm_Monash2013_noCR");

    //    TFile *file = new TFile( "/opt/mygit/PYTHIA_resonances/output10k.root" );
    //    TFile *file = new TFile( "/opt/mygit/PYTHIA_resonances/output10k_ct1cm_Monash2013.root" );
    //    TFile *file = new TFile( "/opt/mygit/PYTHIA_resonances/output100k_ct1cm_Monash2013_newCRmode0.root" );
    //    TFile *file = new TFile( "/opt/mygit/PYTHIA_resonances/output100k_ct1cm_Monash2013_pure.root" );
    //    TFile *file = new TFile( "/opt/mygit/PYTHIA_resonances/output100k_ct1cm_Monash2013_noCR.root" );
    TFile *file = new TFile( Form("/opt/mygit/PYTHIA_resonances/output%s.root", strTuneName.Data()) );

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

    extractPidMap(file);

    //init analysers
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
        analysersArray[iAn]->InitDataMembers();

    TH1D *histEta = new TH1D("histEta", ";#eta;entries", 200, -10, 10);
    TH1D *histPhi = new TH1D("histPhi", ";#phi;entries", 200, -4, 4 );
    TH1D *histPt = new TH1D("histPt", ";p_{#rm T} (GeV/c);entries", 400, 0, 20);

    TH1D *histEtaStrange = new TH1D("histEtaStrange", ";#eta;entries", 200, -10, 10);
    TH1D *histYStrange = new TH1D("histYStrange", ";#eta;entries", 200, -10, 10);


    TH1D *histMult_charged = new TH1D("histMult_charged",";N_{ch}", 400, -0.5, 399.5);
    TH1D *histMult_Lambda0 = new TH1D("histMult_Lambda0",";N_{#Lambda}", 50, -0.5, 49.5);
    TH1D *histMult_K0 = new TH1D("histMult_K0",";N_{K_{0}}", 50, -0.5, 49.5);

    TH1D *histLambda0_daugthers_dEta = new TH1D( "histLambda0_daugthers_dEta", "histLambda0_daugthers_dEta;#delta#eta;entries", 100, 0, 5 );

    //    TH2D *histNN = new TH2D("histNN","histNN", 300, -0.5, 299.5, 200, 0, 10);
    //    TH2D *histPtN = new TH2D("histPtN","histPtN", 300, -0.5, 299.5, 200, 0, 10);



    //for correlations:
    const int nRegularWins = 7;
    const int nSpecWins = 7;
    const int nPairs = nRegularWins+nSpecWins;//9;
    double etaMaxF[nPairs];
    double etaMinF[nPairs];

    double etaMaxB[nPairs];
    double etaMinB[nPairs];

    //arrays to set windows explicitly
    //    double etaMaxF[nPairs] = {  0.8, 0.8,   2.4,  2.4,  2.4,  };
    //    double etaMinF[nPairs] = { -0.8,   0,  -2.4,    0,  1.0,  };

    //    double etaMaxB[nPairs] = {  0.8,    0,  2.4,    0, -1.0,  };
    //    double etaMinB[nPairs] = { -0.8, -0.8, -2.4, -2.4, -2.4,  };
    double arrEtaMaxF[nSpecWins] = {  0.8, 0.8,   2.4,  2.4,  2.4,  2.0,  0.5  };
    double arrEtaMinF[nSpecWins] = { -0.8,   0,  -2.4,    0,  0.5,  1.0, -0.5  };

    double arrEtaMaxB[nSpecWins] = {  0.8,    0,  2.4,    0, -0.5, -1.0,  0.5  };
    double arrEtaMinB[nSpecWins] = { -0.8, -0.8, -2.4, -2.4, -2.4, -2.0, -0.5  };


    // ##### correlation histograms
    // NN
    TH2D *histNN[nPairs];       // N vs N
    TH2D *histNsNs[nPairs];     // Nstrange vs Nstrange
    TH2D *histNsN[nPairs];      // Nstrange vs N
    TH2D *histNskN[nPairs];     // Nstrange+Kaons+/- vs N

    TH2D *histNlN[nPairs];      // Nlambda vs N
    TH2D *histNlNl[nPairs];     // Nlambda vs Nlambda
    TH2D *histNprotonN[nPairs];

    // PtN
    TH2D *histPtN[nPairs];
    TH2D *histPtNs[nPairs];
    TH2D *histPtNsk[nPairs];
    TH2D *histPtNl[nPairs];

    // ratio correlations
    TH2D *histKtoPion_N[nPairs];
    TH2D *histLambdaToPion_N[nPairs];
    TH2D *histProtonToPion_N[nPairs];
    TH2D *histKtoPion_KtoPion[nPairs];
    TH2D *histLambdaToPion_LambdaToPion[nPairs];

    TH2D *histKtoPion_Pt[nPairs];
    TH2D *histLambdaToPion_Pt[nPairs];
    TH2D *histProtonToPion_Pt[nPairs];

    //PtPt
    TH2D *histPtPt[nPairs];
    //EtEt
//    TH2D *histEtEt[nPairs];

    TH2D *histWindowsEtaBounds;       // to remember eta positions of windows
    histWindowsEtaBounds = new TH2D( "histWindowsEtaBounds", ";#eta bounds;n_{win pair}", nPairs, -0.5, nPairs-0.5, 4, -0.5, 3.5);

    const double winEtaSize = 0.5;
    const double winEtaStep = 0.25;
    const double winEtaStartEdge = 0;
    for (int i = 0; i < nPairs; ++i)
    {
        if ( i < nRegularWins ) // symmetric eta wins with increasing gap
        {
            etaMaxF[i] = winEtaStartEdge + i*winEtaStep+winEtaSize;
            etaMinF[i] = winEtaStartEdge + i*winEtaStep;
            etaMaxB[i] = -winEtaStartEdge + -i*winEtaStep;
            etaMinB[i] = -winEtaStartEdge + -i*winEtaStep-winEtaSize;
        }
        else
        {
            etaMaxF[i] = arrEtaMaxF[i-nRegularWins];
            etaMinF[i] = arrEtaMinF[i-nRegularWins];
            etaMaxB[i] = arrEtaMaxB[i-nRegularWins];
            etaMinB[i] = arrEtaMinB[i-nRegularWins];
        }
        histWindowsEtaBounds->SetBinContent(i+1, 3+1, etaMaxF[i]);
        histWindowsEtaBounds->SetBinContent(i+1, 2+1, etaMinF[i]);
        histWindowsEtaBounds->SetBinContent(i+1, 1+1, etaMaxB[i]);
        histWindowsEtaBounds->SetBinContent(i+1, 0+1, etaMinB[i]);

        cout << "bkwd window bounds: " << etaMinB[i] << " " << etaMaxB[i]
                << ", fwd window bounds: " << etaMinF[i] << " " << etaMaxF[i]
                   << endl;

        histNN[i] = new TH2D( Form("histNN%d", i), ";n_{F};n_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);
        histNsNs[i] = new TH2D( Form("histNsNs%d", i), ";ns_{F};ns_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);
        histNsN[i] = new TH2D( Form("histNsN%d", i), ";n_{F};ns_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);
        histNskN[i] = new TH2D( Form("histNskN%d", i), ";n_{F};nsk_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);

        histNlN[i] = new TH2D( Form("histNlN%d", i), ";n_{F};nl_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);
        histNlNl[i] = new TH2D( Form("histNlNl%d", i), ";nl_{F};nl_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);
        histNprotonN[i] = new TH2D( Form("histNprotonN%d", i), ";n_{F};nProton_{B}", 200, -0.5, 199.5, 200, -0.5, 199.5);

        histPtN[i] = new TH2D( Form("histPtN%d", i), ";n_{F};pT_{B}", 200, -0.5, 199.5, 250, 0, 5);
        histPtNs[i] = new TH2D( Form("histPtNs%d", i), ";ns_{F};pT_{B}", 200, -0.5, 199.5, 500, 0, 5);
        histPtNsk[i] = new TH2D( Form("histPtNsk%d", i), ";nsk_{F};pT_{B}", 200, -0.5, 199.5, 500, 0, 5);
        histPtNl[i] = new TH2D( Form("histPtNl%d", i), ";nl_{F};pT_{B}", 200, -0.5, 199.5, 500, 0, 5);

        histKtoPion_N[i] = new TH2D( Form("histKtoPion_N%d", i), ";n_{F};KtoPion_{B}", 200, -0.5, 199.5, 200, 0, 2);
        histLambdaToPion_N[i] = new TH2D( Form("histLambdaToPion_N%d", i), ";n_{F};LambdaToPion_{B}", 200, -0.5, 199.5, 200, 0, 2);
        histProtonToPion_N[i] = new TH2D( Form("histProtonToPion_N%d", i), ";n_{F};ProtonToPion_{B}", 200, -0.5, 199.5, 200, 0, 2);
        histKtoPion_KtoPion[i] = new TH2D( Form("histKtoPion_KtoPion%d", i), ";KtoPion_{F};KtoPion_{B}", 200, 0, 2, 200, 0, 2);
        histLambdaToPion_LambdaToPion[i] = new TH2D( Form("histLambdaToPion_LambdaToPion%d", i), ";LambdaToPion_{F};LambdaToPion_{B}", 200, 0, 2, 200, 0, 2);

        histKtoPion_Pt[i] = new TH2D( Form("histKtoPion_Pt%d", i), ";pT_{F};KtoPion_{B}", 100, 0, 5, 200, 0, 2);
        histLambdaToPion_Pt[i] = new TH2D( Form("histLambdaToPion_Pt%d", i), ";pT_{F};LambdaToPion_{B}", 100, 0, 5, 200, 0, 2);
        histProtonToPion_Pt[i] = new TH2D( Form("histProtonToPion_Pt%d", i), ";pT_{F};ProtonToPion_{B}", 100, 0, 5, 200, 0, 2);

        histPtPt[i] = new TH2D( Form("histPtPt%d", i), ";pT_{F};pT_{B}", 250, 0, 5, 250, 0, 5);
//        histEtEt[i] = new TH2D( Form("histEtEt%d", i), ";ET_{F};ET_{B}", 250, 0, 5, 250, 0, 5);

    }


    SimpleTrack *simpleTrack = new SimpleTrack;


    int prevEvId = -1;
    int nAnalysedEvents = 0;
    int nentries = treeParticles->GetEntries();
    //    cout << "nentries=" << nentries << endl;


    //    const int nEventsToAnalyse = 1000;
    int nCharged_inCuts = 0;
    int nLambda0_inCuts = 0;
    int nK0_inCuts = 0;

    const int nSpecies = 14;
    int arrParticleSpecies[nSpecies] =
    {
        211 , // pi+
        -211 , // pi-

        321 , // K+
        -321 , // K-


        310 , // K_S0
        //        313 , // K*0

        333 , // phi

        2212, // p+
        -2212, // pbar-

        -3122, // Lambdabar0
        3122 , // Lambda0

        3312  , // Xi-
        -3312 , // Xibar+

        3334   , //Omega-
        -3334  , //  Omegabar+

        //        -3322 , // Xibar0

        //        3314  , // Xi*-
        //        3322  , // Xi0
        //        3324  , // Xi*0


        //        213 , // rho+
        //        221 , // eta
        //        223 , // omega

    };
    //    int arrParticleCounterBySpecies[nSpecies];

    TH1D *histPt_species[nSpecies];
    TH1D *histEta_species[nSpecies];
    TH1D *histY_species[nSpecies];

    TH1D *histEta_species_in_windows[nSpecies];
    TH1D *histY_species_in_windows[nSpecies];

    TH1D *histYields = new TH1D("histYields",";;", nSpecies, -0.5, nSpecies-0.5);
    for(int i = 0; i < nSpecies; i++)
    {
        histYields->GetXaxis()->SetBinLabel( i+1, mapPID[arrParticleSpecies[i]].c_str() ) ;
        histPt_species[i] = new TH1D(Form("histPt_%s", mapPID[arrParticleSpecies[i]].c_str())
                , ";p_{#rm T} (GeV/c);entries", 40, 0, 8);
        histEta_species[i] = new TH1D(Form("histEta_%s", mapPID[arrParticleSpecies[i]].c_str())
                , ";#eta;entries", 200, -10, 10);
        histY_species[i] = new TH1D(Form("histY_%s", mapPID[arrParticleSpecies[i]].c_str())
                , ";y;entries", 200, -10, 10);

        histEta_species_in_windows[i] = new TH1D(Form("histEta_in_wins_%s", mapPID[arrParticleSpecies[i]].c_str())
                , ";#eta;entries", 200, -10, 10);
        histY_species_in_windows[i] = new TH1D(Form("histY_in_wins_%s", mapPID[arrParticleSpecies[i]].c_str())
                , ";y;entries", 200, -10, 10);
    }

//    TH1D *histRatios = new TH1D("histRatios",";;", 4, -0.5, 4-0.5);


    //    TH1D *hist_Proton_to_pion_vs_pt = new TH1D("hist_Proton_to_pion_vs_pt",";p_{\rm T};", 0, 0, 12);
    //    TH1D *hist_Kaon_to_pion_vs_pt = new TH1D("hist_Kaon_to_pion_vs_pt",";p_{\rm T};", 0, 0, 12);
    //    TH1D *hist_Lambda0_to_K0s_vs_pt = new TH1D("hist_Lambda0_to_K0s_vs_pt",";p_{\rm T};", 0, 0, 12);




    //    0: -20213 : a_1(1260)-
    //    1: -5122 : Lambda_bbar0
    //    2: -4232 : Xi_cbar-
    //    3: -4122 : Lambda_cbar-
    //    4: -3324 : Xi*bar0
    //    5: -3322 : Xibar0
    //    6: -3314 : Xi*bar+
    //    7: -3312 : Xibar+
    //    8: -3224 : Sigma*bar-
    //    9: -3222 : Sigmabar-
    //    10: -3214 : Sigma*bar0
    //    11: -3212 : Sigmabar0
    //    12: -3122 : Lambdabar0
    //    13: -3114 : Sigma*bar+
    //    14: -3112 : Sigmabar+
    //    15: -2224 : Deltabar--
    //    16: -2214 : Deltabar-
    //    17: -2212 : pbar-
    //    18: -2114 : Deltabar0
    //    19: -2112 : nbar0
    //    20: -1114 : Deltabar+
    //    21: -533 : B*_sbar0
    //    22: -531 : B_sbar0
    //    23: -523 : B*-
    //    24: -521 : B-
    //    25: -433 : D*_s-
    //    26: -431 : D_s-
    //    27: -423 : D*bar0
    //    28: -421 : Dbar0
    //    29: -413 : D*-
    //    30: -411 : D-
    //    31: -323 : K*-
    //    32: -321 : K-
    //    33: -313 : K*bar0
    //    34: -311 : Kbar0
    //    35: -213 : rho-
    //    36: -211 : pi-
    //    37: -16 : nu_taubar
    //    38: -15 : tau+
    //    39: -14 : nu_mubar
    //    40: -13 : mu+
    //    41: -12 : nu_ebar
    //    42: -11 : e+
    //    43: -5 : bbar
    //    44: -4 : cbar
    //    45: -3 : sbar
    //    46: -2 : ubar
    //    47: -1 : dbar
    //    48: 1 : d
    //    49: 2 : u
    //    50: 3 : s
    //    51: 4 : c
    //    52: 5 : b
    //    53: 11 : e-
    //    54: 12 : nu_e
    //    55: 13 : mu-
    //    56: 14 : nu_mu
    //    57: 16 : nu_tau
    //    58: 21 : g
    //    59: 22 : gamma
    //    60: 90 : system
    //    61: 111 : pi0
    //    62: 113 : rho0
    //    63: 130 : K_L0
    //    64: 211 : pi+
    //    65: 213 : rho+
    //    66: 221 : eta
    //    67: 223 : omega
    //    68: 310 : K_S0
    //    69: 311 : K0
    //    70: 313 : K*0
    //    71: 315 : K*_2(1430)0
    //    72: 321 : K+
    //    73: 323 : K*+
    //    74: 331 : eta'
    //    75: 333 : phi
    //    76: 411 : D+
    //    77: 413 : D*+
    //    78: 421 : D0
    //    79: 423 : D*0
    //    80: 431 : D_s+
    //    81: 433 : D*_s+
    //    82: 511 : B0
    //    83: 521 : B+
    //    84: 523 : B*+
    //    85: 990 : Pomeron
    //    86: 1114 : Delta-
    //    87: 2101 : ud_0
    //    88: 2103 : ud_1
    //    89: 2112 : n0
    //    90: 2114 : Delta0
    //    91: 2203 : uu_1
    //    92: 2212 : p+
    //    93: 2214 : Delta+
    //    94: 2224 : Delta++
    //    95: 3112 : Sigma-
    //    96: 3114 : Sigma*-
    //    97: 3122 : Lambda0
    //    98: 3212 : Sigma0
    //    99: 3214 : Sigma*0
    //    100: 3222 : Sigma+
    //    101: 3224 : Sigma*+
    //    102: 3312 : Xi-
    //    103: 3314 : Xi*-
    //    104: 3322 : Xi0
    //    105: 3324 : Xi*0
    //    106: 4114 : Sigma*_c0
    //    107: 4122 : Lambda_c+
    //    108: 5122 : Lambda_b0
    //    109: 5224 : Sigma*_b+
    //    110: 10441 : chi_0c
    //    111: 20213 : a_1(1260)+
    //    112: 20313 : K_1(1400)0
    //    113: 9902210 : p_diffr+


    float nChargedF[nPairs];
    float nChargedB[nPairs];

    float nChargedFs[nPairs];
    float nChargedBs[nPairs];

    float nChargedFsk[nPairs];
    float nChargedBsk[nPairs];

    float nChargedFlambda[nPairs];
    float nChargedBlambda[nPairs];

    float nChargedFproton[nPairs];
    float nChargedBproton[nPairs];

    float nChargedFkaon[nPairs];
    float nChargedBkaon[nPairs];

    float nChargedFpion[nPairs];
    float nChargedBpion[nPairs];



    double meanPt_F[nPairs];
    double meanPt_B[nPairs];




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

        if ( prevEvId != part.eventId ) //i.e. this is the next event
        {
            //note: last event will not go here! i.e. we don't "finish" this event (keep it like this?..)

            if ( part.eventId >= nEventsToAnalyse ) //we analysed enough events
                break;

            // read treeEvent to get info about this new event!
            treeEvent->GetEntry(part.eventId);

            // !!! cut on event using event info from treeEvent
            //            if ( nChargedFromTree < 4 || nChargedFromTree > 12 ) //if this event is not interesting - go to the next particle! (i.e. wait for next event...)
            //                continue;

            //            histMult_charged->Fill(nChargedFromTree);

            //finish prev event if was opened and start new
            if ( prevEvId >= 0 ) //i.e. we have previous event to be closed
            {

                for ( int iAn = 0; iAn < nAnalysers; iAn++ )
                    analysersArray[iAn]->FinishEvent();

                //for correlations: fill hist if any charged particles in eta range
                if ( nCharged_inCuts > 0 )
                {
                    histMult_charged->Fill(nCharged_inCuts);
                    histMult_Lambda0->Fill(nLambda0_inCuts);
                    histMult_K0->Fill(nK0_inCuts);
                    for ( int iPair = 0; iPair < nPairs; ++iPair )
                    {
                        histNN[iPair]->Fill( nChargedF[iPair], nChargedB[iPair] );
                        histNsNs[iPair]->Fill( nChargedFs[iPair], nChargedBs[iPair] );
                        histNsN[iPair]->Fill( nChargedF[iPair], nChargedBs[iPair] );
                        histNskN[iPair]->Fill( nChargedF[iPair], nChargedBsk[iPair] );

                        histNlN[iPair]->Fill( nChargedF[iPair], nChargedBlambda[iPair] );
                        histNlNl[iPair]->Fill( nChargedFlambda[iPair], nChargedBlambda[iPair] );
                        histNprotonN[iPair]->Fill( nChargedF[iPair], nChargedBproton[iPair] );

                        if ( nChargedB[iPair] > 0 )
                        {
                            histPtN[iPair]->Fill( nChargedF[iPair], meanPt_B[iPair]/nChargedB[iPair] );
                            histPtNs[iPair]->Fill( nChargedFs[iPair], meanPt_B[iPair]/nChargedB[iPair] );
                            histPtNsk[iPair]->Fill( nChargedFsk[iPair], meanPt_B[iPair]/nChargedB[iPair] );
                            histPtNl[iPair]->Fill( nChargedFlambda[iPair], meanPt_B[iPair]/nChargedB[iPair] );

                            if ( nChargedF[iPair] > 0 )
                                histPtPt[iPair]->Fill( meanPt_F[iPair]/nChargedF[iPair], meanPt_B[iPair]/nChargedB[iPair] );
                        }


                        if ( nChargedBpion[iPair] > 0 )
                        {
                            histKtoPion_N[iPair]->Fill( nChargedF[iPair], nChargedBkaon[iPair]/nChargedBpion[iPair] );
                            histLambdaToPion_N[iPair]->Fill( nChargedF[iPair], nChargedBlambda[iPair]/nChargedBpion[iPair] );
                            histProtonToPion_N[iPair]->Fill( nChargedF[iPair], nChargedBproton[iPair]/nChargedBpion[iPair] );
                            if ( nChargedFpion[iPair] > 0 )
                            {
                                histKtoPion_KtoPion[iPair]->Fill( nChargedFkaon[iPair]/nChargedFpion[iPair], nChargedBkaon[iPair]/nChargedBpion[iPair] );
                                histLambdaToPion_LambdaToPion[iPair]->Fill( nChargedFlambda[iPair]/nChargedFpion[iPair], nChargedBlambda[iPair]/nChargedBpion[iPair] );
                            }
                        }

                        if ( nChargedBpion[iPair] > 0 && nChargedF[iPair] > 0)
                        {
                            histKtoPion_Pt[iPair]->Fill( meanPt_F[iPair]/nChargedF[iPair], nChargedBkaon[iPair]/nChargedBpion[iPair] );
                            histLambdaToPion_Pt[iPair]->Fill( meanPt_F[iPair]/nChargedF[iPair], nChargedBlambda[iPair]/nChargedBpion[iPair] );
                            histProtonToPion_Pt[iPair]->Fill( meanPt_F[iPair]/nChargedF[iPair], nChargedBproton[iPair]/nChargedBpion[iPair] );
                        }

                    }
                } // if nCh in cuts > 0

            }

            prevEvId = part.eventId;
            nAnalysedEvents++;

            // start first event in analysers...
            for ( int iAn = 0; iAn < nAnalysers; iAn++ )
                analysersArray[iAn]->StartEvent();

            // ##### buffer particles into array
            int iParticleInArray = 0;
            while ( prevEvId == part.eventId && iEntry+iParticleInArray < nentries)
            {
                if ( iParticleInArray >= 5000 )
                {
                    cout << "!!! iParticleInArray exceeds array size !!!" << endl;
                    break;
                }
                //                cout << iParticleInArray << " " << part.eventId << endl;
                treeParticles->GetEntry(iEntry+iParticleInArray); //
                partArr[iParticleInArray] = part;
                iParticleInArray++;
            }
            treeParticles->GetEntry(iEntry); // return to current particle



            //print event counter
            if ( part.eventId % 100 == 0 )
                printf("Processing %d event...\n", part.eventId );

            //refresh
            nCharged_inCuts = 0;
            nLambda0_inCuts = 0;
            nK0_inCuts = 0;

            for ( int iPair = 0; iPair < nPairs; ++iPair )
            {
                nChargedF[iPair] = 0;
                meanPt_F[iPair] = 0;

                nChargedFs[iPair] = 0;
                nChargedFsk[iPair] = 0;
                nChargedFlambda[iPair] = 0;
                nChargedFproton[iPair] = 0;
                nChargedFkaon[iPair] = 0;
                nChargedFpion[iPair] = 0;

                nChargedB[iPair] = 0;
                meanPt_B[iPair] = 0;

                nChargedBs[iPair] = 0;
                nChargedBsk[iPair] = 0;
                nChargedBlambda[iPair] = 0;
                nChargedBproton[iPair] = 0;
                nChargedBkaon[iPair] = 0;
                nChargedBpion[iPair] = 0;
            }
        }

        if (part.eventId < 2000)
            mapTau0_byName.insert ( pair<string,double>( mapPID[part.id], part.tau0 ) );



//        if ( fabs(part.eta) < 2.4 ) // eta cut
        if ( fabs(part.y) < 2.4 ) // eta cut
            //        if ( fabs(part.y) < 2 ) // eta or y cut
        {
            //do we need to consider this particle?
            bool isParticleOfInterest = false;
            int pofIdInArray = -1;
            for ( int i = 0; i < nSpecies; i++ )
            {
                if ( part.id == arrParticleSpecies[i] )
                {
                    isParticleOfInterest = true;
                    pofIdInArray = i;
                    break;
                }
            }

            if (isParticleOfInterest)
            {
                //                if ( part.id == 310 )               getDecayChainById(part);
                bool hasPrimaryMother = hasPrimaryMotherByALICEdefinition(part);
                //                                cout << mapPID[part.id] << " " << hasPrimaryMother << endl;

                if (hasPrimaryMother)
                {
                    histYields->Fill(pofIdInArray);
                    histPt_species[pofIdInArray]->Fill(part.pt);
                    histEta_species[pofIdInArray]->Fill(part.eta);
                    histY_species[pofIdInArray]->Fill(part.y);

                    // fill eta and y for special choice of FB windows
                    if ( part.y > 1 && part.y < 2 )
                    {
                        histEta_species_in_windows[pofIdInArray]->Fill(part.eta);
                        histY_species_in_windows[pofIdInArray]->Fill(part.y);
                    }
                    if ( part.y > -2 && part.y < -1 )
                    {
                        histEta_species_in_windows[pofIdInArray]->Fill(part.eta);
                        histY_species_in_windows[pofIdInArray]->Fill(part.y);
                    }


                    bool isFinalCharged = ( part.isFinal && part.isCharged );
                    bool isStrange = ( fabs(part.id) == 3122 || part.id == 310 );
                    bool isStrangeWithKaonsPM = ( fabs(part.id) == 3122 || part.id == 310 || fabs(part.id) == 321 );
                    bool isLambda = ( fabs(part.id) == 3122 );
                    bool isKaon = ( part.id == 310 || fabs(part.id) == 321 );
                    bool isPion = ( fabs(part.id) == 211 );
                    bool isProton = ( fabs(part.id) == 2212 );

                    if ( isStrangeWithKaonsPM )  //fabs(part.id) == 3122 || part.id == 310 ) //Lambda or K0S
                    {
                        histEtaStrange->Fill(part.eta);
                        histYStrange->Fill(part.y);
                    }

                    //for correlations
                    for ( int iPair = 0; iPair < nPairs; ++iPair )
                    {


                        if ( part.y > etaMinF[iPair] && part.y < etaMaxF[iPair] )
                        {
                            if ( isFinalCharged)            nChargedF[iPair]++;
                            if ( isFinalCharged )           meanPt_F[iPair] += part.pt;

                            if ( isStrange )                nChargedFs[iPair]++;
                            if ( isStrangeWithKaonsPM )     nChargedFsk[iPair]++;
                            if ( isLambda )                 nChargedFlambda[iPair]++;
                            if ( isKaon )                   nChargedFkaon[iPair]++;
                            if ( isPion )                   nChargedFpion[iPair]++;
                            if ( isProton )                 nChargedFproton[iPair]++;
                        }
                        if ( part.y > etaMinB[iPair] && part.y < etaMaxB[iPair] )
                        {
                            if ( isFinalCharged )           nChargedB[iPair]++;
                            if ( isFinalCharged )           meanPt_B[iPair] += part.pt;

                            if ( isStrange )                nChargedBs[iPair]++;
                            if ( isStrangeWithKaonsPM )     nChargedBsk[iPair]++;
                            if ( isLambda )                 nChargedBlambda[iPair]++;
                            if ( isKaon )                   nChargedBkaon[iPair]++;
                            if ( isPion )                   nChargedBpion[iPair]++;
                            if ( isProton )                 nChargedBproton[iPair]++;
                        }

                    }

                } // end of if (hasPrimaryMother)
            } // end of if (isParticleOfInterest)


            if (fabs(part.id) == 3122)
            {
                ++nLambda0_inCuts;

                int d1 = part.daughter1;
                int d2 = part.daughter2;
                histLambda0_daugthers_dEta->Fill( fabs( partArr[d1].eta - partArr[d2].eta ) );

                //                getDecayChainById(part);
                //                bool isPrimary = hasPrimaryMotherByALICEdefinition(part);
                //                cout << mapPID[part.id] << " " << isPrimary << endl;
            }

            if (part.id == 311)
                ++nK0_inCuts;


            if ( part.isFinal && part.isCharged )
            {





                //                getDecayChainById(part);
                //                cout << "particle: " << part.id << ", name: " << mapPID[part.id]
                //                     << "  ---  mother: " << partArr[part.mother1].id << ", name: " << mapPID[partArr[part.mother1].id]
                //                        //<< endl;
                //                //cout
                //                        << "  ---  grandmother: " << partArr[partArr[part.mother1].mother1].id << ", name: " << pidMap[partArr[partArr[part.mother1].mother1].id] << endl;

                //                cout << "decay chain: " <<getDecayChain(part) << endl;

                ++nCharged_inCuts;

                simpleTrack->id       = part.particleId;
                simpleTrack->eta      = part.eta;
                simpleTrack->phi      = part.phi;//    fixPhi(track->phi);
                simpleTrack->pt       = part.pt;
                simpleTrack->charge   = part.charge;
                simpleTrack->pid      = part.id;

                for ( int iAn = 0; iAn < nAnalysers; iAn++ )
                    analysersArray[iAn]->AddTrack( simpleTrack );

            }
        }
        histEta->Fill( part.eta );
        histPhi->Fill( part.phi );
        histPt->Fill( part.pt );

    } // end of particle loop

    //finish last event
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
        analysersArray[iAn]->FinishEvent();

    //last actions in analysers
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
        analysersArray[iAn]->Terminate();

    printf("nAnalysedEvents=%d\n", nAnalysedEvents );



    //print particle lifetimes
    for (auto it = mapTau0_byName.begin(); it != mapTau0_byName.end(); ++it)
    {
        if ((*it).second > 0)
            cout << (*it).first << " : " << (*it).second << endl;
    }



    TFile *fileOutput = new TFile( Form("analysisResults_%s_newWins2_Y_withPtPt.root", strTuneName.Data()),"RECREATE");
    //    TFile *fileOutput = new TFile( "analysisResults.root","RECREATE");
    histMult_charged->Write();
    histMult_Lambda0->Write();
    histMult_K0->Write();
    histLambda0_daugthers_dEta->Write();
    histEta->Write();
    histEtaStrange->Write();
    histYStrange->Write();

    TCanvas *canvChecks = new TCanvas("canvChecks","check plots",10,10,1200,600 );
    canvChecks->Divide(4,2);
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

    canvChecks->cd(4);
    histLambda0_daugthers_dEta->DrawCopy();
    //    histYields->DrawCopy();

    canvChecks->cd(5);
    histEta->DrawCopy();

    canvChecks->cd(6);
    histEtaStrange->DrawCopy();

    canvChecks->cd(7);
    histYStrange->DrawCopy();




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


    if (1)
    {
        TCanvas *canvYields = new TCanvas("canvYields","yields",300,200,1000,600 );
        canvYields->SetTopMargin(0.05);
        canvYields->SetRightMargin(0.15);
        canvYields->SetBottomMargin(0.1);
        canvYields->SetLeftMargin(0.15);
        histYields->SetMarkerColor(kRed);
        histYields->SetMarkerStyle(21);
        histYields->DrawCopy("P");

        histYields->Write();
    }

    if (1) // write pt, eta, y for all species
    {
        for(int i = 0; i < nSpecies; i++)
        {
            histPt_species[i]->Write();
            histEta_species[i]->Write();
            histY_species[i]->Write();
            histEta_species_in_windows[i]->Write();
            histY_species_in_windows[i]->Write();
        }
    }
    histWindowsEtaBounds->Write();

    if (1) //pT for Lambda, K0S
    {
        TCanvas *canvPt = new TCanvas("canvPt","yields",400,250,1000,600 );
        canvPt->SetTopMargin(0.05);
        canvPt->SetRightMargin(0.15);
        canvPt->SetBottomMargin(0.1);
        canvPt->SetLeftMargin(0.15);

        // K0_S
        TH1D *histPt_K0 = histPt_species[ idInArray(arrParticleSpecies, nSpecies, 310) ];
        histPt_K0->SetMarkerColor(kRed);
        histPt_K0->SetMarkerStyle(25);
        histPt_K0->DrawCopy("P");


        // ##### Lambda and Antilambda ration VS pt
        TH1D *histPt_Lambda = histPt_species[ idInArray(arrParticleSpecies, nSpecies, 3122) ];
        TH1D *histPt_AntiLambda = histPt_species[ idInArray(arrParticleSpecies, nSpecies, 3122) ];

        //add L and AL
        TH1D *histPt_Lambda_AntiLambda = new TH1D( *histPt_Lambda );
        // !!! don't add if want only lambda
        histPt_Lambda_AntiLambda->Add( histPt_AntiLambda, 1. );
        histPt_Lambda_AntiLambda->SetName ( "histPt_Lambda_AntiLambda" );

        histPt_Lambda_AntiLambda->SetMarkerColor(kBlue);
        histPt_Lambda_AntiLambda->SetMarkerStyle(21);
        histPt_Lambda_AntiLambda->DrawCopy("P same");


        //canvas with Lambda/K0S ratio
        TCanvas *canvPt_ratio = new TCanvas("canvPt_ratio","yields",480,280,800,600 );
        canvPt_ratio->SetTopMargin(0.05);
        canvPt_ratio->SetRightMargin(0.15);
        canvPt_ratio->SetBottomMargin(0.1);
        canvPt_ratio->SetLeftMargin(0.15);

        TH1D *histPt_divided_Lambda_vs_K0S = new TH1D ( *histPt_Lambda_AntiLambda );
        histPt_divided_Lambda_vs_K0S->Divide ( histPt_K0 );
        histPt_divided_Lambda_vs_K0S->SetName ( "histPt_Lambda_vs_K0S" );
        histPt_divided_Lambda_vs_K0S->DrawCopy("P");

        histPt_divided_Lambda_vs_K0S->Write();
    }

    if (1) //Y_or_eta for Lambda, K0S
    {
        // !!! careful! Y may mean eta-histogram below!
        TCanvas *canvY = new TCanvas("canvY","canvY",400,250,1000,600 );
        canvY->SetTopMargin(0.05);
        canvY->SetRightMargin(0.15);
        canvY->SetBottomMargin(0.1);
        canvY->SetLeftMargin(0.15);

        // K0_S
        TH1D *histEta_K0 = histEta_species[ idInArray(arrParticleSpecies, nSpecies, 310) ];
        histEta_K0->SetMarkerColor(kRed);
        histEta_K0->SetMarkerStyle(25);
        histEta_K0->DrawCopy("P");


        // ##### Lambda and Antilambda ration VS Y
        TH1D *histEta_Lambda = histEta_species[ idInArray(arrParticleSpecies, nSpecies, 3122) ];
        TH1D *histEta_AntiLambda = histEta_species[ idInArray(arrParticleSpecies, nSpecies, 3122) ];

        //add L and AL
        TH1D *histEta_Lambda_AntiLambda = new TH1D( *histEta_Lambda );
        // !!! don't add if want only lambda
        histEta_Lambda_AntiLambda->Add( histEta_AntiLambda, 1. );
        histEta_Lambda_AntiLambda->SetName ( "histEta_Lambda_AntiLambda" );

        histEta_Lambda_AntiLambda->SetMarkerColor(kBlue);
        histEta_Lambda_AntiLambda->SetMarkerStyle(21);
        histEta_Lambda_AntiLambda->DrawCopy("P same");


        //canvas with Lambda/K0S ratio
        TCanvas *canvY_ratio = new TCanvas("canvY_ratio","yields",480,280,800,600 );
        canvY_ratio->SetTopMargin(0.05);
        canvY_ratio->SetRightMargin(0.15);
        canvY_ratio->SetBottomMargin(0.1);
        canvY_ratio->SetLeftMargin(0.15);

        TH1D *histEta_divided_Lambda_vs_K0S = new TH1D ( *histEta_Lambda_AntiLambda );
        histEta_divided_Lambda_vs_K0S->Divide ( histEta_K0 );
        histEta_divided_Lambda_vs_K0S->SetName ( "histEta_Lambda_vs_K0S" );
        histEta_divided_Lambda_vs_K0S->DrawCopy("P");

        histEta_divided_Lambda_vs_K0S->Write();
    }


    if (1)
    {
        for (int i = 0; i < nPairs; ++i)
        {
            histNN[i]->Write();
            histNsNs[i]->Write();
            histNsN[i]->Write();
            histNskN[i]->Write();

            histNlN[i]->Write();
            histNlNl[i]->Write();
            histNprotonN[i]->Write();

            histPtN[i]->Write();
            histPtNs[i]->Write();
            histPtNsk[i]->Write();
            histPtNl[i]->Write();

            histKtoPion_N[i]->Write();
            histLambdaToPion_N[i]->Write();
            histProtonToPion_N[i]->Write();
            histKtoPion_KtoPion[i]->Write();
            histLambdaToPion_LambdaToPion[i]->Write();

            histKtoPion_Pt[i]->Write();
            histLambdaToPion_Pt[i]->Write();
            histProtonToPion_Pt[i]->Write();

            histPtPt[i]->Write();

            // PROFILES
            histNN[i]->ProfileX()->Write();
            histNsNs[i]->ProfileX()->Write();
            histNsN[i]->ProfileX()->Write();
            histNskN[i]->ProfileX()->Write();

            histNlN[i]->ProfileX()->Write();
            histNlNl[i]->ProfileX()->Write();
            histNprotonN[i]->ProfileX()->Write();

            histPtN[i]->ProfileX()->Write();
            histPtNs[i]->ProfileX()->Write();
            histPtNsk[i]->ProfileX()->Write();
            histPtNl[i]->ProfileX()->Write();

            histKtoPion_N[i]->ProfileX()->Write();
            histLambdaToPion_N[i]->ProfileX()->Write();
            histProtonToPion_N[i]->ProfileX()->Write();
            histKtoPion_KtoPion[i]->ProfileX()->Write();
            histLambdaToPion_LambdaToPion[i]->ProfileX()->Write();

            histKtoPion_Pt[i]->ProfileX()->Write();
            histLambdaToPion_Pt[i]->ProfileX()->Write();
            histProtonToPion_Pt[i]->ProfileX()->Write();

            histPtPt[i]->ProfileX()->Write();
        }
    }


    //correlation plots
    if (1)
    {
        TCanvas *canvNN = new TCanvas("canvNN","canvNN",500,250,1200,800 );
        canvNN->Divide(4,3);
        canvNN->SetTopMargin(0.05);
        canvNN->SetRightMargin(0.15);
        canvNN->SetBottomMargin(0.1);
        canvNN->SetLeftMargin(0.15);


        TCanvas *canvNN_profiles = new TCanvas("canvNN_profiles","canvNN_profiles",550,260,1200,800 );
        canvNN_profiles->Divide(4,3);
        canvNN_profiles->SetTopMargin(0.05);
        canvNN_profiles->SetRightMargin(0.15);
        canvNN_profiles->SetBottomMargin(0.1);
        canvNN_profiles->SetLeftMargin(0.15);

        TGraphErrors *graph = new TGraphErrors();
        for (int i = 0; i < nPairs; ++i)
        {
            //2D hists
            canvNN->cd(i+1);
            histPtN[i]->DrawCopy("colz");
            //profiles
            canvNN_profiles->cd(i+1);
            //            histPtN[i]->ProfileX()->DrawCopy();

            TProfile *myProfile = histPtN[i]->ProfileX();
            myProfile->Fit("pol1");
            TF1 *myfit = (TF1*) myProfile->GetFunction("pol1");
            double a = myfit->GetParameter(0);
            double b = myfit->GetParameter(1);
            double bErr = myfit->GetParError(1);
            cout << "corr coeff = " << b << endl;
            graph->SetPoint( i, 2*etaMinF[i], b );
            graph->SetPointError( i, 0, bErr );

        }
        //        histPtN[0]->SetMarkerColor(kRed);
        //        histPtN[0]->SetMarkerStyle(25);

        canvNN_profiles->cd(nPairs+1);
        graph->Draw ( "APL" ) ;

        //graph->SetMarkerStyle(22);
        graph->SetMarkerStyle(kOpenCircle);
        graph->SetMarkerSize(1.5);

        graph->SetTitle( ";#Delta#eta;b" );
        graph->SetMinimum( 0 );

        graph->GetXaxis()->SetTitleSize(0.058);
        graph->GetYaxis()->SetTitleSize(0.058);

        graph->GetXaxis()->SetTitleOffset( 0.75 );
        graph->GetYaxis()->SetTitleOffset( 1. );


        TLatex *tex = new TLatex(0.3,0.84, "#sqrt{s}=7 TeV, p_{T}#in(0.2, 2.0)GeV/c");
        tex->SetNDC(1);
        tex->SetTextFont(42);
        tex->SetTextSize(0.046);
        tex->DrawClone();

    }


    if (0)
    {
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

    }

    fileOutput->Close();

    //    file->Close();

}

