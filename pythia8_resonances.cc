//
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TSystem.h>

#include <TObject.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TString.h>

#include <Pythia8/Pythia.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//#include "../../../commonTools/Tools.cxx"
//#include "../../../commonTools/SimpleTrack.cxx"
//#include "../../../simpleEventAnalyzer/AliSimpleEvent.h"
//#include "../../../simpleEventAnalyzer/AliSimpleEvent.cxx"


using namespace std;

void fixPhi( float &phi )
{
    if ( phi < 0 )
        phi += 2*TMath::Pi();
    if ( phi > 2*TMath::Pi() )
        phi -= 2*TMath::Pi();
}


const double etaRange = 2.0;
const double ptMin = 0.1;

bool checkCuts( Pythia8::Particle &particle )
{
    return ( fabs( particle.eta() ) < etaRange && particle.pT() > ptMin );
}
bool checkEtaCut( Pythia8::Particle &particle )
{
    return ( fabs( particle.eta() ) < etaRange );
}


// prepare variables for the tree:
float myPtToTree, myEtaToTree, myPhiToTree;
int myPid, myCharge;

// the rest of vars is being changed inside main loop
int myEventId = 0;
int myCategory = 0;
int myParticleId = 0;

void fillTree( TTree *tree, Pythia8::Particle &particle )
{
    myPid = particle.id();
    myCharge = particle.charge();
    myPtToTree = particle.pT();
    myEtaToTree = particle.eta();
    myPhiToTree = particle.phi();
    fixPhi( myPhiToTree );

    tree->Fill();
    myParticleId++;
}


//void prepareCanvas(TCanvas *c)
//{
//    c->SetLeftMargin(0.15);
//    c->SetRightMargin(0.05);
//    c->SetTopMargin(0.05);
//    c->SetBottomMargin(0.11);
//}


const int nResonances = 3;
TH1D *fHistResonanceMass[nResonances];
TH1D *fHistResonanceMassRho0ToCrosscheck;
TH1D *fHistResonancePt[nResonances];
TH1D *fHistResonancePtWithoutWeight[nResonances];
TH1D *fHistResonanceEta[nResonances];
TH1D *fHistResonancePhi[nResonances];

TH2D *fHist2DResonancePtWithoutWeight[nResonances];


TH1D *fHistRho0EtaIfTwoPionsDetected;


TH1D *fHistChargedPionsPt;
TH1D *fHistChargedPionsPrimaryPt;
TH1D *fHistChargedPionsFromResonancesPt;
TH1D *fHistChargedPionsFromRho0Pt;
TH1D *fHistChargedPionsFromRhoPlusMinusPt;


TH1D *fHistChargedPionsPrimaryEta;
TH1D *fHistChargedPionsFromResonancesEta;

TH1D *fHistNCharged;
TH1D *fHistNChargedPionsPrimary;
TH1D *fHistNChargedPionsFromResonances;
TH1D *fHistNChargedPionsFromRho0;
TH1D *fHistNChargedPionsFromRhoPlusMinus;
TH1D *fHistNChargedPionsFromRho0BothDetected;


/*
 * Sort the particle ID array and particle name array in such a way that the largest contributions are in the beginning
 *
 *  int particleIDArray[2][100] = particle ID array that needs to sorted
 *  TString particleNames[100] = particle name array that needs to be sorted
 */
void sortArrays(int particleIDArray[2][100], TString particleNames[100])
{
    bool arrayNotSorted = true;
    int temporaryInt = 0;
    TString tempararyString = "";
    while(arrayNotSorted)
    {
        arrayNotSorted = false;
        for( int i = 0; i < 99; i++ )
        {
            // If there is later number is larger, swap the two numbers and the same indices in other array
            if(particleIDArray[1][i] < particleIDArray[1][i+1])
            {
                temporaryInt = particleIDArray[1][i];
                particleIDArray[1][i] = particleIDArray[1][i+1];
                particleIDArray[1][i+1] = temporaryInt;

                temporaryInt = particleIDArray[0][i];
                particleIDArray[0][i] = particleIDArray[0][i+1];
                particleIDArray[0][i+1] = temporaryInt;

                tempararyString = particleNames[i];
                particleNames[i] = particleNames[i+1];
                particleNames[i+1] = tempararyString;

                arrayNotSorted = true;
            }
        }
    }
}

/*
 * Function to extract the particle ID and name from a PYTHIA particle and filling them into arrays.
 * If this is the first appearence of a found particle type, the particle ID and name are given a new spot in the arrays
 * If there is already a similar particle found, then the number of particles found in particleIDArray is incremented by one
 *
 *  Pythia8::Particle &currentParticle = PYTHIA8 particle containing the desired information
 *  int particleIDArray[2][100] = array for particle IDs
 *  TString particleNames[100] = array for particle names
 */
void addParticleToArray(Pythia8::Particle &currentParticle, int particleIDArray[2][100], TString particleNames[100])
{
    // First, find out if the particle is already in the array or not
    int particleID = currentParticle.id();
    int indexInArray = -1;
    bool firstAppearance = false;
    for(int i = 0; i < 100; i++)
    {
        // If we find the searched particle, remember its index
        if(particleIDArray[0][i] == particleID)
        {
            indexInArray = i;
            break;
        }

        // If the number in the first row is zero, we have gone through all the particles without finding a match
        if(particleIDArray[0][i] == 0)
        {
            firstAppearance = true;
            indexInArray = i;
            break;
        }
    }

    // If first appearence, put the particle name and ID to arrays
    if(firstAppearance)
    {
        particleIDArray[0][indexInArray] = particleID;
        particleNames[indexInArray] = Form("%s",currentParticle.name().c_str());
    }

    // Finally, increment the number of found particles of this type
    particleIDArray[1][indexInArray]++;
}

/*
 * Check if certain particle is in the decay chain of the particle considered
 *
 *  TString decayChain = decay chain of a particle
 *  int particleID = ID of the particle that is searched from the decay chain
 */
bool hasParticleInDecayChain(TString decayChain, int particleID)
{

    TObjArray *particlesInDecayChain = decayChain.Tokenize("|");
    const int numberOfParticles = particlesInDecayChain->GetEntries();
    TObjString *currentParticleObject;
    TString currentParticle;
    int currentParticleID = 0;
    for( int i = 0; i < numberOfParticles; i++)
    {
        currentParticleObject = (TObjString *)particlesInDecayChain->At(i);
        currentParticle = currentParticleObject->String();
        currentParticleID = currentParticle.Atoi();
        if(currentParticleID == particleID) return true;
    }

    return false;
}

/*
 * Find out the decay chain of the final state particle after the hadronization
 * This function should be only used for final state particles, since it searches for the decay chain recursively
 * If used for non-final particles, it will find a part of a decay chain for some final state particle
 *
 *  Pythia8::Event &event = PYTHIA8 event information for the current event
 *  Pythia8::Particle &currentParticle = PYTHIA8 particle information for the considered particle
 *  TString decayChain = The decay chain of the particle until this point
 */
TString getDecayChain(Pythia8::Event &event, Pythia8::Particle &currentParticle, TString decayChain = "0")
{
    // The mother of the current particle is born before hadronization. Stop the decay chain here.
    if(event[currentParticle.mother1()].status() > -80)
    {
        return decayChain;
    }

    // The particles that are formed after hadronization can have only one mother. If we are not looking
    // at the particles before the hadronization, mother1 is the only mother we need to care about
    TString newAddition = Form( "|%d",event[currentParticle.mother1()].id() );
    decayChain += newAddition;
    return getDecayChain(event, event[currentParticle.mother1()], decayChain);
}

/*
 * Find out if a decaying particle ends its decay chain to a charged hadron
 */
bool hasChargedDaughter(Pythia8::Event &event, Pythia8::Particle &currentParticle, bool firstIteration)
{
    if (event[currentParticle.mother1()].id() == currentParticle.id() && firstIteration)
        return false;  // The current particle is a recoiled version of its mother. Do not count it twice.

    int daughter1 = currentParticle.daughter1();
    int daughter2 = currentParticle.daughter2();


    // Final state particle. Check if charged hadron.
    if ( daughter1 == 0 && daughter2 == 0)
        //    if( /*currentParticle.isCharged() && currentParticle.isHadron() && */currentParticle.isFinal() )
    {
        if( currentParticle.isCharged() && currentParticle.isHadron() ) //&& currentParticle.isFinal() )
            //            if( currentParticle.isCharged() )//&& currentParticle.isHadron() && currentParticle.isFinal() )
        {
            cout << "final state particle: " << currentParticle.name() << endl;
            return true;
        }
        return false;
    }

    //    const int pdgRho[3] = { 113, 213, -213 };
    if ( currentParticle.id() == 113 || fabs(currentParticle.id()) == 213 )
    {
        if(0)cout << "not final particle: " /*<< currentParticle.id()*/
                  <<  " " << currentParticle.name() << " " << currentParticle.m() << ", first iteration=" << firstIteration
                   <<  " daughter1: " << event[daughter1].name() << " " << event[daughter1].m()
                    <<  ", daughter2: " << event[daughter2].name() << " " << event[daughter2].m() << endl;
        //        fHistRhoMass->Fill( currentParticle.m() );

        int resonanceId = -1000;
        if ( currentParticle.id() == 113 )
            resonanceId = 0;
        else if ( currentParticle.id() == 213 )
            resonanceId = 1;
        else if ( currentParticle.id() == -213 )
            resonanceId = 2;

        if ( resonanceId >= 0 )
        {
            fHistResonanceMass[resonanceId]->Fill( currentParticle.m() );
            fHistResonancePt[resonanceId]->Fill( currentParticle.pT(), 1./currentParticle.pT() );
            if ( checkEtaCut(event[daughter1]) && checkEtaCut(event[daughter2]) )
            {
                fHistResonancePtWithoutWeight[resonanceId]->Fill( currentParticle.pT() );
                fHist2DResonancePtWithoutWeight[resonanceId]->Fill( currentParticle.eta(), currentParticle.pT() );
            }
            fHistResonanceEta[resonanceId]->Fill( currentParticle.eta() );
            fHistResonancePhi[resonanceId]->Fill( currentParticle.phi() );
        }

        //fill eta hist for rho0 if both charged pions are detected
        if ( resonanceId == 0 ) //rho0
        {
            // if decay into two charged pions
            if ( currentParticle.daughterList().size() == 2 && fabs( event[daughter1].id() )  == 211 && fabs( event[daughter2].id() )  == 211 )
            {
                // if both pions are "detected"
                if ( checkCuts(event[daughter1]) && checkCuts(event[daughter2]) )
                    fHistRho0EtaIfTwoPionsDetected->Fill( currentParticle.eta() );

            }
        }


    }

    // Only one decay product. Check if it is final state charged hadron
    if ( (daughter1 == daughter2) || ((daughter1 > daughter2) && (daughter2 == 0)) )
    {
        if( event[daughter1].isFinal() && event[daughter1].isCharged() && event[daughter1].isHadron() )
            return true;
        return hasChargedDaughter(event,event[daughter1],false);
    }

    // Array of decay products from dauhgter1 to daughter2
    if ( (daughter1 < daughter2) && (daughter1 > 0) )
    {
        for(int iDaughter = daughter1; iDaughter <= daughter2; iDaughter++){
            if(event[iDaughter].isFinal() && event[iDaughter].isCharged() && event[iDaughter].isHadron())
                return true;
        }

        bool chargedDaughterExists = false;
        for(int iDaughter = daughter1; iDaughter <= daughter2; iDaughter++)
        {
            chargedDaughterExists == chargedDaughterExists ||
                    (event,event[iDaughter],false);
        }
        return chargedDaughterExists;
    }

    // We are only left with the case daughter 2 < daughter 1. In this case there is two decay particles, daughter 2 and daughter 1
    if( (event[daughter1].isFinal() && event[daughter1].isCharged() && event[daughter1].isHadron()) ||
            (event[daughter2].isFinal() && event[daughter2].isCharged() && event[daughter2].isHadron()) )
    {
        return true;
    }

    return (hasChargedDaughter(event,event[daughter1],false) || hasChargedDaughter(event,event[daughter2],false));
}




/*
 * Main event loop
 *
 * Int_t nEvent = The number of events in PYTHIA run
 * const char *config = cmnd file from which the PYTHIA configuration is read
 * Int_t random_seed = Seed for PYTHIA random number generator
 * const char *outfile = Name for the output file
 */
void pythia8_spectra(Int_t nEvent  = 100, const char *config="configPythia.cmnd", Int_t random_seed=0, const char *outfile="test.root")
{
    // File
    gROOT->cd();
    TFile *fout = new TFile(outfile,"recreate","pythia 8 generation");

    //create the file, the Tree and a few branches
    //    TFile file( "myTree.root", "recreate" );


    // ##### particle tree
    TTree *treeParticles = new TTree( "particleTree","a simple Tree with simple variables");
    treeParticles->Branch("eventId",&myEventId,"ev/I");
    treeParticles->Branch("particleId",&myParticleId,"particleId/I");
    treeParticles->Branch("pt",&myPtToTree,"pt/F");
    treeParticles->Branch("eta",&myEtaToTree,"eta/F");
    treeParticles->Branch("phi",&myPhiToTree,"phi/F");
    treeParticles->Branch("charge",&myCharge,"myCharge/I");
    treeParticles->Branch("pid",&myPid,"myPid/I");
    treeParticles->Branch("particleCategory",&myCategory,"myCategory/I");

    // ##### event info tree
    int nCharged = 0;
    int nPionsPrimary = 0;
    int nPionsFromResonances = 0;
    int nPionsFromRho0 = 0;
    int nPionsFromRhoPlusMinus = 0;
    int pionsFromRho0_bothDetected = 0; //incremented by 2 if both pions within cuts

    TTree *treeEvent = new TTree( "eventMultTree","a simple Tree with simple variables");
    treeEvent->Branch("nCharged",                &nCharged,              "nCharged/I");
    treeEvent->Branch("nPionsPrimary",           &nPionsPrimary,         "nPionsPrimary/I");
    treeEvent->Branch("nPionsFromResonances",    &nPionsFromResonances,  "nPionsFromResonances/I");
    treeEvent->Branch("nPionsFromRho0",          &nPionsFromRho0,        "nPionsFromRho0/I");
    treeEvent->Branch("nPionsFromRhoPlusMinus",  &nPionsFromRhoPlusMinus,"nPionsFromRhoPlusMinus/I");
    treeEvent->Branch("pionsFromRho0_bothDetected",&pionsFromRho0_bothDetected,"pionsFromRho0_bothDetected/I");



    // Histograms
    const int nCJACEK =  74 ;
    double pttJacek[nCJACEK] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
                                10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};

    TH1D *fhChargedPtJacek = new TH1D(Form("hChargedPtJacek%02d",00), Form("hChargedPtJacek%02d",00), nCJACEK-1, pttJacek );
    fhChargedPtJacek->Sumw2();

    TH1D *fhChargedPtJacekOnlyResonance = new TH1D(Form("hChargedPtJacekOnlyResonance%02d",00), Form("hChargedPtJacekOnlyResonance%02d",00), nCJACEK-1, pttJacek );
    fhChargedPtJacekOnlyResonance->Sumw2();

    TH1D *fhChargedPtJacekNoResonance = new TH1D(Form("hChargedPtJacekNoResonance%02d",00), Form("hChargedPtJacekNoResonance%02d",00), nCJACEK-1, pttJacek );
    fhChargedPtJacekNoResonance->Sumw2();

    TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",5,0,5);

    TH1D *fhParticleID = new TH1D("hPID","hPID",100,0,100);

    // ##### histograms for resonance study
    fHistResonanceMass[0] = new TH1D("histRho0Mass","histRho0Mass",1000,0,2);
    fHistResonanceMass[1] = new TH1D("fHistRhoPlusMass","fHistRhoPlusMass",1000,0,2);
    fHistResonanceMass[2] = new TH1D("fHistRhoMinusMass","fHistRhoMinusMass",1000,0,2);

    fHistResonanceMassRho0ToCrosscheck = new TH1D("fHistResonanceMassRho0ToCrosscheck","fHistResonanceMassRho0ToCrosscheck",1000,0,2);


    fHistResonancePt[0] = new TH1D("histRho0Pt","histRho0Pt",200,0,5);
    fHistResonancePt[1] = new TH1D("fHistRhoPlusPt","fHistRhoPlusPt",200,0,5);
    fHistResonancePt[2] = new TH1D("fHistRhoMinusPt","fHistRhoMinusPt",200,0,5);

    fHistResonancePtWithoutWeight[0] = new TH1D("histRho0PtWithoutWeight","histRho0PtWithoutWeight",200,0,5);
    fHistResonancePtWithoutWeight[1] = new TH1D("fHistRhoPlusPtWithoutWeight","fHistRhoPlusPtWithoutWeight",200,0,5);
    fHistResonancePtWithoutWeight[2] = new TH1D("fHistRhoMinusPtWithoutWeight","fHistRhoMinusPtWithoutWeight",200,0,5);

    fHist2DResonancePtWithoutWeight[0] = new TH2D("hist2DRho0PtWithoutWeight","hist2DRho0PtWithoutWeight;#eta;p_{T}",20,-2,2,40,0,5);
    fHist2DResonancePtWithoutWeight[1] = new TH2D("fHist2DRhoPlusPtWithoutWeight","fHist2DRhoPlusPtWithoutWeight;#eta;p_{T}",20,-2,2,40,0,5);
    fHist2DResonancePtWithoutWeight[2] = new TH2D("fHist2DRhoMinusPtWithoutWeight","fHist2DRhoMinusPtWithoutWeight;#eta;p_{T}",20,-2,2,40,0,5);

    fHistResonanceEta[0] = new TH1D("histRho0Eta","histRho0Eta", 200,-10,10 );
    fHistResonanceEta[1] = new TH1D("fHistRhoPlusEta","fHistRhoPlusEta", 200,-10,10 );
    fHistResonanceEta[2] = new TH1D("fHistRhoMinusEta","fHistRhoMinusEta", 200,-10,10 );

    fHistResonancePhi[0] = new TH1D("histRho0Phi","histRho0Phi",200,-TMath::TwoPi(),TMath::TwoPi());
    fHistResonancePhi[1] = new TH1D("fHistRhoPlusPhi","fHistRhoPlusPhi",200,-TMath::TwoPi(),TMath::TwoPi());
    fHistResonancePhi[2] = new TH1D("fHistRhoMinusPhi","fHistRhoMinusPhi",200,-TMath::TwoPi(),TMath::TwoPi());

    fHistRho0EtaIfTwoPionsDetected = new TH1D("fHistRho0EtaIfTwoPionsDetected","fHistRho0EtaIfTwoPionsDetected", 200,-10,10 );

    fHistChargedPionsPt = new TH1D("fHistChargedPionsPt","fHistChargedPionsPt",200,0,5);
    fHistChargedPionsPrimaryPt = new TH1D("fHistChargedPionsPrimaryPt","fHistChargedPionsPrimaryPt",200,0,5);
    fHistChargedPionsFromResonancesPt = new TH1D( "fHistChargedPionsFromResonancesPt","fHistChargedPionsFromResonancesPt",200,0,5);
    fHistChargedPionsFromRho0Pt = new TH1D( "fHistChargedPionsFromRho0Pt","fHistChargedPionsFromRho0Pt",200,0,5);
    fHistChargedPionsFromRhoPlusMinusPt = new TH1D( "fHistChargedPionsFromRhoPlusMinusPt","fHistChargedPionsFromRhoPlusMinusPt",200,0,5);

    fHistChargedPionsPrimaryEta = new TH1D("fHistChargedPionsPrimaryEta","fHistChargedPionsPrimaryEta", 200,-10,10 );
    fHistChargedPionsFromResonancesEta = new TH1D("fHistChargedPionsFromResonancesEta","fHistChargedPionsFromResonancesEta", 200,-10,10 );

    fHistNCharged = new TH1D( "fHistNCharged","fHistNCharged", 201,-0.5,200.5 );
    fHistNChargedPionsPrimary = new TH1D( "fHistNChargedPionsPrimary","fHistNChargedPionsPrimary", 201,-0.5,200.5 );
    fHistNChargedPionsFromResonances = new TH1D("fHistNChargedPionsFromResonances","fHistNChargedPionsFromResonances", 201,-0.5,200.5 );

    fHistNChargedPionsFromRho0 = new TH1D("fHistNChargedPionsFromRho0","fHistNChargedPionsFromRho0", 201,-0.5,200.5 );
    fHistNChargedPionsFromRhoPlusMinus = new TH1D("fHistNChargedPionsFromRhoPlusMinus","fHistNChargedPionsFromRhoPlusMinus", 201,-0.5,200.5 );
    fHistNChargedPionsFromRho0BothDetected = new TH1D("fHistNChargedPionsFromRho0BothDetected","fHistNChargedPionsFromRho0BothDetected", 201,-0.5,200.5 );



    // Array for particle information
    int particleIDArray[2][100] = {0};  // First row = particle ID, second row = number of particles with this ID found
    TString particleNames[100];         // Name for the particle which is in the same column in the particleIDArray

    // Create pythia8 object
    Pythia8::Pythia pythia8;

    // Load Configure file
    pythia8.readFile(config); // set up directly because the interface don't do it.
    pythia8.readString(Form("Random:seed=%d",random_seed));

    TStopwatch timer;
    timer.Start();

    int nPace = TMath::Max(1,nEvent/50);

    // Initialize PYTHIA8
    pythia8.init();
    pythia8.settings.listChanged();

    // Accepted eta range
    //    double fTrackEtaRange = 2.0;//0.8;

    /*-----------------------------------------------------------------------------
     *  Event Loop
     *-----------------------------------------------------------------------------*/

    TClonesArray trackList("TLorentzVector");

    //counters for partices to be within kinematic cuts
    int nBothPionsFromRho0WithinCuts = 0;
    int nOnlyOnePionFromRho0WithinCuts = 0;

    int nChargedInAllEvents = 0;
    for (Int_t iEvent = 0; iEvent < nEvent; iEvent++)
    {
        // Generate Event
        trackList.Clear();
        pythia8.next();
        Pythia8::Event &event = pythia8.event;

        // Print out some information about the first event
        if (iEvent < 1) event.list();

        int numberOfParticles = event.size();

        // Print out information about the progress of the simulation
        if( iEvent > 1000 && ( iEvent%nPace == 0 ) )
        {
            cout << endl << "Event " << iEvent << ", i.e. " << 100.*iEvent/nEvent << " done" << endl;
        }

        //make zero values
        nCharged = 0;
        nPionsPrimary = 0;
        nPionsFromResonances = 0;
        nPionsFromRho0 = 0;
        nPionsFromRhoPlusMinus = 0;
        pionsFromRho0_bothDetected = 0;

        // Particle Loop to fill hist
        myParticleId = 0;
        for (Int_t iParticle = 0; iParticle < numberOfParticles; iParticle++)
        {
            // choose the next particle
            Pythia8::Particle &currentParticle = event[iParticle];

            // Find the particles that decay after hadronization
            if (currentParticle.status() < -80)
            {
                bool thereIsChargedDaughter = hasChargedDaughter(event,currentParticle,true);
                if (thereIsChargedDaughter)
                {
//                    addParticleToArray(currentParticle,particleIDArray,particleNames);
                }
            }

            // IA: ######### find out what is the probability that both pions from rho0 decay are detected within cut conditions
            if ( currentParticle.id() == 113 ) //from rho0
            {
                //check daughterList
                if ( currentParticle.daughterList().size() != 2 )
                {
                    cout << "n rho0 daughters: " << currentParticle.daughterList().size() << endl;
                    for ( int iDau = 0; iDau < currentParticle.daughterList().size(); iDau++ )
                        cout  << "  i=" << iDau << " " << event[currentParticle.daughterList()[iDau]].id() << ", name: " << event[currentParticle.daughterList()[iDau]].name() << endl;
                }



                int daughter1 = currentParticle.daughter1();
                int daughter2 = currentParticle.daughter2();

                //                if( (event[daughter1].isFinal() && event[daughter1].isCharged() && event[daughter1].isHadron()) &&
                //                        (event[daughter2].isFinal() && event[daughter2].isCharged() && event[daughter2].isHadron()) )
                if( currentParticle.daughterList().size() == 2 // two children
                        && fabs( event[daughter1].id() )  == 211  && fabs( event[daughter2].id() )  == 211 ) //both - charged pions
                {
                    if(0) cout << "rho decay: " <<  " daughter1: " << event[daughter1].name() << " " << event[daughter1].m() <<  ", daughter2: " << event[daughter2].name() << " " << event[daughter2].m() << endl;

                    //check whether both pions are within cuts
                    //                    bool pion1_withinCuts = ( fabs( event[daughter1].eta() ) < etaRange && event[daughter1].pT() > ptMin );
                    //                    bool pion2_withinCuts = ( fabs( event[daughter2].eta() ) < etaRange && event[daughter2].pT() > ptMin );

                    bool pion1_withinCuts = checkCuts( event[daughter1] );
                    bool pion2_withinCuts = checkCuts( event[daughter2] );



                    //                    cout << "pion1_withinCuts=" << pion1_withinCuts << ", pion2_withinCuts=" << pion2_withinCuts << endl;
                    if ( pion1_withinCuts && pion2_withinCuts )
                    {
                        nBothPionsFromRho0WithinCuts++;
                        pionsFromRho0_bothDetected += 2;
                    }
                    if ( pion1_withinCuts && !pion2_withinCuts || !pion1_withinCuts && pion2_withinCuts)
                        nOnlyOnePionFromRho0WithinCuts++;

                    //add rho0 and its products to the tree if within cuts
                    // ##### add rho0 if at least one pion within cuts
                    if ( pion1_withinCuts || pion2_withinCuts )
                    {
                        myCategory = 10;
                        fillTree( treeParticles, currentParticle );
                    }

                    // ##### add rho0 if at BOTH pions within cuts
                    if ( pion1_withinCuts && pion2_withinCuts )
                    {
                        myCategory = 11;
                        fillTree( treeParticles, currentParticle );
                    }

                    // ##### add BOTH pions if at BOTH pions are within cuts
                    if ( pion1_withinCuts && pion2_withinCuts )
                    {
                        myCategory = 12;
                        fillTree( treeParticles, event[daughter1] );
                        fillTree( treeParticles, event[daughter2] );
                    }

                    // ##### add first pion
                    if ( pion1_withinCuts )
                    {
                        myCategory = 13;
                        fillTree( treeParticles, event[daughter1] );
                    }

                    // ##### add second pion
                    if ( pion2_withinCuts )
                    {
                        myCategory = 14;
                        fillTree( treeParticles, event[daughter2] );
                    }



                }
                //                else
                //                    cout << "AHTUNG!!!" << endl;
            }



            // IA: ######### take rho+ and rho- to the tree, if charged pion from decay in cut condition
            if ( fabs(currentParticle.id()) == 213 ) //from rho+ or rho-
            {
                int daughter1 = currentParticle.daughter1();
                int daughter2 = currentParticle.daughter2();

                if( currentParticle.daughterList().size() == 2 // two children
                        && event[daughter1].id() == 111 || event[daughter2].id() == 111 ) // if have neutral pions after decay
                {
                    int daughterChargedPion = ( event[daughter1].id() != 111 ) ? daughter1 : daughter2;
                    if(0) cout << "rho decay: " <<  " daughter1: " << event[daughterChargedPion].name() << " " << event[daughterChargedPion].m() << endl;

                    //check whether both pions are within cuts
                    bool chargedPion_withinCuts = checkCuts( event[daughterChargedPion] );

                    //add rho+/- and its product to the tree if charged pion from decay in cut condition
                    if ( chargedPion_withinCuts )
                    {
                        // ##### add rho+/-
                        myCategory = 20;
                        fillTree( treeParticles, currentParticle );

                        // ##### add charged pion (+/-)
                        myCategory = 21;
                        fillTree( treeParticles, event[daughterChargedPion] );
                    }


                }
                //                else
                //                    cout << "AHTUNG!!!" << endl;
            }




            // Select the charged final state hadrons in ALICE acceptance
            if( !currentParticle.isFinal() || !currentParticle.isCharged() || !currentParticle.isHadron() )
                continue;


            // IA: actions for final state stable particles
            if( currentParticle.isCharged() && currentParticle.isFinal() )
            {
                //                if ( fabs( currentParticle.eta() ) < etaRange && currentParticle.pT() > ptMin )
                if ( checkCuts( currentParticle ) )
                {
                    nCharged++;
                    myCategory = 0; //charged final particle
                    //                    cout << "pid = " << currentParticle.id() << ", charge = " << currentParticle.charge() << ", chargeType = " << currentParticle.chargeType() << endl;
                    if (  currentParticle.isHadron() && fabs( currentParticle.id() ) == 211 ) //consider only charged pions
                    {
                        fHistChargedPionsPt->Fill( currentParticle.pT(), 1./currentParticle.pT() );

                        int isFromQuark = ( event[currentParticle.mother1()].status() > -80 ) ? 1 : 0;
                        if(0)cout << "final state particle: " << currentParticle.name() << " " << currentParticle.id()
                                  << ", isPrimary = " << isFromQuark
                                  << ", motherName = " << event[currentParticle.mother1()].name() << endl;


                        // now can distinguish between "primaries" (coming from u,d,s quarks) and resonance products
                        if ( isFromQuark )
                        {
                            nPionsPrimary++;
                            myCategory = 1; //charged final pion from quark ("primary")
                            fHistChargedPionsPrimaryPt->Fill( currentParticle.pT(), 1./currentParticle.pT() );
                            fHistChargedPionsPrimaryEta->Fill( currentParticle.eta() );
                        }
                        else
                        {
                            nPionsFromResonances++;
                            myCategory = 2; //charged final pion from resonance decays
                            fHistChargedPionsFromResonancesPt->Fill( currentParticle.pT(), 1./currentParticle.pT() );
                            fHistChargedPionsFromResonancesEta->Fill( currentParticle.eta() );

                            if ( event[currentParticle.mother1()].id() == 113 ) //pion is from rho0
                            {
                                nPionsFromRho0++;
                                myCategory = 3;
                                fHistChargedPionsFromRho0Pt->Fill( currentParticle.pT(), 1./currentParticle.pT() );
                            }
                            else if ( fabs( event[currentParticle.mother1()].id() ) == 213 ) //pion is from rho+, rho-
                            {
                                nPionsFromRhoPlusMinus++;
                                myCategory = 4;
                                fHistChargedPionsFromRhoPlusMinusPt->Fill( currentParticle.pT(), 1./currentParticle.pT() );
                            }
                            else if ( fabs( event[currentParticle.mother1()].id() ) == 223 ) //pion is from omega
                            {
                                myCategory = 5;
                            }
                        }

                        //increment colomn for the pid-mother of this pion (cound be resonance or quark!)
                        addParticleToArray( event[currentParticle.mother1()], particleIDArray, particleNames );

                        if ( event[currentParticle.mother1()].id() == 113 && currentParticle.id() == 211 ) //take rho0 once (by taking mother of pion+)
                        {
                            fHistResonanceMassRho0ToCrosscheck->Fill( event[currentParticle.mother1()].m() );
                            //                        cout << "rho0 mass = " << event[currentParticle.mother1()].m() << endl;
                            //                        int a;
                            //                        cin >> a;
                        }
                        //check if REALLY no daughters (just a crosscheck)
                        int daughter1 = currentParticle.daughter1();
                        int daughter2 = currentParticle.daughter2();
                        if ( daughter1 != 0 || daughter2 != 0)
                            cout << "!!! have daughters?!." << endl;
                    }
                    else  //print if not charged pions
                    {
                        if(0)cout << "final state particle: " << currentParticle.name() << " " << currentParticle.id()
                                  << ", motherName = " << event[currentParticle.mother1()].name() << endl;

                    }

                    // ##### add particle to the tree (with the category determined above!)
                    fillTree( treeParticles, currentParticle );
                }
            }

            if ( !(TMath::Abs(currentParticle.eta()) < etaRange))
                continue;

            // If the status of the final state particles is less than 90, they have not decayed after hadronization
            // Thus requiring the status to be at least 90 gives only particles that are decay product from decays
            // that have happened after the hadronization
            if (currentParticle.status() > 90)
            {
                fhChargedPtJacekOnlyResonance->Fill(currentParticle.pT(),1./currentParticle.pT());
            }

            // If the status of the final state particles is less than 90, they have not decayed after hadronization
            // Thus they do not have external correlation coming from decays
            if(currentParticle.status() < 90)
            {
                fhChargedPtJacekNoResonance->Fill(currentParticle.pT(),1./currentParticle.pT());
            }

            fhChargedPtJacek->Fill(currentParticle.pT(),1./currentParticle.pT());

        }//particle loop

//        if ( nCharged > 0 )
        {
            treeEvent->Fill();

            fHistNCharged->Fill( nCharged );
            fHistNChargedPionsPrimary->Fill( nPionsPrimary );
            fHistNChargedPionsFromResonances->Fill( nPionsFromResonances );

            fHistNChargedPionsFromRho0->Fill( nPionsFromRho0 );
            fHistNChargedPionsFromRhoPlusMinus->Fill( nPionsFromRhoPlusMinus );

            fHistNChargedPionsFromRho0BothDetected->Fill( pionsFromRho0_bothDetected );

            nChargedInAllEvents += nCharged;
            myEventId++; //important to do it here, inside "if ( nCharged > 0 )" !!!
        }

        //        fSimpleEvent->GetHeader()->SetCentrality(nParticlesWithinEtaPtCuts/2/cutEtaCMS);
        //        fSimpleEvent->FinishEventFilling();

        //        fEventTree->Fill();  //fill the tree

    }//event loop


    // print ratios
    if ( nOnlyOnePionFromRho0WithinCuts+nBothPionsFromRho0WithinCuts > 0 )
    {
        cout << "RATIO: probability for rho0 that both charged pions are captured within acceptance and pt cuts: "
             << (double)nBothPionsFromRho0WithinCuts/(nOnlyOnePionFromRho0WithinCuts+nBothPionsFromRho0WithinCuts) << endl;
        cout << "RATIO: fraction of pions (from rho0), which are captured with their pion-partner together, within acceptance and pt cuts: "
             << (double)2*nBothPionsFromRho0WithinCuts/(nOnlyOnePionFromRho0WithinCuts+2*nBothPionsFromRho0WithinCuts) << endl;
        cout << "RATIO: number of pions (from rho0), which are captured with their pion-partner together, to ALL CHARGED particles, within acceptance and pt cuts: "
             << (double)2*nBothPionsFromRho0WithinCuts/nChargedInAllEvents << endl;
    }


    // Sort the array so that the particles with the most entries are the first
    sortArrays(particleIDArray,particleNames);

    for(int i = 0; i < 100; i++)
    {
        fhParticleID->Fill(i+0.5,particleIDArray[1][i]);
        fhParticleID->GetXaxis()->SetBinLabel(i+1,particleNames[i].Data());
    }


    // Display some statistics
    pythia8.stat();
    long nTried = pythia8.info.nTried();
    long nAccepted = pythia8.info.nAccepted();
    double sigmaGen = pythia8.info.sigmaGen();
    //double sigmaErr = pythia8.info.sigmaErr();
    hCrossSectionInfo->Fill(0.5,nTried);
    hCrossSectionInfo->Fill(1.5,nAccepted);
    hCrossSectionInfo->Fill(2.5,sigmaGen);
    hCrossSectionInfo->Fill(3.5,nEvent);
    fout->Write();

    //    for(int i = 0; i < nResonances; i++)
    //    {
    //        fHistResonanceMass[i]->Write();
    //        fHistResonancePt[i]->Write();
    //        fHistResonanceEta[i]->Write();
    //        fHistResonancePhi[i]->Write();
    //    }

    fout->Close();
    cout <<"Successfully finished."<< endl;

    timer.Print();

}

/*
 * Main program
 */
int main(int argc, char **argv)
{
    TROOT root("pythia","run");
    //    root.LoadMacro( "../../../simpleEventAnalyzer/AliSimpleEvent.cxx+g" );
    //    root.ProcessLine(".L AliSimpleEvent_cxx.so" );

    int i=1;
    TString nEvt = argv[i++];
    TString config = argv[i++];
    TString random = argv[i++];
    TString outfile = argv[i++];
    cout <<nEvt.Atoi()<<"\t"<<random.Atoi()<<"\t"<< outfile.Data()<<endl;
    pythia8_spectra( nEvt.Atoi(), config.Data(), random.Atoi(), outfile.Data());
}  
