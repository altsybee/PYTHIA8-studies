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
#include <TVectorD.h>
#include <TArrayF.h>
#include <TMap.h>

#include <Pythia8/Pythia.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <map>
#include <vector>
#include <string>


//#include "../../../commonTools/Tools.cxx"
//#include "../../../commonTools/SimpleTrack.cxx"
//#include "../../../simpleEventAnalyzer/AliSimpleEvent.h"
//#include "../../../simpleEventAnalyzer/AliSimpleEvent.cxx"


using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class std::map<int,std::string>+;
#endif



void fixPhi( float &phi )
{
    if ( phi < 0 )
        phi += 2*TMath::Pi();
    if ( phi > 2*TMath::Pi() )
        phi -= 2*TMath::Pi();
}


const double etaRange = 2.0;
const double ptMin = 0.1;


map <int,string> pidMap;

void addPidToMap( Pythia8::Particle &particle )
{
    pidMap.insert ( pair<int,string>(particle.id(),particle.name()) );
}


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
 * Main event loop
 *
 * Int_t nEvent = The number of events in PYTHIA run
 * const char *config = cmnd file from which the PYTHIA configuration is read
 * Int_t random_seed = Seed for PYTHIA random number generator
 * const char *outfile = Name for the output file
 */
void pythia8_extractToTree(Int_t nEvent  = 100, const char *config="configPythia.cmnd", Int_t random_seed=0, const char *outfile="test.root")
{
    // File
    gROOT->cd();
    TFile *fout = new TFile(outfile,"recreate","pythia 8 generation");

    //create the file, the Tree and a few branches
    //    TFile file( "myTree.root", "recreate" );

    // Array for particle information
    const int nPIDforArray = 100;
    int particleIDArray[2][nPIDforArray] = {0};  // First row = particle ID, second row = number of particles with this ID found
    TString particleNames[nPIDforArray];         // Name for the particle which is in the same column in the particleIDArray

    TH1D *histPID = new TH1D("histPID","histPID",nPIDforArray,0,nPIDforArray);

    int eventId   ;
    int particleId;
    float pt  ;
    float eta ;
    float phi ;
    int charge    ;
    float y    ;

    int    id       ;
    int    status   ;
    int    mother1  ;
    int    mother2  ;
    int    daughter1;
    int    daughter2;
    float px       ;
    float py       ;
    float pz       ;
    float e        ;
    float m        ;
    bool  hasVertex;
    float xProd    ;
    float yProd    ;
    float zProd    ;
    float tProd    ;
    float tau      ;

    float theta         ;
    int   nDaughters    ;
    int   nMothers      ;
    int   nSisters      ;

    bool isFinal;

    bool   isCharged;
    bool   isNeutral;
    float tau0;
    bool   mayDecay;
    bool   canDecay;
    bool   isResonance;
    bool   isVisible;
    bool   isLepton;
    bool   isQuark;
    bool   isGluon;
    bool   isDiquark;
    bool   isParton;
    bool   isHadron;


    // ##### particle tree
    TTree *treeParticles = new TTree( "particleTree","a simple Tree with simple variables");
    treeParticles->Branch("eventId"     ,&eventId   ,"ev/I"         );
    treeParticles->Branch("particleId"  ,&particleId,"particleId/I" );
    treeParticles->Branch("pt"          ,&pt        ,"pt/F"         );
    treeParticles->Branch("eta"         ,&eta       ,"eta/F"        );
    treeParticles->Branch("phi"         ,&phi       ,"phi/F"        );
    treeParticles->Branch("charge"      ,&charge    ,"charge/I"     );
    treeParticles->Branch("y"           ,&y         ,"y/F"     );

    treeParticles->Branch("id"          ,&id        ,"id/I"         );
    //    treeParticles->Branch("status"      ,&status    ,"status/I"     );
    treeParticles->Branch("mother1"     ,&mother1   ,"mother1/I"    );
    treeParticles->Branch("mother2"     ,&mother2   ,"mother2/I"    );
    treeParticles->Branch("daughter1"   ,&daughter1 ,"daughter1/I"  );
    treeParticles->Branch("daughter2"   ,&daughter2 ,"daughter2/I"  );
    //    treeParticles->Branch("px"          ,&px        ,"px/F"         );
    //    treeParticles->Branch("py"          ,&py        ,"py/F"         );
    //    treeParticles->Branch("pz"          ,&pz        ,"pz/F"         );
    //    treeParticles->Branch("e"           ,&e         ,"e/F"          );
    //    treeParticles->Branch("m"           ,&m         ,"m/F"          );
    //    treeParticles->Branch("hasVertex"   ,&hasVertex ,"hasVertex/O"  );
    //    treeParticles->Branch("xProd"       ,&xProd     ,"xProd/F"      );
    //    treeParticles->Branch("yProd"       ,&yProd     ,"yProd/F"      );
    //    treeParticles->Branch("zProd"       ,&zProd     ,"zProd/F"      );
    //    treeParticles->Branch("tProd"       ,&tProd     ,"tProd/F"      );
    treeParticles->Branch("tau"         ,&tau       ,"tau/F"        );

    //    treeParticles->Branch("theta"       ,&theta     ,"theta/F"      );

    treeParticles->Branch("nDaughters"  ,&nDaughters,"nDaughters/I" );
    //    treeParticles->Branch("nMothers"    ,&nMothers  ,"nMothers/I"   );
    //    treeParticles->Branch("nSisters"    ,&nSisters  ,"nSisters/I"   );

    treeParticles->Branch("isFinal"   ,&isFinal        ,"isFinal/O" );

    treeParticles->Branch("isCharged"   ,&isCharged    ,"isCharged/O"   );
    //    treeParticles->Branch("isNeutral"   ,&isNeutral    ,"isNeutral/O"   );
    treeParticles->Branch("tau0"        ,&tau0         ,"tau0/F"        );
    //    treeParticles->Branch("mayDecay"    ,&mayDecay     ,"mayDecay/O"    );
    //    treeParticles->Branch("canDecay"    ,&canDecay     ,"canDecay/O"    );
    //    treeParticles->Branch("isResonance" ,&isResonance  ,"isResonance/O" );
    //    treeParticles->Branch("isVisible"   ,&isVisible    ,"isVisible/O"   );
    //    treeParticles->Branch("isLepton"    ,&isLepton     ,"isLepton/O"    );
    //    treeParticles->Branch("isQuark"     ,&isQuark      ,"isQuark/O"     );
    //    treeParticles->Branch("isGluon"     ,&isGluon      ,"isGluon/O"     );
    //    treeParticles->Branch("isDiquark"   ,&isDiquark    ,"isDiquark/O"   );
    treeParticles->Branch("isParton"    ,&isParton     ,"isParton/O"    );
    //    treeParticles->Branch("isHadron"    ,&isHadron     ,"isHadron/O"    );


    // ##### event info tree
    int nCharged = 0;
    //    int nPionsPrimary = 0;
    //    int nPionsFromResonances = 0;
    //    int nPionsFromRho0 = 0;
    //    int nPionsFromRhoPlusMinus = 0;
    //    int pionsFromRho0_bothDetected = 0; //incremented by 2 if both pions within cuts

    TTree *treeEvent = new TTree( "eventMultTree","a simple Tree with simple variables");
    treeEvent->Branch("nCharged",                &nCharged,              "nCharged/I");
    //    treeEvent->Branch("nPionsPrimary",           &nPionsPrimary,         "nPionsPrimary/I");
    //    treeEvent->Branch("nPionsFromResonances",    &nPionsFromResonances,  "nPionsFromResonances/I");
    //    treeEvent->Branch("nPionsFromRho0",          &nPionsFromRho0,        "nPionsFromRho0/I");
    //    treeEvent->Branch("nPionsFromRhoPlusMinus",  &nPionsFromRhoPlusMinus,"nPionsFromRhoPlusMinus/I");
    //    treeEvent->Branch("pionsFromRho0_bothDetected",&pionsFromRho0_bothDetected,"pionsFromRho0_bothDetected/I");

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
    for (Int_t iEvent = 0; iEvent < nEvent; iEvent++)
    {
        // Generate Event
        pythia8.next();
        Pythia8::Event &event = pythia8.event;

        // Print out some information about the first event
        if (iEvent < 1) event.list();

        int numberOfParticles = event.size();
        nCharged = 0;

        // Print out information about the progress of the simulation
        if( iEvent % 1000 == 0 )
            cout << "generating event " << iEvent << "..." << endl;

        // Particle Loop to count particles in eta range
        for (Int_t i = 0; i < numberOfParticles; i++)
        {
            if (pythia8.event[i].isFinal() && pythia8.event[i].isCharged()
                    && fabs(pythia8.event[i].eta()) < 2.4 )
                ++nCharged;
        }

        //nCharged cut
        if ( nCharged == 0)
            continue;

        treeEvent->Fill();

        // Particle Loop to fill particle tree
        for (Int_t i = 0; i < numberOfParticles; i++)
        {
            // choose the next particle
            Pythia8::Particle &particle = event[i];

            eventId     = iEvent;
            particleId  = i;
            pt          = particle.pT();
            eta         = particle.eta();
            phi         = particle.phi();            fixPhi( phi );
            charge      = particle.charge();
            y           = particle.y();

            id          = particle.id();
            status      = particle.status();
            mother1     = particle.mother1();
            mother2     = particle.mother2();
            daughter1   = particle.daughter1();
            daughter2   = particle.daughter2();
            px          = particle.px();
            py          = particle.py();
            pz          = particle.pz();
            e           = particle.e ();
            m           = particle.m();
            hasVertex   = particle.hasVertex();
            xProd       = particle.xProd();
            yProd       = particle.yProd();
            zProd       = particle.zProd();
            tProd       = particle.tProd();
            tau         = particle.tau();

            theta         = particle.theta();
            nDaughters    = particle.daughterList().size();
            nMothers      = particle.motherList().size();
            nSisters      = particle.sisterList().size();

            isFinal         = particle.isFinal();

            isCharged       = particle.isCharged();
            isNeutral       = particle.isNeutral();
            tau0            = particle.tau0();
            mayDecay        = particle.mayDecay();
            canDecay        = particle.canDecay();
            isResonance     = particle.isResonance();
            isVisible       = particle.isVisible();
            isLepton        = particle.isLepton();
            isQuark         = particle.isQuark();
            isGluon         = particle.isGluon();
            isDiquark       = particle.isDiquark();
            isParton        = particle.isParton();
            isHadron        = particle.isHadron();

            treeParticles->Fill();

            //add particle pid to map
            if (eventId < 5000)
            {
                addPidToMap(particle);
                //increment colomn for the pid-mother of this pion (cound be resonance or quark!)
                addParticleToArray( particle, particleIDArray, particleNames );
            }


            // Lambda
            if(0)if ( particle.id() == 3122 )
            {
                Pythia8::Particle &d1 = event[particle.daughter1()];
                Pythia8::Particle &d2 = event[particle.daughter2()];

                cout << particle.name() << " mother=" << event[particle.mother1()].name()
                     << " " << d1.name() //particle.daughter1()
                     << " " << d2.name() //particle.daughter2()
                     << endl;
            }

        }//particle loop
    }//event loop


    // Sort the array so that the particles with the most entries are the first
    sortArrays(particleIDArray,particleNames);

    for(int i = 0; i < nPIDforArray; i++)
    {
        histPID->Fill( i+0.5, particleIDArray[1][i] );
        histPID->GetXaxis()->SetBinLabel( i+1, particleNames[i].Data()) ;
    }

    // Display some statistics
    pythia8.stat();
    //    long nTried = pythia8.info.nTried();
    //    long nAccepted = pythia8.info.nAccepted();
    //    double sigmaGen = pythia8.info.sigmaGen();
    //    //double sigmaErr = pythia8.info.sigmaErr();
    //    hCrossSectionInfo->Fill(0.5,nTried);
    //    hCrossSectionInfo->Fill(1.5,nAccepted);
    //    hCrossSectionInfo->Fill(2.5,sigmaGen);
    //    hCrossSectionInfo->Fill(3.5,nEvent);


    //#### write pid-name correspondance
    TArrayF *arrPid = new TArrayF(pidMap.size());
    TH1D *histArrPidNames = new TH1D("histArrPidNames","histArrPidNames",500,0,500);

    cout << "pidMap.size()=" << pidMap.size() << endl;
    // check pid map
    int iPid = 0;
    for (auto it = pidMap.begin(); it != pidMap.end(); ++it)
    {
        cout << iPid << ": " << (*it).first << " : " << (*it).second << endl;
        arrPid->AddAt((*it).first,iPid);
        histArrPidNames->GetXaxis()->SetBinLabel( iPid+1, (*it).second.c_str() ) ;
        iPid++;
    }
    fout->WriteObject(arrPid,"arrPid");

    fout->ls();

    fout->Write();
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
    pythia8_extractToTree( nEvt.Atoi(), config.Data(), random.Atoi(), outfile.Data());
}
