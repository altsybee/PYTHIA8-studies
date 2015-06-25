#ifndef StrangenessAnalysis_H
#define	StrangenessAnalysis_H


class TH1D;
class TRandom3;
class TLorentzVector;
class TF1;
class TH2D;
class TTree;
class TileCorrelations;
//class DiHadronAnalyser;
//struct MiniEvent;

class StrangenessAnalysis
{
public:
    StrangenessAnalysis();
//    void PrepareAnalysisAndRun();
    void RunAnalysis();

private:
//    TH2D* ExtractDihadronHist( DiHadronAnalyser *an );//TList *outList );
////    void UpdateEvent(MiniEvent *event, Double_t Eta , Double_t Phi, Double_t Pt , Int_t Charge);
//    void analyseMiniEvent(MiniEvent *miniEvent );
//    void analyseMixedMiniEvents(MiniEvent *miniEvent1, MiniEvent *miniEvent2 );

//    void decayStringIntoParticles( int nParticlesInString, TLorentzVector *vArr, double fictionRhoPt );
//    int probabilityChargePlusMinusNeutral();
//    void tuneHistogram2D(TH2D *h);
    void setVariablesForTrees();
    void tuneAnalysers();

    TRandom3 *fRand;
//    DiHadronAnalyser **analysersArray;

    TH1D *histEtaRho;
    TH1D *histEtaPion;
    TH1D *histPhiPion;
    TH1D *histPtRho;
    TH1D *histPtRhoWithWeight;
    TH1D *histPtPionFromRho;
    TH1D *histPtPionFromString;
    TH1D *histPtPionAll;
    TH1D *histPtPionWithWeight;

    TH1D *histPtInvM;
    TH1D *histPtMixedInvM;

//    MiniEvent *miniEv1, *miniEv2;

    //string fragmentation
//    StringFragmentation *strFr;

    TileCorrelations **analysersArray;
    int nAnalysers;

    // ##### prepare variables for the trees
    // particle tree vars
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

    // event info tree vars
    int nChargedFromTree;
//    int nPionsPrimary = 0;
//    int nPionsFromResonances = 0;
//    int nPionsFromRho0 = 0;
//    int nPionsFromRhoPlusMinus = 0;
//    int pionsFromRho0_bothDetected = 0; //incremented by 2 if both pions within cuts


    TTree *treeParticles;
    TTree *treeEvent;
};


#endif	/* StrangenessAnalysis_H */
