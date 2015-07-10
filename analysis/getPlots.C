#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLine.h"
#include "TPave.h"

void getYields();
void getGraphWithCoeffs();
void getHistById(int);

TCanvas *canvas;

TFile *fileA;
TFile *fileB0;
TFile *fileB2;
TFile *fileC;
TFile *fileTEST;

TGraphErrors *getGraphFromProfile( TProfile *pr )
{
    TGraphErrors *gr = new TGraphErrors;

    int iPoint = 0;
    for ( int i = 0; i < pr->GetNbinsX(); ++i )
    {
        if (pr->GetBinContent(i+1) > 0)
        {
            gr->SetPoint(iPoint, pr->GetBinCenter(i+1), pr->GetBinContent(i+1));
            //            gr->SetPointError(iPoint, pr->GetBinCenter(i+1), pr->GetBinError(i+1));
            iPoint++;
        }
    }
    return gr;
}


void setGraphPointAsRatio( TGraphErrors *gr, TGraphErrors *gr1, TGraphErrors *gr2 )
{
    double x1, y1, x2, y2;
    double ye1, ye2;
    for ( int i = 0; i < gr->GetN(); ++i )
    {
        gr1->GetPoint(i,x1,y1);
        gr2->GetPoint(i,x2,y2);

        ye1 = gr1->GetErrorY(i);
        ye2 = gr2->GetErrorY(i);
        if (y2>0)
        {
            gr->SetPoint( i, x1, y1/y2 );
            gr->SetPointError( i, 0, sqrt(ye1*ye1+ye2*ye2)/y2 );
        }
    }
}


void fitProfile(TProfile *prof, TGraphErrors *gr)
{
    double meanX = prof->GetMean();
    cout << "meanX=" << meanX << endl;
    prof->Fit("pol1");
    //    prof->Fit("pol1", "", "", 0.05, 1.5);

    TF1 *myfit = (TF1*) prof->GetFunction("pol1");
    double a = myfit->GetParameter(0);
    double b = myfit->GetParameter(1);
    double bErr = myfit->GetParError(1);
    cout << "corr coeff = " << b << endl;
    int iPoint = gr->GetN();
    gr->SetPoint( iPoint, /*2*etaMinF[i]*/iPoint*0.5, b );// /coeffForB[iType] );
    gr->SetPointError( iPoint, 0, bErr );
}


const char *strProfile[] = {
    "histNN"
    , "histNsNs"
    , "histNsN"
    , "histNskN"

    , "histNlN"
    , "histNlNl"
    , "histNprotonN"

    , "histPtN"
    , "histPtNs"
    , "histPtNsk"
    , "histPtNl"

    , "histKtoPion_N"
    , "histLambdaToPion_N"
    , "histProtonToPion_N"
    , "histKtoPion_KtoPion"
    , "histLambdaToPion_LambdaToPion"

    , "histKtoPion_Pt"
    , "histLambdaToPion_Pt"
    , "histProtonToPion_Pt"
};



void getPlots()
{
//    fileA = new TFile("analysisResults_1000k_ct1cm_Monash2013_pure_newWins.root");  //  "analysisResults_100k_Monash2013_pure.root");
    //    //    fileB0 = new TFile("analysisResults_100k_ct1cm_Monash2013_newCR_mode0.root");  //"analysisResults_100k_newCR_mode0.root");
//    fileB2 = new TFile("analysisResults_1000k_ct1cm_Monash2013_newCR_mode2_newWins.root");
    //    fileC = new TFile("analysisResults_1000k_ct1cm_Monash2013_noCR_newWins.root");

        fileA = new TFile("analysisResults_1000k_ct1cm_Monash2013_pure_newWins2_Y.root");  //  "analysisResults_100k_Monash2013_pure.root");
//        fileB2 = new TFile("analysisResults_1000k_ct1cm_Monash2013_newCR_mode2_newWins2_Y.root");
//        fileC = new TFile("analysisResults_1000k_ct1cm_Monash2013_noCR_newWins2_Y.root");


    fileB2 = new TFile("analysisResults_1000k_ct1cm_Monash2013_newCR_mode2_newWins2_Y.root");
    fileC = new TFile("analysisResults_1000k_ct1cm_Monash2013_noCR_newWins2_Y.root");

//        fileA = new TFile("analysisResults_1000k_ct1cm_Monash2013_noCR_newWins2_Y.root");  //  "analysisResults_100k_Monash2013_pure.root");
//        fileC = new TFile("analysisResults_1000k_ct1cm_Monash2013_pure_newWins2_Y.root");

//        fileA = new TFile("analysisResults_1000k_ct1cm_Monash2013_pure_newWins2_Y_withPtPt.root");  //  "analysisResults_100k_Monash2013_pure.root");

        //    fileA = new TFile("analysisResults100k_ct1cm_Monash2013_pure.root");  //  "analysisResults_100k_Monash2013_pure.root");
    //    fileB0 = new TFile("analysisResults10k_ct1cm_Monash2013_newCRmode0.root");  //"analysisResults_100k_newCR_mode0.root");
    //    fileB2 = new TFile("analysisResults_100k_newCR_mode3.root");
    //    fileC = new TFile("analysisResults_100k_Monash2013_no_CR.root");

    //    fileTEST = new TFile("analysisResults10kTEST_ct1cm_Monash2013_newCR_mode0.root");

//        getYields();
        getHistEtaY_in_wins();
        return;

//    getGraphWithCoeffs();
//    return;

    //prepare canvas
    canvas = new TCanvas("canvas","Data Plots",200,10,800,680);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.02);
    canvas->SetTopMargin(0.02);
    canvas->SetBottomMargin(0.12);




    for ( int i = 0; i < 1; ++i )
        getHistById(i);

    //        for ( int i = 0; i < 5; ++i )
    //            getHistById_2(i);


    //    fileA->Close();
    //    fileB0->Close();
    //    fileB2->Close();
    //    fileC->Close();
    //    fileTEST->Close();

}

void getHistById(int histId )
{
    canvas->Clear();

    //        const char *myStr = "histNN"                            ;
    //    const char *myStr = "histNsNs"                          ;
    //    const char *myStr = "histNsN"                           ;
    //    const char *myStr = "histNskN"                          ;
    //    const char *myStr = "                                   ;
    //    const char *myStr = "histNlN"                           ;
    //    const char *myStr = "histNlNl"                          ;
    //        const char *myStr = "histNprotonN"                      ;
    //    const char *myStr = "                                   ;
    const char *myStr = "histPtN"                           ;
    //    const char *myStr = "histPtNs"                          ;
    //    const char *myStr = "histPtNsk"                         ;
    //    const char *myStr = "histPtNl"                          ;
    //    const char *myStr = "                                   ;
    //        const char *myStr = "histKtoPion_N"                     ;
    //        const char *myStr = "histLambdaToPion_N"                ;
    //        const char *myStr = "histProtonToPion_N"                ;
    //    const char *myStr = "histKtoPion_KtoPion"               ;
    //    const char *myStr = "histLambdaToPion_LambdaToPion"     ;
    //    const char *myStr = "                                   ;
    //        const char *myStr = "histKtoPion_Pt"                    ;
    //        const char *myStr = "histLambdaToPion_Pt"               ;
    //        const char *myStr = "histProtonToPion_Pt"               ;

    //    char *myStr = strProfile[histId];

    int winId = 13;//11;
    //    TProfile *profile2 = file->GetObjectChecked( "histLambdaToPion_Pt2_pfx", "TProfile" );
    //    TProfile *profile4 = file->GetObjectChecked( "histLambdaToPion_Pt4_pfx", "TProfile" );
    TProfile *prof_A = (TProfile *)fileA->GetObjectChecked( Form("%s%d_pfx", myStr, winId), "TProfile" );
    //    TProfile *prof_B0 = fileB0->GetObjectChecked( Form("%s%d_pfx", myStr, winId), "TProfile" );
    //    TProfile *prof_B2 = fileB2->GetObjectChecked( Form("%s%d_pfx", myStr, winId), "TProfile" );
    TProfile *prof_C = (TProfile *)fileC->GetObjectChecked( Form("%s%d_pfx", myStr, winId), "TProfile" );

    TProfile *prof_B2 = (TProfile *)fileA->GetObjectChecked( Form("%s%d_pfx", myStr, 12), "TProfile" );

    //    TProfile *prof_TEST = fileTEST->GetObjectChecked( Form("%s%d_pfx", myStr, winId), "TProfile" );


    //    TGraphErrors *gr2 = getGraphFromProfile(prof);
    //    TGraphErrors *gr4 = getGraphFromProfile(profile4);
    //    gr4->SetLineColor(kRed);

    //    gr2->Draw("APL");
    //    gr4->Draw("same PL");

    //    prof_B0->SetLineColor(kRed);
    //    prof_B2->SetLineColor(kGreen+1);
    prof_B2->SetLineColor(kRed);
    prof_C->SetLineColor(kMagenta);
    //    prof_TEST->SetLineColor(kBlack);
    //    prof_B2->SetLineWidth(2);

    prof_A->SetMarkerStyle(21);
    prof_A->SetMarkerSize(0.9);
    prof_A->SetMarkerColor(kBlue-3);

    prof_B2->SetMarkerStyle(20);
    prof_B2->SetMarkerSize(0.9);
    //    prof_B2->SetMarkerColor(kGreen+1);
    prof_B2->SetMarkerColor(kRed);

    prof_C->SetMarkerStyle(24);
    prof_C->SetMarkerSize(0.6);
    prof_C->SetMarkerColor(kMagenta);


    prof_C->SetLineStyle(2);


    // ##### drawing
  gStyle->SetOptStat(false);
    TPad *pad1, *pad2;

    // ##### PAD 1
    pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.017);
    //    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.02);
    pad1->Draw();
    pad1->cd();

    prof_A->Draw( );
    prof_B2->Draw("same");
    //    prof_C->Draw("same");
    //    prof_B0->Draw("same");
    //    prof_TEST->Draw("same");
    //    prof_A->Draw("hist");
    //    prof_B2->Draw("hist same");
    //    prof_C->Draw("hist same");

    prof_A->GetXaxis()->SetLabelSize(0.05);
    prof_A->GetYaxis()->SetLabelSize(0.05);

    prof_A->GetXaxis()->SetTitleSize(0.058);
    prof_A->GetYaxis()->SetTitleSize(0.065);

    prof_A->GetXaxis()->SetTitleOffset( 0.9+1 ); //by hand make in large
    prof_A->GetYaxis()->SetTitleOffset( 0.68 );

    prof_A->GetXaxis()->SetLabelOffset( 0.9+1 ); //by hand make in large

    //    prof_A->GetXaxis()->SetTitle( "#bar{p_{T}}_{F}" );
    //p/pi vs pt
    //    prof_A->GetXaxis()->SetTitle( "#LTp_{T}#GT_{F}, GeV/c" );
    //    prof_A->GetYaxis()->SetTitle( "#LTp/#pi#GT_{B}" );


    TLatex *tex = new TLatex(0.3,0.84, "#sqrt{s}=7 TeV, p_{T}#in(0.2, 2.0)GeV/c");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextSize(0.046);
    //    tex->DrawClone();


    TLegend *leg = new TLegend(0.12,0.7,0.67,0.975);
    //    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    //    leg->AddEntry(prof_A,  "Monash2013",  "pl");
    //    leg->AddEntry(prof_B0, "Monash2013 new CR mode0", "pl");
    //    leg->AddEntry(prof_B2, "Monash2013 new CR", "pl"); // mode2
    //    leg->AddEntry(prof_C,  "Monash2013 no CR",   "pl");

    leg->AddEntry(prof_A,  "|y|<0.5",  "pl");
    leg->AddEntry(prof_B2, "#minus2 < y_{B} < #minus1, 1 < y_{F} < 2", "pl"); // mode2

    leg->Draw();



    // ##### PAD 2
    canvas->cd();
    pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.33);
    //    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.02);
    pad2->Draw();
    pad2->cd();

    TH1D *h1 = prof_A->ProjectionX();
    TH1D *h2 = prof_B2->ProjectionX();
    TH1D *h3 = prof_C->ProjectionX();

    h2->GetXaxis()->SetLabelSize(0.12);
    h2->GetYaxis()->SetLabelSize(0.11);

    h2->GetXaxis()->SetTitleSize(0.15);
    h2->GetYaxis()->SetTitleSize(0.13);

    h2->GetXaxis()->SetTitleOffset( 0.88 );
    h2->GetYaxis()->SetTitleOffset( 0.35 );

    h2->GetYaxis()->SetTitle("ratio");
    h2->GetYaxis()->CenterTitle();


    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(0.9);
    //    h2->SetMarkerColor(kGreen+1);
    //    h2->SetLineColor(kGreen+1);
    h2->SetMarkerColor(kRed);
    h2->SetLineColor(kRed);

    h3->SetMarkerStyle(24);
    h3->SetMarkerSize(0.6);
    h3->SetMarkerColor(kMagenta);
    h3->SetLineColor(kMagenta);

    //    TProfile *histDivided = new TProfile ( *prof_B2 );

    //    histDivided->Sumw2();
    //    prof_A->Sumw2();
    //    prof_B2->SetStats(0);
    h2->Divide(h1);
    h3->Divide(h1);
    //    h1->SetMarkerStyle(21);
    h2->Draw("p");
    //    h3->Draw("same p");



    // ####### SETTING LABELS
    //NlN
    h2->GetXaxis()->SetTitle( "N_{ch}^{F}" );
    //    h2->GetXaxis()->SetTitle( "N_{ch}^{F}" );
    //    h2->GetXaxis()->SetTitle( "#bar{p_{T}^{F}}" );
    //    h2->GetXaxis()->SetTitle( "N_{#Lambda}^{F}" );

    //    prof_A->GetYaxis()->SetTitle( "N_{#Lambda}^{B}" );
    //    prof_A->GetYaxis()->SetTitle( "N_{ch}^{B}" );
    //    prof_A->GetYaxis()->SetTitle( "N_{p}^{B}" );
    //    prof_A->GetYaxis()->SetTitle( "#LTK/#pi#GT_{B}" );
    //    prof_A->GetYaxis()->SetTitle( "#LT#Lambda/#pi#GT_{B}" );
    //    prof_A->GetYaxis()->SetTitle( "#LTp/#pi#GT_{B}" );
    //    prof_A->GetYaxis()->SetTitle( "#LTK/#pi#GT_{B}" );
    //    prof_A->GetYaxis()->SetTitle( "#LT#Lambda/#pi#GT_{B}" );
    //    prof_A->GetYaxis()->SetTitle( "#LTp/#pi#GT_{B}" );
    prof_A->GetYaxis()->SetTitle( "#LT#bar{p_{T}}#GT_{B}, GeV/#it{c}"); //, GeV/c" );


    // ####### SETTING RANGES
    //    h2->GetYaxis()->SetRangeUser(0.5,1.95);
    const double xAxisMax = 65;
    prof_A->GetXaxis()->SetRangeUser(0,xAxisMax);
    h2->GetXaxis()->SetRangeUser(0,xAxisMax);

    //    prof_A->GetYaxis()->SetRangeUser(0.,1.45);


    // draw 1-line
    TLine *line = new TLine(0,1,xAxisMax,1);
    //    TLine *line = new TLine(1,1,80,2);
    //    line->SetNDC();
    line->SetLineColor(kGray+3);
    line->SetLineStyle(9);
    line->Draw();

    canvas->cd();






    //    canvas->Update();
    canvas->Print( Form("plots%d/%s%d.png", winId, myStr, winId) ); //or .png, .pdf
    //    canvas->Update();
}






void getHistById_2(int histId )
{
    const char *strHist[] = {
        "histMult_charged"
        , "histMult_Lambda0"
        , "histMult_K0"
        , "histEta"
        , "histEtaStrange"
    };

    TH1D *histA = (TH1D*)fileA->GetObjectChecked( strHist[histId], "TH1D" );
    //    TH1D *histB0 = fileB0->GetObjectChecked( strHist[histId], "TH1D" );
    TH1D *histB2 = (TH1D*)fileB2->GetObjectChecked( strHist[histId], "TH1D" );
    TH1D *histC = (TH1D*)fileC->GetObjectChecked( strHist[histId], "TH1D" );
    //    TH1D *histTEST = fileTEST->GetObjectChecked( strHist[histId], "TH1D" );

    //    histB0->SetLineColor(kRed);
    histB2->SetLineColor(kGreen);
    histC->SetLineColor(kMagenta);
    //    histTEST->SetLineColor(kBlack);

    histA->DrawNormalized();
    //    histB0->DrawNormalized("same");
    histB2->DrawNormalized("same");
    histC->DrawNormalized("same");
    //    histTEST->DrawNormalized("same");

    if ( histId < 3 )
        gPad->SetLogy();
    else
        gPad->SetLogy(0);


    TLegend *leg = new TLegend(0.45,0.2,0.95,0.345);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(histA,  "Monash2013 old CR",  "l");
    //    leg->AddEntry(histB0, "Monash2013 new CR mode0", "l");
    leg->AddEntry(histB2, "Monash2013 new CR mode2", "l");
    leg->AddEntry(histC,  "Monash2013 no CR",   "l");
    //    leg->AddEntry(histTEST,  "TEST",   "l");
    leg->Draw();


    //    canvas->Update();
    canvas->Print( Form("plots/%s.png", strHist[histId]) ); //or .png, .pdf

}


void getGraphWithCoeffs()
{
    //    TCanvas *canvNN = new TCanvas("canvNN","canvNN",500,250,1200,800 );
    //    canvNN->Divide(3,3);
    //    canvNN->SetTopMargin(0.05);
    //    canvNN->SetRightMargin(0.15);
    //    canvNN->SetBottomMargin(0.1);
    //    canvNN->SetLeftMargin(0.15);


    TCanvas *canv_profiles = new TCanvas("canv_profiles","canv_profiles",150,100,1200,800 );
    canv_profiles->Divide(3,3);
    canv_profiles->SetTopMargin(0.05);
    canv_profiles->SetRightMargin(0.15);
    canv_profiles->SetBottomMargin(0.1);
    canv_profiles->SetLeftMargin(0.15);


    //    const char *myStr = "histNN"                            ;
    //    const char *myStr = "histNsNs"                          ;
    //    const char *myStr = "histNsN"                           ;
    //    const char *myStr = "histNskN"                          ;
    //    const char *myStr = "                                   ;
    //    const char *myStr = "histNlN"                           ;
    //    const char *myStr = "histNlNl"                          ;
    //    const char *myStr = "histNprotonN"                      ;
    //    const char *myStr = "                                   ;
    //    const char *myStr = "histPtN"                           ;
    //    const char *myStr = "histPtNs"                          ;
    //    const char *myStr = "histPtNsk"                         ;
    //    const char *myStr = "histPtNl"                          ;
    //    const char *myStr = "                                   ;
    //    const char *myStr = "histKtoPion_N"                     ;
    //    const char *myStr = "histLambdaToPion_N"                ;
    //    const char *myStr = "histProtonToPion_N"                ;
    //    const char *myStr = "histKtoPion_KtoPion"               ;
    //    const char *myStr = "histLambdaToPion_LambdaToPion"     ;
    //    const char *myStr = "                                   ;
    //    const char *myStr = "histKtoPion_Pt"                    ;
    //    const char *myStr = "histLambdaToPion_Pt"               ;
    //    const char *myStr = "histProtonToPion_Pt"               ;

    //    int winId = 1;

    const int nTypes = 3;
     char *strProfileForGraphs[nTypes] = {
        "histNN"
//        , "histNskN" //
        , "histNsNs"
        , "histNlNl"
    };

    //    const char *strProfileForGraphs[nTypes] = {
    //        "histKtoPion_Pt"
    //        , "histLambdaToPion_Pt"
    //        , "histProtonToPion_Pt"
    //    };


    const int nPairs = 7;

    //    const int nGraphMarkers[nTypes] = { 21, 24, 25};
    const int nGraphMarkers[nTypes] = {  21, 20, 22 };
    const int nGraph2Markers[nTypes] = { 25, 24, 26 };

    const double coeffForB[nTypes] = { 1., 0.15/0.7, 0.07/0.7};

    TGraphErrors *graph[nTypes];
    TGraphErrors *graph2[nTypes];
    for (int iType = 0; iType < nTypes; ++iType)
    {
        graph[iType] = new TGraphErrors();
        graph2[iType] = new TGraphErrors();

        char *myStr = strProfileForGraphs[iType];

        for (int i = 0; i < nPairs; ++i)
        {
            //        TH1D *histA = fileA->GetObjectChecked( strHist[histId], "TH1D" );
            //        TH1D *histB0 = fileB0->GetObjectChecked( strHist[histId], "TH1D" );
            //        TH1D *histB2 = fileB2->GetObjectChecked( strHist[histId], "TH1D" );
            //        TH1D *histC = fileC->GetObjectChecked( strHist[histId], "TH1D" );

            //2D hists
            //        canvNN->cd(i+1);
            //        histB2->DrawCopy("colz");
            //profiles
            canv_profiles->cd(i+1);
            //            histPtN[i]->ProfileX()->DrawCopy();

            //        TProfile *prof_A = fileA->GetObjectChecked( Form("%s%d_pfx", myStr, i), "TProfile" );
            //        TProfile *prof_B0 = fileB0->GetObjectChecked( Form("%s%d_pfx", myStr, i), "TProfile" );
            TProfile *prof_B2 = (TProfile *)fileB2->GetObjectChecked( Form("%s%d_pfx", myStr, i), "TProfile" );
            TProfile *prof_C = (TProfile *)fileC->GetObjectChecked( Form("%s%d_pfx", myStr, i), "TProfile" );

            //            double meanX = prof_B2->GetMean();
            //            cout << "meanX=" << meanX << endl;
            //            prof_B2->Fit("pol1");
            //            TF1 *myfit = (TF1*) prof_B2->GetFunction("pol1");
            //            double a = myfit->GetParameter(0);
            //            double b = myfit->GetParameter(1);
            //            double bErr = myfit->GetParError(1);
            //            cout << "corr coeff = " << b << endl;
            //            graph[iType]->SetPoint( i, /*2*etaMinF[i]*/i, b );// /coeffForB[iType] );
            //            graph[iType]->SetPointError( i, 0, bErr );

            fitProfile( prof_B2, graph[iType] );
            fitProfile( prof_C, graph2[iType] );

        }
    }
    //        histPtN[0]->SetMarkerColor(kRed);
    //        histPtN[0]->SetMarkerStyle(25);


    TCanvas *canvGraph = new TCanvas("canvGraph","canvGraph",110,40,800,600 );
    canvGraph->SetTopMargin(0.05);
    canvGraph->SetRightMargin(0.05);
    canvGraph->SetBottomMargin(0.11);
    canvGraph->SetLeftMargin(0.12);

    // ##### drawing
  gStyle->SetOptStat(false);
    TPad *pad1, *pad2;

    // ##### PAD 1
    pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.017);
    //    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.02);
    pad1->Draw();
    pad1->cd();

    pad1->SetGridy();

    TGraphErrors *graphForFrame = new TGraphErrors();
    graphForFrame->SetPoint(0,-0.1,0);
    graphForFrame->SetPoint(1,3.1,0);
    graphForFrame->SetPoint(2,3.1,1);
    graphForFrame->SetPoint(3,-0.1,1);
    graphForFrame->SetMarkerColor(kWhite);
    graphForFrame->Draw ( "AP" ) ;

    //            graphForFrame->SetTitle( ";#Delta#eta;b" );
    graphForFrame->SetTitle( ";#eta gap;b_{corr}" );
    graphForFrame->SetMinimum( 0 );

    graphForFrame->GetXaxis()->SetTitleSize(0.058);
    graphForFrame->GetYaxis()->SetTitleSize(0.068);

    graphForFrame->GetXaxis()->SetLabelOffset( 2 );
    graphForFrame->GetXaxis()->SetTitleOffset( 0.87 );
    graphForFrame->GetYaxis()->SetTitleOffset( 0.65 );

    graphForFrame->GetXaxis()->SetLabelSize(0.055);
    graphForFrame->GetYaxis()->SetLabelSize(0.055);

    graphForFrame->GetYaxis()->SetRangeUser(0, 0.75);
    graphForFrame->GetYaxis()->CenterTitle(true);

    int kGraphBcorrColors[] = { kBlack, kBlue,  kRed};
//    TLegend *leg = new TLegend(0.57,0.65,0.65,0.97,"#eta " ); //pseudo-rapidity" );// windows");
    TLegend *leg = new TLegend(0.57,0.65,0.65,0.97,"CR " ); //pseudo-rapidity" );// windows");
//    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

//    TLegend *leg2 = new TLegend(0.6,0.65,0.95,0.97,"   y");
    TLegend *leg2 = new TLegend(0.6,0.65,0.95,0.97,"   no CR");
//    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);

//    TString strLegend[] = { "all charged" , "K^{#plus/#minus}, K^{0}_{S}, #Lambda^{0}", "#Lambda^{0}" };
    TString strLegend[] = { "all charged" , "K^{0}_{S}, #Lambda^{0}", "#Lambda^{0}" };

    for (int iType = 0; iType < nTypes; ++iType)
    {
        //    canvGraph->cd(nPairs+1);
        //        if (iType==0)
        //        {
        //            graph[iType]->Draw ( "APL" ) ;

        // //            graph[iType]->SetTitle( ";#Delta#eta;b" );
        //            graph[iType]->SetTitle( ";#eta gap;b_{corr}" );
        //            graph[iType]->SetMinimum( 0 );

        //            graph[iType]->GetXaxis()->SetTitleSize(0.058);
        //            graph[iType]->GetYaxis()->SetTitleSize(0.058);

        //            graph[iType]->GetXaxis()->SetTitleOffset( 0.87 );
        //            graph[iType]->GetYaxis()->SetTitleOffset( 0.9 );

        //            graph[iType]->GetXaxis()->SetLabelSize(0.055);
        //            graph[iType]->GetYaxis()->SetLabelSize(0.055);

        //            graph[iType]->GetYaxis()->SetRangeUser(0, 1);
        //        }
        //        else
        graph[iType]->Draw ( "P" ) ;

        //graph[iType]->SetMarkerStyle(22);
        graph[iType]->SetMarkerStyle(nGraphMarkers[iType]);
        graph[iType]->SetMarkerSize(1.5);
        graph[iType]->SetMarkerColor(kGraphBcorrColors[iType]); //kRed);
        graph[iType]->SetLineColor(kGraphBcorrColors[iType]); //kRed);


        graph2[iType]->Draw ( "P" ) ;
        graph2[iType]->SetMarkerStyle(nGraph2Markers[iType]);
        graph2[iType]->SetMarkerSize(1.5);
        graph2[iType]->SetMarkerColor(kGraphBcorrColors[iType]); //kBlue);
        graph2[iType]->SetLineColor(kGraphBcorrColors[iType]); //kBlue);

        //    leg->SetFillColor(0);
//        leg->AddEntry(graph[iType],  strLegend[iType],  "P");
        leg->AddEntry(graph[iType], "", "P"); // mode2
        leg2->AddEntry(graph2[iType], strLegend[iType], "P"); // mode2
    }
    leg->Draw();
    leg2->Draw();

    TPave *pave = new TPave(1.615998,0.6505508,3.298958,0.9819344,4,"br");

//    ci = 925;
//    color = new TColor(ci, 1, 1, 1, " ", 0);
//    pave->SetFillColor(ci);

//    ci = TColor::GetColor("#000000");
//    pave->SetLineColor(ci);
    pave->Draw();

    TLatex *tex = new TLatex(0.3,0.84, "#sqrt{s}=7 TeV, p_{T}#in(0.2, 2.0)GeV/c");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextSize(0.046);
    //    tex->DrawClone();



    // ##### PAD 2
    canvGraph->cd();
    pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.33);
    //    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.02);
    pad2->Draw();
    pad2->cd();

    TGraphErrors *graphForFrameRatio = new TGraphErrors();
    graphForFrameRatio->SetPoint(0,-0.1,0);
    graphForFrameRatio->SetPoint(1,3.1,0);
    graphForFrameRatio->SetPoint(2,3.1,1);
    graphForFrameRatio->SetPoint(3,-0.1,1);
    graphForFrameRatio->SetMarkerColor(kWhite);
    graphForFrameRatio->Draw ( "AP" ) ;

    pad2->SetGridy();
    graphForFrameRatio->SetTitle( ";#eta (y) gap;ratio" );
    //    graphForFrameRatio->SetMinimum( 0 );

    //    graphForFrameRatio->GetXaxis()->SetTitleSize(0.058);
    //    graphForFrameRatio->GetYaxis()->SetTitleSize(0.058);

    //    graphForFrameRatio->GetXaxis()->SetTitleOffset( 0.87 );
    //    graphForFrameRatio->GetYaxis()->SetTitleOffset( 0.9 );

    //    graphForFrameRatio->GetXaxis()->SetLabelSize(0.055);
    //    graphForFrameRatio->GetYaxis()->SetLabelSize(0.055);

    graphForFrameRatio->GetYaxis()->SetRangeUser(0.63, 1.25);

    graphForFrameRatio->GetXaxis()->SetLabelSize(0.15);
    graphForFrameRatio->GetYaxis()->SetLabelSize(0.13);

    graphForFrameRatio->GetXaxis()->SetTitleSize(0.17);
    graphForFrameRatio->GetYaxis()->SetTitleSize(0.13);

    graphForFrameRatio->GetXaxis()->SetTitleOffset( 0.88 );
    graphForFrameRatio->GetYaxis()->SetTitleOffset( 0.35 );

    graphForFrameRatio->GetYaxis()->SetTitle("ratio");
    graphForFrameRatio->GetYaxis()->CenterTitle();

    //    TH1D *h1 = prof_A->ProjectionX();
    //    TH1D *h2 = prof_B2->ProjectionX();
    //    TH1D *h3 = prof_C->ProjectionX();

    TGraphErrors *graphBcorrRatio[nTypes];
    for (int iType = 0; iType < nTypes; ++iType)
    {
        graphBcorrRatio[iType] = (TGraphErrors*)graph[iType]->Clone( Form("grClone%d",iType ) );
//        setGraphPointAsRatio( graphBcorrRatio[iType], graph[iType], graph2[iType] );
        setGraphPointAsRatio( graphBcorrRatio[iType], graph2[iType], graph[iType] );
        graphBcorrRatio[iType]->Draw ( "P" ) ;

        //graph[iType]->SetMarkerStyle(22);
        graphBcorrRatio[iType]->SetMarkerStyle(nGraphMarkers[iType]);
        graphBcorrRatio[iType]->SetMarkerSize(1.5);
        graphBcorrRatio[iType]->SetMarkerColor(kGraphBcorrColors[iType]); //);
        graphBcorrRatio[iType]->SetLineColor(kGraphBcorrColors[iType]); //);



        //        //    leg->SetFillColor(0);
        //        leg->SetFillStyle(0);
        //        leg->SetBorderSize(0);
        //        leg->AddEntry(graph[iType],  "test",  "pl");
        //        leg->AddEntry(graph2[iType], "test new CR", "pl"); // mode2
    }

    // draw 1-line
    TLine *line = new TLine( -0.4, 1, 3.4, 1 );
    //    TLine *line = new TLine(1,1,80,2);
    //    line->SetNDC();
    line->SetLineColor(kGray+3);
    line->SetLineStyle(9);
    line->Draw();



}

void getYields()
{

    TCanvas *canvasYields = new TCanvas("canvasYields","Data Plots",200,10,800,680);
    canvasYields->SetLeftMargin(0.15);
    canvasYields->SetRightMargin(0.02);
    canvasYields->SetTopMargin(0.02);
    canvasYields->SetBottomMargin(0.12);


    TH1D *histYields_A = (TH1D *)fileA->GetObjectChecked( "histYields", "TH1D" );
    TH1D *histYields_B2 = (TH1D *)fileB2->GetObjectChecked( "histYields", "TH1D" );
    TH1D *histYields_C = (TH1D *)fileC->GetObjectChecked( "histYields", "TH1D" );

    // ##### drawing
  gStyle->SetOptStat(false);
    TPad *pad1, *pad2;

    // ##### PAD 1
    pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.017);
    //    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.02);
    pad1->Draw();
    pad1->cd();



    pad1->SetLogy();

    histYields_A->SetLineColor(kBlue-3);
    histYields_B2->SetLineColor(kGreen+1);
    histYields_C->SetLineColor(kMagenta);

    histYields_A->SetMarkerStyle(21);
    //    histYields_A->SetMarkerSize(0.6);
    histYields_A->SetMarkerColor(kBlue-3);

    histYields_B2->SetMarkerStyle(20);
    //    histYields_B2->SetMarkerSize(0.6);
    histYields_B2->SetMarkerColor(kGreen+1);

    histYields_C->SetMarkerStyle(24);
    //    histYields_C->SetMarkerSize(0.6);
    histYields_C->SetMarkerColor(kMagenta);

    histYields_A->GetYaxis()->SetTitle("yields");


    histYields_A->GetXaxis()->SetLabelSize(0.05);
    histYields_A->GetYaxis()->SetLabelSize(0.05);

    histYields_A->GetXaxis()->SetTitleSize(0.058);
    histYields_A->GetYaxis()->SetTitleSize(0.065);

    histYields_A->GetXaxis()->SetTitleOffset( 0.9+1 ); //by hand make in large
    histYields_A->GetYaxis()->SetTitleOffset( 0.68 );

    histYields_A->GetXaxis()->SetLabelOffset( 0.9+1 ); //by hand make in large

//    histYields_A->DrawNormalized("p");
//    histYields_B2->DrawNormalized("p same");
//    histYields_C->DrawNormalized("p same");
    histYields_A->Draw("p");
    histYields_B2->Draw("p same");
    histYields_C->Draw("p same");

    TLegend *leg = new TLegend(0.15,0.08,0.65,0.42);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(histYields_A,  "Monash2013",  "pl");
    leg->AddEntry(histYields_B2, "Monash2013 new CR", "pl");
    leg->AddEntry(histYields_C,  "Monash2013 no CR",   "pl");
    leg->Draw();

    pad1->SetGridx();



    //    histYields_A->GetXaxis()->SetTitle( "#bar{p_{T}}_{F}" );
    //p/pi vs pt
    //    histYields_A->GetXaxis()->SetTitle( "#LTp_{T}#GT_{F}, GeV/c" );
    //    histYields_A->GetYaxis()->SetTitle( "#LTp/#pi#GT_{B}" );



    // ##### PAD 2
    canvasYields->cd();
    pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.33);
    //    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.02);
    pad2->Draw();
    pad2->cd();

    TH1D *histYields_B2_Copy = new TH1D( *histYields_B2 );
    histYields_B2_Copy->GetXaxis()->SetLabelOffset(0.032);
    histYields_B2_Copy->GetXaxis()->SetLabelSize(0.2);
    histYields_B2_Copy->GetYaxis()->SetLabelSize(0.11);

    histYields_B2_Copy->GetXaxis()->SetTitleSize(0.15);
    histYields_B2_Copy->GetYaxis()->SetTitleSize(0.13);

    histYields_B2_Copy->GetXaxis()->SetTitleOffset( 0.88 );
    histYields_B2_Copy->GetYaxis()->SetTitleOffset( 0.35 );

    histYields_B2_Copy->GetYaxis()->SetTitle("ratio");
    histYields_B2_Copy->GetYaxis()->CenterTitle();


    histYields_B2_Copy->SetMarkerStyle(20);
    //    histYields_B2_Copy->SetMarkerSize(0.6);
    histYields_B2_Copy->SetMarkerColor(kGreen+1);
    histYields_B2_Copy->SetLineColor(kGreen+1);

    TH1D *histYields_C_Copy = new TH1D( *histYields_C );
    histYields_C_Copy->SetMarkerStyle(24);
    //    histYields_C_Copy->SetMarkerSize(0.6);
    histYields_C_Copy->SetMarkerColor(kMagenta);
    histYields_C_Copy->SetLineColor(kMagenta);



    //    TProfile *histDivided = new TProfile ( *prof_B2 );

    //    histDivided->Sumw2();
    //    prof_A->Sumw2();
    //    prof_B2->SetStats(0);
    TH1D *histYields_A_Copy = new TH1D( *histYields_A );
    histYields_B2_Copy->Divide(histYields_A_Copy);
    histYields_C_Copy->Divide(histYields_A_Copy);
    //    h1->SetMarkerStyle(21);

    //set errors for ratios
    for(int i = 0; i < histYields_A_Copy->GetNbinsX(); i++)
    {
        if ( i >= histYields_A_Copy->GetNbinsX()-2 )
        {
            histYields_B2_Copy->SetBinError(i+1, sqrt(2.));
            histYields_C_Copy->SetBinError(i+1, sqrt(2.));
        }
        else
        {
            double binContent = histYields_A_Copy->GetBinContent(i+1);
            double ratioError = sqrt(2.)*sqrt(binContent)/binContent;
            //            cout << "ratio error=" << ratioError << endl;
            histYields_B2_Copy->SetBinError(i+1, ratioError);
            histYields_C_Copy->SetBinError(i+1, ratioError);
        }
    }



    histYields_B2_Copy->Draw("p");
    histYields_C_Copy->Draw("same p");

    pad2->SetGridx();

    TString arrStrPIDnames[] = { "#pi^{#plus}", "#pi^{#minus}", "K", "K", "K_{0}^{S}", "#phi", "p", "#barp", "#bar#Lambda", "#Lambda",
                                 "#Xi^{#minus}", "#bar#Xi^{#plus}", "#Omega^{#minus}", "#Omega^{#plus}" };
    for(int i = 0; i < histYields_B2_Copy->GetNbinsX(); i++)
        histYields_B2_Copy->GetXaxis()->SetBinLabel( i+1, arrStrPIDnames[i] );

    // draw 1-line
    TLine *line = new TLine(-0.5,1,histYields_B2_Copy->GetNbinsX()-0.5,1);
    //    TLine *line = new TLine(1,1,80,2);
    //    line->SetNDC();
    line->SetLineColor(kGray+3);
    line->SetLineStyle(9);
    line->Draw();

}






void getHistEtaY_in_wins()
{
    gStyle->SetOptStat(false);

    TCanvas *canvasYields = new TCanvas("canvasYields","Data Plots",200,10,800,680);
    canvasYields->SetLeftMargin(0.15);
    canvasYields->SetRightMargin(0.02);
    canvasYields->SetTopMargin(0.06);
    canvasYields->SetBottomMargin(0.12);


//    TH1D *histEtaInWins = (TH1D *)fileA->GetObjectChecked( "histEta_in_wins_pi+", "TH1D" );
//    TH1D *histYInWins = (TH1D *)fileA->GetObjectChecked( "histY_in_wins_pi+", "TH1D" );
//    TH1D *histEtaInWins = (TH1D *)fileA->GetObjectChecked( "histEta_in_wins_p+", "TH1D" );
//    TH1D *histYInWins = (TH1D *)fileA->GetObjectChecked( "histY_in_wins_pbar-", "TH1D" );
//    TH1D *histEtaInWins = (TH1D *)fileA->GetObjectChecked( "histEta_in_wins_Lambda0", "TH1D" );
//    TH1D *histYInWins = (TH1D *)fileA->GetObjectChecked( "histY_in_wins_Lambda0", "TH1D" );

//    TH1D *histEtaInWins = (TH1D *)fileA->GetObjectChecked( "histEta_in_wins_K+", "TH1D" );
//    TH1D *histYInWins = (TH1D *)fileA->GetObjectChecked( "histY_in_wins_K-", "TH1D" );

    TH1D *histEtaInWins = (TH1D *)fileA->GetObjectChecked( "histMult_charged", "TH1D" );
    TH1D *histYInWins = (TH1D *)fileB2->GetObjectChecked( "histMult_charged", "TH1D" );
    TH1D *hist3 = (TH1D *)fileC->GetObjectChecked( "histMult_charged", "TH1D" );



//    histEtaInWins->SetLineColor(kRed);
//    histYInWins->SetLineColor(kGreen+1);

    histEtaInWins->SetLineColor(kBlue-3);
    histYInWins->SetLineColor(kGreen+1);
    hist3->SetLineColor(kMagenta);


    histEtaInWins->SetLineWidth(2);
    histYInWins->SetLineWidth(2);
    hist3->SetLineWidth(2);

//    histYInWins->SetFillStyle(3004);
//    histYInWins->SetFillColor(kGreen+1);


    histEtaInWins->GetYaxis()->SetTitle("entries");


    histEtaInWins->GetXaxis()->SetLabelSize(0.05);
    histEtaInWins->GetYaxis()->SetLabelSize(0.05);

    histEtaInWins->GetXaxis()->SetTitleSize(0.058);
    histEtaInWins->GetYaxis()->SetTitleSize(0.065);

    histEtaInWins->GetXaxis()->SetTitleOffset( 0.92 ); //by hand make in large
    histEtaInWins->GetYaxis()->SetTitleOffset( 1.01 );

//    histEtaInWins->GetXaxis()->SetLabelOffset( 0.9 ); //by hand make in large

//    histEtaInWins->DrawNormalized("p");
//    histYInWins->DrawNormalized("p same");
//    histYields_C->DrawNormalized("p same");
    histEtaInWins->Draw();
    histYInWins->Draw("same");

    hist3->Draw("same");


    TLegend *leg = new TLegend(0.8,0.75,0.98,0.95);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
//    leg->AddEntry(histEtaInWins,  "#eta",  "l");
//    leg->AddEntry(histYInWins, "y", "l");
    leg->AddEntry(histEtaInWins,  "Monash2013",  "pl");
    leg->AddEntry(histYInWins, "Monash2013 new CR", "pl");
    leg->AddEntry(hist3,  "Monash2013 no CR",   "pl");
    leg->Draw();




//    histEtaInWins->GetXaxis()->SetRangeUser( -5, 5 );
//    histEtaInWins->GetYaxis()->SetRangeUser( 0, 255000 );

        histEtaInWins->GetXaxis()->SetTitle( "#eta, y" );
    //p/pi vs pt
    //    histEtaInWins->GetXaxis()->SetTitle( "#LTp_{T}#GT_{F}, GeV/c" );
    //    histEtaInWins->GetYaxis()->SetTitle( "#LTp/#pi#GT_{B}" );




}
