#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "iostream"
#include "vector"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <math.h>  
#include "TLegend.h"
using namespace std;

#define SIZE_Res 9
#define SIZE_hist 12

//int Mom[SIZE_M] = {1, 3, 5, 10, 15, 25, 50, 100};
int Res[SIZE_Res] = {2.5, 5, 10, 15, 20, 25, 50, 100, 200};
float LimAxis;
int color, marker;

char a1[][50] = { 
  "  Measurement TPC only (with smearing)         ",
  //  "MomentumResolutions no Si layer              ",
  " Measurement only 2nd Si layer                ",
  " Measurement only 3rd Si layer                ",
  " Measurement only 4th Si layer                ",
  " Measurement only 1st Si layer                ",
  " Measurement only 2nd and 3rd  Si layer       ",
  " Measurement only 1st and 4th Si layer        ",
  " Measurement only 1st, 2nd and 3rd Si layer   ",
  " Measurement only 2nd, 3rd and 4th Si layer   ",
  " Measurement all 4 Si layers                  ",
  " Measurement only 4 Si layers - no TPC        ",
  " Measurement only 2nd and 3rd layer Perfect Si"
};


void Compare_diff_reso(){
  
  
  //  TGraphErrors* angle_res[SIZE_angle][SIZE_res] = new TGraphErrors[SIZE_angle][SIZE_res]();
  TH1D* single_angle[SIZE_hist][SIZE_Res];
  
  for ( int j=1; j<SIZE_Res; j++ ) {
    

    TFile *f = new TFile(Form("PResolution_%dum.root",Res[j]), "read");
    //    TFile *f = new TFile(Form("PResolution_"+Res[j]+"um.root", "read");
    
    
    // 85 degrees
    single_angle[0][j] = (TH1D*)f->FindObjectAny("MomentumResolutionSmearing;1");
    //    single_angle[0][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;1");

    single_angle[1][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;1");
    single_angle[2][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;2");
    single_angle[3][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;3");
    single_angle[4][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;4");
    single_angle[5][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;5");
    single_angle[6][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;6");
    single_angle[7][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;7");
    single_angle[8][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;8");
    single_angle[9][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;9");
    single_angle[10][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;10");
    single_angle[11][j] = (TH1D*)f->FindObjectAny("MomentumResolutionReference;11");

    std::cout<<f->GetTitle()<<Res[j]<<" um"<<std::endl;
    for ( int i=0; i<SIZE_hist; i++ ) {
      if ( i==0 ) TF1 *func = (TF1*)single_angle[i][j]->GetFunction("CrystalBallSmear");
      else  TF1 *func = (TF1*)single_angle[i][j]->GetFunction("CrystalBallReso");
      Int_t npar = func->GetNpar();
      std::cout<<single_angle[i][j]->GetTitle()<<a1[i]<<" Fit sigma parameter "<<func->GetParName(3)<<
	" "<< func->GetParameter(3)<<" +/- "<< func->GetParError(3)<<std::endl;
    }
    
  }
}
/*
    for ( int i=0; i<SIZE_angle; i++ ) {
      
      single_angle[i][j]->SetTitle("Momentum Resolution");
      single_angle[i][j]->SetMarkerColor(1+j);
      if ( j>3 )  single_angle[i][j]->SetMarkerColor(2+j);
      single_angle[i][j]->SetMarkerStyle(20+j);
      if ( j==4 ) single_angle[i][j]->SetMarkerStyle(29);
      if ( j==5 ) single_angle[i][j]->SetMarkerStyle(33);
      single_angle[i][j]->SetMarkerSize(1);
      single_angle[i][j]->GetXaxis() -> SetTitle("P (GeV)");
      single_angle[i][j]->GetYaxis() -> SetTitle("#sigma_{1/p_{T}}(GeV^{-1})");
      //    single_angle[i][j].SetMinimum( pow(10,-5) );
      //    single_angle[i][j].SetMaximum( 2*pow(10, -1) );
      
    }
     
  }
  
  
  TFile* outfile_plots = TFile::Open("Comparison_DiffReso.root","recreate");
  
  TLegend *leg = new TLegend(0.55,0.6,0.85,0.85);
  leg->SetHeader("Sensor Spatial resolution"); //name of the legend
  leg->SetFillColor(kWhite);
  leg->AddEntry(single_angle[0][0],"5 #mum","ep");
  leg->AddEntry(single_angle[0][1],"7 #mum","ep");
  leg->AddEntry(single_angle[0][2],"10 #mum","ep");
  leg->AddEntry(single_angle[0][3],"15 #mum","ep");
  leg->AddEntry(single_angle[0][4],"20 #mum","ep");
  leg->AddEntry(single_angle[0][5],"200 #mum","ep");
  
  

  TCanvas c1;
  single_angle[0][0]->SetTitle("Momentum Resolution at 85 degrees");
  single_angle[0][0]->SetMinimum(0.00001);
  single_angle[0][0]->SetMaximum(0.002);
  single_angle[0][0]->Draw("AP");
  single_angle[0][1]->Draw("P");
  single_angle[0][2]->Draw("P");
  single_angle[0][3]->Draw("P");
  single_angle[0][4]->Draw("P");
  single_angle[0][5]->Draw("P");
  leg->Draw();


  TCanvas c2;
  single_angle[1][0]->SetTitle("Momentum Resolution at 40 degrees");
  single_angle[1][0]->SetMinimum(0.00002);
  single_angle[1][0]->SetMaximum(0.005);
  single_angle[1][0]->Draw("AP");
  single_angle[1][1]->Draw("P");
  single_angle[1][2]->Draw("P");
  single_angle[1][3]->Draw("P");
  single_angle[1][4]->Draw("P");
  single_angle[1][5]->Draw("P");
  leg->Draw();

  TCanvas c3;
  single_angle[2][0]->SetTitle("Momentum Resolution at 20 degrees");
  single_angle[2][0]->SetMinimum(0.0001);
  single_angle[2][0]->SetMaximum(0.05);
  single_angle[2][0]->Draw("AP");
  single_angle[2][1]->Draw("P");
  single_angle[2][2]->Draw("P");
  single_angle[2][3]->Draw("P");
  single_angle[2][4]->Draw("P");
  single_angle[2][5]->Draw("P");
  leg->Draw();

  TCanvas c4;
  single_angle[3][0]->SetTitle("Momentum Resolution at 10 degrees");
  single_angle[3][0]->SetMinimum(0.001);
  single_angle[3][0]->SetMaximum(0.2);
  single_angle[3][0]->Draw("AP");
  single_angle[3][1]->Draw("P");
  single_angle[3][2]->Draw("P");
  single_angle[3][3]->Draw("P");
  single_angle[3][4]->Draw("P");
  single_angle[3][5]->Draw("P");
  leg->Draw();
  
  c1.SetLogx();
  c1.SetLogy();
  c1.SetTickx(1);
  c1.SetTicky(1);
  
  c2.SetLogx();
  c2.SetLogy();
  c2.SetTickx(1);
  c2.SetTicky(1);
  
  c3.SetLogx();
  c3.SetLogy();
  c3.SetTickx(1);
  c3.SetTicky(1);
  
  c4.SetLogx();
  c4.SetLogy();
  c4.SetTickx(1);
  c4.SetTicky(1);
  
  
  
  c1.Print("angle85.png");
  c2.Print("angle40.png");
  c3.Print("angle20.png");
  c4.Print("angle10.png");
  
  c1.Write();
  c2.Write();
  c3.Write();
  c4.Write();
  outfile_plots->Close();

}
*/
