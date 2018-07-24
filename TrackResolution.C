
#include <vector>
#include <utility>
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TCanvas.h"

std::pair<double,double> linearRegression(std::vector<double>* x,std::vector<double>* y){
	double * xArray = new double [x->size()];
	double * yArray = new double [x->size()];
	double * x2 = new double [x->size()];
	double * xy = new double [x->size()];
	for(int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
	{
		xArray[vectorIterator]=x->at(vectorIterator);
		yArray[vectorIterator]=y->at(vectorIterator);
		x2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator);
		xy[vectorIterator]=x->at(vectorIterator)*y->at(vectorIterator);

	}

	double xMean=TMath::Mean(x->size(),xArray);
	double yMean=TMath::Mean(x->size(),yArray);
	double xyMean=TMath::Mean(x->size(),xy);
	double x2Mean=TMath::Mean(x->size(),x2);
	double inclination = (xyMean-xMean*yMean)/(x2Mean-xMean*xMean);
	double offset = yMean-inclination*xMean;

	return (std::pair<double,double>(inclination,offset));
}

double calcResidual(double px,double py, std::pair <double,double> track){
	double pcax= px+(-py+track.second+track.first*px)/(track.first*track.first+1)*(-track.first);
	double pcay= py+(-py+track.second+track.first*px)/(track.first*track.first+1);
	return (TMath::Sign(TMath::Sqrt((px-pcax)*(px-pcax)+(py-pcay)*(py-pcay)),py-pcay));
}



//void TrackResolutionAnalysis(){
void TrackResolution(){
  //	TFile *f1  = TFile::Open("~/workspace/tpcsimulation/run_B0_E0.root","Update");
  //	TFile *f1  = TFile::Open("run_0.root","Update");
  TFile *f1  = TFile::Open("run0_BeamSource.root","Update");
  //	TFile *f1  = TFile::Open("run0_BeamSource_mv0.root","Update");
  //  TFile *f1  = TFile::Open("SiDistance_5mm/run0_BeamSource_mv3.5cm.root");
  //  TFile *f1  = TFile::Open("SiDistance_1cm/run0_BeamSource_mv3cm.root");
  //  TFile *f1  = TFile::Open("SiDistance_2cm/run0_BeamSource_mv2cm.root");
  //  TFile *f1  = TFile::Open("SiDistance_3cm/run0_BeamSource_mv1cm.root");
	TTree *t  = (TTree*)f1->Get("tree");

	std::vector<double>* x =0;
	std::vector<double>* y =0;
	std::vector<double>* ysmearing =new std::vector<double>;
	std::vector<double>* refx =new std::vector<double>;
	std::vector<double>* refy =new std::vector<double>;

	TVector3* entryPosition;
	TVector3* exitPosition;
	TVector3* entryPositionGas;
	TVector3* exitPositionGas;
	TVector3* entryPositionMagnet;
	TVector3* exitPositionMagnet;
	double entryEnergy;
	TBranch *brxHit=0;
	TBranch *bryHit=0;
	TBranch *brentryPosition;
	TBranch *brexitPosition;
	TBranch *brentryPositionGas;
	TBranch *brexitPositionGas;
	TBranch *brentryPositionMagnet;
	TBranch *brexitPositionMagnet;
	TBranch *brentryEnergy;

	t->SetBranchAddress("hitXPositions", &x,&brxHit);
	t->SetBranchAddress("hitYPositions", &y,&bryHit);
	t->SetBranchAddress("entryPosition", &entryPosition, &brentryPosition);
	t->SetBranchAddress("exitPosition", &exitPosition, &brexitPosition);
	t->SetBranchAddress("entryPositionGas", &entryPositionGas, &brentryPositionGas);
	t->SetBranchAddress("exitPositionGas", &exitPositionGas, &brexitPositionGas);
	t->SetBranchAddress("entryPositionMagnet", &entryPositionMagnet, &brentryPositionMagnet);
	t->SetBranchAddress("exitPositionMagnet", &exitPositionMagnet, &brexitPositionMagnet);
	t->SetBranchAddress("entryEnergy", &entryEnergy,&brentryEnergy);

	int nEvents=t->GetEntries();

	TRandom3* smearing = new TRandom3();
	TCanvas *c1=new TCanvas("c1","c1",0,0,1200,900);
	c1->Divide(2,2);
	TH1D *residualRefTPC=new TH1D("residualRefTPC","residualRefTPC",1000, -30., 30.);
	TH1D *residualRefGas=new TH1D("residualRefGas","residualRefGas",1000, -30., 30.);
	TH1D *residualRefMagnet=new TH1D("residualRefMagnet","residualRefMagnet",1000, -30., 30.);
	TH1D *residual=new TH1D("residuals","residuals",1000, -30., 30.);
	TCanvas *c2=new TCanvas("c2","c2",0,0,1200,900);
	c2->Divide(2,2);
	TH1D *residualSmearingRefTPC=new TH1D("residualSmearingRefTPC","residualSmearingRefTPC",1000, -30., 30.);
	TH1D *residualSmearingRefGas=new TH1D("residualSmearingRefGas","residualSmearingRefGas",1000, -30., 30.);
	TH1D *residualSmearingRefMagnet=new TH1D("residualSmearingRefMagnet","residualSmearingRefMagnet",1000, -30., 30.);
	TH1D *residualSmearing=new TH1D("residualSmearing","residualSmearing",1000, -30., 30.);
	TCanvas *c3=new TCanvas("c3","c3",0,0,1200,900);
	c3->Divide(2,2);
	TH1D *residualSmearingRefSmearingTPC=new TH1D("residualSmearingRefSmearingTPC","residualSmearingRefSmearingTPC",1000, -30., 30.);
	TH1D *residualSmearingRefSmearingGas=new TH1D("residualSmearingRefSmearingGas","residualSmearingRefSmearingGas",1000, -30., 30.);
	TH1D *residualSmearingRefSmearingMagnet=new TH1D("residualSmearingRefSmearingMagnet","residualSmearingRefSmearingMagnet",1000, -30., 30.);
	nEvents=500;
	for (int iterator=0;iterator<nEvents;iterator++){
		if(iterator%100==0)
		{
			std::cout<<"Event "<<iterator<<" out of "<<nEvents<<std::endl;
		}

		t->GetEvent(iterator);
		ysmearing->clear();
		refx->clear();
		refy->clear();
		refx->push_back(entryPosition->x());
		refy->push_back(entryPosition->y());
		refx->push_back(exitPosition->x());
		refy->push_back(exitPosition->y());

		// this one takes the Si position and makes a fit (least square)
		// then it compares the Si track with every TPC hit in each event
		std::pair<double,double> refTrack1 = linearRegression(refx,refy);
		for(int i=0;i<y->size();i++){
			ysmearing->push_back(y->at(i)+smearing->Gaus(0,0.1));
			residualRefTPC->Fill(calcResidual(x->at(i),y->at(i),refTrack1));
		}
		refx->clear();
		refy->clear();
		refx->push_back(entryPositionGas->x());
		refy->push_back(entryPositionGas->y());
		refx->push_back(exitPositionGas->x());
		refy->push_back(exitPositionGas->y());

		std::pair<double,double> refTrack2 = linearRegression(refx,refy);
		for(int i=0;i<x->size();i++){
			residualRefGas->Fill(calcResidual(x->at(i),y->at(i),refTrack2));
		}
		refx->clear();
		refy->clear();
		refx->push_back(entryPositionMagnet->x());
		refy->push_back(entryPositionMagnet->y());
		refx->push_back(exitPositionMagnet->x());
		refy->push_back(exitPositionMagnet->y());

		std::pair<double,double> refTrack3 = linearRegression(refx,refy);
		for(int i=0;i<x->size();i++){
			residualRefMagnet->Fill(calcResidual(x->at(i),y->at(i),refTrack3));
		}

		// looks at the distribution of TPC hits around the TPC track
		std::pair<double,double> refTrack4 = linearRegression(x,y);
		for(int i=0;i<x->size();i++){
			residual->Fill(calcResidual(x->at(i),y->at(i),refTrack4));
		}

		for(int i=0;i<x->size();i++){
			residualSmearingRefTPC->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack1));
		}

		for(int i=0;i<x->size();i++){
			residualSmearingRefGas->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack2));
		}

		for(int i=0;i<x->size();i++){
			residualSmearingRefMagnet->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack3));
		}

		std::pair<double,double> refTrack5 = linearRegression(x,ysmearing);
		for(int i=0;i<x->size();i++){
			residualSmearing->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack5));
		}

		refx->clear();
		refy->clear();
		refx->push_back(entryPosition->x());
		refy->push_back(entryPosition->y()+smearing->Gaus(0,0.01));
		refx->push_back(exitPosition->x());
		refy->push_back(exitPosition->y()+smearing->Gaus(0,0.01));

		std::pair<double,double> refTrack6 = linearRegression(refx,refy);
		for(int i=0;i<x->size();i++){
			residualSmearingRefSmearingTPC->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack6));
		}
		refx->clear();
		refy->clear();
		refx->push_back(entryPositionGas->x());
		refy->push_back(entryPositionGas->y()+smearing->Gaus(0,0.01));
		refx->push_back(exitPositionGas->x());
		refy->push_back(exitPositionGas->y()+smearing->Gaus(0,0.01));

		std::pair<double,double> refTrack7 = linearRegression(refx,refy);
		for(int i=0;i<x->size();i++){
			residualSmearingRefSmearingGas->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack7));
		}
		refx->clear();
		refy->clear();
		refx->push_back(entryPositionMagnet->x());
		refy->push_back(entryPositionMagnet->y()+smearing->Gaus(0,0.01));
		refx->push_back(exitPositionMagnet->x());
		refy->push_back(exitPositionMagnet->y()+smearing->Gaus(0,0.01));

		std::pair<double,double> refTrack8 = linearRegression(refx,refy);
		for(int i=0;i<x->size();i++){
			residualSmearingRefSmearingMagnet->Fill(calcResidual(x->at(i),ysmearing->at(i),refTrack8));
		}


	}
	c1->cd(1);
	residualRefTPC->Draw();
	residualRefTPC->Fit("gaus");
	c1->cd(2);
	residualRefGas->Draw();
	residualRefGas->Fit("gaus");
	c1->cd(3);
	residualRefMagnet->Draw();
	residualRefMagnet->Fit("gaus");
	c1->cd(4);
	residual->Draw();
	residual->Fit("gaus");

	c2->cd(1);
	residualSmearingRefTPC->Draw();
	residualSmearingRefTPC->Fit("gaus");
	c2->cd(2);
	residualSmearingRefGas->Draw();
	residualSmearingRefGas->Fit("gaus");
	c2->cd(3);
	residualSmearingRefMagnet->Draw();
	residualSmearingRefMagnet->Fit("gaus");
	c2->cd(4);
	residualSmearing->Draw();
	residualSmearing->Fit("gaus");

	c3->cd(1);
	residualSmearingRefSmearingTPC->Draw();
	residualSmearingRefSmearingTPC->Fit("gaus");
	c3->cd(2);
	residualSmearingRefSmearingGas->Draw();
	residualSmearingRefSmearingGas->Fit("gaus");
	c3->cd(3);
	residualSmearingRefSmearingMagnet->Draw();
	residualSmearingRefSmearingMagnet->Fit("gaus");
	c3->cd(4);
	residualSmearing->Draw();
	residualSmearing->Fit("gaus");

	TFile* outfile = TFile::Open("TrackResolution.root","recreate");
	c1->Write();
	c2->Write();
	c3->Write();

	outfile->Close();
}
