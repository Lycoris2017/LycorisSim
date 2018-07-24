
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TStyle.h"

//function to calculate the momentum of a track. you just need to insert the x and y positions as std::vectors in units of mm to get the momentum in MeV
double calculateMomentum(std::vector<double>* x,std::vector<double>* y)
{
	double * xArray = new double [x->size()];
	double * yArray = new double [x->size()];
	double * x2 = new double [x->size()];
	double * y2 = new double [x->size()];
	double * xy = new double [x->size()];
	double * xr2 = new double [x->size()];
	double * yr2 = new double [x->size()];
	double * r2 = new double [x->size()];
	double * r4 = new double [x->size()];
	for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
	{
		xArray[vectorIterator]=x->at(vectorIterator);
		yArray[vectorIterator]=y->at(vectorIterator);
		x2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator);
		y2[vectorIterator]=y->at(vectorIterator)*y->at(vectorIterator);
		xy[vectorIterator]=x->at(vectorIterator)*y->at(vectorIterator);
		r2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator)+y->at(vectorIterator)*y->at(vectorIterator);
		xr2[vectorIterator]=x->at(vectorIterator)*r2[vectorIterator];
		yr2[vectorIterator]=y->at(vectorIterator)*r2[vectorIterator];
		r4[vectorIterator]=r2[vectorIterator]*r2[vectorIterator];
	}
	double Cxx,Cxy,Cyy,Cxr2,Cyr2,Cr2r2,phi,kappa,delta,rho,momentum;
	double xMean,yMean,x2Mean,y2Mean,xyMean,xr2Mean,yr2Mean,r2Mean,r4Mean;
	xMean=TMath::Mean(x->size(),xArray);

	x2Mean=TMath::Mean(x->size(),x2);
	yMean=TMath::Mean(x->size(),yArray);
	y2Mean=TMath::Mean(x->size(),y2);
	xyMean=TMath::Mean(x->size(),xy);
	xr2Mean=TMath::Mean(x->size(),xr2);
	yr2Mean=TMath::Mean(x->size(),yr2);
	r2Mean=TMath::Mean(x->size(),r2);
	r4Mean=TMath::Mean(x->size(),r4);
	Cxx=x2Mean-(xMean*xMean);
	Cxy=xyMean-xMean*yMean;
	Cyy=y2Mean-yMean*yMean;
	Cxr2=xr2Mean-xMean*r2Mean;
	Cyr2=yr2Mean-yMean*r2Mean;
	Cr2r2=r4Mean-r2Mean*r2Mean;

	phi=0.5*TMath::ATan(2*(Cr2r2*Cxy-Cxr2*Cyr2)/(Cr2r2*(Cxx-Cyy)-Cxr2*Cxr2+Cyr2*Cyr2));
	kappa=(TMath::Sin(phi)*Cxr2-TMath::Cos(phi)*Cyr2)/Cr2r2;
	delta=-kappa*r2Mean+TMath::Sin(phi)*xMean-TMath::Cos(phi)*yMean;

	rho=2*kappa/(TMath::Sqrt(1-4*delta*kappa));
	momentum=0.3*1/rho;
	delete [] xArray;
	delete [] yArray;
	delete [] x2;
	delete [] y2;
	delete [] xy;
	delete [] xr2;
	delete [] yr2;
	delete [] r2;
	delete [] r4;
	return (momentum);
};

//function to calculate the momentum of a track with weighted points. You just need to insert the x and y positions and x and y error as std::vectors in units of mm to get the momentum in MeV
double calculateWeightedMomentum(std::vector<double>* x,std::vector<double>* y,std::vector<double>* xerr,std::vector<double>* yerr)
{
	double * xArray = new double [x->size()];
	double * yArray = new double [x->size()];
	double * x2 = new double [x->size()];
	double * y2 = new double [x->size()];
	double * xy = new double [x->size()];
	double * r2 = new double [x->size()];
	double * xr2 = new double [x->size()];
	double * yr2 = new double [x->size()];
	double * r4 = new double [x->size()];
	double * weight = new double [x->size()];
	for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
	{
		//std::cout<<"test0"<<std::endl;
		weight[vectorIterator]=1/(yerr->at(vectorIterator)*yerr->at(vectorIterator));
		xArray[vectorIterator]=x->at(vectorIterator);
		yArray[vectorIterator]=y->at(vectorIterator);
		x2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator);
		y2[vectorIterator]=y->at(vectorIterator)*y->at(vectorIterator);
		xy[vectorIterator]=x->at(vectorIterator)*y->at(vectorIterator);
		r2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator)+y->at(vectorIterator)*y->at(vectorIterator);
		xr2[vectorIterator]=x->at(vectorIterator)*r2[vectorIterator];
		yr2[vectorIterator]=y->at(vectorIterator)*r2[vectorIterator];
		r4[vectorIterator]=r2[vectorIterator]*r2[vectorIterator];
	}
	double Cxx,Cxy,Cyy,Cxr2,Cyr2,Cr2r2,phi,kappa,delta,rho,momentum;
	double xMean,yMean,x2Mean,y2Mean,xyMean,xr2Mean,yr2Mean,r2Mean,r4Mean;
	xMean=TMath::Mean(x->size(),xArray,weight);
	x2Mean=TMath::Mean(x->size(),x2,weight);
	yMean=TMath::Mean(x->size(),yArray,weight);
	y2Mean=TMath::Mean(x->size(),y2,weight);
	xyMean=TMath::Mean(x->size(),xy,weight);
	xr2Mean=TMath::Mean(x->size(),xr2,weight);
	yr2Mean=TMath::Mean(x->size(),yr2,weight);
	r2Mean=TMath::Mean(x->size(),r2,weight);
	r4Mean=TMath::Mean(x->size(),r4,weight);

	Cxx=x2Mean-(xMean*xMean);
	Cxy=xyMean-xMean*yMean;
	Cyy=y2Mean-yMean*yMean;
	Cxr2=xr2Mean-xMean*r2Mean;
	Cyr2=yr2Mean-yMean*r2Mean;
	Cr2r2=r4Mean-r2Mean*r2Mean;

	phi=0.5*TMath::ATan(2*(Cr2r2*Cxy-Cxr2*Cyr2)/(Cr2r2*(Cxx-Cyy)-Cxr2*Cxr2+Cyr2*Cyr2));
	kappa=(TMath::Sin(phi)*Cxr2-TMath::Cos(phi)*Cyr2)/Cr2r2;
	delta=-kappa*r2Mean+TMath::Sin(phi)*xMean-TMath::Cos(phi)*yMean;

	rho=2*kappa/(TMath::Sqrt(1-4*delta*kappa));
	momentum=0.3*1/rho;
	delete [] xArray;
	delete [] weight;
	delete [] yArray;
	delete [] x2;
	delete [] y2;
	delete [] xy;
	delete [] xr2;
	delete [] yr2;
	delete [] r2;
	delete [] r4;

	return (momentum);
};

double CrystalBall(double* x, double* par){
double xcur = x[0];
double alpha = par[0];
double n = par[1];
double mu = par[2];
double sigma = par[3];
double N = par[4];
TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);
double A; double B;
if (alpha < 0){
A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2);
B = n/(-1*alpha) + alpha;}
else {
A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2);
B = n/alpha - alpha;}
double f;
if ((xcur-mu)/sigma > (-1)*alpha)
f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/
(2*sigma*sigma));
else
f = N*A*pow((B- (xcur-mu)/sigma),(-1*n));
delete exp;
return (f);
};

double si_reso = 0.01;

//void MomentumResolutionAnalysis(){
void MomentumResolution_MyPlots(){
  //	TFile *f1  = TFile::Open("run_0.root","Update");
  
  gStyle->SetOptFit(1);

  //  TFile *f1  = TFile::Open("run0_BeamSource.root");
  TFile *f1  = TFile::Open("run0_BeamSource_mv0.root");
  //  TFile *f1  = TFile::Open("run0_BeamSourceNoMSc.root");
  //  TFile *f1  = TFile::Open("run0_PointSourceNoMSc.root");
  //  TFile *f1  = TFile::Open("run0_PointSource.root" );
  
	TTree *t  = (TTree*)f1->Get("tree");
	double momentum;
	TBranch *bMomentum=t->Branch("momentum", &momentum, "momentum/D");
	std::vector<double>* x =0;
	std::vector<double>* xerr = new std::vector<double>;
	std::vector<double>* y =0;
	std::vector<double>* yerr = new std::vector<double>;
	//	std::vector<double>* refx1 =new std::vector<double>;
	//	std::vector<double>* refy1 =new std::vector<double>;

	TVector3* entryPosition=0;
	TVector3* exitPosition=0;
	TVector3* entryPositionGas=0;
	TVector3* exitPositionGas=0;
	TVector3* entryPositionMagnet=0;
	TVector3* exitPositionMagnet=0;
	double entryEnergy,entryEnergyGas,entryEnergyMagnet;
	TBranch *brxHit=0;
	TBranch *bryHit=0;
	TBranch *brentryPosition=0;
	TBranch *brexitPosition=0;
	TBranch *brentryPositionGas=0;
	TBranch *brexitPositionGas=0;
	TBranch *brentryPositionMagnet=0;
	TBranch *brexitPositionMagnet=0;
	TBranch *brentryEnergy=0;
	TBranch *brentryEnergyGas=0;
	TBranch *brentryEnergyMagnet=0;

	t->SetBranchAddress("hitXPositions", &x,&brxHit);
	t->SetBranchAddress("hitYPositions", &y,&bryHit);
	t->SetBranchAddress("entryPosition", &entryPosition, &brentryPosition);
	t->SetBranchAddress("exitPosition", &exitPosition, &brexitPosition);
	t->SetBranchAddress("entryPositionGas", &entryPositionGas, &brentryPositionGas);
	t->SetBranchAddress("exitPositionGas", &exitPositionGas, &brexitPositionGas);
	t->SetBranchAddress("entryPositionMagnet", &entryPositionMagnet, &brentryPositionMagnet);
	t->SetBranchAddress("exitPositionMagnet", &exitPositionMagnet, &brexitPositionMagnet);
	t->SetBranchAddress("entryEnergy", &entryEnergy,&brentryEnergy);
	t->SetBranchAddress("entryEnergyGas", &entryEnergyGas,&brentryEnergyGas);
	t->SetBranchAddress("entryEnergyMagnet", &entryEnergyMagnet,&brentryEnergyMagnet);


	int nEvents=t->GetEntries();
	TCanvas *c1=new TCanvas("MomentumResolutions no Si layer","Momentum Resolutions",0,0,1200,900);
	c1->Divide(2,2);
	TH1D *MomentumResolution=new TH1D("MomentumResolution","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionSmearing=new TH1D("MomentumResolutionSmearing","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionSmearing->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearing->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionSmearingEnergy=new TH1D("MomentumResolutionSmearingEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionSmearingEnergy->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearingEnergy->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionSmearing2=new TH1D("MomentumResolutionSmearing2","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionSmearing2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearing2->GetYaxis()->SetTitle("number of events");

	//	TH1D *MomentumResolutionRef=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	//	MomentumResolutionRef->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	//	MomentumResolutionRef->GetYaxis()->SetTitle("number of events");

	TCanvas *c2=new TCanvas("Measurement only 2nd Si layer","Measurement only 2nd Si layer",0,0,1200,900);
	c2->Divide(2,2);
	TH1D *MomentumResolution2=new TH1D("MomentumResolution2","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution2->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight2=new TH1D("MomentumResolutionWeight2","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight2->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy2=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy2->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef2=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef2->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE2=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE2->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE2=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE2->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE2->GetYaxis()->SetTitle("number of events");



	TCanvas *c3=new TCanvas("Measurement only 3rd Si layer","Measurement only 3rd Si layer",0,0,1200,900);
	c3->Divide(2,2);
	TH1D *MomentumResolution3=new TH1D("MomentumResolution3","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution3->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution3->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight3=new TH1D("MomentumResolutionWeight3","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight3->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight3->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy3=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy3->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy3->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef3=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef3->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef3->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE3=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE3->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE3->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE3=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE3->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE3->GetYaxis()->SetTitle("number of events");

	TCanvas *c4=new TCanvas("Measurement only 4th Si layer","Measurement only 4th Si layer",0,0,1200,900);
	c4->Divide(2,2);
	TH1D *MomentumResolution4=new TH1D("MomentumResolution4","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution4->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution4->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight4=new TH1D("MomentumResolutionWeight4","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight4->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight4->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy4=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy4->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy4->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef4=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef4->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef4->GetYaxis()->SetTitle("number of events");


	TH1D *MomentumResolutionSmearE4=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE4->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE4->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE4=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE4->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE4->GetYaxis()->SetTitle("number of events");

	TCanvas *c5=new TCanvas("Measurement only 1st Si layer","Measurement only 1st Si layer",0,0,1200,900);
	c5->Divide(2,2);
	TH1D *MomentumResolution5=new TH1D("MomentumResolution5","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution5->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution5->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight5=new TH1D("MomentumResolutionWeight5","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight5->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight5->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy5=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy5->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy5->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef5=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef5->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef5->GetYaxis()->SetTitle("number of events");


	TH1D *MomentumResolutionSmearE5=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE5->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE5->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE5=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE5->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE5->GetYaxis()->SetTitle("number of events");
 
	TCanvas *c6=new TCanvas("Measurement only 2nd and 3rd  Si layer","Measurement only 2nd and 3rd  Si layer",0,0,1200,900);
	c6->Divide(2,2);
	TH1D *MomentumResolution6=new TH1D("MomentumResolution6","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution6->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution6->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight6=new TH1D("MomentumResolutionWeight6","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight6->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight6->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy6=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy6->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy6->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef6=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef6->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef6->GetYaxis()->SetTitle("number of events");


	TH1D *MomentumResolutionSmearE6=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE6->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE6->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE6=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE6->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE6->GetYaxis()->SetTitle("number of events");


	TCanvas *c7=new TCanvas("Measurement only 1st and 4th Si layer","Measurement only 1st and 4th Si layer",0,0,1200,900);
	c7->Divide(2,2);
	TH1D *MomentumResolution7=new TH1D("MomentumResolution7","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution7->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution7->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight7=new TH1D("MomentumResolutionWeight7","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight7->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight7->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy7=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy7->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy7->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef7=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef7->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef7->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE7=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE7->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE7->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE7=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE7->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE7->GetYaxis()->SetTitle("number of events");

	TCanvas *c8=new TCanvas("Measurement only 1st, 2nd and 3rd Si layer","Measurement only 1st, 2nd and 3rd Si layer",0,0,1200,900);
	c8->Divide(2,2);
	TH1D *MomentumResolution8=new TH1D("MomentumResolution8","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution8->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution8->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight8=new TH1D("MomentumResolutionWeight8","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight8->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight8->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy8=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy8->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy8->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef8=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef8->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef8->GetYaxis()->SetTitle("number of events");


	TH1D *MomentumResolutionSmearE8=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE8->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE8->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE8=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE8->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE8->GetYaxis()->SetTitle("number of events");


	TCanvas *c9=new TCanvas("Measurement only 2nd, 3rd and 4th Si layer","Measurement only 2nd, 3rd and 4th Si layer",0,0,1200,900);
	c9->Divide(2,2);
	TH1D *MomentumResolution9=new TH1D("MomentumResolution9","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution9->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution9->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight9=new TH1D("MomentumResolutionWeight9","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight9->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight9->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy9=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy9->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy9->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef9=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef9->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef9->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE9=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE9->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE9->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE9=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE9->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE9->GetYaxis()->SetTitle("number of events");

	TCanvas *c10=new TCanvas("Measurement all 4 Si layers","Measurement all 4 Si layers",0,0,1200,900);
	c10->Divide(2,2);
	TH1D *MomentumResolution10=new TH1D("MomentumResolution10","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution10->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution10->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight10=new TH1D("MomentumResolutionWeight10","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight10->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight10->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy10=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy10->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy10->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef10=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef10->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef10->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE10=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE10->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE10->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE10=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE10->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE10->GetYaxis()->SetTitle("number of events");

	TCanvas *c11=new TCanvas("simulatedMomentum","simulated momentum",0,0,1200,900);
	c11->Divide(2,2);
	TH1D *tpcMomentum=new TH1D("Momentum","Simulated TPC Momentum",512, 0, 6000);
	tpcMomentum->GetXaxis()->SetTitle("p_{t,sim}/MeV");
	tpcMomentum->GetXaxis()->SetNdivisions(505);
	tpcMomentum->GetYaxis()->SetTitle("number of events");
	TH1D *magnetMomentum=new TH1D("magnetMomentum","Simulated Momentum after Magnet",512, 0, 6000);
	magnetMomentum->GetXaxis()->SetTitle("p_{t,sim}/MeV");
	magnetMomentum->GetXaxis()->SetNdivisions(505);
	magnetMomentum->GetYaxis()->SetTitle("number of events");
	TH1D *beamMomentum=new TH1D("beamMomentum","Beam Momentum",512, 0, 6000);
	beamMomentum->GetXaxis()->SetTitle("p_{t,sim}/MeV");
	beamMomentum->GetXaxis()->SetNdivisions(505);
	beamMomentum->GetYaxis()->SetTitle("number of events");

	TCanvas *c12=new TCanvas("reconstructedMomentum","reconstructed momentum",0,0,1200,900);
	c12->Divide(2,2);
	TH1D *recoMomentum=new TH1D("recoMomentum","Reconstructed Momentum",512, 0, 6000);
	recoMomentum->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentum->GetXaxis()->SetNdivisions(505);
	recoMomentum->GetYaxis()->SetTitle("number of events");
	TH1D *recoMomentumSmearing=new TH1D("recoMomentumSmearing","Reconstructed Momentum with Smearing",512, 0, 6000);
	recoMomentumSmearing->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumSmearing->GetXaxis()->SetNdivisions(505);
	recoMomentumSmearing->GetYaxis()->SetTitle("number of events");
	TH1D *recoMomentumRef=new TH1D("recoMomentumRef","Reconstructed Momentum with Reference",512, 0, 6000);
	recoMomentumRef->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumRef->GetXaxis()->SetNdivisions(505);
	recoMomentumRef->GetYaxis()->SetTitle("number of events");
	TH1D *recoMomentumWeight=new TH1D("recoMomentumWeight","Reconstructed Momentum with Weighted Track Points",512, 0, 6000);
	recoMomentumWeight->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumWeight->GetXaxis()->SetNdivisions(505);
	recoMomentumWeight->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE12=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE12->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE12->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE12=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE12->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE12->GetYaxis()->SetTitle("number of events");

	TCanvas *c13=new TCanvas("Measurement","Measurement",0,0,1200,900);
	c13->Divide(2,2);
	TH1D *MomentumResolutionWeight=new TH1D("MomentumResolutionWeight","Momentum Resolution with Weighted Track Points", 512, -0.00006, 0.00006);
	MomentumResolutionWeight->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy=new TH1D("MomentumResolutionWeightEenrgy","Momentum Resolution with Weighted Track Points and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE13=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE13->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE13->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE13=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE13->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE13->GetYaxis()->SetTitle("number of events");

	TCanvas *c14=new TCanvas("Measurement only 4 Si layers - no TPC","Measurement only 4 Si layers - no TPC",0,0,1200,900);
	c14->Divide(2,2);
	TH1D *MomentumResolution14=new TH1D("MomentumResolution14","Momentum Resolution without smearing", 512, -0.0003, 0.0003);
	MomentumResolution14->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution14->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight14=new TH1D("MomentumResolutionWeight14","Momentum Resolution with Smearing", 512, -0.0003, 0.0003);
	MomentumResolutionWeight14->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight14->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy14=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0003, 0.0003);
	MomentumResolutionWeightEnergy14->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy14->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef14=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.0003, 0.0003);
	MomentumResolutionRef14->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef14->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE14=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE14->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE14->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE14=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE14->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE14->GetYaxis()->SetTitle("number of events");

	TCanvas *c15=new TCanvas("Measurement only 2nd and 3rd layer Perfect Si", "Measurement only 2nd and 3rd layer Perfect Si",0,0,1200,900);
	c15->Divide(2,2);
	TH1D *MomentumResolution15=new TH1D("MomentumResolution15","Momentum Resolution without smearing", 512, -0.00006, 0.00006);
	MomentumResolution15->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution15->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight15=new TH1D("MomentumResolutionWeight15","Momentum Resolution with Smearing", 512, -0.00006, 0.00006);
	MomentumResolutionWeight15->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight15->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy15=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy15->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy15->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef15=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef15->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef15->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE15=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE15->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE15->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE15=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE15->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE15->GetYaxis()->SetTitle("number of events");


	TCanvas *c16=new TCanvas("Measurement only 2 outer Si layers - no TPC","Measurement only 2 outer layers - no TPC",0,0,1200,900);
	c16->Divide(2,2);
	TH1D *MomentumResolution16=new TH1D("MomentumResolution16","Momentum Resolution without smearing", 512, -0.0003, 0.0003);
	MomentumResolution16->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution16->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight16=new TH1D("MomentumResolutionWeight16","Momentum Resolution with Smearing", 512, -0.0003, 0.0003);
	MomentumResolutionWeight16->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight16->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy16=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0003, 0.0003);
	MomentumResolutionWeightEnergy16->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy16->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef16=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.0003, 0.0003);
	MomentumResolutionRef16->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef16->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE16=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE16->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE16->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE16=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE16->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE16->GetYaxis()->SetTitle("number of events");

	TCanvas *c17=new TCanvas("Measurement only 2 inner Si layers - no TPC","Measurement only 2 inner Si layers - no TPC",0,0,1200,900);
	c17->Divide(2,2);
	TH1D *MomentumResolution17=new TH1D("MomentumResolution17","Momentum Resolution without smearing", 512, -0.0003, 0.0003);
	MomentumResolution17->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolution17->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeight17=new TH1D("MomentumResolutionWeight17","Momentum Resolution with Smearing", 512, -0.0003, 0.0003);
	MomentumResolutionWeight17->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight17->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy17=new TH1D("MomentumResolutionWeightEnergy","Momentum Resolution with Smearing and unknown Energy", 512, -0.0003, 0.0003);
	MomentumResolutionWeightEnergy17->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy17->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionRef17=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.0003, 0.0003);
	MomentumResolutionRef17->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef17->GetYaxis()->SetTitle("number of events");

	TH1D *MomentumResolutionSmearE17=new TH1D("MomentumResolution","Momentum Resolution", 512, -0.0001, 0.00006);
	MomentumResolutionSmearE17->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionSmearE17->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightE17=new TH1D("MomentumResolution","Momentum Resolution with reference", 512, -0.0001, 0.00006);
	MomentumResolutionWeightE17->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightE17->GetYaxis()->SetTitle("number of events");

	TRandom2* smearing = new TRandom2();
	//iterate over tree
//	nEvents=10000;
	for (int iterator=0;iterator<nEvents;iterator++){
		if(iterator%5000==0)
		{
			std::cout<<"Event "<<iterator<<" out of "<<nEvents<<std::endl;
		}
		t->GetEvent(iterator);
		if (x->size()<49)
			continue;

		if ( iterator<2 ) {
		  for (uint i=0; i<x->size(); i++ ) {
		    std::cout<<" x "<<x->at(i)<<" y "<<y->at(i)<<std::endl;
		  }
		}

		// why are you doing this? Are the first and last 24 outside the detector?
		// A: we don't have GEMs in all the gas volume
		x->erase(x->begin(),x->begin()+24);
		y->erase(y->begin(),y->begin()+24);
		x->erase(x->end()-24,x->end());
		y->erase(y->end()-24,y->end());
		beamMomentum->Fill(entryEnergyMagnet);
		magnetMomentum->Fill(entryEnergy);
		tpcMomentum->Fill(entryEnergyGas);

		if ( iterator<2 ) {
		  std::cout<<"After"<<std::endl;
		  for (uint i=0; i<x->size(); i++ ) {
		    std::cout<<" x "<<x->at(i)<<" y "<<y->at(i)<<std::endl;
		  }
		}

		momentum=calculateMomentum(x,y);
		recoMomentum->Fill(momentum);
		MomentumResolution->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));

		if(iterator<100) std::cout<<"x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		xerr->resize(y->size());
		yerr->resize(y->size());
		for (uint i=0;i<y->size();i++)
		{
			y->at(i)=y->at(i)+smearing->Gaus(0,0.1);
			xerr->at(i)=1;
			yerr->at(i)=0.1;
		}
		/*refx1->push_back(entryPosition->x());
		refy1->push_back(entryPosition->y());
		refx1->push_back(exitPosition->x());
		refy1->push_back(exitPosition->y());
		refx1->push_back(entryPositionMagnet->x());
		refy1->push_back(entryPositionMagnet->y());
		refx1->push_back(exitPositionMagnet->x());
		refy1->push_back(exitPositionMagnet->y());
		 */
		momentum=calculateMomentum(x,y);
		recoMomentumSmearing->Fill(momentum);
		MomentumResolutionSmearing->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearingEnergy->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionSmearing2->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));

		if(iterator<100) std::cout<<"Smearing: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;



		// here add only 2nd layer
		x->push_back(entryPosition->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef2->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE2->Fill((momentum-5000)/(5000*5000));
		//		y->pop_back(); y->push_back(entryPosition->y()+smearing->Gaus(0,0.01));
		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		//momentum=calculateMomentum(x,y);
		MomentumResolutionWeight2->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionWeightE2->Fill((momentum-5000)/(5000*5000));

		// here add only 3rd layer
		x->pop_back(); xerr->pop_back(); 
		y->pop_back(); yerr->pop_back(); 

		x->push_back(exitPosition->x());
		xerr->push_back(1); //default 1
		y->push_back(exitPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef3->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE3->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight3->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy3->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE3->Fill((momentum-5000)/(5000*5000));

		// here add only fourth layer
		x->pop_back(); xerr->pop_back(); 
		y->pop_back(); yerr->pop_back(); 

		x->push_back(exitPositionGas->x());
		xerr->push_back(1); //default 1
		y->push_back(exitPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef4->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE4->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight4->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy4->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE4->Fill((momentum-5000)/(5000*5000));

		// here add only first layer
		x->pop_back(); xerr->pop_back(); 
		y->pop_back(); yerr->pop_back(); 

		x->push_back(entryPositionGas->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef5->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE5->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight5->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy5->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE5->Fill((momentum-5000)/(5000*5000));

		// here add 2nd and 3rd
		x->pop_back(); xerr->pop_back();
		y->pop_back(); yerr->pop_back(); 

		x->push_back(entryPosition->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPosition->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01
		y->push_back(exitPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef6->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE6->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight6->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy6->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE6->Fill((momentum-5000)/(5000*5000));

		// here add 1st and 4th
		x->pop_back(); xerr->pop_back();
		x->pop_back(); xerr->pop_back();
		y->pop_back(); yerr->pop_back(); 
		y->pop_back(); yerr->pop_back(); 

		x->push_back(entryPositionGas->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPositionGas->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01
		y->push_back(exitPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef7->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE7->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight7->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy7->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE7->Fill((momentum-5000)/(5000*5000));

		// here add 1st, 2nd and 3rd
		x->pop_back(); xerr->pop_back();
		y->pop_back(); yerr->pop_back(); 

		x->push_back(entryPosition->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPosition->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01
		y->push_back(exitPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef8->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE8->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight8->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy8->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE8->Fill((momentum-5000)/(5000*5000));

		// here add 2nd, 3rd and 4th
		x->pop_back(); xerr->pop_back();
		x->pop_back(); xerr->pop_back();
		x->pop_back(); xerr->pop_back();
		y->pop_back(); yerr->pop_back(); 
		y->pop_back(); yerr->pop_back(); 
		y->pop_back(); yerr->pop_back(); 

		x->push_back(entryPosition->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPosition->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPositionGas->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01
		y->push_back(exitPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01
		y->push_back(exitPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef9->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE9->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight9->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy9->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE9->Fill((momentum-5000)/(5000*5000));

		// here add all 4
		x->push_back(entryPositionGas->x());
		xerr->push_back(1); //default 1
		y->push_back(entryPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		yerr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(x,y);
		MomentumResolutionRef10->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE10->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight10->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy10->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE10->Fill((momentum-5000)/(5000*5000));

		// here add only 2nd and 3rd but perfect TPC
		x->pop_back(); x->pop_back(); x->pop_back(); x->pop_back();
		xerr->pop_back(); xerr->pop_back(); xerr->pop_back(); xerr->pop_back();
		y->pop_back(); y->pop_back(); y->pop_back(); y->pop_back();
		yerr->pop_back(); yerr->pop_back(); yerr->pop_back(); yerr->pop_back();

		x->push_back(entryPosition->x());
		xerr->push_back(0);
		x->push_back(exitPosition->x());
		xerr->push_back(0);
		y->push_back(entryPosition->y());
		yerr->push_back(0.); //default si_reso
		y->push_back(exitPosition->y());
		yerr->push_back(0.); //default 0.01


		momentum=calculateMomentum(x,y);
		MomentumResolutionRef15->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE15->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		MomentumResolutionWeight15->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy15->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE15->Fill((momentum-5000)/(5000*5000));

		// here add all 4 but without the TPC
		std::vector<double>* xsi = new std::vector<double>;
		std::vector<double>* xsierr = new std::vector<double>;
		std::vector<double>* ysi =new std::vector<double>;
		std::vector<double>* ysierr = new std::vector<double>;

		xsi->push_back(entryPositionGas->x());
		xsi->push_back(entryPosition->x());
		xsi->push_back(exitPosition->x());
		xsi->push_back(exitPositionGas->x());
		xsierr->push_back(1);
		xsierr->push_back(1);
		xsierr->push_back(1);
		xsierr->push_back(1);

		ysi->push_back(entryPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysierr->push_back(si_reso); //default 0.01
		ysi->push_back(entryPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysierr->push_back(si_reso); //default 0.01
		ysi->push_back(exitPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysierr->push_back(si_reso); //default 0.01
		ysi->push_back(exitPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysierr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(xsi,ysi);
		MomentumResolutionRef14->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE14->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(xsi,ysi,xsierr,ysierr);
		MomentumResolutionWeight14->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy14->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE14->Fill((momentum-5000)/(5000*5000));


		// 2 outer si
		xsi->clear();
		xsierr->clear();
		ysi->clear();
		ysierr->clear();

		xsi->push_back(entryPositionGas->x());
		xsi->push_back(exitPositionGas->x());
		xsierr->push_back(1);
		xsierr->push_back(1);

		ysi->push_back(entryPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysi->push_back(exitPositionGas->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysierr->push_back(si_reso); //default 0.01
		ysierr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(xsi,ysi);
		MomentumResolutionRef16->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE16->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(xsi,ysi,xsierr,ysierr);
		MomentumResolutionWeight16->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy16->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE16->Fill((momentum-5000)/(5000*5000));

		// 2 inner si
		xsi->clear();
		xsierr->clear();
		ysi->clear();
		ysierr->clear();

		xsi->push_back(entryPosition->x());
		xsi->push_back(exitPosition->x());
		xsierr->push_back(1);
		xsierr->push_back(1);

		ysi->push_back(entryPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysi->push_back(exitPosition->y()+smearing->Gaus(0,si_reso)); //default 0.01
		ysierr->push_back(si_reso); //default 0.01
		ysierr->push_back(si_reso); //default 0.01

		momentum=calculateMomentum(xsi,ysi);
		MomentumResolutionRef17->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionSmearE17->Fill((momentum-5000)/(5000*5000));

		momentum=calculateWeightedMomentum(xsi,ysi,xsierr,ysierr);
		MomentumResolutionWeight17->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		//		MomentumResolutionWeightEnergy17->Fill((momentum-5000)/(5000*5000));
		MomentumResolutionWeightE17->Fill((momentum-5000)/(5000*5000));


		if (iterator > 99998) { std::cout<<"Hello "<<iterator<<std::endl; }
//Weight: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		x->clear();
		xerr->clear();
		y->clear();
		yerr->clear();


		if (iterator > 99998) { std::cout<<"Hello2 "<<std::endl; }
		//		refx1->clear();
		//		refy1->clear();

	}

	TF1 *CrystalBallReso=new TF1("CrystalBallReso",CrystalBall,-0.00002,0.00001,5);
	CrystalBallReso->SetParameters(1,1,0.000001,0.000001,3000);
	CrystalBallReso->SetParNames("alpha","n","mu","sigma","N");
	
	TF1 *CrystalBallResoHigh=new TF1("CrystalBallResoHigh",CrystalBall,-0.0003,0.0003,5);
	CrystalBallResoHigh->SetParameters(1,1,0.000001,0.000005,3000);
	CrystalBallResoHigh->SetParNames("alpha","n","mu","sigma","N");


		TF1 *CrystalBallSmear=new TF1("CrystalBallSmear",CrystalBall,-0.00008,0.00002,5);
		CrystalBallSmear->SetParameters(1,1,0.000001,0.000001,3000);
		CrystalBallSmear->SetParNames("alpha","n","mu","sigma","N");

		TF1 *CrystalBallFnc=new TF1("CrystalBallFit",CrystalBall,1000,6000,5);
		CrystalBallFnc->SetParameters(1,1,5000,25,3000);
		CrystalBallFnc->SetParNames("alpha","n","mu","sigma","N");


		c1->cd(1);
		MomentumResolution->Fit("CrystalBallReso","R");
		MomentumResolution->Draw();
		c1->cd(2);
		MomentumResolutionSmearing->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearing->Draw();
		c1->cd(3);
		MomentumResolutionSmearing2->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearing2->Draw();
		//		MomentumResolutionRef->Fit("CrystalBallReso","R");
		//		MomentumResolutionRef->Draw();
		c1->cd(4);
		MomentumResolutionWeight->Fit("CrystalBallReso","R");
		MomentumResolutionWeight->Draw();

	std::cout<<"Hello5 "<<std::endl; 

		c2->cd(1);
		MomentumResolutionRef2->Fit("CrystalBallReso","R");
		MomentumResolutionRef2->Draw();
		c2->cd(2);
		MomentumResolutionWeight2->Fit("CrystalBallReso","R");
		MomentumResolutionWeight2->Draw();
		c2->cd(3);
		MomentumResolutionSmearE2->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE2->Draw();
		c2->cd(4);
		MomentumResolutionWeightE2->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE2->Draw();

		c3->cd(1);
		MomentumResolutionRef3->Fit("CrystalBallReso","R");
		MomentumResolutionRef3->Draw();
		c3->cd(2);
		MomentumResolutionWeight3->Fit("CrystalBallReso","R");
		MomentumResolutionWeight3->Draw();
		c3->cd(3);
		MomentumResolutionSmearE3->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE3->Draw();
		c3->cd(4);
		MomentumResolutionWeightE3->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE3->Draw();

		c4->cd(1);
		MomentumResolutionRef4->Fit("CrystalBallReso","R");
		MomentumResolutionRef4->Draw();
		c4->cd(2);
		MomentumResolutionWeight4->Fit("CrystalBallReso","R");
		MomentumResolutionWeight4->Draw();
		c4->cd(3);
		MomentumResolutionSmearE4->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE4->Draw();
		c4->cd(4);
		MomentumResolutionWeightE4->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE4->Draw();

		c5->cd(1);
		MomentumResolutionRef5->Fit("CrystalBallReso","R");
		MomentumResolutionRef5->Draw();
		c5->cd(2);
		MomentumResolutionWeight5->Fit("CrystalBallReso","R");
		MomentumResolutionWeight5->Draw();
		c5->cd(3);
		MomentumResolutionSmearE5->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE5->Draw();
		c5->cd(4);
		MomentumResolutionWeightE5->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE5->Draw();

		c6->cd(1);
		MomentumResolutionRef6->Fit("CrystalBallReso","R");
		MomentumResolutionRef6->Draw();
		c6->cd(2);
		MomentumResolutionWeight6->Fit("CrystalBallReso","R");
		MomentumResolutionWeight6->Draw();
		c6->cd(3);
		MomentumResolutionSmearE6->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE6->Draw();
		c6->cd(4);
		MomentumResolutionWeightE6->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE6->Draw();

		c7->cd(1);
		MomentumResolutionRef7->Fit("CrystalBallReso","R");
		MomentumResolutionRef7->Draw();
		c7->cd(2);
		MomentumResolutionWeight7->Fit("CrystalBallReso","R");
		MomentumResolutionWeight7->Draw();
		c7->cd(3);
		MomentumResolutionSmearE7->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE7->Draw();
		c7->cd(4);
		MomentumResolutionWeightE7->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE7->Draw();

		c8->cd(1);
		MomentumResolutionRef8->Fit("CrystalBallReso","R");
		MomentumResolutionRef8->Draw();
		c8->cd(2);
		MomentumResolutionWeight8->Fit("CrystalBallReso","R");
		MomentumResolutionWeight8->Draw();
		c8->cd(3);
		MomentumResolutionSmearE8->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE8->Draw();
		c8->cd(4);
		MomentumResolutionWeightE8->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE8->Draw();

		c9->cd(1);
		MomentumResolutionRef9->Fit("CrystalBallReso","R");
		MomentumResolutionRef9->Draw();
		c9->cd(2);
		MomentumResolutionWeight9->Fit("CrystalBallReso","R");
		MomentumResolutionWeight9->Draw();
		c9->cd(3);
		MomentumResolutionSmearE9->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE9->Draw();
		c9->cd(4);
		MomentumResolutionWeightE9->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE9->Draw();

		c10->cd(1);
		MomentumResolutionRef10->Fit("CrystalBallReso","R");
		MomentumResolutionRef10->Draw();
		c10->cd(2);
		MomentumResolutionWeight10->Fit("CrystalBallReso","R");
		MomentumResolutionWeight10->Draw();
		c10->cd(3);
		MomentumResolutionSmearE10->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE10->Draw();
		c10->cd(4);
		MomentumResolutionWeightE10->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE10->Draw();


		if ( si_reso > 0.15 ) {
		  c14->cd(1);
		  MomentumResolutionRef14->Fit("CrystalBallResoHigh","R");
		  MomentumResolutionRef14->Draw();
		  c14->cd(2);
		  MomentumResolutionWeight14->Fit("CrystalBallResoHigh","R");
		  MomentumResolutionWeight14->Draw();
		}
		else {
		  c14->cd(1);
		  MomentumResolutionRef14->Fit("CrystalBallReso","R");
		  MomentumResolutionRef14->Draw();
		  c14->cd(2);
		  MomentumResolutionWeight14->Fit("CrystalBallReso","R");
		  MomentumResolutionWeight14->Draw();
		}
		c14->cd(3);
		MomentumResolutionSmearE14->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE14->Draw();
		c14->cd(4);
		MomentumResolutionWeightE14->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE14->Draw();

		c15->cd(1);
		MomentumResolutionRef15->Fit("CrystalBallReso","R");
		MomentumResolutionRef15->Draw();
		c15->cd(2);
		MomentumResolutionWeight15->Fit("CrystalBallReso","R");
		MomentumResolutionWeight15->Draw();
		c15->cd(3);
		MomentumResolutionSmearE15->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE15->Draw();
		c15->cd(4);
		MomentumResolutionWeightE15->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE15->Draw();

		if ( si_reso > 0.15 ) {
		  c16->cd(1);
		  MomentumResolutionRef16->Fit("CrystalBallResoHigh","R");
		  MomentumResolutionRef16->Draw();
		  c16->cd(2);
		  MomentumResolutionWeight16->Fit("CrystalBallResoHigh","R");
		  MomentumResolutionWeight16->Draw();
		}
		else {
		  c16->cd(1);
		  MomentumResolutionRef16->Fit("CrystalBallReso","R");
		  MomentumResolutionRef16->Draw();
		  c16->cd(2);
		  MomentumResolutionWeight16->Fit("CrystalBallReso","R");
		  MomentumResolutionWeight16->Draw();
		}
		c16->cd(3);
		MomentumResolutionSmearE16->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE16->Draw();
		c16->cd(4);
		MomentumResolutionWeightE16->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE16->Draw();

		if ( si_reso > 0.15 ) {
		  c17->cd(1);
		  MomentumResolutionRef17->Fit("CrystalBallResoHigh","R");
		  MomentumResolutionRef17->Draw();
		  c17->cd(2);
		  MomentumResolutionWeight17->Fit("CrystalBallResoHigh","R");
		  MomentumResolutionWeight17->Draw();
		}
		else {
		  c17->cd(1);
		  MomentumResolutionRef17->Fit("CrystalBallReso","R");
		  MomentumResolutionRef17->Draw();
		  c17->cd(2);
		  MomentumResolutionWeight17->Fit("CrystalBallReso","R");
		  MomentumResolutionWeight17->Draw();
		}
		c17->cd(3);
		MomentumResolutionSmearE17->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearE17->Draw();
		c17->cd(4);
		MomentumResolutionWeightE17->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightE17->Draw();

		TFile* outfile = TFile::Open("MomentumResolution_MyPlots.root","recreate");
		c1->Write();
		c2->Write();
		c3->Write();
		c4->Write();
		c5->Write();
		c6->Write();
		c7->Write();
		c8->Write();
		c9->Write();
		c10->Write();
		c14->Write();
		c15->Write();
		c16->Write();
		c17->Write();
		MomentumResolution->Write();
		MomentumResolutionSmearing->Write();
		MomentumResolutionSmearing2->Write();
		//		MomentumResolutionRef->Write();
		MomentumResolutionWeight->Write();
		MomentumResolutionSmearingEnergy->Write();
		MomentumResolutionWeightEnergy->Write();
		beamMomentum->Write();
		magnetMomentum->Write();
		tpcMomentum->Write();
		recoMomentum->Write();
		recoMomentumSmearing->Write();
		recoMomentumRef->Write();
		recoMomentumWeight->Write();

		MomentumResolutionRef2->Write();
		MomentumResolutionWeight2->Write();
		MomentumResolutionRef3->Write();
		MomentumResolutionWeight3->Write();
		MomentumResolutionRef4->Write();
		MomentumResolutionWeight4->Write();
		MomentumResolutionRef5->Write();
		MomentumResolutionWeight5->Write();
		MomentumResolutionRef6->Write();
		MomentumResolutionWeight6->Write();
		MomentumResolutionRef7->Write();
		MomentumResolutionWeight7->Write();
		MomentumResolutionRef8->Write();
		MomentumResolutionWeight8->Write();
		MomentumResolutionRef9->Write();
		MomentumResolutionWeight9->Write();
		MomentumResolutionRef10->Write();
		MomentumResolutionWeight10->Write();
		MomentumResolutionRef14->Write();
		MomentumResolutionWeight14->Write();	
		MomentumResolutionRef15->Write();
		MomentumResolutionWeight15->Write();
		MomentumResolutionRef16->Write();
		MomentumResolutionWeight16->Write();
		MomentumResolutionRef17->Write();
		MomentumResolutionWeight17->Write();
		outfile->Close();

}

