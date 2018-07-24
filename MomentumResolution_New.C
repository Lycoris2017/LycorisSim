
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


//void MomentumResolutionAnalysis(){
void MomentumResolution_New(){
  //	TFile *f1  = TFile::Open("run_0.root","Update");
  
  gStyle->SetOptFit(1);

  TFile *f1  = TFile::Open("run0_BeamSource.root");
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
	TCanvas *c1=new TCanvas("MomentumResolutions","Momentum Resolutions",0,0,1200,900);
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
	TH1D *MomentumResolutionRef=new TH1D("MomentumResolutionReference","Momentum Resolution with reference", 512, -0.00006, 0.00006);
	MomentumResolutionRef->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRef->GetYaxis()->SetTitle("number of events");
	TCanvas *c2=new TCanvas("Measurement","Measurement",0,0,1200,900);
	c2->Divide(2,2);
	TH1D *MomentumResolutionWeight=new TH1D("MomentumResolutionWeight","Momentum Resolution with Weighted Track Points", 512, -0.00006, 0.00006);
	MomentumResolutionWeight->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeight->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergy=new TH1D("MomentumResolutionWeightEenrgy","Momentum Resolution with Weighted Track Points and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergy->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergy->GetYaxis()->SetTitle("number of events");
	TCanvas *c3=new TCanvas("simulatedMomentum","simulated momentum",0,0,1200,900);
	c3->Divide(2,2);
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
	TCanvas *c4=new TCanvas("reconstructedMomentum","reconstructed momentum",0,0,1200,900);
	c4->Divide(2,2);
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


	TCanvas *c5=new TCanvas("MomentumResolutionsLast2","Momentum Resolutions Last 2",0,0,1200,900);
	c5->Divide(2,2);
	TH1D *MomentumResolutionRefLastTwoSiLayers=new TH1D("MomentumResolutionReference","Momentum Resolution with reference (Last 2 Si)", 512, -0.00006, 0.00006);
	MomentumResolutionRefLastTwoSiLayers->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRefLastTwoSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightLastTwoSiLayers=new TH1D("MomentumResolutionWeight","Momentum Resolution with Weighted Track Points", 512, -0.00006, 0.00006);
	MomentumResolutionWeightLastTwoSiLayers->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightLastTwoSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergyLastTwoSiLayers=new TH1D("MomentumResolutionWeightEenrgy","Momentum Resolution with Weighted Track Points and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergyLastTwoSiLayers->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergyLastTwoSiLayers->GetYaxis()->SetTitle("number of events");


	TCanvas *c6=new TCanvas("MomentumResolutions All SI","Momentum Resolutions All Si",0,0,1200,900);
	c6->Divide(2,2);
	TH1D *MomentumResolutionRefAllSiLayers=new TH1D("MomentumResolutionReference","Momentum Resolution with reference (All 4 Si)", 512, -0.00006, 0.00006);
	MomentumResolutionRefAllSiLayers->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionRefAllSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightAllSiLayers=new TH1D("MomentumResolutionWeight","Momentum Resolution with Weighted Track Points", 512, -0.00006, 0.00006);
	MomentumResolutionWeightAllSiLayers->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightAllSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *MomentumResolutionWeightEnergyAllSiLayers=new TH1D("MomentumResolutionWeightEenrgy","Momentum Resolution with Weighted Track Points and unknown Energy", 512, -0.0001, 0.00006);
	MomentumResolutionWeightEnergyAllSiLayers->GetXaxis()->SetTitle("#Delta p_{t,reco}/p_{t}^{2} [1/MeV]");
	MomentumResolutionWeightEnergyAllSiLayers->GetYaxis()->SetTitle("number of events");

	TCanvas *c7=new TCanvas("MomentumResolutions Tests","Momentum Resolutions Tests",0,0,1200,900);
	c7->Divide(2,2);
	TH1D *recoMomentumRefLastTwoSiLayers=new TH1D("recoMomentumRef","Reconstructed Momentum with RefLastTwoSiLayerserence",512, 0, 6000);
	recoMomentumRefLastTwoSiLayers->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumRefLastTwoSiLayers->GetXaxis()->SetNdivisions(505);
	recoMomentumRefLastTwoSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *recoMomentumWeightLastTwoSiLayers=new TH1D("recoMomentumWeight","Reconstructed Momentum with Weighted Track Points",512, 0, 6000);
	recoMomentumWeightLastTwoSiLayers->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumWeightLastTwoSiLayers->GetXaxis()->SetNdivisions(505);
	recoMomentumWeightLastTwoSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *recoMomentumRefAllSiLayers=new TH1D("recoMomentumRef","Reconstructed Momentum with RefAllSiLayerserence",512, 0, 6000);
	recoMomentumRefAllSiLayers->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumRefAllSiLayers->GetXaxis()->SetNdivisions(505);
	recoMomentumRefAllSiLayers->GetYaxis()->SetTitle("number of events");
	TH1D *recoMomentumWeightAllSiLayers=new TH1D("recoMomentumWeight","Reconstructed Momentum with Weighted Track Points",512, 0, 6000);
	recoMomentumWeightAllSiLayers->GetXaxis()->SetTitle("p_{t,reco}/MeV");
	recoMomentumWeightAllSiLayers->GetXaxis()->SetNdivisions(505);
	recoMomentumWeightAllSiLayers->GetYaxis()->SetTitle("number of events");




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

		if(iterator<100) std::cout<<"Smearing: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		// here only the first 2 Si layers are added
		x->push_back(entryPosition->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPosition->x());
		xerr->push_back(1);

		if(iterator<100) std::cout<<" x entry "<<entryPosition->x()<<" x exit "<<exitPosition->x()<<std::endl;
		y->push_back(entryPosition->y()+smearing->Gaus(0,0.01)); //default 0.01
		yerr->push_back(0.01); //default 0.01
		y->push_back(exitPosition->y()+smearing->Gaus(0,0.01));
		yerr->push_back(0.01);

		momentum=calculateMomentum(x,y);
		recoMomentum->Fill(momentum);
		MomentumResolutionRef->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));

		if(iterator<100) std::cout<<"Ref: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		recoMomentumWeight->Fill(momentum);
		MomentumResolutionWeight->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionWeightEnergy->Fill((momentum-5000)/(5000*5000));

		if(iterator<100) std::cout<<"Weight: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;


		// here only the last 2 Si layers are added

		// ** delete the first 2 Si layers
		x->pop_back();
		x->pop_back();
		xerr->pop_back();
		xerr->pop_back();
		y->pop_back();  y->pop_back();
		yerr->pop_back();  yerr->pop_back();
		// ** add the last 2 Si layers
		x->push_back(entryPositionGas->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPositionGas->x());
		xerr->push_back(1);

		if(iterator<100) std::cout<<" Last two: x entry "<<entryPositionGas->x()<<" x exit "<<exitPositionGas->x()<<std::endl;
		y->push_back(entryPositionGas->y()+smearing->Gaus(0,0.01)); //default 0.01
		yerr->push_back(0.01); //default 0.01
		y->push_back(exitPositionGas->y()+smearing->Gaus(0,0.01));
		yerr->push_back(0.01);

		momentum=calculateMomentum(x,y);
		recoMomentumRefLastTwoSiLayers->Fill(momentum);
		MomentumResolutionRefLastTwoSiLayers->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));

		if(iterator<100) std::cout<<"Ref: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		recoMomentumWeightLastTwoSiLayers->Fill(momentum);
		MomentumResolutionWeightLastTwoSiLayers->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionWeightEnergyLastTwoSiLayers->Fill((momentum-5000)/(5000*5000));

		if(iterator<100) std::cout<<"Weight: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		// here all the 4 2 Si layers are added
		// ** delete the last 2 Si layers
		x->pop_back();
		xerr->pop_back();
		x->pop_back();
		xerr->pop_back();
		y->pop_back();  y->pop_back();
		yerr->pop_back();  yerr->pop_back();
		// ** add the first 2 Si layers
		x->push_back(entryPosition->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPosition->x());
		xerr->push_back(1);
		// ** add the last 2 Si layers
		x->push_back(entryPositionGas->x());
		xerr->push_back(1); //default 1
		x->push_back(exitPositionGas->x());
		xerr->push_back(1);

		if(iterator<100) std::cout<<" x entry "<<entryPosition->x()<<" x exit "<<exitPosition->x()<<std::endl;
		y->push_back(entryPosition->y()+smearing->Gaus(0,0.01)); //default 0.01
		yerr->push_back(0.01); //default 0.01
		y->push_back(exitPosition->y()+smearing->Gaus(0,0.01));
		yerr->push_back(0.01);
		y->push_back(entryPositionGas->y()+smearing->Gaus(0,0.01)); //default 0.01
		yerr->push_back(0.01); //default 0.01
		y->push_back(exitPositionGas->y()+smearing->Gaus(0,0.01));
		yerr->push_back(0.01);

		momentum=calculateMomentum(x,y);
		recoMomentumRefAllSiLayers->Fill(momentum);
		MomentumResolutionRefAllSiLayers->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));

		if(iterator<100) std::cout<<"Ref: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		momentum=calculateWeightedMomentum(x,y,xerr,yerr);
		recoMomentumWeightAllSiLayers->Fill(momentum);
		MomentumResolutionWeightAllSiLayers->Fill((momentum-entryEnergy)/(entryEnergy*entryEnergy));
		MomentumResolutionWeightEnergyAllSiLayers->Fill((momentum-5000)/(5000*5000));

		if(iterator<100) std::cout<<"Weight: x "<<x->size() <<" y "<<y->size()<<" momentum "<<momentum<<" reso "<<(momentum-entryEnergy)/(entryEnergy*entryEnergy)<<std::endl;

		// here only the 4 Si layers are used without the TPC

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
		MomentumResolutionRef->Fit("CrystalBallReso","R");
		MomentumResolutionRef->Draw();
		c1->cd(4);
		MomentumResolutionWeight->Fit("CrystalBallReso","R");
		MomentumResolutionWeight->Draw();

	std::cout<<"Hello5 "<<std::endl; 

		c2->cd(1);
		MomentumResolutionSmearingEnergy->Fit("CrystalBallSmear","R");
		MomentumResolutionSmearingEnergy->Draw();
		c2->cd(2);
		MomentumResolutionWeightEnergy->Fit("CrystalBallSmear","R");
		MomentumResolutionWeightEnergy->Draw();
		c3->cd(1);
//		beamMomentum->Fit("CrystalBallFit","R");
		beamMomentum->Draw();
		c3->cd(2);
//		magnetMomentum->Fit("CrystalBallFit","R");
		magnetMomentum->Draw();
		c3->cd(3);
//		tpcMomentum->Fit("CrystalBallFit","R");
		tpcMomentum->Draw();
		c4->cd(1);
		//		recoMomentum->Fit("CrystalBallFit","R");
		recoMomentum->Draw();
		c4->cd(2);
		//		recoMomentumSmearing->Fit("CrystalBallFit","R");
		recoMomentumSmearing->Draw();
		c4->cd(3);
		//		recoMomentumRef->Fit("CrystalBallFit","R");
		recoMomentumRef->Draw();
		c4->cd(4);
		//		recoMomentumWeight->Fit("CrystalBallFit","R");
		recoMomentumWeight->Draw();



		c5->cd(1);
		MomentumResolutionRefLastTwoSiLayers->Fit("CrystalBallReso","R");
		MomentumResolutionRefLastTwoSiLayers->Draw();
		c5->cd(2);
		MomentumResolutionWeightLastTwoSiLayers->Fit("CrystalBallReso","R");
		MomentumResolutionWeightLastTwoSiLayers->Draw();
		c5->cd(4);
		MomentumResolutionWeightEnergyLastTwoSiLayers->Fit("CrystalBallSmear","R");

		c6->cd(1);
		MomentumResolutionRefAllSiLayers->Fit("CrystalBallReso","R");
		MomentumResolutionRefAllSiLayers->Draw();
		c6->cd(2);
		MomentumResolutionWeightAllSiLayers->Fit("CrystalBallReso","R");
		MomentumResolutionWeightAllSiLayers->Draw();
		c6->cd(4);
		MomentumResolutionWeightEnergyAllSiLayers->Fit("CrystalBallSmear","R");

		TFile* outfile = TFile::Open("MomentumResolution_New.root","recreate");
		c1->Write();
		c2->Write();
		c3->Write();
		c4->Write();
		c5->Write();
		c6->Write();
		MomentumResolution->Write();
		MomentumResolutionSmearing->Write();
		MomentumResolutionRef->Write();
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

		MomentumResolutionRefLastTwoSiLayers->Write();
		MomentumResolutionWeightLastTwoSiLayers->Write();
		MomentumResolutionWeightEnergyLastTwoSiLayers->Write();
		recoMomentumRefLastTwoSiLayers->Write();
		recoMomentumWeightLastTwoSiLayers->Write();

		MomentumResolutionRefAllSiLayers->Write();
		MomentumResolutionWeightAllSiLayers->Write();
		MomentumResolutionWeightEnergyAllSiLayers->Write();
		recoMomentumRefAllSiLayers->Write();
		recoMomentumWeightAllSiLayers->Write();

		outfile->Close();

}

