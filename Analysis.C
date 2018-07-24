void Analysis(){
	TCut cutSelection("");
	//TCut cutSelection("entryEnergy<5500 && trackLength>764 && noOfChargedParticleEnteringDriftVolume==0");
	TFile *f1  = TFile::Open("run_0.root");
	TTree *t  = (TTree*)f1->Get("tree");
	TCanvas  *c1=new TCanvas("c1","c1",0,0,1200,900);
	c1->Divide(2,2);
	c1->cd(1);
	TH1F* hNoOfSecPart = new TH1F("NumberOfSecondaryParticles","number of secondary particles",1000,-0.5,999.5);
	t->Draw("noOfSecondaryParticle>>NumberOfSecondaryParticles",cutSelection);
	hNoOfSecPart->GetXaxis()->SetTitle("# secondary particles");
	hNoOfSecPart->Draw();

	//TCanvas c2("c2","c2",0,0,1200,900);
	c1->cd(2);
	TH1F* hNoOfSecChargedPart = new TH1F("NumberOfSecondaryChargedParticles","number of secondary charged particles",1000,-0.5,999.5);
	t->Draw("noOfSecondaryChargedParticle>>NumberOfSecondaryChargedParticles",cutSelection);
	hNoOfSecChargedPart->GetXaxis()->SetTitle("# secondary charged particles");
	hNoOfSecChargedPart->Draw();


	//TCanvas c3("c3","c3",0,0,1200,900);
	c1->cd(3);
	TH1F* hNoOfSecPartEnteringVolume = new TH1F("NumberOfSecondaryParticlesEnteringVolume","number of secondary charged particles entering the drift volume",16,-0.5,15.5);
	t->Draw("noOfChargedParticleEnteringDriftVolume>>NumberOfSecondaryParticlesEnteringVolume",cutSelection);
	hNoOfSecPartEnteringVolume->GetXaxis()->SetTitle("# secondary particles");
	hNoOfSecPartEnteringVolume->Draw();

	TCanvas  *c2=new TCanvas("c2","c2",0,0,1200,900);
	c2->Divide(2,2);
	c2->cd(2);
	TH1F* hEntryE = new TH1F("EntryEnergyOfPrimaryParticle","energy of primary particle during entry",5500,0,5500);
	t->Draw("entryEnergy>>EntryEnergyOfPrimaryParticle",cutSelection);
	hEntryE->GetXaxis()->SetTitle("E_{prim} (MeV) ");
	hEntryE->Draw();

	//TCanvas c5("c5","c5",0,0,1200,900);
	c2->cd(3);
	TH1F* hExitE = new TH1F("ExitEnergyOfPrimaryParticle","energy of primary particle during exit",5500,0,5500);
	t->Draw("exitEnergy>>ExitEnergyOfPrimaryParticle",cutSelection);
	hExitE->GetXaxis()->SetTitle("E_{prim} (MeV) ");
	hExitE->Draw();

	c2->cd(4);
	TH1F* hExitEMagnet = new TH1F("ExitEnergyofPrimaryParticleMagnet","energy of primary particle during exit in the magnet",5500,0,5500);
		t->Draw("exitEnergyMagnet>>ExitEnergyofPrimaryParticleMagnet",cutSelection);
		hExitEMagnet->GetXaxis()->SetTitle("E_{prim} (MeV) ");
		hExitEMagnet->Draw();

	//TCanvas c6("c6","c6",0,0,1200,900);
	c2->cd(1);
	TH1F* hEntryEMagnet = new TH1F("EntryEnergyofPrimaryParticleMagnet","energy of primary particle during entry in the magnet",5500,0,5500);
	t->Draw("entryEnergyMagnet>>EntryEnergyofPrimaryParticleMagnet",cutSelection);
	hEntryEMagnet->GetXaxis()->SetTitle("E_{prim} (MeV) ");	
	hEntryEMagnet->Draw();

	TCanvas  *c3=new TCanvas("c3","c3",0,0,1200,900);
	c3->Divide(2,2);
	c3->cd(1);
	TH1F* hTotE = new TH1F("TotalEnergyInDriftVolume","deposited energy in drift volume by primary particle",256,0,1);
	t->Draw("eDepTPCGas>>TotalEnergyInDriftVolume",cutSelection);
	hTotE->GetXaxis()->SetTitle("E_{dep} (MeV) ");
	hTotE->Draw();

	//TCanvas c8("c8","c8",0,0,1200,900);
	c3->cd(2);
	TH1F* henergyOfSecChargedPart = new TH1F("EnergyOfSecChargedPart","energy of secondary charged part",256,0,0.1);
	t->Draw("energyOfSecondaryChargedPart>>EnergyOfSecChargedPart",cutSelection);
	henergyOfSecChargedPart->GetXaxis()->SetTitle("E (MeV) ");
	henergyOfSecChargedPart->Draw();

	TCanvas  *c4=new TCanvas("c4","c4",0,0,1200,900);
	c4->Divide(2,2);
	c4->cd(1);
	TH2F* hEntryPosMagnet = new TH2F("EntryPositionOfPrimaryParticleMagnet","entry position of primary particle at the magnet",256,-20,20,256,-20,20);
	t->Draw("entryPositionMagnet.fZ:entryPositionMagnet.fY>>EntryPositionOfPrimaryParticleMagnet",cutSelection); 
	hEntryPosMagnet->GetXaxis()->SetTitle("y_{entry} ");
	hEntryPosMagnet->GetYaxis()->SetTitle("z_{entry} ");
	hEntryPosMagnet->Draw("colz");

	//TCanvas c10("c10","c10",0,0,1200,900);
	c4->cd(2);
	TH2F* hEntryPos = new TH2F("EntryPositionOfPrimaryParticle","entry position of primary particle",256,-20,20,256,-20,20);
	t->Draw("entryPosition.fZ:entryPosition.fY>>EntryPositionOfPrimaryParticle",cutSelection); 
	hEntryPos->GetXaxis()->SetTitle("y_{entry} ");
	hEntryPos->GetYaxis()->SetTitle("z_{entry} ");
	hEntryPos->Draw("colz");

	//TCanvas c11("c11","c11",0,0,1200,900);
	c4->cd(3);
	TH2F* hExitPos = new TH2F("ExitPositionOfPrimaryParticle","exit position of primary particle",256,-40,0,256,-20,20);
	hExitPos->GetDirectory()->cd();
	t->Draw("exitPosition.fZ:exitPosition.fY>>ExitPositionOfPrimaryParticle",cutSelection);
	hExitPos->GetXaxis()->SetTitle("y_{exit} ");
	hExitPos->GetYaxis()->SetTitle("z_{exit} ");
	hExitPos->Draw("colz");
	

	c4->cd(4);
	TH2F* hExitPosMagnet = new TH2F("ExitPositionOfPrimaryParticleMagnet","exit position of primary particle ath the magnet",256,-40,0,256,-20,20);
	hExitPosMagnet->GetDirectory()->cd();
	t->Draw("exitPositionMagnet.fZ:exitPositionMagnet.fY>>ExitPositionOfPrimaryParticleMagnet",cutSelection);
	hExitPosMagnet->GetXaxis()->SetTitle("y_{exit} ");
	hExitPosMagnet->GetYaxis()->SetTitle("z_{exit} ");
	hExitPosMagnet->Draw("colz");

	TCanvas  *c5=new TCanvas("c5","c5",0,0,1200,900);
	c5->Divide(2,2);
	c5->cd(1);
	TH2F* hDEDX = new TH2F("dEdx","dEdx",256,0, 5500,128,0,0.0008);
	t->Draw("eDepTPCGas/trackLength:entryEnergy>>dEdx",cutSelection);
	hDEDX->GetYaxis()->SetTitle("dE/dx");
	hDEDX->GetXaxis()->SetTitle("E");
	hDEDX->Draw("colz");
	
	//TCanvas c13("c13","c13",0,0,1200,900);
	c5->cd(2);
	TH1F* hTrackLength = new TH1F("TrackLength","track length",1000,0,1000);
	t->Draw("trackLength>>TrackLength",cutSelection);
	hTrackLength->GetXaxis()->SetTitle("track length");
	hTrackLength->Draw();
	TFile* outfile = TFile::Open("Analysis.root","recreate");
	c1->Write();
	c2->Write();
	c3->Write();
	c4->Write();
	c5->Write();
	outfile->Close();

}
