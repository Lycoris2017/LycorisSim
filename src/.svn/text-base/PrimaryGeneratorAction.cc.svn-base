// $Id: PrimaryGeneratorAction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file
 * @brief implements mandatory user class PrimaryGeneratorAction
 */

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
  : outfile(0)
{
  gun = InitializeGPS();
  
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  gun->GeneratePrimaryVertex(anEvent);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete gun;
}

G4VPrimaryGenerator* PrimaryGeneratorAction::InitializeGPS()
{
  G4GeneralParticleSource * gps = new G4GeneralParticleSource();
  
  // setup details easier via UI commands see gps.mac

  // particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* electron = particleTable->FindParticle("e-");
   gps->GetCurrentSource()->SetParticleDefinition(electron);



  // set energy distribution
  G4SPSEneDistribution *eneDist = gps->GetCurrentSource()->GetEneDist() ;
  eneDist->SetEnergyDisType("Gauss"); // or gauss
  eneDist->SetMonoEnergy(5000*MeV);
  eneDist->SetBeamSigmaInE(100*MeV);

  // set position distribution
  G4SPSPosDistribution* posDist = gps->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Beam");  // or Point,Plane,Volume,Beam
  posDist->SetPosDisShape("Circle");
  posDist->SetCentreCoords(G4ThreeVector(100.0*cm,0*cm,0*cm));
  posDist->SetPosRot1(G4ThreeVector(0.,1.,0.));
  posDist->SetPosRot2(G4ThreeVector(0.,0.,1.));
  //posDist->SetRadius(1*mm);
  //posDist->SetBeamSigmaInX(0.001*mm);
  //posDist->SetBeamSigmaInY(0.001*mm);
  posDist->SetBeamSigmaInR(1*mm);

  // set angular distribution
   G4SPSAngDistribution *angDist = gps->GetCurrentSource()->GetAngDist();
   angDist->SetParticleMomentumDirection( G4ThreeVector(-1., 0., 0.) );
//   angDist->SetAngDistType("beam2d");
//   angDist->SetBeamSigmaInAngX(0.1*mrad);
//   angDist->SetBeamSigmaInAngY(0.1*mrad);
   //angDist->DefineAngRefAxes("angref1",G4ThreeVector(-1.,0.,0.));

  return gps;
}

