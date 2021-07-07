#include <TROOT.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <TChainElement.h>
#include <TMath.h>
#include <TH2Poly.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <string>
#include <TString.h>
#include <TCanvas.h>
#include <vector>
#include <algorithm>
#include <TLorentzVector.h>
// #include <TIter.h>

// #include "G4SBSRunData.hh"
// #include "G4SBSTextFile.hh"
// #include "g4sbs_tree_mtpc.C"
#include "g4sbstree_mtpc.C"

#include "GetPadLoaded.cpp"
#include "SampaShapingAndSample_ReturnVals.cpp"

#include <chrono>

using namespace std;


//==============================================================================
// FUNCTIONS
vector<Double_t> CalcCharge(Double_t Edep, Double_t PathLength);
Double_t ZtoGEM(Double_t zpos);


//==============================================================================
// GLOBAL VARIABLES
const Int_t GEMGain = 1000;//1000; //is typical gem gain 1k?


//==============================================================================
// START MAIN FUNCTION
void Analyse_mTPC_Trk(const char *infile = "/scratch/g4sbs_out_protons_60-400MeV_30-70degrees_1kevents"){
  // infile is the name of the g4sbs output file like "pathtofile/filename"
  // remove the .root from infile

  // to execute macro 
  // root -l
  // .L Analyse_mTPC_Trk.cpp+
  // Analyse_mTPC_Trk("path/filename")
  
  // varaibles for chrono
  std::chrono::time_point<std::chrono::steady_clock> start_all = 
    std::chrono::steady_clock::now();
  double t0 = 0.0, t1 = 0.0, t2 = 0.0;
  
  //============================================================================
  // FLAGS AND VARIABLES
  Int_t EventCounter = 0;
  Int_t nevent = 0;
  TRandom3 RandomVar(0);

  //
  // values for diffusion from 07/04/10 update tdis daq
  //
  // he:ch4 90:10 at 1.2kV/cm:
  // Double_t DL = 0.022;//0.015;//0.03;//0.015;//sqrt(cm) up to 0.03sqrt(cm)
  // Double_t DT = 0.025;
  // vd = 19um/ns = max 2.6us
  //
  // he:ch4 70:30 at 1.2kV/cm:
  Double_t DL = 0.015;//0.015;//0.03;//0.015;//sqrt(cm) up to 0.03sqrt(cm)
  Double_t DT = 0.015;
  Double_t MaxDriftTime = 1700.0;// convert 1.7e-6 to ns
  Double_t DriftVelocity = 32.1;//30um/ns
  // vd = 30um/ns = max 1.7us

  // he:ch4 70:30 at 3kV/cm for gem regions
  // Double_t DL_GEM = 0.017;////sqrt(cm)
  // Double_t DT_GEM = 0.025;//sqrt(cm)
  // Double_t DriftVelocity_GEM = 32.0; // um/ns

  // he:ch4 70:30 at 3.4kV/cm for gem regions
  Double_t DL_GEM = 0.014468; //sqrt(cm)
  Double_t DT_GEM = 0.0263103; //sqrt(cm)
  Double_t DriftVelocity_GEM = 33.04; // um/ns


  // gem dimensions
  const Double_t GEMGap = 0.1; // take 1mm, units cm to match diffusion units
  // smer the transverse drift of each electron to get charge cloud area
  Double_t sigmatransverse_GEM = DT_GEM * TMath::Sqrt(GEMGap);
  // smear the longitudinal drift of each electron to get time window 
  Double_t sigmalongitudinal_GEM = DL_GEM * TMath::Sqrt(GEMGap); 

  //============================================================================
  
  // get histo of readout pads
  // TFile *fPads =  TFile::Open("hPads.root"); // my original pads with staggering by ring id
  TFile *fPads =  TFile::Open("hPads2.root"); //hPads2 is steves suggested pad phis
  TCanvas *c1 = (TCanvas*)fPads->Get("c1");
  TH2Poly *hReadoutPads = (TH2Poly*)c1->GetPrimitive("hPolyBins");
  // TCanvas *cPads = new TCanvas();
  // cPads->cd();
  // hReadoutPads->Draw();
  fPads->Close();
  // get bin list
  TList *bin_list=hReadoutPads->GetBins();
  if( !bin_list ){
    cout << "Problem in getting bin list from hReadoutPads, exitting." << endl;
    exit (EXIT_FAILURE);
  }
  TIter next(bin_list);
  TObject *obj;
  TH2PolyBin *bin;
  double area;
  vector<Double_t> vbinareas;
  while((obj=next())){
    bin = (TH2PolyBin*) obj;
    vbinareas.push_back(bin->GetArea());
  }

  //  READ IN SIM FILE
  TChain *C = new TChain("T");
  TString simfile = infile;
  simfile += ".root";
  // C->Add("/scratch/g4sbs_out_protons_60-400MeV_30-70degrees_1kevents.root");
  C->Add(simfile);
  
  // ASSIGN TREE
  g4sbstree_mtpc *T;
  T = new g4sbstree_mtpc(C);
  // g4sbs_tree_mtpc *T;
  // T = new g4sbs_tree_mtpc(C); 
  
  Int_t NEventsChain = C->GetEntries();
  cout << "N events in chain is " << NEventsChain << endl << endl;
  
  
  // Text file for hit outputs
  TString txtfile = infile;
  txtfile += ".txt";
  ofstream TxtFileOut;
  TxtFileOut.open (txtfile);
  
  //============================================================================
  // HISTOGRAM DEFINITIONS
  TH2D *hEdepVRHit = new TH2D("hEdepHitvR","hEdepHitvR",30,0.0,15.0,100, 0.0, 50.0);
  TH1D *hEDepHit = new TH1D("hEDepHit","hEDepHit",100,0.0,100.0);
  TH1D *hEDepTrack = new TH1D("hEDepTrack","hEDepTrack",400,0.0,800.0);
  TH1D *hNHitsPerEvent = new TH1D("hNHitsPerEvent","hNHitsPerEvent",100,0.0,100.0);
  TH1D *hPathLengthPerHit = new TH1D("hPathLengthPerHit","hPathLengthPerHit",100,0.0,16.0);
  TH1D *hChargePerHit = new TH1D("hChargePerHit","hChargePerHit",500,0.0,500.0);
  TH1D *hnEPerHit = new TH1D("hnEPerHit","hnEPerHit",100,0.0,1000.0);
  TH1D *hDriftTime = new TH1D("hDriftTime","hDriftTime",100,0.0,10.0);
  TH1D *hTimeWindow = new TH1D("hTimeWindow","hTimeWindow",50,0.0,100.0);
  TH1D *hStartTime = new TH1D("hStartTime","hStartTime",200,0.0,2000.0);
  TH1D *hChargeCloudDimemsion = new TH1D("hChargeCloudDimemsion","hChargeCloudDimemsion",100,0.0,10.0);
  TH1D *hChargeCloudArea = new TH1D("hChargeCloudArea","hChargeCloudArea",100,0.0,10.0);
  TH1D *hPadMultiplicityPerHit = new TH1D("hPadMultiplicityPerHit","hPadMultiplicityPerHit",100,0.0,100.0);
  TH1D *hHitBinArea = new TH1D("hHitBinArea","hHitBinArea",500,0.0,50.0);
  TH1D *hMomentum = new TH1D("hMomentum","hMomentum",250,0.0,500.0);
  TH1D *hTheta = new TH1D("hTheta","hTheta",360,0.0,360.0);
  TH1D *hPhi = new TH1D("hPhi","hPhi",360,-180.0,180.0);
  TH1D *hDriftTimeGEM1 = new TH1D("hDriftTimeGEM1","hDriftTimeGEM1",100,0.0,10.0);
  TH1D *hDriftTimeGEM2 = new TH1D("hDriftTimeGEM2","hDriftTimeGEM2",100,0.0,10.0);
  TH1D *hTimeWindowGEM1 = new TH1D("hTimeWindowGEM1","hTimeWindowGEM1",50,0.0,100.0);
  TH1D *hStartTimeGEM1 = new TH1D("hStartTimeGEM1","hStartTimeGEM1",200,0.0,2000.0);
  TH1D *hTimeWindowGEM2 = new TH1D("hTimeWindowGEM2","hTimeWindowGEM2",50,0.0,100.0);
  TH1D *hStartTimeGEM2 = new TH1D("hStartTimeGEM2","hStartTimeGEM2",200,0.0,2000.0);
  TH1D *hChargeCloudDimemsionGEM1 = new TH1D("hChargeCloudDimemsionGEM1","hChargeCloudDimemsionGEM1",100,0.0,10.0);
  TH1D *hChargeCloudAreaGEM1 = new TH1D("hChargeCloudAreaGEM1","hChargeCloudAreaGEM1",100,0.0,10.0);
  TH1D *hChargeCloudDimemsionGEM2 = new TH1D("hChargeCloudDimemsionGEM2","hChargeCloudDimemsionGEM2",100,0.0,10.0);
  TH1D *hChargeCloudAreaGEM2 = new TH1D("hChargeCloudAreaGEM2","hChargeCloudAreaGEM2",100,0.0,10.0);

  // sampa histos
  TH1D *hSampleLength = new TH1D("hSampleLength","hSampleLength",50,0.0,1000.0);
  TH1D *hNSamples = new TH1D("hNSamples","hNSamples",100,0.0,100.0);
  TH1D *hNSamplesBefore = new TH1D("hNSamplesBefore","hNSamplesBefore",10,0.0,10.0);
  TH1D *hNSamplesAfter = new TH1D("hNSamplesAfter","hNSamplesAfter",10,0.0,10.0);
  TH2D *hNSamplesBeforeVStartTime = new TH2D("hNSamplesBeforeVStartTime","hNSamplesBeforeVStartTime",200,0.0,2000.0,10,0.0,10.0);
  TH2D *hStartTimeVFracZDrift = new TH2D("hStartTimeVFracZDrift","hStartTimeVFracZDrift",100,0.0,1.0,100,0.0,2000.0);

  
  //============================================================================
  // EVENT LOOP
  
  // Variables outside loop
  Int_t RateInc = -1;
  vector<Double_t> vSummedSampaPulseAmps;
  vector<Double_t> vSummedSampaPulseTimes;
  Double_t Period = 1.0/3.0e6;
  //------------------------------------------------------------ start event loop
  // while ( T->GetEntry(nevent++) ){
  //for(Int_t i=5; i<20;i++){//7; i++){//11
  for(Int_t i=0; i<NEventsChain; i++){
    cout << endl;
    cout << "On event " << i << " / " << NEventsChain << endl;
    T->GetEntry(i);
    if (EventCounter % 1000000 == 0) cout << EventCounter << "/" << NEventsChain << endl;

    TxtFileOut << "Event " << i << "\n";
    
    // variables re set every event
    Double_t EDepSumHits = 0.0;
    Int_t protonhitcounter = 0;
    Int_t zplane = 0;

    // mTPC hits
    Int_t mTPC_nHitsPerEvent = T->SBS_mTPC_hit_nhits;
    // cout << "mTPC_nHitsPerEvent " << mTPC_nHitsPerEvent << endl;
    hNHitsPerEvent->Fill(mTPC_nHitsPerEvent);
    
    
    //-------------------------------------------- hitnumber of particles loop
    Int_t mTPCnPart = T->SBS_mTPC_npart_mTPC;
    Int_t mTPCPartnPart = T->SBS_mTPC_part_npart;
    // cout << "mTPCnPart " << mTPCnPart << endl;//" mTPCPartnPart " << mTPCPartnPart << endl;
    bool printmom = true;
    for(Int_t particle=0; particle<mTPCnPart; particle++){
    // for(Int_t particle=0; particle<1; particle++){
      // cout << endl << "i " << particle << endl;
      
      Int_t ParticleHitID = (*(T->SBS_mTPC_ihit))[particle];
      Int_t mTPCpid =  (*(T->SBS_mTPC_pid))[particle];
      Int_t mTPCtid =  (*(T->SBS_mTPC_trid))[particle];
      Int_t mTPCmid =  (*(T->SBS_mTPC_mid))[particle];
      // cout << "pid " << mTPCpid << " tid " << mTPCtid << " mid " << mTPCmid << endl;
      Double_t steplength = (*(T->SBS_mTPC_L))[particle];
      Double_t HitTime = (*(T->SBS_mTPC_t))[particle];
      Double_t Px = (*(T->SBS_mTPC_px))[particle];
      Double_t Py = (*(T->SBS_mTPC_py))[particle];
      Double_t Pz = (*(T->SBS_mTPC_pz))[particle];
      Double_t mom = TMath::Sqrt(Px*Px + Py*Py + Pz*Pz);

      // cout << "mTPCpid " << mTPCpid << endl;
      // cout << "mTPCmid " << mTPCmid << endl;
      if(mTPCpid==2212 && mTPCmid==0){
	RateInc++;
	// cout << "in proton loop" << endl;
	Int_t hitinc = (*(T->SBS_mTPC_ihit))[particle];
	Double_t energydepositedbyhit;
	Double_t hitx;
	Double_t hity;
	Double_t hitz;
	Double_t hitr;

	// cout << "Momentum of track is " << mom*1.0e3 << " MeV/c" << endl;
	TLorentzVector vProton;
	vProton.SetPxPyPzE(Px,Py,Pz,(*(T->SBS_mTPC_E))[particle]);
	Double_t Theta = vProton.Theta()*TMath::RadToDeg();
	Double_t Phi = vProton.Phi()*TMath::RadToDeg();
	// cout << "Theta " << Theta << " phi " << Phi << endl;
	if(protonhitcounter==0){
	  hMomentum->Fill(mom*1.0e3);
	  hTheta->Fill(Theta);
	  hPhi->Fill(Phi);
	}// if first proton

	if(printmom){
	  TxtFileOut << " " << mom
		     << " " << Theta
		     << " " << Phi
		     << " " << T->ev_vz
		     << "\n";
	  printmom = false;
	}

	energydepositedbyhit = (*(T->SBS_mTPC_hit_sumedep))[hitinc];
	EDepSumHits += energydepositedbyhit;
	hitx =  (*(T->SBS_mTPC_hit_xhit))[hitinc];
	hity =  (*(T->SBS_mTPC_hit_yhit))[hitinc];
	hitz =  (*(T->SBS_mTPC_hit_zhit))[hitinc];
	hitr = TMath::Sqrt(hitx*hitx + hity*hity);

	// Get charge from edep and number of electrons
	// g4 energies are in gev
	vector<Double_t> vCharge = CalcCharge(energydepositedbyhit, steplength);
	Double_t hitcharge = vCharge[0]; // charge in fC
	Double_t nDeltaE = vCharge[1]; // number of electrons before gem
	Double_t nEAfterGem = vCharge[2]; // number after gem gain

	// Adding diffusion - loop over electrons and smear each one
	// ztogem must be in cm since diff coeffs in cm
	Double_t ztogem =  ZtoGEM(hitz*1000.0); // convert argument, from g4, to mm for fncn, return is mm
	ztogem = ztogem/10.0; // convert to cm
	// smear the transverse drift of each electron to get the charge cloud size
	Double_t sigmatransverse = DT * TMath::Sqrt(ztogem);
	// smear the longitudinal drift of each electron to get time window 
	Double_t sigmalongitudinal = DL * TMath::Sqrt(ztogem); 
	Int_t nelectrons = floor(nDeltaE);
	// Int_t nelectrons = floor(nEAfterGem);
	// cout << "Number of electrons is " << nelectrons << endl;
	vector<Double_t> vSmearedDriftTimes;
	vector<Double_t> vSmearedZDrift;
	vector<Double_t> vSmearedZDriftFrac;
	vector<Double_t> vSmearedXHit;
	vector<Double_t> vSmearedYHit;
	vector<Double_t> vSmearedRHit;
	vector<Int_t> vSmearedHitPads;
	vector<Int_t> vSmearedHitPlanes;

	for(Int_t electronInc=0; electronInc<nelectrons; electronInc++){
	  // cout << "In electron loop" << endl;

	  // time smearing
	  Double_t TimeSmearDist = RandomVar.Gaus(0.0,sigmalongitudinal); //cm
	  Double_t SmearedDistToGEM = ztogem + TimeSmearDist; //cm
	  vSmearedZDrift.push_back(SmearedDistToGEM);
	  vSmearedZDriftFrac.push_back(SmearedDistToGEM/5.0);
	  if(SmearedDistToGEM<0.0) SmearedDistToGEM = 0.0;
	  Double_t SmearedTimeDriftToGEM = HitTime + (SmearedDistToGEM*10000.0)/DriftVelocity;
	  vSmearedDriftTimes.push_back(SmearedTimeDriftToGEM);
	  hDriftTime->Fill(SmearedTimeDriftToGEM/1000.0); // convert to us

	  // position smearing, all should be in cm, must convert g4 from m
	  Double_t xe_smeared = (hitx*100.0) + RandomVar.Gaus(0.0,sigmatransverse);
	  vSmearedXHit.push_back(xe_smeared);
	  Double_t ye_smeared = (hity*100.0) + RandomVar.Gaus(0.0,sigmatransverse);
	  vSmearedYHit.push_back(ye_smeared);
	  vSmearedRHit.push_back(TMath::Sqrt(xe_smeared*xe_smeared + ye_smeared*ye_smeared));

	  // get hit pads
	  // g4 out is in m, smearing is cm, so convert cm to mm in end
	  // vector<int> smearedhitpad;
	  // smearedhitpad = GetPadLoaded(xe_smeared*10.0,
	  // 			       ye_smeared*10.0,
	  // 			       hitz*1000.0,
	  // 			       hReadoutPads);
	  // vSmearedHitPads.push_back(smearedhitpad[0]);
	  // vSmearedHitPlanes.push_back(smearedhitpad[3]);

	  // text file for electrons
	  // myfile << xe_smeared-(hitx*100.0) << "\t"
	  // 	 << ye_smeared-(hity*100.0) << "\t"
	  // 	 << SmearedDistToGEM << "\t"
	  // 	 << SmearedTimeDriftToGEM << "\n";

	}// end electron loop
	// myfile.close();

	// if we had electrons

	// variables for smearing by gems
	// gem 1
	vector<Double_t> vSmearedZGem1;
	vector<Double_t> vSmearedZFracGem1;
	vector<Double_t> vSmearedTimeGem1;
	vector<Double_t> vSmearedXGem1;
	vector<Double_t> vSmearedYGem1;
	vector<Double_t> vSmearedRGem1;
	// gem 2
	vector<Double_t> vSmearedZGem2;
	vector<Double_t> vSmearedZFracGem2;
	vector<Double_t> vSmearedTimeGem2;
	vector<Double_t> vSmearedXGem2;
	vector<Double_t> vSmearedYGem2;
	vector<Double_t> vSmearedRGem2;

	if(nelectrons>0){
	  // loop through the electrons for the hit
	  for(Int_t HitElectronInc=0; HitElectronInc<nelectrons; HitElectronInc++){
	    // loop through the edep hit electrons and add gem gain, smearing for each
	    // gem layer and each electron after gain factor
	    for(Int_t GemGainInc=0; GemGainInc<GEMGain; GemGainInc++){

	      // GEM LAYER 1
	      // time smearing, in ns
	      Double_t TimeSmearDistGEM1 = RandomVar.Gaus(0.0,sigmalongitudinal_GEM); //cm
	      Double_t SmearedDistToGEM1 = GEMGap + TimeSmearDistGEM1; //cm
	      vSmearedZGem1.push_back(SmearedDistToGEM1);
	      Double_t SmearedFracDistToGEM1 = SmearedDistToGEM1/GEMGap;
	      vSmearedZFracGem1.push_back(SmearedFracDistToGEM1);
	      Double_t SmearedTimeToGEM1 = vSmearedDriftTimes[HitElectronInc] + (SmearedDistToGEM1*10000.0)/DriftVelocity;
	      vSmearedTimeGem1.push_back(SmearedTimeToGEM1);
	      // position smearing in cm
	      Double_t xSmearedGEM1 = vSmearedXHit[HitElectronInc] + RandomVar.Gaus(0.0,sigmatransverse_GEM);
	      Double_t ySmearedGEM1 = vSmearedYHit[HitElectronInc] + RandomVar.Gaus(0.0,sigmatransverse_GEM);
	      vSmearedXGem1.push_back(xSmearedGEM1);
	      vSmearedYGem1.push_back(ySmearedGEM1);
	      vSmearedRGem1.push_back(TMath::Sqrt(xSmearedGEM1*xSmearedGEM1 + ySmearedGEM1*ySmearedGEM1));
	      // histos
	      hDriftTimeGEM1->Fill(SmearedTimeToGEM1/1000.0); // convert to us

	      // GEM LAYER 2
	      // time smearing, in ns
	      Double_t TimeSmearDistGEM2 = RandomVar.Gaus(0.0,sigmalongitudinal_GEM); //cm
	      Double_t SmearedDistToGEM2 = GEMGap + TimeSmearDistGEM2; //cm
	      vSmearedZGem2.push_back(SmearedDistToGEM2);
	      Double_t SmearedFracDistToGEM2 = SmearedDistToGEM2/GEMGap;
	      vSmearedZFracGem2.push_back(SmearedFracDistToGEM2);
	      Double_t SmearedTimeToGEM2 = SmearedTimeToGEM1 + (SmearedDistToGEM2*10000.0)/DriftVelocity;
	      vSmearedTimeGem2.push_back(SmearedTimeToGEM2);
	      // position smearing in cm
	      Double_t xSmearedGEM2 = xSmearedGEM1 + RandomVar.Gaus(0.0,sigmatransverse_GEM);
	      Double_t ySmearedGEM2 = ySmearedGEM1 + RandomVar.Gaus(0.0,sigmatransverse_GEM);
	      vSmearedXGem2.push_back(xSmearedGEM2);
	      vSmearedYGem2.push_back(ySmearedGEM2);
	      vSmearedRGem2.push_back(TMath::Sqrt(xSmearedGEM2*xSmearedGEM2 + ySmearedGEM2*ySmearedGEM2));
	      // histos
	      hDriftTimeGEM2->Fill(SmearedTimeToGEM2/1000.0); // convert to us

	    } // gem gain loop
	  } // original electrons loop

	  
	  // for evaluating electrons direct from hits and no electrons from gem gain
	  // get the start time from the shortest drift time
	  Int_t SizevSmearedDriftTimes = vSmearedDriftTimes.size();
	  vector<Double_t> vSmearedDriftTimesSorted(SizevSmearedDriftTimes);
	  copy(vSmearedDriftTimes.begin(),vSmearedDriftTimes.end(), vSmearedDriftTimesSorted.begin());
	  sort(vSmearedDriftTimesSorted.begin(),vSmearedDriftTimesSorted.end());
	  Double_t StartTime = vSmearedDriftTimesSorted[0];
	  Double_t TimeWindow = vSmearedDriftTimesSorted[SizevSmearedDriftTimes-1] - StartTime;
	  
	  // get z distance fraction for shortest drift time
	  Int_t SizevSmearedZDriftFrac = vSmearedZDriftFrac.size();
	  vector<Double_t> vSmearedZDriftFracSorted(SizevSmearedZDriftFrac);
	  copy(vSmearedZDriftFrac.begin(),vSmearedZDriftFrac.end(), vSmearedZDriftFracSorted.begin());
	  sort(vSmearedZDriftFracSorted.begin(),vSmearedZDriftFracSorted.end());
	  Double_t ZFracStart = vSmearedZDriftFracSorted[0];
	  hStartTimeVFracZDrift->Fill(ZFracStart,StartTime);
	  
	  // get charge cloud size
	  //smeared is in cm
	  Int_t SizevSmearedRHit = vSmearedRHit.size();
	  vector<Double_t> vSmearedRHitSorted(SizevSmearedRHit);
	  copy(vSmearedRHit.begin(),vSmearedRHit.end(), vSmearedRHitSorted.begin());
	  sort(vSmearedRHitSorted.begin(),vSmearedRHitSorted.end());
	  Double_t ChargeDimension = vSmearedRHitSorted[SizevSmearedRHit-1] - vSmearedRHitSorted[0];
	  hChargeCloudDimemsion->Fill(ChargeDimension*10.0);// put in mm
	  // to get area: divide charge cloud dimension by 2 and convert to mm
	  Double_t ChargeCloudRadius = (ChargeDimension*10.0)/2.0;
	  Double_t ChargeCloudArea = TMath::Pi()*ChargeCloudRadius*ChargeCloudRadius;
	  hChargeCloudArea->Fill(ChargeCloudArea);
	  
	  // after gem 1
	  
	  // get the start time from the shortest drift time
	  Int_t SizevSmearedTimeGem1 = vSmearedTimeGem1.size();
	  vector<Double_t> vSmearedTimeGem1Sorted(SizevSmearedTimeGem1);
	  copy(vSmearedTimeGem1.begin(),vSmearedTimeGem1.end(), vSmearedTimeGem1Sorted.begin());
	  sort(vSmearedTimeGem1Sorted.begin(),vSmearedTimeGem1Sorted.end());
	  Double_t StartTimeGem1 = vSmearedTimeGem1Sorted[0];
	  Double_t TimeWindowGem1 = vSmearedTimeGem1Sorted[SizevSmearedTimeGem1-1] - StartTimeGem1;
	  hTimeWindowGEM1->Fill(TimeWindowGem1);
	  hStartTimeGEM1->Fill(StartTimeGem1);
	  
	  // get z distance fraction for shortest drift time
	  Int_t SizevSmearedZFracGem1 = vSmearedZFracGem1.size();
	  vector<Double_t> vSmearedZFracGem1Sorted(SizevSmearedZFracGem1);
	  copy(vSmearedZFracGem1.begin(),vSmearedZFracGem1.end(), vSmearedZFracGem1Sorted.begin());
	  sort(vSmearedZFracGem1Sorted.begin(),vSmearedZFracGem1Sorted.end());
	  Double_t ZFracStartGem1 = vSmearedZFracGem1Sorted[0];
	  
	  // get charge cloud size
	  //smeared is in cm
	  Int_t SizevSmearedRHitGem1 = vSmearedRGem1.size();
	  vector<Double_t> vSmearedRGem1Sorted(SizevSmearedRHitGem1);
	  copy(vSmearedRGem1.begin(),vSmearedRGem1.end(), vSmearedRGem1Sorted.begin());
	  sort(vSmearedRGem1Sorted.begin(),vSmearedRGem1Sorted.end());
	  Double_t ChargeDimensionGem1 = vSmearedRGem1Sorted[SizevSmearedRHitGem1-1] - vSmearedRGem1Sorted[0];
	  hChargeCloudDimemsionGEM1->Fill(ChargeDimensionGem1*10.0);// put in mm
	  // to get area: divide charge cloud dimension by 2 and convert to mm
	  Double_t ChargeCloudRadiusGem1 = (ChargeDimensionGem1*10.0)/2.0;
	  Double_t ChargeCloudAreaGem1 = TMath::Pi()*ChargeCloudRadiusGem1*ChargeCloudRadiusGem1;
	  hChargeCloudAreaGEM1->Fill(ChargeCloudAreaGem1);
	  
	  
	  // after gem 2
	  
	  // get the start time from the shortest drift time
	  Int_t SizevSmearedTimeGem2 = vSmearedTimeGem2.size();
	  vector<Double_t> vSmearedTimeGem2Sorted(SizevSmearedTimeGem2);
	  copy(vSmearedTimeGem2.begin(),vSmearedTimeGem2.end(), vSmearedTimeGem2Sorted.begin());
	  sort(vSmearedTimeGem2Sorted.begin(),vSmearedTimeGem2Sorted.end());
	  Double_t StartTimeGem2 = vSmearedTimeGem2Sorted[0];
	  Double_t TimeWindowGem2 = vSmearedTimeGem2Sorted[SizevSmearedTimeGem2-1] - StartTimeGem2;
	  hTimeWindowGEM2->Fill(TimeWindowGem2);
	  hStartTimeGEM2->Fill(StartTimeGem2);
	  
	  // get z distance fraction for shortest drift time
	  Int_t SizevSmearedZFracGem2 = vSmearedZFracGem2.size();
	  vector<Double_t> vSmearedZFracGem2Sorted(SizevSmearedZFracGem2);
	  copy(vSmearedZFracGem2.begin(),vSmearedZFracGem2.end(), vSmearedZFracGem2Sorted.begin());
	  sort(vSmearedZFracGem2Sorted.begin(),vSmearedZFracGem2Sorted.end());
	  Double_t ZFracStartGem2 = vSmearedZFracGem2Sorted[0];
	  
	  // get charge cloud size
	  //smeared is in cm
	  Int_t SizevSmearedRHitGem2 = vSmearedRGem2.size();
	  vector<Double_t> vSmearedRGem2Sorted(SizevSmearedRHitGem2);
	  copy(vSmearedRGem2.begin(),vSmearedRGem2.end(), vSmearedRGem2Sorted.begin());
	  sort(vSmearedRGem2Sorted.begin(),vSmearedRGem2Sorted.end());
	  Double_t ChargeDimensionGem2 = vSmearedRGem2Sorted[SizevSmearedRHitGem2-1] - vSmearedRGem2Sorted[0];
	  hChargeCloudDimemsionGEM2->Fill(ChargeDimensionGem2*10.0);// put in mm
	  // to get area: divide charge cloud dimension by 2 and convert to mm
	  Double_t ChargeCloudRadiusGem2 = (ChargeDimensionGem2*10.0)/2.0;
	  Double_t ChargeCloudAreaGem2 = TMath::Pi()*ChargeCloudRadiusGem2*ChargeCloudRadiusGem2;
	  hChargeCloudAreaGEM2->Fill(ChargeCloudAreaGem2);

	  // // get pads hit after gem 2
	  // // cout << "enterring loop to study hit pads" << endl;
	  // vector <int> vHitBin;
	  // vector <int> vRingID;
	  // vector <int> vPadID;
	  // vector <int> vPlane;
	  // vector<int> smearedhitpadsaftergem2;
	  // // cout << " got " << vSmearedRGem2.size() << " electrons" << endl;
	  // // *************************************************************************
	  // // this method of getting the hit pad for every electron is
	  // // really slowing things down!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  // for(int gem2inc=0; gem2inc<vSmearedRGem2.size(); gem2inc++){
	  //   // cout << "on electron " << gem2inc << endl;
	  //   smearedhitpadsaftergem2.clear();
	  //   smearedhitpadsaftergem2 = GetPadLoaded(vSmearedXGem2[gem2inc]*10.0,
	  // 					   vSmearedYGem2[gem2inc]*10.0,
	  // 					   hitz*1000.0,
	  // 					   hReadoutPads);
	  //   if(smearedhitpadsaftergem2.size()==4)
	  //     {
	  // 	vHitBin.push_back(smearedhitpadsaftergem2[0]);
	  // 	vRingID.push_back(smearedhitpadsaftergem2[1]);
	  // 	vPadID.push_back(smearedhitpadsaftergem2[2]);
	  // 	vPlane.push_back(smearedhitpadsaftergem2[3]);
	  //   }
	  // }
	  // // sort and find out how many pads are hit
	  // cout << "vHitBin size before " << vHitBin.size() << endl;
	  // Int_t SizevHitBin = vHitBin.size();
	  // vector<Int_t> vHitBinCopy(SizevHitBin);
	  // copy(vHitBin.begin(),vHitBin.end(), vHitBinCopy.begin());
	  // sort(vHitBinCopy.begin(),vHitBinCopy.end());
	  // vHitBinCopy.erase(unique(vHitBinCopy.begin(),vHitBinCopy.end()),vHitBinCopy.end());
	  // cout << "size vHitBin after removing duplicates " << vHitBinCopy.size();
	  // for(int&x:vHitBinCopy)
	  //   cout << x << " " << endl;
	  // // for(int y=0; y<vHitBinCopy.size(); y++){
	  // //   cout << "vHitBinCopy[y] " << vHitBinCopy[y] << endl;
	  // // }
	  
	  // get pad hit
	  // get mean pad hit
	  // original mtpc elecs
	  // Double_t xhitsmearedmeanstart = 0.0;
	  // Double_t xhitsmearedmeansum = accumulate(vSmearedXHit.begin(), vSmearedXHit.end(), xhitsmearedmeanstart);
	  // Double_t xhitsmearedmean = xhitsmearedmeansum/(vSmearedXHit.size()*1.0); // in cm
	  // Double_t yhitsmearedmeanstart = 0.0;
	  // Double_t yhitsmearedmeansum = accumulate(vSmearedYHit.begin(), vSmearedYHit.end(), yhitsmearedmeanstart);
	  // Double_t yhitsmearedmean = yhitsmearedmeansum/(vSmearedYHit.size()*1.0); // in cm

	  // after gem 2
	  Double_t xhitsmearedmeanstart = 0.0;
	  Double_t xhitsmearedmeansum = accumulate(vSmearedXGem2.begin(), vSmearedXGem2.end(), xhitsmearedmeanstart);
	  Double_t xhitsmearedmean = xhitsmearedmeansum/(vSmearedXGem2.size()*1.0); // in cm
	  Double_t yhitsmearedmeanstart = 0.0;
	  Double_t yhitsmearedmeansum = accumulate(vSmearedYGem2.begin(), vSmearedYGem2.end(), yhitsmearedmeanstart);
	  Double_t yhitsmearedmean = yhitsmearedmeansum/(vSmearedYGem2.size()*1.0); // in cm
	  
	  // get pad loaded script deals in mm
	  // g4 out is in m, smearing is cm, so convert cm to mm in end
	  vector<int> meansmearedhitpad;
	  std::chrono::time_point<std::chrono::steady_clock> start_t0 = 
	    std::chrono::steady_clock::now();
	  meansmearedhitpad = GetPadLoaded(xhitsmearedmean*10.0,
					   yhitsmearedmean*10.0,
					   hitz*1000.0,
					   hReadoutPads, 
					   t1, t2);
	  std::chrono::time_point<std::chrono::steady_clock> end_t0 = 
	    std::chrono::steady_clock::now();
	  std::chrono::duration<double> diff_t0 = end_t0-start_t0;
	  t0+= diff_t0.count();
	  // Get the area of the hit bin in mm
	  Double_t BinArea = vbinareas[meansmearedhitpad[0]-1];
	  hHitBinArea->Fill(BinArea);
	  // cout << "Hit bin " << meansmearedhitpad[0]
	  //      << " hit ring " << meansmearedhitpad[1]
	  //      << " hit pad " << meansmearedhitpad[2]
	  //      << " hit plane " << meansmearedhitpad[3]
	  //      << " hit bin area " << BinArea << endl;
	  
	  // does charge cloud fit in one pad?
	  // Int_t NPads = ceil(ChargeCloudArea/BinArea); // using only mtpc smear b4 gem
	  Int_t NPads = ceil(ChargeCloudAreaGem2/BinArea); // using only mtpc smear b4 gem
	  hPadMultiplicityPerHit->Fill(NPads);
	  
	  // if(protonhitcounter==0) zplane = meansmearedhitpad[3];//hitpad[3];
	  // // if(zplane != hitpad[3]) hReadoutPads->Fill(hitx*1000.0, hity*1000.0, Smeared_TimeDriftToGEM);
	  // if(zplane == meansmearedhitpad[3])
	  // hReadoutPads->Fill(xhitsmearedmean*10.0, yhitsmearedmean*10.0, StartTimeGem2);
	  
	  // charge in fC and time in ns for sampa analysis
	  vector<Double_t> vSampaPulseValues;
	  // 2 vals per event in return vector. first is time, second is value/charge
	  // vSampaPulseValues = SampaShapingAndSample_ReturnVals(80.0, // ShapingTime
	  // 						       11.86, // Gain
	  // 						       80.0, // baseline
	  // 						       hitcharge, // input charge
	  // 						       TimeWindow, // time window of charge
	  // 						       StartTime/50.0, // charge start time
	  // 						       90.0, // threshold
	  // 						       3, // samples before
	  // 						       7); // samples after
	  vSampaPulseValues = SampaShapingAndSample_ReturnVals(80.0, // ShapingTime
							       11.86, // Gain
							       80.0, // baseline
							       hitcharge, // input charge
							       TimeWindowGem2, // time window of charge
							       StartTimeGem2/50.0, // charge start time
							       90.0, // threshold
							       3, // samples before
							       7); // samples after
	  // cout << "Sampa Pulse N samples " << vSampaPulseValues.size()/2.0 << endl;
	  

	  // get time stamp and amplitude of sampa peak
	  Double_t PeakSampleTime = 0.0;
	  Double_t PeakAmplitude = 0.0;
	  Int_t PeakFlag = -1;
	  Int_t PeakInc = -1;
	  for(Int_t i=0; i<vSampaPulseValues.size(); i++){
	    if(i%2==1){// this is second in a pair so amp
	      if(vSampaPulseValues[i]>PeakAmplitude){
		PeakAmplitude = vSampaPulseValues[i];
		PeakSampleTime = vSampaPulseValues[i-1];
	      }
	    }
	  }
	  // convert peak sample time for 50ns binning;
	  PeakSampleTime = PeakSampleTime*50.0;//e-9;
	  if(meansmearedhitpad[1]>=0 && meansmearedhitpad[2]>=0 && meansmearedhitpad[3]>=0)
	  TxtFileOut << PeakSampleTime << "\t"
		     << PeakAmplitude  << "\t"
		     << ((meansmearedhitpad[3] << 18) | (meansmearedhitpad[0]-1)) << "\t"
		     << meansmearedhitpad[1] << "\t"
		     << meansmearedhitpad[2] << "\t"
		     << meansmearedhitpad[3] << "\t"
		     << endl;

	  // // Pile up study
	  // // Double_t PileUpStartTime;
	  // // Double_t Period = (1.0/3.0e6)/1.0e-9;// convert to ns
	  // Double_t Period_50ns = (1.0/6.0e6)/50.0e-9;// convert to 50ns binning

	  // // Get the sample information
	  // Double_t StartSampleTime = 0.0;
	  // Double_t StopSampleTime = 0.0;
	  // Int_t AboveThreshFlag = -1;
	  // Int_t AboveThreshInc = -1;
	  // Double_t AboveThreshVal = 0.0;
	  // Int_t BelowThreshFlag = -1;
	  // Int_t BelowThreshInc = -1;
	  // Double_t BelowThreshVal = 0.0;
	  // cout << endl;
	  // for(Int_t i=0; i<vSampaPulseValues.size(); i++){
	  //   if(i%2==0){// this is first in a pair so time
	  //     // cout << "RateInc " << RateInc << endl;
	  //     // cout << "vSampaPulseValues[i] " << vSampaPulseValues[i] << endl;
	  //     // cout << "(RateInc*Period_50ns) + (vSampaPulseValues[i])" << (RateInc*Period_50ns) + (vSampaPulseValues[i]) << endl;
	  //     vSummedSampaPulseTimes.push_back((RateInc*Period_50ns) + (vSampaPulseValues[i]));
	  //     if(i==0){
	  // 	StartSampleTime = vSampaPulseValues[i];
	  // 	// cout << "on first hit" << endl;
	  //     }
	  //     if(i==(vSampaPulseValues.size()-2)){
	  // 	// cout << "on last hit" << endl;
	  // 	StopSampleTime = vSampaPulseValues[i];
	  //     }
	  //   }
	  //   if(i%2==1){// this is second in a pair so amp
	  //     vSummedSampaPulseAmps.push_back(vSampaPulseValues[i]);
	  //     if(vSampaPulseValues[i]>90.0 && AboveThreshFlag<0){
	  // 	AboveThreshInc = i;
	  // 	AboveThreshFlag = 999;
	  // 	AboveThreshVal = vSampaPulseValues[i];
	  //     }
	  //     if(vSampaPulseValues[i]<90.0 && AboveThreshFlag>0 && BelowThreshFlag<0){
	  // 	BelowThreshInc = i;
	  // 	BelowThreshFlag = 999;
	  // 	BelowThreshVal = vSampaPulseValues[i];
	  //     }
	  //   }
	  // }// loop over sample points
	  // //calculate time period for sample points, multiply by 50ns as that is binning
	  // // in sampa function/macro
	  // Double_t nSamplesLength = (StopSampleTime*50.0) - (StartSampleTime*50.0);
	  
	  
	  //fill histos
	  hEdepVRHit->Fill(hitr*100.0, energydepositedbyhit*1000000.0); //convert to keV
	  hEDepHit->Fill(energydepositedbyhit*1000000.0);// convert to keV
	  hPathLengthPerHit->Fill(steplength*100.0); // convert to cm
	  hChargePerHit->Fill(hitcharge); // in fC
	  hnEPerHit->Fill(nelectrons); // number
	  hTimeWindow->Fill(TimeWindow); // ns
	  hStartTime->Fill(StartTime); // ns
	  // hSampleLength->Fill(nSamplesLength); // ns
	  // hNSamplesBefore->Fill(AboveThreshInc/2); // number
	  // hNSamplesAfter->Fill((vSampaPulseValues.size()/2.0)-(BelowThreshInc/2)); // number
	  hNSamples->Fill(vSampaPulseValues.size()/2.0);
	  
	  // count protons
	  protonhitcounter++;
	  
	  // clear vectors
	  vCharge.clear();
	  vSmearedDriftTimes.clear();
	  vSmearedZDrift.clear();
	  vSmearedZDriftFrac.clear();
	  vSmearedXHit.clear();
	  vSmearedYHit.clear();
	  vSmearedRHit.clear();
	  vSmearedHitPads.clear();
	  vSmearedHitPlanes.clear();
	  vSmearedDriftTimesSorted.clear();
	  vSmearedRHitSorted.clear();
	  vSmearedZDriftFracSorted.clear();
	  vSmearedZGem1.clear();
	  vSmearedZFracGem1.clear();
	  vSmearedTimeGem1.clear();
	  vSmearedXGem1.clear();
	  vSmearedYGem1.clear();
	  vSmearedRGem1.clear();
	  vSmearedZGem2.clear();
	  vSmearedZFracGem2.clear();
	  vSmearedTimeGem2.clear();
	  vSmearedXGem2.clear();
	  vSmearedYGem2.clear();
	  vSmearedRGem2.clear();
	  vSmearedTimeGem1Sorted.clear();
	  vSmearedZFracGem1Sorted.clear();
	  vSmearedRGem1Sorted.clear();
	  vSmearedTimeGem2Sorted.clear();
	  vSmearedZFracGem2Sorted.clear();
	  vSmearedRGem2Sorted.clear();
	  
	}// if we had electrons
	
      }// if we are looking at the original TDIS protons
      
    }//particle loop
    
    hEDepTrack->Fill(EDepSumHits*1000000.0); //convert to keV


    EventCounter++;
  } //------------------------------------------------------------ end event loop
  
  // // pile up plot
  // Int_t nPoints = vSummedSampaPulseTimes.size();
  // TGraph *gPileUp = new TGraph(nPoints,
  // 			       &vSummedSampaPulseTimes[0],
  // 			       &vSummedSampaPulseAmps[0]);
  // TCanvas *cPileUp = new TCanvas();
  // cPileUp->cd();
  // gPileUp->SetTitle("");
  // gPileUp->GetXaxis()->SetTitle("Time [50ns binning]");
  // // gSampledPoints->GetXaxis()->SetLimits(StartViewWindow,
  // // 					    EndViewWindow);
  // // gSampledPoints->GetYaxis()->SetLimits(0.0,
  // // 					1100.0);
  // gPileUp->GetYaxis()->SetTitle("ADC Counts");
  // gPileUp->SetMarkerStyle(20);
  // gPileUp->Draw("ALP");

  cout << endl;

  //============================================================================

  // WRITE HISTOS
  // Create file in which to store histos
  TString soutfile = infile;
  soutfile += "_analysed.root";
  // TFile * f = new TFile("/scratch/Analysed_g4sbs_out_30-70Theta_AllPhi_Field_AllZ_60-400MeV_1kEvents.root","RECREATE");
  TFile * f = new TFile(soutfile,"RECREATE");

  hMomentum->GetXaxis()->SetTitle("Track Momentum [MeV/c]");
  hMomentum->Write();

  hTheta->GetXaxis()->SetTitle("Track Theta [#circ]");
  hTheta->Write();

  hPhi->GetXaxis()->SetTitle("Track Phi [#circ]");
  hPhi->Write();

  hEdepVRHit->GetXaxis()->SetTitle("Hit R [cm]");
  hEdepVRHit->GetYaxis()->SetTitle("Hit EDep [keV]");
  hEdepVRHit->Write();

  hEDepHit->GetXaxis()->SetTitle("Hit EDep [keV]");
  hEDepHit->Write();

  hEDepTrack->GetXaxis()->SetTitle("EDep Summed Over Track [keV]");
  hEDepTrack->Write();

  hNHitsPerEvent->GetXaxis()->SetTitle("Number Hits Per Track");
  hNHitsPerEvent->Write();

  hPathLengthPerHit->GetXaxis()->SetTitle("Path Length of Hit [cm]");
  hPathLengthPerHit->Write();

  hChargePerHit->GetXaxis()->SetTitle("Charge/Hit [fC]");
  hChargePerHit->Write();

  hnEPerHit->GetXaxis()->SetTitle("Number of Electrons/Hit (before GEM)");
  hnEPerHit->Write();

  hDriftTime->GetXaxis()->SetTitle("Electron Drift Times [microseconds]");
  hDriftTime->Write();

  hDriftTimeGEM1->GetXaxis()->SetTitle("Electron Times After 1st GEM [microseconds]");
  hDriftTimeGEM1->Write();

  hDriftTimeGEM2->GetXaxis()->SetTitle("Electron Times After 2nd GEM [microseconds]");
  hDriftTimeGEM2->Write();

  hTimeWindow->GetXaxis()->SetTitle("Time Window of Charge [ns]");
  hTimeWindow->Write();

  hTimeWindowGEM1->GetXaxis()->SetTitle("Time Window of Charge after 1st GEM [ns]");
  hTimeWindowGEM1->Write();

  hTimeWindowGEM2->GetXaxis()->SetTitle("Time Window of Charge after 2nd GEM [ns]");
  hTimeWindowGEM2->Write();

  hStartTime->GetXaxis()->SetTitle("Start Time of Charge [ns]");
  hStartTime->Write();

  hStartTimeGEM1->GetXaxis()->SetTitle("Start Time of Charge after 1st GEM [ns]");
  hStartTimeGEM1->Write();

  hStartTimeGEM2->GetXaxis()->SetTitle("Start Time of Charge after 2nd GEM [ns]");
  hStartTimeGEM2->Write();

  hChargeCloudDimemsion->GetXaxis()->SetTitle("Charge Cloud Diameter [mm]");
  hChargeCloudDimemsion->Write();

  hChargeCloudArea->GetXaxis()->SetTitle("Charge Cloud Area [mm^{2}]");
  hChargeCloudArea->Write();

  hChargeCloudDimemsionGEM1->GetXaxis()->SetTitle("Charge Cloud Diameter After 1st GEM [mm]");
  hChargeCloudDimemsionGEM1->Write();

  hChargeCloudAreaGEM1->GetXaxis()->SetTitle("Charge Cloud Area After 1st GEM [mm^{2}]");
  hChargeCloudAreaGEM1->Write();

  hChargeCloudDimemsionGEM2->GetXaxis()->SetTitle("Charge Cloud Diameter After 2nd GEM [mm]");
  hChargeCloudDimemsionGEM2->Write();

  hChargeCloudAreaGEM2->GetXaxis()->SetTitle("Charge Cloud Area After 2nd GEM [mm^{2}]");
  hChargeCloudAreaGEM2->Write();

  hPadMultiplicityPerHit->GetXaxis()->SetTitle("Pad Multiplicity Per Hit");
  hPadMultiplicityPerHit->Write();

  hHitBinArea->GetXaxis()->SetTitle("Hit Bin Area [mm^{2}]");
  hHitBinArea->Write();

  // sampa histos
  hSampleLength->GetXaxis()->SetTitle("Length of sampled pulse [ns]");
  hSampleLength->Write();

  hNSamplesBefore->GetXaxis()->SetTitle("Number of sampled points before threshold");
  hNSamplesBefore->Write();
  
  hNSamplesAfter->GetXaxis()->SetTitle("Number of sampled points after threshold");
  hNSamplesAfter->Write();

  hNSamples->GetXaxis()->SetTitle("Number of sampled points");
  hNSamples->Write();

  hStartTimeVFracZDrift->GetXaxis()->SetTitle("Fraction of cell Z-length drifted");
  hStartTimeVFracZDrift->GetYaxis()->SetTitle("Start Time [ns]");
  hStartTimeVFracZDrift->Write();

  // // readout pad histos
  // hReadoutPads->GetXaxis()->SetTitle("X [mm]");
  // hReadoutPads->GetXaxis()->SetTitle("Y [mm]");
  // hReadoutPads->Write();

  // CLOSE FILES
  f->Close();
  TxtFileOut.close();
  
  cout << " retrieve pad time " << std::setprecision(9) << t0 << " s "<< endl;
  cout << " tlist declaration time " << std::setprecision(9) << t1 << " s "<< endl;
  cout << " isinsidebin time " << std::setprecision(9) << t2 << " s "<< endl;
  
  std::chrono::time_point<std::chrono::steady_clock> end_all = 
    std::chrono::steady_clock::now();
  
  std::chrono::duration<double> diff_all = end_all-start_all;
  cout << " Total time " << std::setprecision(9) << diff_all.count() << " s "<< endl;
  
} //------------------------------------------------------------ end main function

vector<Double_t> CalcCharge(Double_t Edep, Double_t PathLength){
  // helium values from pdg table 33.5
  Double_t He_density = 0.179; // mg/cm^3
  Double_t He_Ex = 19.8; // eV
  Double_t He_Ei = 24.6; // eV
  Double_t He_Wi = 41.3; // eV
  Double_t He_dEdx_MIP = 0.32; // 0.32 keV/cm
  Double_t He_Np = 3.5; // cm^-1
  Double_t He_Nt = 8.0; // cm^-1
  // methane values from pdg table 33.5
  Double_t CH4_density = 0.667; // mg/cm^3
  Double_t CH4_Ex = 8.8; // eV
  Double_t CH4_Ei = 12.6; // eV
  Double_t CH4_Wi = 30.0; // eV
  Double_t CH4_dEdx_MIP = 1.61; // 0.32 keV/cm
  Double_t CH4_Np = 28.0; // cm^-1
  Double_t CH4_Nt = 54.0; // cm^-1

  Double_t He_Frac = 0.7;//0.9;
  Double_t CH4_Frac = 0.3;//0.1;
  Double_t Mix_Nt = (He_Frac*He_Nt) + (CH4_Frac*CH4_Nt); // cm^-1
  Double_t Mix_Wi = (He_Frac*He_Wi) + (CH4_Frac*CH4_Wi); // cm^-1
  // cout << "Wi of mix " << Mix_Wi << endl;

  Double_t PathLengthIncm = PathLength*100.0; // assume g4 in [m]
  // cout << "PathLengthincm " << PathLengthIncm << endl;

  Double_t EDepInkeV = Edep*1000000.0;
  // cout << "EDep in keV " << EDepInkeV << endl;

  Double_t EDep_keVPercm = EDepInkeV/PathLengthIncm;
  Double_t EDep_eVPercm = EDep_keVPercm*1000.0;
  Double_t NDeltaElecPercm = EDep_eVPercm/Mix_Wi; // this is per cm
  // cout << "NDeltaElec per cm " << NDeltaElecPercm << endl;

  Double_t NDeltaElec = NDeltaElecPercm * PathLengthIncm;
  // cout << "NDeltaElec per path length " << NDeltaElec << endl;

  Double_t GainDoubleGem = GEMGain;
  Double_t TotalElec = GainDoubleGem * NDeltaElec;

  Double_t ChargeElectronInfC = 1.60217662e-19/1.0e-15;
  // Double_t Charge = NDeltaElec * ChargeElectronInfC;
  Double_t Charge = TotalElec * ChargeElectronInfC;

  vector<Double_t> vreturn;
  vreturn.push_back(Charge);
  vreturn.push_back(NDeltaElec);
  vreturn.push_back(TotalElec);

  return vreturn;
}


Double_t ZtoGEM(Double_t zpos){ // it is input in m, from g4

  Double_t DistanceToGEM = -999.0;
  if(zpos>=-250.0 && zpos<=-200.0)
    DistanceToGEM = fabs(zpos - (-250.0));
  if(zpos>-200.0 && zpos<=-150.0)
    DistanceToGEM = fabs(zpos - (-150.0));
  if(zpos>-150.0 && zpos<=-100.0)
    DistanceToGEM = fabs(zpos - (-150.0));
  if(zpos>-100.0 && zpos<=-50.0)
    DistanceToGEM = fabs(zpos - (-50.0));
  if(zpos>-50.0 && zpos<=0.0)
    DistanceToGEM = fabs(zpos - (-50.0));
  if(zpos>0.0 && zpos<=50.0)
    DistanceToGEM = fabs(zpos - (50.0));
  if(zpos>50.0 && zpos<=100.0)
    DistanceToGEM = fabs(zpos - (50.0));
  if(zpos>100.0 && zpos<=150.0)
    DistanceToGEM = fabs(zpos - (150.0));
  if(zpos>150.0 && zpos<=200.0)
    DistanceToGEM = fabs(zpos - (150.0));
  if(zpos>200.0 && zpos<=250.0)
    DistanceToGEM = fabs(zpos - (250.0));
  
  return DistanceToGEM; // it is in mm
}

// incCell 0 Z of readout is -0.24987;
// incCell 1 Z of readout is -0.15013;
// incCell 2 Z of readout is -0.14987;
// incCell 3 Z of readout is -0.0501305;
// incCell 4 Z of readout is -0.0498695;
// incCell 5 Z of readout is 0.0498695;
// incCell 6 Z of readout is 0.0501305;
// incCell 7 Z of readout is 0.14987;
// incCell 8 Z of readout is 0.15013;
// incCell 9 Z of readout is 0.24987;
