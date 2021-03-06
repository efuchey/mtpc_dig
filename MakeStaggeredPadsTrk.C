
void Arc(int n, double a0, double a, double r, double *px, double *py);

TH2Poly* MakeStaggeredPadsTrk(){

  Double_t OuterRadius = 150.00;//mm
  Double_t InnerRadius = 50.00;
  cout << endl << "Inner and Outer radii are " << InnerRadius
       << " and " << OuterRadius << " mm " << endl;
  // Double_t RingWidth = 5.00;
  // cout << endl << "We have input Inner Radius " << InnerRadius
  //      << "mm, outer radius " << OuterRadius << "mm, and RingWidth "
  //      << RingWidth << "mm" << endl;

  // const Int_t nRings = ceil( (OuterRadius-InnerRadius) / RingWidth);
  // RingWidth = (OuterRadius-InnerRadius)/nRings;
  // cout << "The resulting number of rings is " << nRings <<
  //   " and the new ring width is " << RingWidth << endl;

  const Int_t nRings = 21;
  Double_t RingWidth = (OuterRadius-InnerRadius)/nRings;
  cout << "The number of rings is " << nRings <<
    " and the ring width is " << RingWidth << endl;

  Int_t nPads = 122;//120;
  cout << "We have input " << nPads << " pads per ring" << endl
       << endl;

  const Int_t NP = 4000;
  const Int_t nPforB = 2.0 * NP;
  TH2Poly *hPolyBins = new TH2Poly();
  TH2Poly *hPolyBinAreas = new TH2Poly();

  for(Int_t rInc=0; rInc<nRings; rInc++){
  // for(Int_t rInc=0; rInc<1; rInc++){
    // cout << "We are on radius loop inc " << rInc << endl << endl;

    Double_t RIn = InnerRadius + ((rInc*1.00)*RingWidth);
    Double_t CircumIn = TMath::Pi() * 2.00 * RIn;
    Double_t CPerPadIn = CircumIn/(nPads*1.000);
    // nPads = (CircumIn/CPerPadIn);
    // CPerPadIn = CircumIn/(nPads*1.00000);
    Double_t AngleSubtendedIn = CPerPadIn / RIn;
    // To get the angle for each rotation CPerPadIn = theta * RIn
    // cout << "Inner radius is " << RIn << "mm, with circumference "
    // 	 << CircumIn << "mm, and " << CPerPadIn << "mm length per pad"
    // 	 << " and angle subtended per segment " << AngleSubtendedIn
    // 	 << endl;

    Double_t ROut = InnerRadius + (((rInc+1)*1.00)*RingWidth);
    Double_t CircumOut = TMath::Pi() * 2.0 * ROut;
    Double_t CPerPadOut = CircumOut/(nPads*1.000);
    // nPads = CircumOut/CPerPadOut;
    // CPerPadOut = CircumOut/(nPads*1.00000);
    Double_t AngleSubtendedOut = CPerPadOut / ROut;
    // cout << "Outer radius is " << ROut << "mm, with circumference "
    // 	 << CircumOut << "mm, and " << CPerPadOut << "mm length per pad"
    // 	 << " and angle subtended per segment " << AngleSubtendedOut
    // 	 << endl << endl;

    for(Int_t pInc=0; pInc<nPads; pInc++){
    // for(Int_t pInc=0; pInc<2; pInc++){
      // cout << "On pad " << pInc << endl << endl;

      Double_t FirstAngleIn = -1.0*AngleSubtendedIn/2.0 + (pInc*1.00)*AngleSubtendedIn + (rInc*1.0)*AngleSubtendedIn/2.0;
      // cout << "First angle in " << FirstAngleIn << endl;

      Double_t FirstXIn = 0.0 + RIn*TMath::Sin(AngleSubtendedIn);
      Double_t FirstYIn = RIn - RIn*(1.0-TMath::Cos(AngleSubtendedIn));;
      // cout << "FirstXIn " << FirstXIn << " FirstYIn " << FirstYIn << endl;


      // FirstAngleIn = (pInc*1.00) * AngleSubtendedIn;
      Double_t pxIn[NP];
      Double_t pyIn[NP];
      Arc(NP, FirstAngleIn, AngleSubtendedIn, RIn, pxIn, pyIn);

      Double_t FirstAngleOut = -1.0*AngleSubtendedOut/2.0 + (pInc*1.00)*AngleSubtendedOut + (rInc*1.0)*AngleSubtendedOut/2.0;
      Double_t FirstXOut = 0.0 + ROut*TMath::Sin(AngleSubtendedOut);
      Double_t FirstYOut = ROut - ROut*(1.0-TMath::Cos(AngleSubtendedOut));

      Double_t pxOut[NP];
      Double_t pyOut[NP];
      Arc(NP, FirstAngleOut, AngleSubtendedOut, ROut, pxOut, pyOut);

      Double_t pxOutReverse[NP];
      Double_t pyOutReverse[NP];
      for(int pointInc=0; pointInc<NP; pointInc++){
	pxOutReverse[pointInc] = pxOut[NP-1-pointInc];
	pyOutReverse[pointInc] = pyOut[NP-1-pointInc];
      }

      Double_t AddPointX1 = pxIn[NP-1] + (pxOutReverse[0] - pxIn[NP-1])/2.0;
      Double_t AddPointY1 = pyIn[NP-1] + (pyOutReverse[0] - pyIn[NP-1])/2.0;
      Double_t AddPointX2 = pxIn[0] + (pxOutReverse[NP-1] - pxIn[0])/2.0;
      Double_t AddPointY2 = pyIn[0] + (pyOutReverse[NP-1] - pyIn[0])/2.0;

 
      Double_t pBinCoordsX[nPforB];//+2];
      Double_t pBinCoordsY[nPforB];//+2];  
      for(int pointInc2=0; pointInc2<NP; pointInc2++){
	pBinCoordsX[pointInc2] = pxIn[pointInc2];
	pBinCoordsY[pointInc2] = pyIn[pointInc2];
      }
      for(int pointInc3=NP; pointInc3<2*NP; pointInc3++){
	pBinCoordsX[pointInc3] = pxOutReverse[pointInc3-NP];
	pBinCoordsY[pointInc3] = pyOutReverse[pointInc3-NP];
      }


      hPolyBins->AddBin(nPforB,pBinCoordsX,pBinCoordsY);

    }//pad loop

  }//ring loop



  // gStyle->SetPalette(kGreyScale);
  // TColor::InvertPalette();
  TCanvas *cPolyBins = new TCanvas();
  cPolyBins->cd();
  hPolyBins->SetName("hPolyBins");
  hPolyBins->SetTitle("hPolyBins");
  hPolyBins->Draw();
  cout << "Number of pads is " << hPolyBins->GetNumberOfBins() << endl;

  // // draw bin areas
  // TList *bin_list=hPolyBins->GetBins();
  // TIter next(bin_list);
  // TObject *obj;
  // TH2PolyBin *bin;
  // double area;
  // vector<Double_t> vbinareas;
  // while((obj=next())){
  //   bin = (TH2PolyBin*) obj;
  //   vbinareas.push_back(bin->GetArea());
  //   bin->Fill(bin->GetArea());
  // }
  // // TCanvas *cPolyBinAreas = new TCanvas();
  // // cPolyBinAreas->cd();
  // hPolyBins->Draw("COLZ");



  return hPolyBins;

}//make pads

void Arc(int n, double a0, double a, double r, double *px, double *py) {
  // cout << endl;
  // cout << "n " << n << " a0 " << a0 << " a " << a << " r " << r << endl;
  // cout << "px[0] " << px[0] << " py[0] " << py[0] << endl;
  Double_t AngleStep = a/(1.00*n);
  // cout << "AngleStep " << AngleStep << endl;
  Double_t Angle = a0;
  // cout << "Angle " << Angle << endl;
  for (int i=0; i<n; i++) {//i started from 0
    Angle = a0 + ((i*1.00)*AngleStep);
    px[i] = r*TMath::Cos(Angle);//+ px[0];
    py[i] = r*TMath::Sin(Angle);//+ py[0];
    // cout << "i " << i << " Angle " << Angle << " px[i] " << px[i] << " py[i] " << py[i] << endl;
   }
}








