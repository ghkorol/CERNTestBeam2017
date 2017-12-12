//root
#include <TLine.h>
//#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>

//C, C++
#include <stdio.h>
#include <assert.h>
//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
//#include <string>
//#include <iomanip>

using namespace std;

float SP = 0.3125;
float pe = 47.46;//mV*ns
int wavesPrintRate = 1000000;
int ch0PrintRate = 1000000;
int trigPrintRate = 1000000;//100
int signalPrintRate = 100000;//100
double coef = 2.5 / (4096 * 10);

//Geometry
//functions
std::vector<float> getStartPos(float horisontal, float vertical, float angle);//return start position of track in the box (x,y,z,angle) [mm,mm,mm,rad]
std::vector<float> solidAngleFactor(const std::vector<float> &startPos,const std::vector<float> &pmtPos);//return solid angle correction factor and path length in the box (factor,length) [,mm]
float solidAngleABH(float A,float B, float H);//return solid angle of pyramid with rectangular base (a,b) and hight of H
int getZone(float x,float y,float z,float pmtX, float pmtY);
float getSolidAngle(float x,float y,float z,float pmtX, float pmtY);


//box dimensions
float boxL = 500;
float boxH = 250;
//wom dimensions
float womL = 200;//length
float womDout = 70;//outer diameter [mm]
float womDin = 42;//inner diameter [mm]
float womR = 30;//(inner radius = 25)//wom radius [mm]
float womRin = 26;
//wom position
std::vector<float> pmt2Pos = {410,410};



  int runNr = -999;
  float horizontal=-999;
  float vertical=-999;
  float angle=-999;
  int pdgID=-999;
  float energy = -999; // [GeV]
  int isSP = -999;
  int mp = -999;
  int safPMT2 = -999;//solid angle factor of pmt 2
  int safPMT1 = -999;//solid angle factor of pmt 1
  int safSiPM = -999;//solid angle factor of SiPM
  int trackL = -999;//track length
  
void cleanEventMemory(std::vector<TObject*>& trash);
float CDF(TH1F* hWave,TF1* fTrigFit,float thr);
float CDFinvert(TH1F* hWave,float thr);
float iCFD(TH1F* hWave,float t,float thrNpe,float BL);
float integral(TH1F* hWave,float t1,float t2,float BL);

float* getBL(TH1F* hWave, float* BL, float t1, float t2); 
TString getRunName(TString inDataFolder);
void read(TString inFileList, TString inDataFolder, TString outFile);

int main(int argc, char *argv[]){
  TString inFileList;
  TString inDataFolder;
  TString outFile;

  if(argc == 4){
    inFileList = argv[1];
    inDataFolder = argv[2];
    outFile = argv[3];
    cout<<"In data file list : "<<inFileList<<endl
        <<"In data path      : "<<inDataFolder<<endl
        <<"Out root file     : "<<outFile<<endl;
    read(inFileList, inDataFolder, outFile);
  }
  else if(argc == 13){
    inFileList = argv[1];
    inDataFolder = argv[2];
    outFile = argv[3];
    runNr=atoi(argv[4]);
    horizontal = atof(argv[5])/1000; //units: [cm]
    vertical = atof(argv[6])/1000; //units: [cm]
    angle = atof(argv[7]);
    //argv[8] - nevents from table
    pdgID=atoi(argv[9]);// - pdgID
    energy=atof(argv[10]);// - energy [GeV]
    isSP=atoi(argv[11]);// - isSteelPlate
    mp=atoi(argv[12]);//  measure-point


    cout<<"In data file list : "<<inFileList<<endl
        <<"In data path      : "<<inDataFolder<<endl
        <<"Out root file     : "<<outFile<<endl
	<<"Run number         : "<<runNr<<endl;
	printf("hor,ver,ang: %4.2f %4.2f %4.2f\n",horizontal,vertical,angle);
	std::vector<float> starPos = getStartPos(horizontal,vertical,angle);
	//printf("start pos: %4.2f %4.2f %4.2f %4.2f\n",starPos[0],starPos[1],starPos[2],starPos[3]);
	std::vector<float> saf = solidAngleFactor(starPos,pmt2Pos);
	//printf("length: %4.2f %4.2f \n", saf[1],10*calculateDistance(horizontal,angle*TMath::Pi()/180));
	printf("length: %4.2f \n", saf[1]);
	printf("saf: %4.2f \n", saf[0]);
	safPMT2 = saf[0];
	trackL = saf[1];
	
    read(inFileList, inDataFolder, outFile);

  }
  else{
    cout<<" ERROR --->  in input arguments "<<endl
        <<"        [1] - in data file list"<<endl
        <<"        [2] - in data path"<<endl
        <<"        [3] - out root file"<<endl;
  }
  return 0;
}

void read(TString _inFileList, TString _inDataFolder, TString _outFile){

  TF1* fTrigFit = new TF1("fTrigFit","gaus");
  fTrigFit->SetParameter(0,800);
  fTrigFit->SetParameter(2,1);
  fTrigFit->SetLineWidth(1);
  
  
  
  ///////////////////Root file with data/////////////////
  TFile *rootFile = new TFile( _outFile, "RECREATE");
  if (rootFile->IsZombie()) {
    cout << "PROBLEM with the initialization of the output ROOT ntuple "
	 << _outFile << ": check that the path is correct!!!"
	 << endl;
    exit(-1);
  }
  TTree *tree = new TTree("T", "USBWC Data Tree");
  //rootFile->SetCompressionLevel(2);
  //tree->SetAutoSave(1000000);
  // Create new event
  TTree::SetBranchStyle(0);

  Int_t EventNumber=-999;
  Int_t LastEventNumber=-999;
  Float_t SamplingPeriod = -999;
  Double_t EpochTime = -999;
  Int_t Year = -999;
  Int_t Month = -999;
  Int_t Day = -999;
  Int_t Hour = -999;
  Int_t Minute = -999;
  Int_t Second = -999;
  Int_t Millisecond = -999;
  Float_t trigT = -999;//t_trig = (t0+t1+t2+t3)/4
  Float_t tPMT1 = -999;
  Float_t tPMT2 = -999;
  Float_t tPMT2i = -999;
  Float_t tSUMp = -999;
  Float_t tSUMm = -999;
  Float_t tSiPM = -999;
  Float_t trigTp = -999;//t_trig' = [(t0+t1)-(t2+t3)]/4
  Float_t t0t1 = -999;//t0t1 = [(t0-t1)]
  Float_t t2t3 = -999;//t2t3 = [(t2-t3)]
  Int_t isVeto = -999; //variable to define veto, 1 if veto, 0 if not, -999 if undefined
  Int_t isTrig = -999;
  Int_t isLastEvt = -999;
  Int_t isGoodSignal_5 = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  std::vector<float> amp(16,-999);
  std::vector<float> max(16,-999);
  std::vector<float> min(16,-999);
  Float_t t[16];
  Float_t BL[16];//store baseline for 16 channels
  Float_t BL_RMS[16];//store rms of baseline for 16 channels
  float BL_output[2];//array used for output getBL-function
  float Integral_0_300[16];//array used to store Integral of signal from 0 to 300ns
  float Integral[16];
  int NumberOfBins;
  Int_t EventIDsamIndex[16];
  Int_t FirstCellToPlotsamIndex[16];
  
  std::vector<TH1F*> hChSum;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChSum_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChSum.push_back(h);
  }
  std::vector<TH1F*> hChShift;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChShift_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChShift.push_back(h);
  }
  std::vector<TH1F> hChtemp;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChtemp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChtemp.push_back(h);
  }
  std::vector<TH1F> hChShift_temp;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChShift_temp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChShift_temp.push_back(h);
  }
  
  Short_t amplValues[16][1024];
  TH1F hCh("hCh","dummy;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
  TString plotSaveFolder  = _outFile;
  plotSaveFolder.ReplaceAll("out.root","");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(4,4);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);

  tree->Branch("EventNumber",&EventNumber, "EventNumber/I");
  tree->Branch("SamplingPeriod", &SamplingPeriod,  "SamplingPeriod/F");
  tree->Branch("EpochTime",&EpochTime, "EpochTime/D");
  tree->Branch("Year",&Year, "Year/I");
  tree->Branch("Month",&Month, "Month/I");
  tree->Branch("Day",&Day, "Day/I");
  tree->Branch("Hour",&Hour, "Hour/I");
  tree->Branch("Minute",&Minute, "Minute/I");
  tree->Branch("Second",&Second, "Second_/I");
  tree->Branch("Millisecond",&Millisecond, "Millisecond/I");
  tree->Branch("trigT",&trigT, "trigT/F");
  tree->Branch("tPMT1",&tPMT1, "tPMT1/F");
  tree->Branch("tPMT2",&tPMT2, "tPMT2/F");
  tree->Branch("tPMT2i",&tPMT2i, "tPMT2i/F");
  tree->Branch("tSUMp",&tSUMp, "tSUMp/F");
  tree->Branch("tSUMm",&tSUMm, "tSUMm/F");
  tree->Branch("tSiPM",&tSiPM, "tSiPM/F");

  tree->Branch("runNr",&runNr, "runNr/I");//run number in google table
  tree->Branch("horiz",&horizontal,"horiz/F");// horizontal position of the box units: [cm]
  tree->Branch("vert",&vertical,"vert/F");//vertical position of the box, units: [cm]
  tree->Branch("angle",&angle,"angle/F");
  tree->Branch("pdgID",&pdgID,"pdgID/I");
  tree->Branch("energy",&energy,"energy/F");
  tree->Branch("isSP",&isSP,"isSP/I");
  tree->Branch("mp",&mp,"mp/I");
  tree->Branch("safPMT2",&safPMT2,"safPMT2/I");//solid angle factor
  tree->Branch("safPMT1",&safPMT1,"safPMT1/I");//solid angle factor
  tree->Branch("safSiPM",&safSiPM,"safSiPM/I");//solid angle factor
  tree->Branch("trackL",&trackL,"trackL/I");//track length
  tree->Branch("isLastEvt",&isLastEvt,"isLastEvt/I");
  
  tree->Branch("trigGate",&trigGate,"trigGate/F");
  tree->Branch("trigTp",&trigTp, "trigTp/F");
  tree->Branch("t0t1",&t0t1, "t0t1/F");//t0t1 = [(t0-t1)]
  tree->Branch("t2t3",&t2t3, "t2t3/F");
  tree->Branch("isVeto",&isVeto,"isVeto/I");
  tree->Branch("isTrig",&isTrig,"isTrig/I");
  tree->Branch("isGoodSignal_5",&isGoodSignal_5,"isGoodSignal_5/I");
  
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  tree->Branch("amp",amp.data(), "amp[nCh]/F");
  tree->Branch("max",max.data(), "max[nCh]/F");
  tree->Branch("min",min.data(), "min[nCh]/F");
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("BL", BL, "BL[nCh]/F");
  tree->Branch("BL_RMS", BL_RMS, "BL_RMS[nCh]/F");
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

// // //   tree->Branch("MeasuredBaseline_usbwc", MeasuredBaseline_usbwc, baseline_ss.Data());
// // //   tree->Branch("AmplitudeValue_usbwc", AmplitudeValue_usbwc, amplitude_ss.Data());
// // //   tree->Branch("ComputedCharge_usbwc", ComputedCharge_usbwc, charge_ss.Data());
// // //   tree->Branch("RiseTimeInstant_usbwc", RiseTimeInstant_usbwc, leadingEdgeTime_ss.Data());
// // //   tree->Branch("FallTimeInstant_usbwc", FallTimeInstant_usbwc, trailingEdgeTime_ss.Data());
// // //   tree->Branch("RawTriggerRate_usbwc", RawTriggerRate_usbwc, rateCounter_ss.Data());
  //tree->Branch("amplValues", amplValues, "amplValues[nCh][1024]/S");
 // tree->Branch("hCh","TH1F",&hCh,128000,1);
  ///////////////////////////////////////////////////////

    int nitem = 1;
    ifstream inList;
    TString fileName;
    inList.open(_inFileList);
    assert(inList.is_open());

    int wavePrintStatus=-1;
    int ch0PrintStatus=-1;
    int trigPrintStatus=-1;
    int signalPrintStatus=-1;
    while(inList >> fileName){
      fileName = _inDataFolder + fileName;
      cout << endl;
      cout << fileName << endl;
      FILE* pFILE = fopen(fileName.Data(),"rb");
      if (pFILE==NULL) {fputs ("File error",stderr); assert(0);}
      //cout<<" ---> File to convert : " << fileName << endl;
      fseek (pFILE , 0 , SEEK_END);
      int totFileSizeByte = ftell (pFILE);
      rewind (pFILE);
      cout<<"totFileSizeByte = "<<totFileSizeByte<<endl;
      char header[328];
      nitem=fread(header,1,328,pFILE);
      cout << "Header:\n" << header << endl;

      char* word;
      word = strtok(header," \n");
      while(word != NULL){
	  if(strcmp("ACQUIRED:",word) == 0){
	    word = strtok(NULL, " \n");
	    nActiveCh = atoi(word);
	    break;
	  }
	 //printf ("%s\n",word);
	 word = strtok(NULL, " \n");
      }

      if(nActiveCh>9){
	cout << endl;
	char dummy;
	nitem=fread(&dummy,1,1,pFILE);
      }

      int whileCounter = 0;
      while(nitem>0){ //event loop
      std::vector<TObject*> eventTrash;
      
      whileCounter++;
	nitem = fread (&EventNumber	,sizeof(int), 1,pFILE);
	nitem = fread (&EpochTime	,sizeof(double)      , 1,pFILE);
	nitem = fread (&Year		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Month		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Day		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Hour		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Minute		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Second		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Millisecond	,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&nCh	,sizeof(unsigned int),1,pFILE); // since V2.8.14 the number of stored channels is written for each event


	if(EventNumber%100==0)printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber,Year,Month,Day,Hour,Minute,Second,Millisecond,nActiveCh,nCh);

	float	MeasuredBaseline[16];
	float	AmplitudeValue[16];
	float	ComputedCharge[16];
	float	RiseTimeInstant[16];
	float	FallTimeInstant[16];
	float	RawTriggerRate[16];
	float floatR=-1;
        for(int i = 0;i<nCh;i++){
	  //printf("i, currentPositionByte %d %ld\n",i,ftell(pFILE));
	  nitem = fread (&ChannelNr[i]	       ,sizeof(int),1,pFILE);
          nitem = fread (&EventIDsamIndex[i]        ,sizeof(int),1,pFILE);
	  nitem = fread (&FirstCellToPlotsamIndex[i],sizeof(int),1,pFILE);
	  nitem = fread (&floatR,1,4,pFILE); MeasuredBaseline[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); AmplitudeValue[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); ComputedCharge[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); RiseTimeInstant[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); FallTimeInstant[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); RawTriggerRate[i] = floatR;
	  ChannelNr[i]=i;

	  TString title("");
	  title.Form("ch %d, ev %d",i,EventNumber);
	  hCh.Reset();
	  hCh.SetTitle(title);
	  
	  if(i<=5||i==15){
	    for(int j = 0;j<1024;j++){
	      nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
	      hCh.SetBinContent(j+1,-(amplValues[i][j]*coef*1000));
	    }//for 1024
	  }
	  else{
	    for(int j = 0;j<1024;j++){
	      nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
	      hCh.SetBinContent(j+1,(amplValues[i][j]*coef*1000));
	    }//for 1024
	  }
	  
	  //for(int t=0;t<nActiveCh-nCh;t++){
	  //  int dummy;
	  //  nitem = fread(&dummy,sizeof(int),1,pFILE);
	  //  printf("trigger channel number: %d\n",dummy);
	  //}

	  max[i]=hCh.GetMaximum();
	  min[i]=hCh.GetMinimum();
	  amp[i]=hCh.GetMaximum();
	  
	  getBL(&hCh, BL_output,0,30);
	  BL[i] = BL_output[0];
	  BL_RMS[i] = BL_output[1];
	  
	  if(EventNumber%wavesPrintRate==0){
	    cWaves.cd(1+4*(i%4)+(i)/4);
	    hCh.DrawCopy();
	  }
	  for(int j=1;j<=hCh.GetXaxis()->GetNbins();j++){
	    hCh.SetBinError(j,BL_RMS[i]);
	  }
	  hChtemp.at(i) = hCh;
	  t[i]=CDF(&hCh,fTrigFit,0.1);

	  if(i<=5||i==15)Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-BL[i]*1024*SP)/pe;//Calculating Integral of histogram from 0 to 300ns; starting from bin 1 (0 is the overflow bin) to bin corresponding to 300ns. Option "width" multiplies by bin-width such that the integral is independant of the binning
	  else Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-BL[i]*1024*SP);
	  
          if(EventNumber%ch0PrintRate==0&&i==0){
	    //cCh0.cd(1);
	    cCh0.cd();
	    hCh.DrawCopy("hist");
	    TLine ln(t[0],-2000,t[0],2000);
	    ln.SetLineColor(2);
	    //hCh.GetYaxis()->SetRangeUser(-10,1200);
	    ln.Draw("same");
	    if(ch0PrintStatus<0){cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf("),"pdf");ch0PrintStatus=0;}
	    else cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf"),"pdf");
	  }

	 if(EventNumber%trigPrintRate==0&&(i<4)){
	    cTrig.cd(i+1);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000);
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    //fTrigFit->DrawCopy("same");
	    eventTrash.push_back(ln);
	  }

	 if(EventNumber%signalPrintRate==0&&(i>=4&&i<=7)){
	    cSignal.cd(i+1-4);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000); 
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    eventTrash.push_back(ln);
	  }

       }//for nCh

      trigT = (t[0]+t[1]+t[2]+t[3])/4;
      trigTp = (t[0]+t[1]-t[2]-t[3])/4;
      t0t1 = (t[0]-t[1]);
      t2t3 = (t[2]-t[3]);
      tPMT1 = t[4]-trigT;
      tPMT2 = t[5]-trigT;
      if(tPMT2<-52){
	t[5]=CDFinvert(&hChtemp.at(5),0.1);
	tPMT2 = t[5]-trigT;
      }
      //tPMT2i = iCFD(&hChtemp.at(5),trigT-55,2,BL[5])-trigT;
      Integral[5] = integral(&hChtemp.at(5),t[5]-5,t[5]+65,BL[5])/pe;
      
      
      hChtemp.at(6).Add(&hChtemp.at(6),&hChtemp.at(7),-1,1);
       if(EventNumber%wavesPrintRate==0){
	    cWaves.cd(1+4*(6%4)+(6)/4);
	    hChtemp.at(6).DrawCopy();
	  }
      
      Integral_0_300[6] = -Integral_0_300[6] + Integral_0_300[7];
      Integral_0_300[7] = Integral_0_300[6]-Integral_0_300[8]-Integral_0_300[9]-Integral_0_300[10]-Integral_0_300[11]-Integral_0_300[12]-Integral_0_300[13]-Integral_0_300[14];
// //       hChtemp.at(7).Add(&hChtemp.at(6),&hChtemp.at(8),1,-1);
// //       hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(9),1,-1);
// //       hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(10),1,-1);
// //       hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(11),1,-1);
// //       hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(12),1,-1);
// //       hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(13),1,-1);
// //       hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(14),1,-1);
// //       if(EventNumber%wavesPrintRate==0){
// // 	    cWaves.cd(1+4*(7%4)+(7)/4);
// // 	    hChtemp.at(7).DrawCopy();
// // 	  }
// //       
      
      tSUMp = t[6] - trigT;
      tSUMm = t[7] - trigT;
      t[6] = CDF(&hChtemp.at(6),fTrigFit,0.1);
      tSiPM = t[6]-trigT;
      if(tSiPM<-56){
	t[6]=CDFinvert(&hChtemp.at(6),0.1);
	tSiPM = t[6]-trigT;
      }
      
      
      
      
      if(max[15]>5||min[15]<-5)isVeto=1;
      else isVeto = 0;
      trigGate = abs(*(std::max_element(t,t+4))-*(std::min_element(t,t+4)));  
      
      if(max[0]<1240&&max[1]<1240&&max[2]<1240&&max[3]<1240&&isVeto==0){
	isTrig=1;
	if(isTrig&&BL[0]<1.1&&BL[1]<1.1&&BL[2]<1.1&&BL[3]<1.1){
	  isTrig=1;
	  if(trigT<140&&trigT>90&&trigGate<10){
	    isTrig=1;
	  }
	  else isTrig=0;
	}
	else isTrig=0;
      }
      else isTrig=0;
      if(isTrig==1){
	int shift = (int)((140-trigT)/SP);
	for(int j=0;j<(int)hChtemp.size();j++){
	  hChSum.at(j)->Add(&hChtemp.at(j),1);
	  hChShift_temp.at(j).Reset();
	  for(int bin=1;bin<=hCh.GetXaxis()->GetNbins()-shift;bin++){
	   hChShift_temp.at(j).SetBinContent(shift+bin,hChtemp.at(j).GetBinContent(bin));
	  }
	  hChShift.at(j)->Add(&hChShift_temp.at(j),1);
	}
      }
      
      
      if(isTrig==1&&max[5]<1240)isGoodSignal_5=1;
      else isGoodSignal_5=0;
      
      if(EventNumber%wavesPrintRate==0){
	    //TString plotSaveName("");
	    //plotSaveName.Form("%s/wave-%d.png",plotSaveFolder.Data(),EventNumber);
	    if(wavePrintStatus<0){cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");wavePrintStatus=0;}
	    else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
      }
      if(EventNumber%trigPrintRate==0){
	    if(trigPrintStatus<0){cTrig.Print((TString)(plotSaveFolder+"/trig.pdf("),"pdf");trigPrintStatus=0;}
	    else cTrig.Print((TString)(plotSaveFolder+"/trig.pdf"),"pdf");
      }
      if(EventNumber%signalPrintRate==0){
	    if(signalPrintStatus<0){cSignal.Print((TString)(plotSaveFolder+"/signal.pdf("),"pdf");signalPrintStatus=0;}
	    else cSignal.Print((TString)(plotSaveFolder+"/signal.pdf"),"pdf");
      }

      tree->Fill();
      //cleanEventMemory(eventTrash);
      }//while events

    }//while
    inList.close();
    cWaves.Clear();
    cWaves.Print((TString)(plotSaveFolder+"/waves.pdf)"),"pdf");
    cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf)"),"pdf");
    cTrig.Print((TString)(plotSaveFolder+"/trig.pdf)"),"pdf");
    cSignal.Print((TString)(plotSaveFolder+"/signal.pdf)"),"pdf");

  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}

float CDF(TH1F* hWave, TF1* fTrigFit,float thr){
  float peak=hWave->GetMaximum();
  int timePos=1;
  float val = 0;
  while(abs(val)<thr*peak){
    timePos+=1;
    val = hWave->GetBinContent(timePos);
  }
  
  double x1 = SP*(timePos-1);
  double x2 = SP*(timePos);
  double y1 = hWave->GetBinContent(timePos-1);
  double y2 = hWave->GetBinContent(timePos);
  double k = (x2-x1)/(y2-y1);
  return  x1+k*(thr*peak-y1);
  
  //fit procedure
  //fTrigFit->SetParameter(1,SP*(timePos)+1);
  //fTrigFit->SetRange(SP*(timePos-6),SP*(timePos+2));
  //hWave->Fit(fTrigFit,"RNQ");
  //double p0=fTrigFit->GetParameter(0);
  //double p1=fTrigFit->GetParameter(1);
  //double p2=fTrigFit->GetParameter(2);
  //return p1-sqrt(2*p2*p2*log(p0/(0.1*abs(peak))));
}

float CDFinvert(TH1F* hWave,float thr){
  float peak=hWave->GetMaximum();
  int timePos=hWave->GetMaximumBin();
  float val = peak;
  while(val>thr*peak){
    val = hWave->GetBinContent(timePos);
    timePos-=1;
  }
  
  double x1 = SP*(timePos);
  double x2 = SP*(timePos+1);
  double y1 = hWave->GetBinContent(timePos);
  double y2 = hWave->GetBinContent(timePos+1);
  double k = (x2-x1)/(y2-y1);
  return  x1+k*(thr*peak-y1);
  
}


float iCFD(TH1F* hWave,float t,float thrNpe,float BL){
  int bin1 = hWave->GetXaxis()->FindBin(t);
  float bin1_UpEdge = hWave->GetXaxis()->GetBinUpEdge(bin1);
  float sum = (bin1_UpEdge-t)*(hWave->GetBinContent(bin1)-BL);
  int bin=bin1+1;
  while(sum<thrNpe){
      sum+=(hWave->GetBinContent(bin)-BL)*SP;
      bin+=1;
      if(bin>1024)return -999;
  }
  float bin_LowEdge = hWave->GetXaxis()->GetBinUpEdge(bin);
  return bin_LowEdge + (sum-thrNpe)/(hWave->GetBinContent(bin)-BL); 
}


float integral(TH1F* hWave,float t1,float t2,float BL){
  float BW = hWave->GetXaxis()->GetBinWidth(1);
  int bin1 = hWave->FindBin(t1);
  int bin2 = hWave->FindBin(t2);
  float c1 = hWave->GetBinContent(bin1);
  float c2 = hWave->GetBinContent(bin2);
  return hWave->Integral(bin1,bin2,"width")-BL*(t2-t1)-c1*(t1-hWave->GetXaxis()->GetBinLowEdge(bin1))-c2*(hWave->GetXaxis()->GetBinUpEdge(bin2)-t2);
}



float* getBL(TH1F* hWave, float* BL, float t1, float t2){
  /*
  Function to calculate the baseline of the given TH1F-Object.
  Input: TH1F-Object to calculate baseline for; bool isNegative; float-array BL for the output
  Output: baseline and rms of baseline written to 1st and 2nd component of BL-array
  The baseline is calculated as the mean of the values in range of (t1,t2) of the TH1F-Object. The rms
  value is also calculated from the same values.
  
  The float-array BL that is used for the output must be declared before the function call using 'float BL[2];'.
  This is to insure that the output is stored on the heap and not deleted when the memory on the stack is freed up.
   
  Dependencies: function uses C++ vector-class and needs the TMath-header
  */
  
  vector<float> amp;
  for (int i = int(t1/SP); i < int(t2/SP); i++){
    
    amp.push_back(hWave->GetBinContent(i+1));
  }
  BL[0] = TMath::Mean(amp.begin(), amp.end());
  BL[1] = TMath::RMS(amp.begin(), amp.end());
  return BL;
}


TString getRunName(TString inDataFolder){
      char* word;
      char* lastWord;
      word = strtok((char*)inDataFolder.Data(),"/");
      while(word != NULL){
	 //printf ("%s\n",word);
	 lastWord = word;
	 word = strtok(NULL, "/");
      }

      return (TString)lastWord;
}

void cleanEventMemory(std::vector<TObject*>& trash){
  for(int i=0;i<(int)trash.size();i++){
    trash.at(i)->Delete();
  }
  trash.clear();
}
 

float solidAngleABH(float A,float B, float H){//return solid angle of pyramid with rectangular base (a,b) and hight of H
  return 4*TMath::ASin((A*B)/sqrt((A*A+H*H)*(B*B+H*H)));
}

std::vector<float> getStartPos(float horisontal, float vertical, float angle){//return start position of track in the box (x,y,z,angle) [mm,mm,mm,rad]
    std::vector<float> pos(4,0);
    pos[3]=angle*TMath::Pi()/180;
    pos[1]=(boxL/2-vertical*10);
    pos[0]=boxL/2+horisontal*10*TMath::Cos(pos[3])-(boxH/2-horisontal*10*TMath::Sin(pos[3]))*TMath::Tan(pos[3]);
    if(pos[0]<0){
      pos[2]=fabs(pos[0])/TMath::Tan(fabs(pos[3]));
      pos[0]=0;
    }
    else if(pos[0]>boxL){
      pos[2]=(pos[0]-boxL)/TMath::Tan(fabs(pos[3]));
      pos[0]=boxL;
    }
    else {
      pos[2]=0;
    }
    return pos;
}

int getZone(float x,float y,float z,float pmtX, float pmtY){
  float r2 = (x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY);
  if(z>=(boxH-womL)){
    if(r2>(womDout*womDout/4))return 1;
    else if(r2>(womDout*womDout/4))return 2;
    else return 3;
  }
  else{
    if(r2>=(womR*womR))return 4;
    else return 5;
  }
}

std::vector<float> solidAngleFactor(const std::vector<float> &startPos,const std::vector<float> &pmtPos){
  std::vector<float> result(2,0);
  float x0 = startPos[0];
  float y0 = startPos[1];
  float z0 = startPos[2];
  float angle = startPos[3];
  
  
  float length = 0;
  float dLength = 1; //mm
  float integratedSolidAngleFactor = 0;
  float x = x0;
  float y = y0;
  float z = z0;  
  while(z<=boxH && x<=boxL && z>=0 && x>=0){
    x=x+TMath::Sin(angle);
    z=z+TMath::Cos(angle);
    length=length+dLength;
    integratedSolidAngleFactor += dLength*getSolidAngle(x,y,z,pmtPos[0],pmtPos[1])/(4*TMath::Pi());
    //printf("solidAngleFactor: %4.2f %4.2f %4.2f %4.2f %d\n",x,y,z,length,getZone(x,y,z,pmt2Pos[0],pmt2Pos[1]));
    //printf("%d ",getZone(x,y,z,pmtPos[0],pmtPos[1]));
  } 
  //printf("integratedSolidAngleFactor: %4.2f\n",integratedSolidAngleFactor);
  result[0]=integratedSolidAngleFactor;
  result[1]=length;
  return result;
}

float getSolidAngle(float x,float y,float z,float pmtX, float pmtY){//return solid angle per one step
  float solidAngle = 0;
  float x1=-999,y1=-999,z1=-999;
  float x2=-999,y2=-999,z2=-999;
  float dr = (4-TMath::Pi())*womR/4;
  float r = TMath::Pi()*womR/4;
  int zone = getZone(x,y,z,pmtX,pmtY);
  if(zone==1){
    float h = sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY)) - r;
    float l1 = z+womL-boxH;
    float l2 = boxH-z;
    float a1 = sqrt(h*h+l1*l1);
    float a2 = sqrt(h*h+l2*l2);
    float A = a1*sqrt(2*(1-(h*h-l1*l2)/(a1*a2)));
    float B = 2*womR;
    float H = sqrt(a1*a1-A*A/4);
    return solidAngleABH(A,B,H);
  }
  else if(zone==4){
    float h1 = sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY)) + r;
    float h2 = sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY)) - r;
    float l1 = boxH-womL-z;
    float l2 = boxH-z;
    float a1 = sqrt(h1*h1+l1*l1);
    float a2 = sqrt(h2*h2+l2*l2);
    float cosA = (a1*a1+a2*a2-4*r*r-womL*womL)/(2*a1*a2);
    float A = a1*sqrt(2*(1-cosA));
    float B = 2*womR;
    float H = sqrt(a1*a1-A*A/4);
    return solidAngleABH(A,B,H);
  }
  else if(zone==2){
    float l1 =womRin - sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float l2 =womRin + sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float h = boxH - z;
    float H = womL - h;
    float a1 = sqrt(l1*l1+h*h);
    float a2 = sqrt(l2*l2+h*h);
    float A1 = sqrt(l1*l1+H*H);
    float A2 = sqrt(l2*l2+H*H);
    return 2*TMath::Pi()*(sqrt((1+(h*h-l1*l2)/(a1*a2))/2)+sqrt((1+(H*H-l1*l2)/(A1*A2))/2));
  }
  else if(zone==5){
    float l1 =womR - sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float l2 =womR + sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float h = boxH - z -womL;
    float H = boxH - z;
    float a1 = sqrt(l1*l1+h*h);
    float a2 = sqrt(l2*l2+h*h);
    float A1 = sqrt(l1*l1+H*H);
    float A2 = sqrt(l2*l2+H*H);
    return 2*TMath::Pi()*(sqrt((1+(H*H-l1*l2)/(A1*A2))/2)-sqrt((1+(h*h-l1*l2)/(a1*a2))/2));
  }
  else return 0;
}

