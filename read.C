//root
#include <TLine.h>
//#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

//C, C++
#include <stdio.h>
#include <assert.h>
//#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <string>
//#include <iomanip>

using namespace std;

float SP = 0.3125;
int wavesPrintRate = 1000000;
int ch0PrintRate = 1000000;
int trigPrintRate = 100;
int signalPrintRate = 100;
double coef = 2.5 / (4096 * 10);

  int runNr = -999;
  float horizontal=-999;
  float vertical=-999;
  float angle=-999;
  int Ptype=-999; // 1=hadrons, 2=electron, 3=muons

float CFD(TH1F* hWave,bool isNegative = 1);
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
  else if(argc == 12){
    inFileList = argv[1];
    inDataFolder = argv[2];
    outFile = argv[3];
    runNr=atoi(argv[4]);
    horizontal = atof(argv[5]);
    vertical = atof(argv[6]);
    angle = atof(argv[7]);
    //argv[8] - nevents


    cout<<"In data file list : "<<inFileList<<endl
        <<"In data path      : "<<inDataFolder<<endl
        <<"Out root file     : "<<outFile<<endl
	<<"Run nuber         : "<<runNr<<endl;
	printf("hor,ver,ang: %4.2f %4.2f %4.2f\n",horizontal,vertical,angle);
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
  Float_t tSUMp = -999;
  Float_t tSUMm = -999;
  Float_t trigTp = -999;//t_trig' = [(t0+t1)-(t2+t3)]/4
  Float_t t0t1 = -999;//t0t1 = [(t0-t1)]
  Float_t t2t3 = -999;//t2t3 = [(t2-t3)]
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  Float_t amp[16];
  Float_t t[16];
  Int_t EventIDsamIndex[16];
  Int_t FirstCellToPlotsamIndex[16];

  Short_t amplValues[16][1024];
  TH1F hCh("hCh","dummy;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
  TString plotSaveFolder  = getRunName(_inFileList);
  plotSaveFolder.ReplaceAll(".list","");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(4,4);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);

   //Event USBWC
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
  tree->Branch("tSUMp",&tSUMp, "tSUMp/F");
  tree->Branch("tSUMm",&tSUMm, "tSUMm/F");

  tree->Branch("trigTp",&trigTp, "trigTp/F");
  tree->Branch("t0t1",&t0t1, "t0t1/F");
  tree->Branch("t2t3",&t2t3, "t2t3/F");
  tree->Branch("runNr",&runNr, "runNr/I");
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  tree->Branch("amp",amp, "amp[nCh]/F");
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

  ///tree->Branch();
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
      while(nitem>0){
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


	printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber,Year,Month,Day,Hour,Minute,Second,Millisecond,nActiveCh,nCh);

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
	  for(int j = 0;j<1024;j++){
	    nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
	    hCh.SetBinContent(j+1,fabs(amplValues[i][j]*coef*1000));
	  }//for 1024
//        for(int t=0;t<nActiveCh-nCh;t++){
// 	    int dummy;
// 	    nitem = fread(&dummy,sizeof(int),1,pFILE);
// 	    printf("trigger channel number: %d\n",dummy);
// 	}


	  if(EventNumber%wavesPrintRate==0){
	    cWaves.cd(1+4*(i%4)+(i)/4);
	    hCh.DrawCopy();
	  }
	  amp[i]=hCh.GetMaximum();
	  t[i] = CFD(&hCh);

          if(EventNumber%ch0PrintRate==0&&i==0){
	    cCh0.cd(1);
	    hCh.DrawCopy();
	    TLine ln(t[0],-2000,t[0],2000);
	    ln.SetLineColor(2);
	    //hCh.GetYaxis()->SetRangeUser(-10,1200);
	    ln.Draw("same");
	    if(ch0PrintStatus<0){cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf("),"pdf");ch0PrintStatus=0;}
	    else cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf"),"pdf");
	  }

	 if(EventNumber%trigPrintRate==0&&(i>=0&&i<=3)){
	    cTrig.cd(i+1);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000); //stop point
	    ln->SetLineColor(2);
	    ln->Draw("same");
	  }

	 if(EventNumber%signalPrintRate==0&&(i>=4&&i<=7)){
	    cSignal.cd(i+1-4);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000); //stop point
	    ln->SetLineColor(2);
	    ln->Draw("same");
	  }



       }//for nCh


       if(EventNumber%wavesPrintRate==0){
	    //TString plotSaveName("");
	    //plotSaveName.Form("%s/wave-%d.png",plotSaveFolder.Data(),EventNumber);
	    if(wavePrintStatus<0){cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");wavePrintStatus=0;}
	    else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
       }

      trigT = (t[0]+t[1]+t[2]+t[3])/4;
      tPMT1 = t[4] - trigT;
      tPMT2 = t[5] - trigT;
      tSUMp = t[6] - trigT;
      tSUMm = t[7] - trigT;
      trigTp = (t[0]+t[1]-t[2]-t[3])/4;
      t0t1 = (t[0]-t[1]);
      t2t3 = (t[2]-t[3]);

       if(EventNumber%trigPrintRate==0){
	    if(trigPrintStatus<0){cTrig.Print((TString)(plotSaveFolder+"/trig.pdf("),"pdf");trigPrintStatus=0;}
	    else cTrig.Print((TString)(plotSaveFolder+"/trig.pdf"),"pdf");
	  }
      if(EventNumber%signalPrintRate==0){
	    if(signalPrintStatus<0){cSignal.Print((TString)(plotSaveFolder+"/signal.pdf("),"pdf");signalPrintStatus=0;}
	    else cSignal.Print((TString)(plotSaveFolder+"/signal.pdf"),"pdf");
	  }



      tree->Fill();
      }//while ivents

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

float CFD(TH1F* hWave,bool isNegative){
  float time = -999;
  int timePos=-999;
  float peak=-999;
  int peakPos = -999;
  peak = hWave->GetMaximum();
  peakPos = hWave->GetMaximumBin();
  float val = peak;
  timePos=peakPos;
  while(abs(val)>0.1*abs(peak)){
    val = hWave->GetBinContent(timePos);
    timePos-=1;
  }
  time = SP*(timePos);
  //time = SP*(peakPos);
  return time;
}

TString getRunName(TString inDataFolder){
      char* word;
      char* lastWord;
      word = strtok((char*)inDataFolder.Data(),"/");
      while(word != NULL){
	 printf ("%s\n",word);
	 lastWord = word;
	 word = strtok(NULL, "/");
      }

      return (TString)lastWord;
}
