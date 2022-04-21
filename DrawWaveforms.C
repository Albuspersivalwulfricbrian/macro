#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TString.h>
#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "ChannelEntry.h"

class MyMainFrame : public TGMainFrame {

private:
   static const Int_t total_channels = 35;
   TRootEmbeddedCanvas  *fEcan;
   TGStatusBar          *fStatusBar;
   TGNumberEntry       *fNumber;
   Int_t entrynum;

   TTree *t1;
   ChannelEntry channel_info[total_channels];
   Int_t n;
   TGNumberEntry       *fChannel;
   Int_t ChNum;
public:
   MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString s);
   virtual ~MyMainFrame();
   void DoExit();
   void DoDraw();
   void SetStatusText(const char *txt, Int_t pi);
   void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);
   void SetEntry();
   void SetChannel();
   void ReadFile(TString s);

   ClassDef(MyMainFrame, 0)
};

void MyMainFrame::DoDraw()
{
   n = channel_info[ChNum].wf_size;
   Printf("Slot DoDraw()");

   TCanvas *c1 = fEcan->GetCanvas();
   c1->SetFillColor(42);
   c1->SetGrid();
    cout << n << endl;
   Double_t x[2048] = {0.}; Double_t y[2048] = {0.};
   
   for (Int_t i=0;i<n;i++) {
     x[i] = i*16;
     y[i] = channel_info[ChNum].wf[i];
   }

   ///////////
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(1);
   gr->SetMarkerStyle(21);
   gr->SetMarkerSize(1);

   gr->SetTitle("ADC Waveform");
   gr->GetXaxis()->SetTitle("Pulse time [ns]");
   gr->GetYaxis()->SetTitle("ADC channels");
   gr->Draw("APL");

   // TCanvas::Update() draws the frame, after which it can be changed
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
   c1->Update();
}

void MyMainFrame::DoExit()
{
   printf("Exit application...");
   gApplication->Terminate(0);
}

void MyMainFrame::ReadFile(TString s)
{
   TFile *f = new TFile(s);
   t1 = (TTree*)f->Get("adc64_data");
   for (Int_t channel = 0; channel < total_channels; channel++) 
   {
      t1->SetBranchAddress(Form("channel_%i",channel), &channel_info[channel]);

   }

}


void MyMainFrame::SetStatusText(const char *txt, Int_t pi)
{
   fStatusBar->SetText(txt,pi);
}

void MyMainFrame::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString s) :
   TGMainFrame(p, w, h)
{
   
    MyMainFrame::ReadFile(s);
   ChNum = 3500;
   // Create the embedded canvas
   fEcan = new TRootEmbeddedCanvas(0,this,1920,1080);
   Int_t wid = fEcan->GetCanvasWindowId();
   TCanvas *myc = new TCanvas("MyCanvas", 10,10,wid);
   fEcan->AdoptCanvas(myc);
   myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,
               "EventInfo(Int_t,Int_t,Int_t,TObject*)");

   AddFrame(fEcan, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                                     kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
   TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 200, 40);

////////chnum

   TGLabel* fLchannel = new TGLabel(this, "Channel Number");
   AddFrame(fLchannel,  new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   fChannel = new TGNumberEntry(this, 0, 9,999, TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
                                               0,34);
   fChannel->Connect("ValueSet(Long_t)", "MyMainFrame", this, "SetChannel()");
   //(fChannel->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "SetChannel()");
   AddFrame(fChannel, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   //AddFrame(fChannel, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
////////////




   TGLabel* fLevent = new TGLabel(this, "Event Number");
   AddFrame(fLevent,  new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   fNumber = new TGNumberEntry(this, 0, 9,999, TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
                                               0,1000000000000);
   fNumber->Connect("ValueSet(Long_t)", "MyMainFrame", this, "SetEntry()");
   (fNumber->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this,"SetEntry()");

   AddFrame(fNumber, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   //AddFrame(fNumber, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));


   TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
   exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
   hframe->AddFrame(exit, new TGLayoutHints(kLHintsRight, 5, 5, 3, 4));
   AddFrame(hframe, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   
   // Set a name to the main frame
   SetWindowName("Embedded Canvas Status Info");
   MapSubwindows();

   // Initialize the layout algorithm via Resize()
   Resize(GetDefaultSize());

   // Map main frame
   MapWindow();
}

MyMainFrame::~MyMainFrame()
{
   // Clean up main frame...
   Cleanup();
   delete fEcan;
}

void MyMainFrame::SetEntry()
{
    entrynum = fNumber->GetNumberEntry()->GetIntNumber();
   t1->GetEntry(entrynum);

    //cout << ChNum << " " << n << endl;
    DoDraw();
}
void MyMainFrame::SetChannel()
{
   //channel_info.Initialize();
   //if (CHNum == 32 || CHNum == 33 || CHNum == 34);
   ChNum = fChannel->GetNumberEntry()->GetIntNumber();
   cout << ChNum << endl;
    DoDraw();

}






























void DrawWaveforms(    
        TString file_path = "/home/doc/Downloads",
        TString file_name = "07a8de9a_20220418_100507.root"
)

{
   // Popup the GUI...
   new MyMainFrame(gClient->GetRoot(), 200, 200, file_path+"/"+file_name);
}
