
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"

// Draw a Bidimensional Histogram in many ways
// together with its profiles and projections


void macro1(){
   //make array of pallette 

   Int_t MyPalette[100];
   Double_t Red[]    = {0., 0.0, 1.0, 1.0, 1.0};
   Double_t Green[]  = {0., 0.0, 0.0, 1.0, 1.0};
   Double_t Blue[]   = {0., 1.0, 0.0, 0.0, 1.0};
   Double_t Length[] = {0., .25, .50, .75, 1.0};
   Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 100);
   for (int i=0;i<100;i++) MyPalette[i] = FI+i;
   
   
   gStyle->SetPalette(100, MyPalette);
   gStyle->SetOptStat(0); //gets rid of legend 


   //change title and barriers etc 
   auto bidi_h = new TH2F("bidi_h","vol_layer;z;r",
                          10000,-3000,3000,  // X axis
                          10000,-3000,3000); // Y axis


   //change to filling with my data 
   TRandom3 rgen;
   for (int i=0;i<500000;i++) {
      bidi_h->Fill(rgen.Gaus(0,2),10-rgen.Exp(4),.1);
   }

//this all works 
   auto c=new TCanvas("Canvas","Canvas",800,800);
   c->cd(1); bidi_h->Draw("Colz");

   c->Print("graph_with_law.pdf");
}
int main(){
    macro1();
    }