// #include <iostream>
// #include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TH2F.h>
// #include <TObjArray.h>
#include <TString.h>
#include <TGraphErrors.h>
#include "ChannelEntry.h"
#include "CHSH_calculator.h"
#include "CHSH_calculator_class.h"
#include <string.h>
#include "Check_entanglement_cuts.h"
#include "like_ivashkin_wants_it.h"
using namespace std;
using namespace CHSH;
using namespace CUTS;

#define UseDecoherentPhotons 1
#define Analyze_Monte_Carlo 0
#define UseNoCutG 1
#define UseTimeCut 1


void MiniDST_analysis_1(TString source_path = "/home/doc/entanglement/with_spline/entangled/")
{
    const Int_t module_angles = 8;
    const Int_t all_angles = module_angles*2;

#if UseDecoherentPhotons
    TString middle_path  = source_path + "decoh_mini_tree";
#else
    TString middle_path = source_path + "entangled_mini_tree";
#endif
    MiniTree *short_tree = new MiniTree;
    
    mini_tree_time *time_tree = new mini_tree_time;
    TFile *f = TFile::Open(middle_path+".root");
    TTree *MiniDST_tree = (TTree*)f->Get("Signals");
    
    TString result_path = middle_path + "_analyzed";
    TFile *result_root = new TFile(result_path+".root","RECREATE");
    //short_tree->SetBranch(MiniDST_tree);
    MiniDST_tree->SetBranchAddress("MiniTree",short_tree);
    time_tree->SetBranch(MiniDST_tree);

    Float_t total_events_for_angle_diff[16] = {0};
    Float_t total_events_for_angle_diff_err[16] = {0};    
    Int_t NumEvents[16][16] = {0};     
    Float_t angle_arr[16]  = {0.};
    Float_t angle_arr_err[16] = {0.};    
    Float_t high_time_cut[33] = {0};
    Float_t low_time_cut[33] = {0};

    for (Int_t i = 0; i <16; i++) angle_arr[i] = (float)i*22.5;
    /////////////////////////////
    //////////////////////////////

        Float_t low_det0_cut = 180;
        Float_t high_det0_cut = 300;
        Float_t low_det1_cut = 180;
        Float_t high_det1_cut = 300;

        Float_t low_scat0_cut = 150;
        Float_t high_scat0_cut = 300; 
        Float_t low_scat1_cut = 180;
        Float_t high_scat1_cut  = 300;

        Float_t high_Intermediate_cut = 60;
        Float_t low_Intermediate_cut = 2;
/////////////////////////////////
////////////////////////////////

//////////////Setting TCUTG boarders
    TCanvas *canv_0 = new TCanvas("canv_0,canv_0");
    canv_0->cd();

    #if UseDecoherentPhotons
    TH2F *h2 = new TH2F("h2","h2", 75,0,400,75,0,150); h2->GetZaxis()->SetRangeUser(5,70);
    MiniDST_tree->Draw("MiniTree.EdepIntermediate:MiniTree.EdepDet0 >> h2","TimeTree.TimeIntermediate-TimeScat0 < 75 && TimeTree.TimeIntermediate-TimeScat0 > 50", "colz");
        //MiniDST_tree->Draw("MiniTree.EdepDet0:MiniTree.EdepScat0 >> (75,0,500,75,0,500)","", "colz");

    #else
    MiniDST_tree->Draw("MiniTree.EdepDet0:MiniTree.EdepScat0 >> (300,0,500,300,0,500)","", "colz");
    #endif
    
   TCutG *cutg = new TCutG("CUTG",16);
   cutg->SetVarX("EdepDet0 ");
   cutg->SetVarY("MiniTree.EdepIntermediate");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,181.822,40.0);
   cutg->SetPoint(1,185.536,30.935);
   cutg->SetPoint(2,185.536,30.935);
   cutg->SetPoint(3,198.045,16.5658);
   cutg->SetPoint(4,211.532,10.1295);
   cutg->SetPoint(5,211.532,10.1295);
   cutg->SetPoint(6,235.379,5.82626);
   cutg->SetPoint(7,259.812,4.59141);
   cutg->SetPoint(8,259.812,4.59141);
   cutg->SetPoint(9,288.937,5.67659);
   cutg->SetPoint(10,312.002,11.2521);
   cutg->SetPoint(11,321.579,19.8587);
   cutg->SetPoint(12,325.684,28.2033);
   cutg->SetPoint(13,330.571,36.5854);
   cutg->SetPoint(14,331.353,40.0);
   cutg->SetPoint(15,181.822,40.0);
   cutg->Draw("same");
    canv_0->SaveAs(result_path+".pdf(",".pdf");
    canv_0->Write();

//////////////////////////////////////////
        Int_t nbins = 125;
        Int_t left_bin = 0;
        Int_t right_bin = 450;
        TH1F *det0 = new TH1F("det0","det0", nbins, left_bin, right_bin);
        TH1F *det1 = new TH1F("det1","det1", nbins, left_bin, right_bin);
        TH1F *scat0 = new TH1F("scat0","scat0", nbins, left_bin, right_bin);
        TH1F *scat1 = new TH1F("scat1","scat1", nbins, left_bin, right_bin);
        TCanvas *new_canv = new TCanvas("new", "new");
        MiniDST_tree->Draw("MiniTree.EdepDet0 >> det0"

#if UseDecoherentPhotons

        ,Form("MiniTree.EdepIntermediate < %f && MiniTree.EdepIntermediate > %f", high_Intermediate_cut, low_Intermediate_cut)
#endif        
        ); 
        double_gauss_fit(det0, low_det0_cut, high_det0_cut);
        det0->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();

        MiniDST_tree->Draw("MiniTree.EdepDet1 >> det1"
#if UseDecoherentPhotons

        ,Form("MiniTree.EdepIntermediate < %f && MiniTree.EdepIntermediate > %f", high_Intermediate_cut, low_Intermediate_cut)
#endif        
        );
        double_gauss_fit(det1, low_det1_cut, high_det1_cut);
        det1->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();

        MiniDST_tree->Draw("MiniTree.EdepScat1 >> scat0", 
        Form("MiniTree.EdepDet1 < %f && MiniTree.EdepDet1 > %f", high_det1_cut, low_det1_cut));
        double_gauss_fit(scat0, low_scat1_cut, high_scat1_cut);
        scat0->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();


        MiniDST_tree->Draw("MiniTree.EdepScat0 >> scat1", 
        Form("MiniTree.EdepDet0 < %f && MiniTree.EdepDet0 > %f" 
#if UseDecoherentPhotons
        "&& MiniTree.EdepIntermediate < %f && MiniTree.EdepIntermediate > %f"
#endif
        , high_det0_cut, low_det0_cut 
#if UseDecoherentPhotons
        ,high_Intermediate_cut, low_Intermediate_cut
#endif
        ));
        double_gauss_fit(scat1, low_scat0_cut, high_scat0_cut);
        scat1->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();

    delete det0;
    delete det1;
    delete scat0;
    delete scat1;
#if UseTimeCut
   ////Drawing_time_spectra
        TH1F *hist_time = new TH1F("hist_time","hist_time", 1201,-600,600);
        //TCanvas *time = new TCanvas("time", "time");
    for (Int_t chnum = 0; chnum < 16; chnum++)
    {
        TH1F *hist_time = new TH1F("hist_time","hist_time", 1201,-600,600);
        MiniDST_tree->Draw("TimeTree.TimeDet0-TimeScat0 >> hist_time", 
        Form("MiniTree.DetNum0 == %i",chnum));
        double_gauss_fit(hist_time, low_time_cut[chnum], high_time_cut[chnum],1.5,1.5);
    }    
    for (Int_t chnum = 16; chnum < 32; chnum++)
    {
        TH1F *hist_time = new TH1F("hist_time","hist_time", 1201,-600,600);
        MiniDST_tree->Draw("TimeTree.TimeDet1-TimeScat1 >> hist_time", 
        Form("MiniTree.DetNum1 == %i",chnum));
        double_gauss_fit(hist_time, low_time_cut[chnum], high_time_cut[chnum],1.5,1.5);
    }  
        hist_time = new TH1F("hist_time","hist_time", 1201,-600,600);
        MiniDST_tree->Draw("TimeTree.TimeIntermediate-TimeScat0 >> hist_time","MiniTree.EdepIntermediate > 2 && EdepIntermediate < 40");

        double_gauss_fit(hist_time, low_time_cut[32], high_time_cut[32],1.5,1.5);
#endif
/////////////////////Selecting coicidences for all counters
    for (Int_t NumEvent = 0; NumEvent < MiniDST_tree->GetEntries(); NumEvent++)
    {
        MiniDST_tree->GetEntry(NumEvent);     
        Short_t num0 = short_tree->DetNum0;
        Short_t num1 = short_tree->DetNum1;
        //cout << num0 << " "<<num1<<endl;

        #if Analyze_Monte_Carlo
        num0-- ; num1--;
        //&& short_tree->EdepIntermediate == 0

        #endif
        if(short_tree->EdepDet1 > low_det1_cut
        && short_tree->EdepDet1 < high_det1_cut
        && short_tree->EdepScat1 > low_scat1_cut
        && short_tree->EdepScat1 < high_scat1_cut 
#if UseNoCutG
        && short_tree->EdepDet0 > low_det0_cut
        && short_tree->EdepDet0 < high_det0_cut        
        && short_tree->EdepScat0 > low_scat0_cut
        && short_tree->EdepScat0 < high_scat0_cut
#else
        // && short_tree->EdepDet0 < 340
        // && short_tree->EdepDet0 > 200 
        // && short_tree->EdepIntermediate > 0
        // && short_tree->EdepIntermediate < 40

        && cutg->IsInside(short_tree->EdepDet0,short_tree->EdepIntermediate)
        
#endif 

        #if UseTimeCut
        && time_tree->TimeScat1-time_tree->TimeScat0 > -3
        && time_tree->TimeScat1-time_tree->TimeScat0 < 6        
        && time_tree->TimeDet0-time_tree->TimeScat0 > low_time_cut[num0]
        && time_tree->TimeDet1-time_tree->TimeScat1 < high_time_cut[num1]
        #endif

#if UseDecoherentPhotons
        && short_tree->EdepIntermediate < high_Intermediate_cut
        && short_tree->EdepIntermediate > low_Intermediate_cut  
        #if UseTimeCut         
        && time_tree->TimeIntermediate - time_tree->TimeScat0 < high_time_cut[32]
        && time_tree->TimeIntermediate - time_tree->TimeScat0 > low_time_cut[32]
        #endif
#endif
        && num0 > -1 && num1 > -1)
            NumEvents[num0][num1-16] +=1;
    }
    /////////////////////Calculate coincidences for ratio

        for (Int_t channel_number = 0; channel_number < 16; channel_number++)
        {
            for(Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)
            {
                Int_t angle = channel_number_2-channel_number-16;
                if (angle < 0) angle += 16;    
                total_events_for_angle_diff[angle]+=NumEvents[channel_number][channel_number_2-16];
                // Int_t i = true_number(channel_number);
                // Int_t j = true_number(channel_number_2);
                // Int_t i1 = true_number(channel_number-2);
                // Int_t j1 = true_number(channel_number_2-2);
                // total_events_for_angle_diff[angle]+=NumEvents[i][j]+NumEvents[j][i]+NumEvents[i][j1]+NumEvents[j][i1];
            
            }
        }

/////////////////////
/////////////////////CHSH inequality and Correlation coefficients
        TF1 *my_fit = new TF1 ("my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",-190,190);
        TF1 *e_fit = new TF1("e_fit","[0]*(cos(2*3.141593/180*x))",-190,190);


            Float_t phi_angle[all_angles] = {0.};
            for (Int_t i = 0; i < module_angles; i++) phi_angle[i] = (float)(i+1)*22.5;
            for (Int_t i = 0; i < module_angles; i++) phi_angle[module_angles+i] = -(float)(i+1)*22.5;
            CHSH_class CHSH_cl;
            //CHSH_cl.Initialize();      
            CHSH_cl.SetNumEvents(NumEvents);
            CHSH_cl.Create_average_CHSH_arrays();
            CHSH_cl.Create_global_CHSH_arrays();
///////////////////////////////////////
///////////////////////////////////////
    Float_t average_E = 0; Float_t sigma_E_average = 0;

    for (Int_t a = 0; a < 16; a++)
    {
        Int_t b = true_number(a+4); Int_t b_1 = true_number(a-4);
        average_E += (E_coeff(NumEvents,a,b) + E_coeff(NumEvents,a,b_1))/32.;
        sigma_E_average += (sqr_E_error(NumEvents,a,b)+sqr_E_error(NumEvents,a,b_1));
    }
    sigma_E_average = sqrt(sigma_E_average)/32;
    TCanvas *canvas_CHSH = new TCanvas ( "canvas_CHSH", "canvas_CHSH");
    canvas_CHSH->Divide(2);

    TGraphErrors * gr_CHSH = new TGraphErrors(module_angles,phi_angle,CHSH_cl.Get_CHSH(), angle_arr_err, CHSH_cl.Get_CHSH_error());
	gr_CHSH->Fit("my_fit","R");
    Float_t a_par_CHSH = abs(my_fit->GetParameter(0));
    Float_t a_error_CHSH = abs(my_fit->GetParError(0));
    canvas_CHSH->cd(1);
    graph_like_ivashkin_wants_it(gr_CHSH,"angle [degrees]","S", Form("Average E <> E(90) = %4.3f+-%4.3f",average_E,sigma_E_average),1);
    gr_CHSH->Draw("AP");

    TGraphErrors * gr_CHSH_all = new TGraphErrors(all_angles-1,phi_angle,CHSH_cl.Get_CHSH_all(), angle_arr_err, CHSH_cl.Get_CHSH_error_all());
	gr_CHSH_all->Fit("my_fit","R");
    a_par_CHSH = abs(my_fit->GetParameter(0));
    a_error_CHSH = abs(my_fit->GetParError(0));
    canvas_CHSH->cd(2);
    graph_like_ivashkin_wants_it(gr_CHSH_all,"angle [degrees]","S", Form("Average E <> E(90) = %4.3f+-%4.3f",average_E,sigma_E_average),1);
    gr_CHSH_all->Draw("AP");
    gr_CHSH_all->Write();    
    canvas_CHSH->SaveAs(result_path+".pdf",".pdf");

    TCanvas *global_canvas_CHSH = new TCanvas ( "global_canvas_CHSH", "global_canvas_CHSH");
    global_canvas_CHSH->Divide(2);

    TGraphErrors * global_gr_CHSH = new TGraphErrors(module_angles,phi_angle,CHSH_cl.Get_global_CHSH(), angle_arr_err, CHSH_cl.Get_global_CHSH_error());
	global_gr_CHSH->Fit("my_fit","R");
    a_par_CHSH = abs(my_fit->GetParameter(0));
    a_error_CHSH = abs(my_fit->GetParError(0));
    global_canvas_CHSH->cd(1);
    graph_like_ivashkin_wants_it(global_gr_CHSH,"angle [degrees]","S", Form("Sum N <> E(90) = %4.3f+-%4.3f",global_E_coeff(4,NumEvents,"both"), sqrt(global_sqr_E_error(4,NumEvents, "both"))),1);
    global_gr_CHSH->Draw("AP");
    
    TGraphErrors * global_gr_CHSH_all = new TGraphErrors(all_angles-1,phi_angle,CHSH_cl.Get_global_CHSH_all(), angle_arr_err, CHSH_cl.Get_global_CHSH_error_all());
	
    global_gr_CHSH_all->Fit("my_fit","R");
    a_par_CHSH = abs(my_fit->GetParameter(0));
    a_error_CHSH = abs(my_fit->GetParError(0));
    global_canvas_CHSH->cd(2);
    graph_like_ivashkin_wants_it(global_gr_CHSH_all,"angle [degrees]","S", Form("Sum N <> E(90) = %4.3f+-%4.3f",global_E_coeff(4,NumEvents,"both"), sqrt(global_sqr_E_error(4,NumEvents, "both"))),1);
    
    global_gr_CHSH_all->Draw("AP");
    global_gr_CHSH_all->Write();
    global_canvas_CHSH->SaveAs(result_path+".pdf",".pdf");
  
   //////////////////////// Draw E_coeffs

    TCanvas *E_canv = new TCanvas ( "global_E", "global_E");
    E_canv->Divide(2);
    TGraphErrors * E_gr = new TGraphErrors(module_angles,phi_angle,CHSH_cl.Get_global_Correlation_E(), angle_arr_err, CHSH_cl.Get_global_Corr_error());
    E_gr->Fit("e_fit","R");
    E_canv->cd(1);
    graph_like_ivashkin_wants_it(E_gr,"angle [degrees]","E", "Average E for mirroring angles", 1);
    E_gr->Draw("AP");

    TGraphErrors * E_gr_all = new TGraphErrors(all_angles-1,phi_angle,CHSH_cl.Get_global_Correlation_E_all(), angle_arr_err, CHSH_cl.Get_global_Corr_error_all());
    E_gr_all->Fit("e_fit","R");
    E_canv->cd(2);
    graph_like_ivashkin_wants_it(E_gr_all,"angle [degrees]","E", "Average E for all angles", 1);
    E_gr_all->Draw("AP");
    E_gr_all->Write();    
    E_canv->SaveAs(result_path+".pdf",".pdf");

    TCanvas *E_total = new TCanvas ( "total_E", "total_E");
    E_total->Divide(2);
    TGraphErrors * E_gr_total = new TGraphErrors(module_angles,phi_angle, CHSH_cl.Get_total_Correlation_E(), angle_arr_err, CHSH_cl.Get_total_Corr_error());
    E_gr_total->Fit("e_fit","R");
    E_total->cd(1);
    graph_like_ivashkin_wants_it(E_gr_total,"angle [degrees]","E", "Summ E for mirroring angles", 1);
    E_gr_total->Draw("AP");

    TGraphErrors * E_gr_all_total = new TGraphErrors(all_angles-1,phi_angle, CHSH_cl.Get_total_Correlation_E_all(), angle_arr_err, CHSH_cl.Get_total_Corr_error_all());
    E_gr_all_total->Fit("e_fit","R");
    E_total->cd(2);
    graph_like_ivashkin_wants_it(E_gr_all_total,"angle [degrees]","E", "Summ E for all angles", 1);
    E_gr_all_total->Draw("AP");
    E_gr_all_total->Write();
    E_total->SaveAs(result_path+".pdf",".pdf");


   ///////////////////////Draw Asymmetry

    for (Int_t i = 0; i < 16; i++) total_events_for_angle_diff_err[i] = pow(total_events_for_angle_diff[i],0.5);

    TGraphErrors * gr = new TGraphErrors(16,angle_arr,total_events_for_angle_diff, angle_arr_err, total_events_for_angle_diff_err);
	TF1 *sin_fit = new TF1 ("sin_fit","[0]+[1]*cos(3.14159265359/180*2*x)",-10,360);
	gr->Fit("sin_fit","R");
    Float_t a_par = abs(sin_fit->GetParameter(0));
    Float_t b_par = abs(sin_fit->GetParameter(1));
    Float_t a_error = abs(sin_fit->GetParError(0));
    Float_t b_error = abs(sin_fit->GetParError(1));
    Float_t diff = (a_par+b_par)/(a_par-b_par);
    Float_t diff_err = 2*diff*sqrt(pow((a_error/a_par),2)+pow((b_error/b_par),2))/(a_par/b_par-b_par/a_par);
    //Float_t diff_err = 2*a_par/(pow(a_par-b_par,2))*sqrt(pow((a_error*b_par/a_par),2)+pow((b_error),2));

    TCanvas *canvas = new TCanvas ( "canvas", "canvas");
    canvas->cd();
    graph_like_ivashkin_wants_it(gr,"angle [degrees]","Counts", Form("%5.1f-%5.1f*cos(2x) ratio = %4.3f +- %4.3f",a_par,b_par,diff,diff_err),1);
    gr->Draw("AP");
    gr->Write();
    canvas->SaveAs(result_path+".pdf)",".pdf");    
    result_root->Close();    
    }