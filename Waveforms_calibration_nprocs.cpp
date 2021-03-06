#include <iostream>
#include <fstream>
#include "like_ivashkin_wants_it.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TF1.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TObjString.h>
#include "ChannelEntry.h"
#include "Coeffs.h"
#include <stdlib.h>

#define UsedScatterer 1

const int GATE_BEG = 50;
const int GATE_END = 150;
const Int_t N_zl_intervals = 2;
const Int_t Gate_end_for_main_scatterers = 71;
const Int_t Intermediate_gate = 82;

// const int GATE_BEG = 150;
// const int GATE_END = 300;
// const Int_t N_zl_intervals = 2;
// const Int_t Gate_end_for_main_scatterers = 300;
////////////////////////
///////////////////////Define fit constants range
Int_t l_integr_NaI = 30000;
Int_t r_integr_NaI = 300000;

///////////////////////
///////////////////////

void Waveforms_calibration_nprocs(
    TString source_path = "/home/doc/entanglement/new_files_data/",
    TString file_name = "",TString temp_DIR = "", TString postfix = "_0")
{

    #if UsedScatterer
    const Int_t total_channels = 35;
    #else
    const Int_t total_channels = 34;
    #endif

    ChannelEntry channel_info[total_channels];
    short_energy_ChannelEntry short_channel_info[total_channels];
    diff_short_energy_ChannelEntry diff_short_channel_info[total_channels];
    Coeffs_struct cal_coeff[total_channels];
    
    //PMT_tree->SetBranchAddress(Form("channel_%d", channel_number), &channel_info[channel_number]);    
    TFile *combined_root = new TFile (source_path+temp_DIR+"calibrated_time" + postfix+".root", "RECREATE");
    TTree *combined_tree = new TTree ("adc64_data","adc64_data");
    TTree *coefficients_tree = new TTree ("coefficients","coefficients");
    TTree *diff_tree = new TTree ("diff_tree","diff_tree");

    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
        (short_channel_info[channel_number]).CreateBranch(combined_tree, channel_number);    
    for (Int_t channel_number = 34; channel_number < 35; channel_number++)
        (diff_short_channel_info[channel_number]).CreateBranch(combined_tree, channel_number);

    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
        (cal_coeff[channel_number]).CreateBranch(coefficients_tree, channel_number);
        
    Float_t interval_width = 1.;


        Int_t gr_point_counter = 0;


                TChain *PMT_tree = new TChain;
				PMT_tree->AddFile( (file_name + "/adc64_data").Data() );
                for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                    (channel_info[channel_number]).SetBranch(PMT_tree,channel_number);
                    Int_t total_entries = PMT_tree->GetEntries()/1;
                TLeaf *tf = PMT_tree->GetLeaf("time_in_seconds");
///////////////////////
///////////////////////
                Int_t intervals = 1;
                Int_t int_width = total_entries;
                for (intervals = 1; intervals < 1000; intervals++) 
                {
                    int_width=total_entries/intervals;
                    if (int_width < 2000000) break;
                }
                
                for (Int_t num_int = 0; num_int < intervals; num_int++)
                {
                    
                    Int_t low_cut[total_channels] = {0}; Int_t high_cut[total_channels] = {0};
                    Float_t sigma_i[total_channels] = {0}; Float_t mean_int[total_channels] = {0};
                    Int_t Events_for_cuts_left = num_int*int_width;
                    Int_t Events_for_cuts_right = (num_int+1)*int_width;
                    if (num_int==intervals-1) Events_for_cuts_right = total_entries;
                    Float_t time_of_run = 0; tf->GetBranch()->GetEntry((Int_t)(Events_for_cuts_left+Events_for_cuts_right)/2); time_of_run=tf->GetValue()/3600;

                    TH1F *peak_histo[total_channels];
                    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                    {
                        TString hernya = Form(postfix+"_integral_in_gate_ADC_ch_%i_int_%i",channel_number,num_int);
                        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,l_integr_NaI,r_integr_NaI);
                    }

                    for (Long_t iEvent = Events_for_cuts_left; iEvent < Events_for_cuts_right; iEvent++)
                    {

                        //for (int ch = 0; ch < total_channels; ch++) channel_info[ch].Initialize();
                        PMT_tree->GetEntry(iEvent);
                        //cout << " " << endl;
                        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                        {
                            Int_t zero_level = channel_info[channel_number].Get_Zero_Level(GATE_BEG);
                            Int_t gateIntegral = channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END);
                            if (gateIntegral > 4000
                            && channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, GATE_END) > 1000
                            && channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, GATE_END) < 60000

                        // #if UsedScatterer
                        //         && abs(channel_info[34].Get_peak_position(channel_info[34].Get_Zero_Level(GATE_BEG))-channel_info[32].Get_peak_position(channel_info[32].Get_Zero_Level(GATE_BEG))) > 10
                        //         && channel_info[34].Get_Amplitude(channel_info[34].Get_Zero_Level(GATE_BEG), GATE_BEG, GATE_END) < 400
                        // #endif
                            ) 
                            {
                                cout<<"ch="<<channel_number<<"Amp = "<<channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, GATE_END)<<" zl="<<zero_level<<endl;
                                peak_histo[channel_number]->Fill(gateIntegral);
                            }
                        }
                        for (int ch = 0; ch < total_channels; ch++) channel_info[ch].Initialize();
                    }

                    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                    {
                        //PMT_tree->Draw(Form("channel_%i.integral_in_gate >> %s",channel_number, hernya.Data()), Form("channel_%i.integral_in_gate > 4000",channel_number));
                        TF1 *gaus_fit = new TF1 ("gaus_func","gaus",
                        peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())-20000,
                        peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())+20000);
                        peak_histo[channel_number]->Fit("gaus_func","R");

                        TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
                        (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
                        (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
                        peak_histo[channel_number]->Fit("gaus_func_2","R");

                        mean_int[channel_number] = gaus_2_fit->GetParameter(1);
                        sigma_i[channel_number] = gaus_2_fit->GetParameter(2);
                        low_cut[channel_number] = gaus_2_fit->GetParameter(1)-interval_width*gaus_2_fit->GetParameter(2);
                        high_cut[channel_number] = gaus_2_fit->GetParameter(1)+interval_width*gaus_2_fit->GetParameter(2);
                        peak_histo[channel_number]->Write();

                        cal_coeff[channel_number].peak_for_calibration = mean_int[channel_number];
                        cal_coeff[channel_number].time_in_hours = time_of_run;
                        delete peak_histo[channel_number];

                    }
                //////////////////////////////////
                //Calibrating intermediate scatterer
                    for (Long_t iEvent = Events_for_cuts_left; iEvent < Events_for_cuts_right; iEvent++)
                    {
                    #if UsedScatterer
                        PMT_tree->GetEntry(iEvent);
                        Int_t zero_level = channel_info[total_channels-1].Get_Zero_Level(GATE_BEG);
                        Int_t gateIntegral = channel_info[total_channels-1].Get_Integral_in_Gate(zero_level, GATE_BEG, Intermediate_gate);
                        if (gateIntegral > 4000
                        && channel_info[total_channels-1].Get_Amplitude(zero_level, GATE_BEG, Intermediate_gate) > 1000
                        && channel_info[total_channels-1].Get_Amplitude(zero_level, GATE_BEG, Intermediate_gate) < 60000
                        // && abs(channel_info[34].Get_peak_position(channel_info[34].Get_Zero_Level(GATE_BEG))-channel_info[32].Get_peak_position(channel_info[32].Get_Zero_Level(GATE_BEG))) > 10
                        // && channel_info[34].Get_Amplitude(channel_info[34].Get_Zero_Level(GATE_BEG), GATE_BEG, GATE_END) < 400
                        ) 
                        {
                            //cout<<"ch="<<total_channels-1<<"Amp = "<<channel_info[total_channels-1].Get_Amplitude(zero_level, GATE_BEG, GATE_END)<<" zl="<<zero_level<<endl;
                            peak_histo[total_channels-1]->Fill(gateIntegral);
                        }
                        channel_info[total_channels-1].Initialize();

                    }

                        //PMT_tree->Draw(Form("channel_%i.integral_in_gate >> %s",total_channels-1, hernya.Data()), Form("channel_%i.integral_in_gate > 4000",total_channels-1));
                        TF1 *gaus_fit = new TF1 ("gaus_func","gaus",
                        peak_histo[total_channels-1]->GetBinCenter(peak_histo[total_channels-1]->GetMaximumBin())-20000,
                        peak_histo[total_channels-1]->GetBinCenter(peak_histo[total_channels-1]->GetMaximumBin())+20000);
                        peak_histo[total_channels-1]->Fit("gaus_func","R");

                        TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
                        (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
                        (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
                        peak_histo[total_channels-1]->Fit("gaus_func_2","R");

                        mean_int[total_channels-1] = gaus_2_fit->GetParameter(1);
                        sigma_i[total_channels-1] = gaus_2_fit->GetParameter(2);
                        low_cut[total_channels-1] = gaus_2_fit->GetParameter(1)-interval_width*gaus_2_fit->GetParameter(2);
                        high_cut[total_channels-1] = gaus_2_fit->GetParameter(1)+interval_width*gaus_2_fit->GetParameter(2);
                        peak_histo[total_channels-1]->Write();
                        cal_coeff[total_channels-1].peak_for_calibration = mean_int[total_channels-1];
                        cal_coeff[total_channels-1].time_in_hours = time_of_run;
                        delete peak_histo[total_channels-1];

                    #endif            
                    //
                ////////////////////////////////////
                    TString left_name = postfix+Form("_integral_in_gate_left_arm_scatterer_interval_number_%i",num_int);
                    TString right_name = postfix+Form("_integral_in_gate_right_arm_scatterer_interval_number_%i",num_int);
                    Int_t scatterer_left_boarder_peak_search = 0000;


                    peak_histo[32] = new TH1F (left_name.Data(),left_name.Data(),200,scatterer_left_boarder_peak_search,500000);
                    peak_histo[33] = new TH1F (right_name.Data(),right_name.Data(),200,scatterer_left_boarder_peak_search,500000);

                        //TH1F *peak_histo[32] = nullptr;
                    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                    {
                        TString hernya = Form(postfix+"_integral_in_gate_scatterer_ch_cut_%i_interval_%i",channel_number,num_int);
                        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),200,scatterer_left_boarder_peak_search,500000);
                    }


                    for (Long_t iEvent = Events_for_cuts_left; iEvent < Events_for_cuts_right; iEvent++)
                    {
                        //for (int ch = 0; ch < total_channels; ch++) channel_info[ch].Initialize();

                        PMT_tree->GetEntry(iEvent);
                        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                        {
                            Int_t zero_level = channel_info[channel_number].Get_Zero_Level(GATE_BEG);

                            Int_t sc_number = 0;
                            if (channel_number < 16) sc_number = 32;
                            if (channel_number >= 16) sc_number = 33;
                            
                            if (channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END) > low_cut[channel_number]
                            && channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END) < high_cut[channel_number]
                            //&& abs(channel_info[32].Get_peak_position(zero_level, GATE_BEG, GATE_END)-channel_info[33].Get_peak_position(zero_level, GATE_BEG, GATE_END)) < 2
                            && (channel_info[channel_number].Get_peak_position(zero_level)-channel_info[sc_number].Get_peak_position(channel_info[sc_number].Get_Zero_Level(GATE_BEG))) > 0
                            && (channel_info[channel_number].Get_peak_position(zero_level)-channel_info[sc_number].Get_peak_position(channel_info[sc_number].Get_Zero_Level(GATE_BEG))) < 14

                            && channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, GATE_END) < 60000
                            && channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, GATE_END) > 400
                            
                            //&& channel_info[sc_number].Get_peak_position(zero_level) < 31000
                            

                #if UsedScatterer
                            && abs(channel_info[34].Get_peak_position(channel_info[34].Get_Zero_Level(GATE_BEG))-channel_info[32].Get_peak_position(channel_info[32].Get_Zero_Level(GATE_BEG))) > 10
                            && channel_info[34].Get_Amplitude(channel_info[34].Get_Zero_Level(GATE_BEG), GATE_BEG, GATE_END) < 400
                #endif
                            )
                            peak_histo[channel_number]->Fill(channel_info[sc_number].Get_Integral_in_Gate(channel_info[sc_number].Get_Zero_Level(GATE_BEG), GATE_BEG, Gate_end_for_main_scatterers));
                        }
                    }
                    for (int ch = 0; ch < total_channels; ch++) channel_info[ch].Initialize();


                    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                    {

                        //peak_histo[channel_number]->Write();
                        if (channel_number < 16) peak_histo[32]->Add(peak_histo[channel_number]);

                        if (channel_number >= 16) peak_histo[33]->Add(peak_histo[channel_number]);

                        delete peak_histo[channel_number];
                    }

                    for (Int_t channel_number = 32; channel_number < 34; channel_number++)
                    {
                        TF1 *gaus_fit = new TF1 ("gaus_func","gaus",
                        peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())-20000,
                        peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())+20000);
                        peak_histo[channel_number]->Fit("gaus_func","R");

                        TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
                        (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
                        (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
                        peak_histo[channel_number]->Fit("gaus_func_2","R");
                        mean_int[channel_number] = gaus_2_fit->GetParameter(1);
                        peak_histo[channel_number]->Write();
                        cal_coeff[channel_number].peak_for_calibration = mean_int[channel_number];
                        cal_coeff[channel_number].time_in_hours = time_of_run;
                        delete peak_histo[channel_number];

                    }
                    gr_point_counter++;

    //////////////////////
    //////////////////////Fill_clibrated_tree
                    for (Long_t entry_num = Events_for_cuts_left;  entry_num< Events_for_cuts_right; entry_num++)
                    {
                        //for (int ch = 0; ch < total_channels; ch++) channel_info[ch].Initialize();

                        PMT_tree->GetEntry(entry_num);
                        for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                        {
                            short_channel_info[channel_number].Initialize();
                            Int_t zero_level = channel_info[channel_number].Get_Zero_Level(GATE_BEG);
                            short_channel_info[channel_number].zl_rms = channel_info[channel_number].Get_Zero_Level_RMS(GATE_BEG);
                            Int_t CH_GATE_END = GATE_END; if (channel_number == 32 || channel_number ==33) CH_GATE_END = Gate_end_for_main_scatterers; if (channel_number == 34) CH_GATE_END = Intermediate_gate;

                            Int_t energy_in_peak = 255;
                            if (channel_number == 32 || channel_number == 33) energy_in_peak = 255;
                            if (channel_number < 34) 
                            {
                                short_channel_info[channel_number].integral_in_gate = 
                                (float)energy_in_peak/mean_int[channel_number]*
                                channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, CH_GATE_END);
                                if ((entry_num + 1) % 10000 == 0) cout << channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, CH_GATE_END) << " " << short_channel_info[channel_number].integral_in_gate << endl;
                            }
                            else 
                            {   if (channel_number == 34) 
                                {
                                    channel_info[channel_number].SplineWf();
                                    channel_info[channel_number].SplineWf();

                                    short_channel_info[channel_number].integral_in_gate = 
                                    (Float_t)170.3/mean_int[34]*channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, CH_GATE_END);
                                }
                            }
                        //Int_t inv_amp = channel_info[channel_number].Get_Amplitude(zero_level);
                            //short &wf_size = channel_info[channel_number].wf_size;
                            //Short_t *wf = channel_info[channel_number].wf;

                            if ((entry_num + 1) % 100000 == 0) cout << channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, CH_GATE_END) << endl<<endl<<endl;
                            short_channel_info[channel_number].amp = channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, CH_GATE_END);//inv_amp;

                            short_channel_info[channel_number].peak_pos = 
                            channel_info[channel_number].Get_peak_position_gauss(zero_level, short_channel_info[channel_number].amp, CH_GATE_END);
                            if (channel_info[channel_number].wf_size == 350 && (entry_num + 1) % 100000 == 0) cout << short_channel_info[channel_number].peak_pos << " " << channel_info[channel_number].Get_peak_position_gauss(zero_level, short_channel_info[channel_number].amp, CH_GATE_END) << "; channel_number = " << channel_number<< endl;
                                
                            diff_short_channel_info[channel_number].Initialize();
                            if (channel_number == 34 && channel_info[channel_number].wf_size != 0) 
                            {   Short_t max = 0; Short_t min = 0; Short_t maxt = 0; Short_t mint = 0;
                                channel_info[channel_number].DiffWf();
                                channel_info[channel_number].FindDiffWfPars(min,mint,max,maxt);
                                diff_short_channel_info[channel_number].min_diff = min;
                                diff_short_channel_info[channel_number].min_diff_time = mint; 
                                diff_short_channel_info[channel_number].max_diff = max; 
                                diff_short_channel_info[channel_number].max_diff_time  = maxt;                          
                                if (max > 1000) cout << "min_diff="<<min<<";min_time="<<mint<<";max_diff="<<max<<";max_time="<<maxt<<endl;      
                            }          
                        }
                        combined_tree->Fill();
                        diff_tree->Fill();
                        if ((entry_num + 1) % 100000 == 0) cout << (float)entry_num/total_entries*100<<"%"<<endl;
                        for (int ch = 0; ch < total_channels; ch++) channel_info[ch].Initialize();
                    }

                coefficients_tree->Fill();
                
                }
                delete PMT_tree;
    ////////////////////////

        //for (Int_t ch = 0; ch < total_channels; ch++) gr[ch]->Write("AP");
    combined_tree->Write();
    diff_tree->Write();
    coefficients_tree->Write();
    combined_root->Close();
}
        



