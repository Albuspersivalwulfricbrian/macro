#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>

#include "like_ivashkin_wants_it.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

#include <TChain.h>
#include <TObjArray.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TGraph.h>
#include <TObjString.h>
#include "../calo_analysis/Readme_reader.h"
#include "ChannelEntry.h"
#define UsedScatterer 1

void Draw_waveforms()
{

    TString source_path = "/home/doc/entanglement/new_files_data/new_root/temp_directory/new_scat/";
    TChain *PMT_tree = new TChain;

	TObjArray files_names;
	TSystemDirectory dir(source_path, source_path);
	TList *files = dir.GetListOfFiles();
	if(files) 
	{
		TSystemFile *file; 
		TString fname; 
		TIter next(files); 
		while((file=(TSystemFile*)next())) 
		{ 
			fname = file->GetName(); 
			if(!file->IsDirectory() && fname.EndsWith("07a8de9a_20220118_152151.root")) 
			{
				files_names.Add(new TObjString(fname));
				PMT_tree->AddFile( (source_path + fname + "/adc64_data").Data() );
			}
		}
	}
    Int_t total_entries = PMT_tree->GetEntries();

#if UsedScatterer
    const Int_t total_channels = 35;
#else 
    const Int_t total_channels = 34;
#endif

     Int_t looK_channel = 5;
    const Int_t tot_canv = 300;

    Int_t global_counter = 0;
    Int_t samples[2050];
    ChannelEntry *channel_info = new ChannelEntry[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
	    (channel_info+ch)->SetBranch(PMT_tree, ch);

    for (Int_t i = 0; i < 2050; i++) samples[i] = i;
        for (Int_t entryNum = 0; entryNum < total_entries; entryNum++)
        {
            PMT_tree->GetEntry(entryNum);
            for (Int_t channel_number = 34; channel_number < 35; channel_number++)
            {
            const Int_t waveform_size = channel_info[channel_number].wf_size;
            Short_t flag = 0;
            if (global_counter < tot_canv)
            {
                if (waveform_size > 10)
                {
                    const Int_t sz = waveform_size;
                    Int_t wf[1000] = {0};
                    TCanvas *canv = new TCanvas("canv","canv");
                    Float_t rn = (Float_t)std::rand()/numeric_limits<Int_t>::max();

                    if (global_counter < tot_canv && rn < 0.1) 
                    {


                    for (Int_t i = 0; i < waveform_size; i++) wf[i] = channel_info[channel_number].wf[i] - channel_info[channel_number].Get_Zero_Level(50);

                        TGraph *gr = new TGraph(waveform_size,samples,wf);
                        graph_like_ivashkin_wants_it(gr,"Time [a.u.]","ADC channels", Form("low_ampl_waveform_%i_ch_%i",global_counter,channel_number),2);
                        gr->GetXaxis()->SetRangeUser(0,150);
                        gr->Draw("AL");
                        channel_info[channel_number].SplineWf();
                        channel_info[channel_number].SplineWf();   
                    for (Int_t i = 0; i < waveform_size; i++) wf[i] = channel_info[channel_number].wf[i] - channel_info[channel_number].Get_Zero_Level(50);
                     
                        TGraph *gr1 = new TGraph(waveform_size,samples,wf);
                        graph_like_ivashkin_wants_it(gr1,"Time [a.u.]","ADC channels", Form("low_ampl_waveform_%i_ch_%i",global_counter,channel_number),2);
                        //gr->GetXaxis()->SetRangeUser(0,150);
                        gr1->SetLineColor(2);

                        gr1->Draw("LS");

                        channel_info[channel_number].DiffWf();
                    for (Int_t i = 0; i < waveform_size; i++) wf[i] = channel_info[channel_number].wf[i];

                        TGraph *gr2 = new TGraph(waveform_size,samples,wf);
                        graph_like_ivashkin_wants_it(gr2,"Time [a.u.]","ADC channels", Form("low_ampl_waveform_%i_ch_%i",global_counter,channel_number),2);
                        //gr->GetXaxis()->SetRangeUser(0,150);
                        gr2->SetLineColor(4);

                        gr2->Draw("LS");
                        global_counter++;
                        flag++;
                    }    
                    if (flag > 0) canv->SaveAs(source_path+"waveforms.pdf(","pdf");
                }
            }
            else
            {entryNum = total_entries;}
            }
        }
                    TCanvas *canv = new TCanvas("canv","canv");
                    canv->SaveAs(source_path+"waveforms.pdf)","pdf");

}