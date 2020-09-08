// this script uses output of ana.C as input. Note, it has two inputs, one from raw data, the other from dataprep (caldata).
// It plots the normalized signal and noise rms, and signal-to-noise ratio in the same canvase.


double fit_gaus(double *x, double *par) {
  double m=par[0];
  double s=par[1];
  double n=par[2];

  double g = n*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
  // double g = n/(s*sqrt(2*3.14159))*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
  
  return g;
}

TF1* fit_rms(TH1F* h, int color) {
  // pre-fit parameters
  float pre_mean=h->GetBinCenter(h->GetMaximumBin());
  float pre_max=h->GetMaximum();
  float pre_rms=h->GetRMS();
  cout<<"mean: "<<pre_mean<<endl;
  cout<<"rms: "<<pre_rms<<endl;
  cout<<"max: "<<pre_max<<endl;
  cout<<""<<endl;

  // 1st fitting
  TF1 *gg=new TF1("gg",fit_gaus,pre_mean-2.*pre_rms,pre_mean+2.*pre_rms,3);
  gg->SetParameter(0,pre_mean);
  gg->SetParameter(1,pre_rms);
  gg->SetParameter(2,pre_max);
  gg->SetLineColor(color);
  gg->SetLineStyle(2);
  h->Fit("gg","remn");

  // 2nd fitting
  TF1 *g=new TF1("g",fit_gaus,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
  g->SetParameter(0,gg->GetParameter(0));
  g->SetParameter(1,gg->GetParameter(1));
  g->SetParameter(2,gg->GetParameter(2));
  g->SetLineColor(color);
  g->SetLineStyle(1);
  g->SetLineWidth(1);
  h->Fit("g","remn");

  return g;
}

TF1* fit_s2n(TH1F* h, int color) {
  // pre-fit parameters
  float pre_mean=h->GetBinCenter(h->GetMaximumBin());
  //float pre_max=h->GetMaximum();
  float pre_max=h->GetBinContent(h->GetMaximumBin());
  float pre_rms=h->GetRMS();
  cout<<"mean: "<<pre_mean<<endl;
  cout<<"rms: "<<pre_rms<<endl;
  cout<<"max: "<<pre_max<<endl;
  cout<<""<<endl;
  double lfit = pre_mean-1.*pre_rms;
  double hfit = pre_mean+1.*pre_rms;
  if (lfit<0) {
    lfit = 0.5;
    hfit = pre_mean + pre_mean - lfit;
  }

  pre_rms = pre_rms/1.1;

  // 1st fitting
  //TF1 *gg=new TF1("gg",fit_gaus,pre_mean-1.*pre_rms,pre_mean+1.*pre_rms,3);
  TF1 *gg=new TF1("gg",fit_gaus, lfit, hfit, 3);
  gg->SetParameter(0,pre_mean);
  gg->SetParameter(1,pre_rms);
  gg->SetParameter(2,pre_max);
  gg->SetLineColor(color);
  gg->SetLineStyle(2);
  h->Fit("gg","remn");
  cout<<"1st fitting"<<endl;
  cout<<"mean: "<<gg->GetParameter(0)<<endl;
  cout<<"rms: "<<gg->GetParameter(1)<<endl;
  cout<<"max: "<<gg->GetParameter(2)<<endl;
  cout<<""<<endl;

  // 2nd fitting
  TF1 *g=new TF1("g",fit_gaus,gg->GetParameter(0)-.5*gg->GetParameter(1),gg->GetParameter(0)+.5*gg->GetParameter(1),3);
  //TF1 *g=new TF1("g",fit_gaus,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
  g->SetParameter(0,gg->GetParameter(0));
  g->SetParameter(1,gg->GetParameter(1));
  g->SetParameter(2,gg->GetParameter(2));
  g->SetLineColor(color);
  g->SetLineStyle(1);
  g->SetLineWidth(1);

  h->Fit("g","remn");

  cout<<"2st fitting"<<endl;
  cout<<"mean: "<<g->GetParameter(0)<<endl;
  cout<<"rms: "<<g->GetParameter(1)<<endl;
  cout<<"max: "<<g->GetParameter(2)<<endl;
  cout<<""<<endl;

  return g;
}





void compare_plot() {

  // plot options
  const char * plane_options[3] = {"u","v","y"};
  const char * draw_options[10] = {"same hist", "SAME hist", "SAME hist", "SAME hist", "SAME", "SAME", "SAME", "SAME", "SAME", "SAME"};
  const int color_options[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 20};


  // input
  TFile *input_1 = TFile::Open("../ana_1/output_ana.root"); // raw
  TFile *input_2 = TFile::Open("../ana_2/output_ana.root"); // caldata

  TH1F* h1_s2n_raw[3];
  TH1F* h1_rms_raw[3];
  TH1F* h1_rmsfit_raw[3];
  TH1F* h1_amp_raw[3];

  TH1F* h1_s2n_dataprep[3];
  TH1F* h1_rms_dataprep[3];
  TH1F* h1_rmsfit_dataprep[3];
  TH1F* h1_amp_dataprep[3];

  for (int p=0; p<3; p++) {
    h1_s2n_raw[p] = (TH1F*)input_1->Get(TString::Format("h1_s2n_plane_%d", p));
    h1_s2n_raw[p]->Rebin(4);
    h1_s2n_raw[p]->Scale(1./h1_s2n_raw[p]->GetMaximum());
    h1_rms_raw[p] = (TH1F*)input_1->Get(TString::Format("h1_rms_plane_%d", p)); 
    h1_rmsfit_raw[p] = (TH1F*)input_1->Get(TString::Format("h1_rmsfit_plane_%d", p)); 
    h1_amp_raw[p] = (TH1F*)input_1->Get(TString::Format("h1_amp_plane_%d", p)); 
    h1_amp_raw[p]->Rebin(2);
    h1_amp_raw[p]->Scale(1./h1_amp_raw[p]->GetMaximum());

    h1_s2n_dataprep[p] = (TH1F*)input_2->Get(TString::Format("h1_s2n_plane_%d", p));
    h1_s2n_dataprep[p]->Rebin(4);
    h1_s2n_dataprep[p]->Scale(1./h1_s2n_dataprep[p]->GetMaximum());
    h1_rms_dataprep[p] = (TH1F*)input_2->Get(TString::Format("h1_rms_plane_%d", p)); 
    h1_rmsfit_dataprep[p] = (TH1F*)input_2->Get(TString::Format("h1_rmsfit_plane_%d", p)); 
    h1_amp_dataprep[p] = (TH1F*)input_2->Get(TString::Format("h1_amp_plane_%d", p)); 
    h1_amp_dataprep[p]->Rebin(2);
    h1_amp_dataprep[p]->Scale(1./h1_amp_dataprep[p]->GetMaximum());
  }

  // rms and signal
  TCanvas *c1_sn[3];

  for (int p=0; p<3; p++) {
    c1_sn[p] = new TCanvas(TString::Format("c1_sn_%d",p), TString::Format("c1_sn_%d",p), 800, 600);
    
    c1_sn[p]->SetLogx();

    TH2F* frame = new TH2F("frame", "", 400, 0.5, 400.5, 125, 0, 1.25); 
    frame->SetStats(0);
    frame->GetXaxis()->CenterTitle();
    frame->GetXaxis()->SetTitle("ADC counts");
    frame->GetYaxis()->CenterTitle();
    frame->GetYaxis()->SetTitle("Arbitrary units");
    frame->Draw();
    
    auto lg_sn = new TLegend(0.28, 0.45, 0.48, 0.55);
    auto lg_sn_dataprep = new TLegend(0.01, 0.75, 0.7, 0.85);

    h1_amp_raw[p]->SetStats(0);
    h1_amp_raw[p]->SetLineColor(1);
    h1_amp_raw[p]->SetLineWidth(4);
    h1_amp_raw[p]->SetLineStyle(9);
    h1_amp_raw[p]->Draw("hist same");

    h1_amp_dataprep[p]->SetStats(0);
    h1_amp_dataprep[p]->SetLineColor(kRed);
    h1_amp_dataprep[p]->SetLineWidth(4);
    h1_amp_dataprep[p]->SetLineStyle(1);
    h1_amp_dataprep[p]->Draw("hist same");

    h1_rms_raw[p]->SetStats(0);
    h1_rms_raw[p]->Scale(1./h1_rms_raw[p]->GetMaximum());
    h1_rms_raw[p]->SetLineColor(1);
    h1_rms_raw[p]->SetLineWidth(4);
    h1_rms_raw[p]->SetLineStyle(9);
    h1_rms_raw[p]->Draw("HIST SAME");

    h1_rms_dataprep[p]->SetStats(0);
    h1_rms_dataprep[p]->Scale(1./h1_rms_dataprep[p]->GetMaximum());
    h1_rms_dataprep[p]->SetLineColor(kRed);
    h1_rms_dataprep[p]->SetLineWidth(4);
    h1_rms_dataprep[p]->SetLineStyle(1);
    h1_rms_dataprep[p]->Draw("HIST SAME");

    TLegendEntry* leg_raw = lg_sn->AddEntry(h1_rms_raw[p], "Raw", "l");
    leg_raw->SetTextColor(1);
    TLegendEntry* leg_dataprep = lg_sn->AddEntry(h1_rms_dataprep[p], "Noise-filtered", "l");
    leg_dataprep->SetTextColor(kRed);

    lg_sn->SetTextSize(0.05);
    lg_sn->SetFillStyle(0);
    lg_sn->SetLineWidth(0);
    lg_sn->Draw();
    
    TString txt_sn(Form("Noise sigma                      Signal size"));
    lg_sn_dataprep->AddEntry(txt_sn, txt_sn.Data(), "");
    lg_sn_dataprep->SetTextSize(0.05);
    lg_sn_dataprep->SetFillStyle(0);
    lg_sn_dataprep->SetLineWidth(0);
    lg_sn_dataprep->Draw();

    TLatex **txt_iceberg = new TLatex*[1];
    txt_iceberg[0] = new TLatex(0.5, 1.27, "IceBerg");
    txt_iceberg[0]->SetTextColor(1);
    txt_iceberg[0]->Draw();

    TLatex **txt_plane = new TLatex*[1];
    txt_plane[0] = new TLatex(80, 1.27, Form("#bf{Data, plane %s}",plane_options[p]));
    txt_plane[0]->SetTextColor(1);
    txt_plane[0]->Draw();

    c1_sn[p]->SaveAs(TString::Format("sn_plane_%d.png",p));

    delete frame;
  }


  // s2n
  TCanvas *c1_s2n  = new TCanvas("c1_sn3", "c1_sn", 800, 600);

  TH2F* frame_s2n = new TH2F("frame_s2n", "", 160, 0.0, 160., 105, 0, 1.05); 
  frame_s2n->SetStats(0);
  frame_s2n->GetXaxis()->CenterTitle();
  frame_s2n->GetXaxis()->SetTitle("Angle-Corrected Signal-to-Noise Ratio");
  frame_s2n->GetYaxis()->CenterTitle();
  frame_s2n->GetYaxis()->SetTitle("Arbitrary units");
  frame_s2n->Draw();

  auto lg_s2n = new TLegend(0.5, 0.4, 0.9, 0.85);

  for (int p=0; p<3; p++) {

    h1_s2n_raw[p]->SetStats(0);
    h1_s2n_raw[p]->SetLineColor(5-color_options[p]);
    h1_s2n_raw[p]->SetLineWidth(4);
    h1_s2n_raw[p]->SetLineStyle(9);
    h1_s2n_raw[p]->Draw("hist same");

    h1_s2n_dataprep[p]->SetStats(0);
    h1_s2n_dataprep[p]->SetLineColor(5-color_options[p]);
    h1_s2n_dataprep[p]->SetLineWidth(4);
    h1_s2n_dataprep[p]->SetLineStyle(1);
    h1_s2n_dataprep[p]->Draw("hist same");

    TLegendEntry* leg_raw = lg_s2n->AddEntry(h1_s2n_raw[p], Form("Plane %s (Raw)", plane_options[p]), "l");
    leg_raw->SetTextColor(5-color_options[p]);

    TLegendEntry* leg_dataprep = lg_s2n->AddEntry(h1_s2n_dataprep[p], Form("Plane %s (Noise-filtered)", plane_options[p]), "l");
    leg_dataprep->SetTextColor(5-color_options[p]);
  }

  lg_s2n->SetTextSize(0.04);
  lg_s2n->SetFillStyle(0);
  lg_s2n->SetLineWidth(0);
  lg_s2n->Draw();

  TLatex **txt_iceberg = new TLatex*[1];
  txt_iceberg[0] = new TLatex(0.0, 1.07, "IceBerg");
  txt_iceberg[0]->SetTextColor(1);
  txt_iceberg[0]->Draw();

  TLatex **txt_plane = new TLatex*[1];
  txt_plane[0] = new TLatex(120, 1.07, "#bf{Cosmics Data}");
  txt_plane[0]->SetTextColor(1);
  txt_plane[0]->Draw();

  c1_s2n->SaveAs("s2n.png");
  delete frame_s2n;
  

  /*
  // rmsfit
  TCanvas *c1_rmsfit = new TCanvas("c1_rmsfit","c1_rmsfit", 800, 600);
  
  auto lg_rmsfit = new TLegend(0.35, 0.6, 0.9, 0.9);
  float mean_rmsfit[3];
  double ymax_rmsfit=-99.;
  for (int p=0; p<3; p++) {
    h1_rmsfit[p]->SetStats(0);
    h1_rmsfit[p]->GetXaxis()->SetTitle("Noise #sigma [ADC counts]");
    h1_rmsfit[p]->GetYaxis()->SetTitle("");
    h1_rmsfit[p]->SetLineColor(5-color_options[p]);
    h1_rmsfit[p]->SetLineWidth(3);
    if (ymax_rmsfit<h1_rmsfit[p]->GetMaximum()) ymax_rmsfit = h1_rmsfit[p]->GetMaximum();
    //h1_rmsfit[p]->Draw(draw_options[p]);

    TLegendEntry* leg = lg_rmsfit->AddEntry(h1_rmsfit[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    TF1 *f1 = fit_rms(h1_rmsfit[p], 1);
    mean_rmsfit[p] = f1->GetParameter(0);
    printf("rms[%d]: %.10f\n", p, mean_rms[p]);
    printf("sigma of rms[%d]: %.10f\n", p, f1->GetParameter(1));
  }
  for (int p=0; p<3; p++) {
    h1_rmsfit[p]->GetXaxis()->SetRangeUser(0,10);
    h1_rmsfit[p]->GetYaxis()->SetRangeUser(0,ymax_rmsfit*1.05);
    h1_rmsfit[p]->Draw(draw_options[p]);
  }
  
  TString rmsfit_txt(Form("Noise #sigma (u/v/y): %.1f/%.1f/%.1f", mean_rmsfit[0], mean_rmsfit[1], mean_rmsfit[2]));
  lg_rmsfit->AddEntry(rms_txt, rms_txt.Data(), "");
  lg_rmsfit->SetTextColor(color_options[0]);
  lg_rmsfit->SetTextSize(0.05);
  lg_rmsfit->SetFillStyle(0);
  lg_rmsfit->SetLineWidth(0);

  lg_rmsfit->Draw();
  c1_rmsfit->SaveAs("noise_rmsfit.png");

  // signal amp
  TCanvas *c1_amp = new TCanvas("c1_amp","c1_amp", 800, 600);
  
  auto lg_amp = new TLegend(0.3, 0.6, 0.9, 0.9);
  float mean_amp[3];
  double ymax_amp=-99.;
  for (int p=0; p<3; p++) {
    h1_amp[p]->SetStats(0);
    h1_amp[p]->GetXaxis()->SetTitle("Signal pulse height [ADC counts]");
    h1_amp[p]->GetYaxis()->SetTitle("");
    h1_amp[p]->SetLineColor(5-color_options[p]);
    h1_amp[p]->SetLineWidth(3);
    if (ymax_amp<h1_amp[p]->GetMaximum()) ymax_amp = h1_amp[p]->GetMaximum();

    TLegendEntry* leg = lg_amp->AddEntry(h1_amp[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    TF1 *f1 = fit_rms(h1_amp[p], 1);
    mean_amp[p] = f1->GetParameter(0);
    printf("rms[%d]: %.10f\n", p, mean_amp[p]);
    printf("sigma of rms[%d]: %.10f\n", p, f1->GetParameter(1));
  }
  for (int p=0; p<3; p++) {
    h1_amp[p]->GetXaxis()->SetRangeUser(0,200);
    h1_amp[p]->GetYaxis()->SetRangeUser(0,ymax_amp*1.05);
    h1_amp[p]->Draw(draw_options[p]);
  }
  
  TString amp_txt(Form("Signal (u/v/y): %.1f/%.1f/%.1f", mean_amp[0], mean_amp[1], mean_amp[2]));
  lg_amp->AddEntry(amp_txt, amp_txt.Data(), "");
  lg_amp->SetTextColor(color_options[0]);
  lg_amp->SetTextSize(0.05);
  lg_amp->SetFillStyle(0);
  lg_amp->SetLineWidth(0);

  lg_amp->Draw();
  c1_amp->SaveAs("signal_amp.png");

  // raw waveform
  TCanvas *c1_raw_wf = new TCanvas("c1_raw_wf", "c1_raw_wf", 800, 600);

  auto lg_wf = new TLegend(0.15, 0.6, 0.9, 0.9);
  float mean_s2n[3];
  float average_amax[3];
  double ymax_s2n = -99.;
  for (int p=0; p<3; p++) {
    h1_s2n[p]->SetStats(0);
    h1_s2n[p]->GetXaxis()->SetTitle("Signal-to-Noise ratio [a.u.]");
    h1_s2n[p]->GetYaxis()->SetTitle("Counts [#]");
    if (ymax_s2n < h1_s2n[p]->GetMaximum()) ymax_s2n = h1_s2n[p]->GetMaximum();
    h1_s2n[p]->SetLineColor(5-color_options[p]);
    h1_s2n[p]->SetLineWidth(3);
    h1_s2n[p]->SetLineStyle(2);
    //h1_s2n[p]->Draw(draw_options[p]);

    TLegendEntry* leg = lg_wf->AddEntry(h1_s2n[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    average_amax[p] = h1_s2n[p]->GetMean();
    
    cout << "plane: " << p << endl;
    TF1 *f1 = fit_s2n(h1_s2n[p], 1);
    mean_s2n[p] = f1->GetParameter(0);
    printf("peak position[%d]: %.10f\n", p, mean_s2n[p]);
    printf("average [%d]: %.10f\n", p, average_amax[p]);
  }

  for (int p=0; p<3; p++) {
    h1_s2n[p]->GetYaxis()->SetRangeUser(0, 1.05*ymax_s2n);
    h1_s2n[p]->Draw(draw_options[p]);
  }

  TString s2n_txt(Form("Signal-to-Noise (u/v/y): %.1f/%.1f/%.1f", mean_s2n[0], mean_s2n[1], mean_s2n[2]));
  lg_wf->AddEntry(s2n_txt, s2n_txt.Data(), "");
  lg_wf->SetTextColor(color_options[0]);
  lg_wf->SetTextSize(0.05);
  lg_wf->SetFillStyle(0);
  lg_wf->SetLineWidth(0);

  lg_wf->Draw();
  c1_raw_wf->SaveAs("signal_to_noise_peak.png");

  // signal2noise mean
  TCanvas *c1_s2n = new TCanvas("c1_s2n", "c1_s2n", 800, 600);

  auto lg_s2n = new TLegend(0.15, 0.6, 0.9, 0.9);
  float average_s2n[3];
  double ymax_s2nmean = -99.;
  for (int p=0; p<3; p++) {
    h1_s2n[p]->SetStats(0);
    h1_s2n[p]->GetXaxis()->SetTitle("Signal-to-Noise ratio [a.u.]");
    h1_s2n[p]->GetYaxis()->SetTitle("Counts [#]");
    if (ymax_s2nmean < h1_s2n[p]->GetMaximum()) ymax_s2nmean = h1_s2n[p]->GetMaximum();
    h1_s2n[p]->SetLineColor(5-color_options[p]);
    h1_s2n[p]->SetLineWidth(3);
    h1_s2n[p]->SetLineStyle(2);
    //h1_s2n[p]->Draw(draw_options[p]);

    TLegendEntry* leg = lg_s2n->AddEntry(h1_s2n[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    average_amax[p] = h1_s2n[p]->GetMean();
    
    cout << "plane: " << p << endl;
    printf("average [%d]: %.10f\n", p, average_amax[p]);
  }

  for (int p=0; p<3; p++) {
    //h1_s2n[p]->GetYaxis()->SetRangeUser(0, 1.05*ymax_s2nmean);
    h1_s2n[p]->Draw(draw_options[p]);
  }

  TString s2nmean_txt(Form("Signal-to-Noise (u/v/y): %.1f/%.1f/%.1f", average_amax[0], average_amax[1], average_amax[2]));
  lg_s2n->AddEntry(s2nmean_txt, s2nmean_txt.Data(), "");
  lg_s2n->SetTextColor(color_options[0]);
  lg_s2n->SetTextSize(0.05);
  lg_s2n->SetFillStyle(0);
  lg_s2n->SetLineWidth(0);

  lg_s2n->Draw();
  c1_s2n->SaveAs("signal_to_noise_mean.png");
  // dqdx
  TCanvas *c1_dqdx = new TCanvas("c1_dqdx", "c1_dqdx", 800, 600);
  auto lg_dqdx = new TLegend(0.7,0.7,0.9,0.9);
  for (int p=0; p<3; p++) {
    
    h1_dqdx[p]->SetStats(0);
    h1_dqdx[p]->GetXaxis()->SetRangeUser(0, 600.);
    h1_dqdx[p]->SetLineColor(5-color_options[p]);
    h1_dqdx[p]->SetLineWidth(2);
    h1_dqdx[p]->Draw(draw_options[p]);
    lg_dqdx->AddEntry(h1_dqdx[p], plane_options[p], "l");
  }
  lg_dqdx->Draw();
  c1_dqdx->SaveAs("h1_dqdx.png");


  // dqdx vs time
  TCanvas *c1_dqdxtime[3];
  for (int p=0; p<3; p++) {
    
    c1_dqdxtime[p] = new TCanvas(TString::Format("c1_dqdxtime_plane_%d",p), TString::Format("c1_dqdxtime_plane_%d",p), 800, 600);
    h2_dqdxtime[p]->SetStats(0);
    h2_dqdxtime[p]->GetYaxis()->SetRangeUser(0, 600.);
    h2_dqdxtime[p]->Draw("colz");

    c1_dqdxtime[p]->SaveAs(TString::Format("h2_dqdxtime_plane_%d.png",p));
  }

  // angle dqdx
  TCanvas *c1_angle_dqdx[3];

  for (int p=0; p<3; p++) {
    
    c1_angle_dqdx[p] = new TCanvas(TString::Format("c1_angle_dqdx_plane_%d",p), TString::Format("c1_angle_dqdx_plane_%d",p), 800, 600);
    hp2d_angle_dqdx[p]->SetStats(0);
    hp2d_angle_dqdx[p]->Draw("colz");

    c1_angle_dqdx[p]->SaveAs(TString::Format("hp2d_angle_dqdx_plane_%d.png",p));
  }

*/

}

