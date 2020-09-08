// this script uses the output of ana.C as input and  plots noise rms, signal size, signal-to-noise ratio, and track angle, length, start/end position.


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












void plot() {

  // plot options
  const char * plane_options[3] = {"u","v","y"};
  const char * draw_options[10] = {"same", "SAME", "SAME", "SAME", "SAME", "SAME", "SAME", "SAME", "SAME", "SAME"};
  const int color_options[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 20};


  // input
  TFile *input = TFile::Open("output_ana.root");

  TH1F *h1_trk_startx = (TH1F*)input->Get("h1_trk_startx");
  TH1F *h1_trk_starty = (TH1F*)input->Get("h1_trk_starty");
  TH1F *h1_trk_startz = (TH1F*)input->Get("h1_trk_startz");

  TH1F *h1_trk_endx = (TH1F*)input->Get("h1_trk_endx");
  TH1F *h1_trk_endy = (TH1F*)input->Get("h1_trk_endy");
  TH1F *h1_trk_endz = (TH1F*)input->Get("h1_trk_endz");

  TH1F *h1_trk_len = (TH1F*)input->Get("h1_trk_len");

  TH1F *h1_trk_thetaxz = (TH1F*)input->Get("h1_trk_thetaxz");
  TH1F *h1_trk_thetayz = (TH1F*)input->Get("h1_trk_thetayz");
  TH2F *h2_trk_angle = (TH2F*)input->Get("h2_trk_angle");

  TH1F *h1_trk_thetaxz_1 = (TH1F*)input->Get("h1_trk_thetaxz_1");
  TH1F *h1_trk_thetayz_1 = (TH1F*)input->Get("h1_trk_thetayz_1");
  TH2F *h2_trk_angle_1 = (TH2F*)input->Get("h2_trk_angle_1");

  TH1F *h1_trk_thetaxz_2 = (TH1F*)input->Get("h1_trk_thetaxz_2");
  TH1F *h1_trk_thetayz_2 = (TH1F*)input->Get("h1_trk_thetayz_2");
  TH2F *h2_trk_angle_2 = (TH2F*)input->Get("h2_trk_angle_2");

  TH1F* h1_s2n[3];
  TH1F* h1_rms[3];
  TH1F* h1_rmsfit[3];
  TH1F* h1_amp[3];
  TH1F* h1_t_maxpulseheight[3];

  TProfile2D *hp2d_angle[3];
  TProfile2D *hp2d_angle_1[3];
  TProfile2D *hp2d_angle_2[3];

  //TH1F* h1_dqdx[3];
  //TH2F* h2_dqdxtime[3];
  //TProfile2D *hp2d_angle_dqdx[3];


  for (int p=0; p<3; p++) {
    h1_s2n[p] = (TH1F*)input->Get(TString::Format("h1_s2n_plane_%d", p));
    h1_s2n[p]->Rebin(4);
    h1_rms[p] = (TH1F*)input->Get(TString::Format("h1_rms_plane_%d", p)); 
    h1_rmsfit[p] = (TH1F*)input->Get(TString::Format("h1_rmsfit_plane_%d", p)); 
    h1_amp[p] = (TH1F*)input->Get(TString::Format("h1_amp_plane_%d", p)); 
    h1_amp[p]->Rebin(2);
    h1_t_maxpulseheight[p] = (TH1F*)input->Get(TString::Format("h1_t_maxpulseheight_plane_%d", p));

    //h1_s2n_aftercut[p] = (TH1F*)input->Get(TString::Format("h1_s2n_aftercut_plane_%d", p));

    hp2d_angle[p] = (TProfile2D*)input->Get(TString::Format("hp2d_angle_plane_%d",p));
    hp2d_angle_1[p] = (TProfile2D*)input->Get(TString::Format("hp2d_angle_1_plane_%d",p));
    hp2d_angle_2[p] = (TProfile2D*)input->Get(TString::Format("hp2d_angle_2_plane_%d",p));
  
    //h1_dqdx[p] = (TH1F*)input->Get(TString::Format("h1_dqdx_plane_%d",p));
    //h2_dqdxtime[p] = (TH2F*)input->Get(TString::Format("h2_dqdxtime_plane_%d",p));
    
    //hp2d_angle_dqdx[p] = (TProfile2D*)input->Get(TString::Format("hp2d_angle_dqdx_plane_%d",p));
  }

  // track start and end
  TCanvas *c1_start_end = new TCanvas("c1_start_end", "c1_start_end", 1200, 800);
  c1_start_end->Divide(3,2);
  
  c1_start_end->cd(1);
  h1_trk_startx->SetTitle("Start");
  h1_trk_startx->Draw();
  c1_start_end->cd(2);
  h1_trk_starty->SetTitle("Start");
  h1_trk_starty->Draw();
  c1_start_end->cd(3);
  h1_trk_startz->SetTitle("Start");
  h1_trk_startz->Draw();

  c1_start_end->cd(4);
  h1_trk_endx->SetTitle("End");
  h1_trk_endx->Draw();
  c1_start_end->cd(5);
  h1_trk_endy->SetTitle("End");
  h1_trk_endy->Draw();
  c1_start_end->cd(6);
  h1_trk_endz->SetTitle("End");
  h1_trk_endz->Draw();
  
  c1_start_end->SaveAs("track_start_end.png");

  // track length
  TCanvas *c1_len = new TCanvas("c1_len", "c1_len", 800, 600);
  
  gPad->SetLogy();

  h1_trk_len->Draw();

  c1_len->SaveAs("track_length.png");

  // track anagle
  TCanvas *c1_theta = new TCanvas("c1_theta", "c1_theta", 1200, 600);
  c1_theta->Divide(2,1);

  c1_theta->cd(1);
  h1_trk_thetaxz->SetStats(0);
  h1_trk_thetaxz->Draw();

  c1_theta->cd(2);
  h1_trk_thetayz->SetStats(0);
  h1_trk_thetayz->Draw();

  c1_theta->SaveAs("track_angle_1d.png");
  // h2
  TCanvas *c1_theta_2d = new TCanvas("c1_theta_2d", "c1_theta_2d", 800, 600);

  h2_trk_angle->SetStats(0);
  h2_trk_angle->Draw("colz");
  
  c1_theta_2d->SaveAs("track_angle_2d.png");

  // track anagle 1 (x<0)
  TCanvas *c1_theta_1 = new TCanvas("c1_theta_1", "c1_theta_1", 1200, 600);
  c1_theta_1->Divide(2,1);

  c1_theta_1->cd(1);
  h1_trk_thetaxz_1->SetStats(0);
  h1_trk_thetaxz_1->Draw();

  c1_theta_1->cd(2);
  h1_trk_thetaxz_2->SetStats(0);
  h1_trk_thetaxz_2->Draw();

  c1_theta_1->SaveAs("track_angle_1d_1.png");
  // h2
  TCanvas *c1_theta_2d_1 = new TCanvas("c1_theta_2d_1", "c1_theta_2d_1", 800, 600);

  h2_trk_angle_1->SetStats(0);
  h2_trk_angle_1->Draw("colz");
  
  c1_theta_2d_1->SaveAs("track_angle_2d_1.png");
  
  // track anagle 2 (x>0)
  TCanvas *c1_theta_2 = new TCanvas("c1_theta_2", "c1_theta_2", 1200, 600);
  c1_theta_2->Divide(2,1);

  c1_theta_2->cd(1);
  h1_trk_thetayz_1->SetStats(0);
  h1_trk_thetayz_1->Draw();

  c1_theta_2->cd(2);
  h1_trk_thetayz_2->SetStats(0);
  h1_trk_thetayz_2->Draw();

  c1_theta_2->SaveAs("track_angle_1d_2.png");
  // h2
  TCanvas *c1_theta_2d_2 = new TCanvas("c1_theta_2d_2", "c1_theta_2d_2", 800, 600);

  h2_trk_angle_2->SetStats(0);
  h2_trk_angle_2->Draw("colz");
  
  c1_theta_2d_2->SaveAs("track_angle_2d_2.png");

  //hp2d
  TCanvas *c1_theta_hp2d[3];

  for (int p=0; p<3; p++) {
    c1_theta_hp2d[p] = new TCanvas(TString::Format("c1_theta_hp2d_plane_%d",p), TString::Format("c1_theta_hp2d_plane_%d",p), 1000, 400);
    c1_theta_hp2d[p]->Divide(2,1);
    
    c1_theta_hp2d[p]->cd(1);
    hp2d_angle_1[p]->SetStats(0);
    hp2d_angle_1[p]->GetZaxis()->SetRangeUser(0, 100);
    hp2d_angle_1[p]->Draw("colz");
    
    c1_theta_hp2d[p]->cd(2);
    hp2d_angle_2[p]->SetStats(0);
    hp2d_angle_2[p]->GetZaxis()->SetRangeUser(0, 100);
    hp2d_angle_2[p]->Draw("colz");

    c1_theta_hp2d[p]->SaveAs(TString::Format("track_angle_sn2d_plane_%d.png",p));
  }

  // rms
  TCanvas *c1_rms = new TCanvas("c1_rms","c1_rms", 800, 600);
  
  auto lg_rms = new TLegend(0.35, 0.6, 0.9, 0.9);
  float mean_rms[3];
  double ymax_rms=-99.;
  for (int p=0; p<3; p++) {
    h1_rms[p]->SetStats(0);
    h1_rms[p]->GetXaxis()->SetTitle("Noise #sigma [ADC counts]");
    h1_rms[p]->GetYaxis()->SetTitle("");
    h1_rms[p]->SetLineColor(5-color_options[p]);
    h1_rms[p]->SetLineWidth(3);
    if (ymax_rms<h1_rms[p]->GetMaximum()) ymax_rms = h1_rms[p]->GetMaximum();
    //h1_rms[p]->Draw(draw_options[p]);

    TLegendEntry* leg = lg_rms->AddEntry(h1_rms[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    TF1 *f1 = fit_rms(h1_rms[p], 1);
    mean_rms[p] = f1->GetParameter(0);
    printf("rms[%d]: %.10f\n", p, mean_rms[p]);
    printf("sigma of rms[%d]: %.10f\n", p, f1->GetParameter(1));
  }
  for (int p=0; p<3; p++) {
    h1_rms[p]->GetXaxis()->SetRangeUser(0,10);
    h1_rms[p]->GetYaxis()->SetRangeUser(0,ymax_rms*1.05);
    h1_rms[p]->Draw(draw_options[p]);
  }
  
  TString rms_txt(Form("Noise #sigma (u/v/y): %.1f/%.1f/%.1f", mean_rms[0], mean_rms[1], mean_rms[2]));
  lg_rms->AddEntry(rms_txt, rms_txt.Data(), "");
  lg_rms->SetTextColor(color_options[0]);
  lg_rms->SetTextSize(0.05);
  lg_rms->SetFillStyle(0);
  lg_rms->SetLineWidth(0);

  lg_rms->Draw();
  c1_rms->SaveAs("noise_rms.png");

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
 /* 
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

