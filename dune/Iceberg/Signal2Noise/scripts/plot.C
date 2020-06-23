


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
  TF1 *gg=new TF1("gg",fit_gaus,pre_mean-5.*pre_rms,pre_mean+5.*pre_rms,3);
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

  TH1F *h1_trk_startx_selected = (TH1F*)input->Get("h1_trk_startx_selected");
  TH1F *h1_trk_starty_selected = (TH1F*)input->Get("h1_trk_starty_selected");
  TH1F *h1_trk_startz_selected = (TH1F*)input->Get("h1_trk_startz_selected");

  TH1F *h1_trk_endx_selected = (TH1F*)input->Get("h1_trk_endx_selected");
  TH1F *h1_trk_endy_selected = (TH1F*)input->Get("h1_trk_endy_selected");
  TH1F *h1_trk_endz_selected = (TH1F*)input->Get("h1_trk_endz_selected");

  

  TH1F *h1_trk_len = (TH1F*)input->Get("h1_trk_len");

  TH1F *h1_trk_thetaxz = (TH1F*)input->Get("h1_trk_thetaxz");
  TH1F *h1_trk_thetayz = (TH1F*)input->Get("h1_trk_thetayz");
  TH2F *h2_trk_angle = (TH2F*)input->Get("h2_trk_angle");

  TH1F* h1_s2n[3];
  TH1F* h1_s2n_aftercut[3];
  TH1F* h1_rms[3];
  TH1F* h1_t_maxpulseheight[3];

  TProfile2D *hp2d_angle[3];

  TH1F* h1_dqdx[3];
  TH2F* h2_dqdxtime[3];
  TProfile2D *hp2d_angle_dqdx[3];


  for (int p=0; p<3; p++) {
    h1_s2n[p] = (TH1F*)input->Get(TString::Format("h1_s2n_plane_%d", p));
    h1_rms[p] = (TH1F*)input->Get(TString::Format("h1_rms_plane_%d", p)); 
    h1_t_maxpulseheight[p] = (TH1F*)input->Get(TString::Format("h1_t_maxpulseheight_plane_%d", p));

    h1_s2n_aftercut[p] = (TH1F*)input->Get(TString::Format("h1_s2n_aftercut_plane_%d", p));

    hp2d_angle[p] = (TProfile2D*)input->Get(TString::Format("hp2d_angle_plane_%d",p));
  
    h1_dqdx[p] = (TH1F*)input->Get(TString::Format("h1_dqdx_plane_%d",p));
    h2_dqdxtime[p] = (TH2F*)input->Get(TString::Format("h2_dqdxtime_plane_%d",p));
    
     hp2d_angle_dqdx[p] = (TProfile2D*)input->Get(TString::Format("hp2d_angle_dqdx_plane_%d",p));
  }

  // track start and end
  TCanvas *c1_start_end = new TCanvas("c1_start_end", "c1_start_end", 1200, 800);
  c1_start_end->Divide(3,2);
  
  c1_start_end->cd(1);
  h1_trk_startx_selected->SetTitle("Start");
  h1_trk_startx_selected->Draw();
  c1_start_end->cd(2);
  h1_trk_starty_selected->SetTitle("Start");
  h1_trk_starty_selected->Draw();
  c1_start_end->cd(3);
  h1_trk_startz_selected->SetTitle("Start");
  h1_trk_startz_selected->Draw();

  c1_start_end->cd(4);
  h1_trk_endx_selected->SetTitle("End");
  h1_trk_endx_selected->Draw();
  c1_start_end->cd(5);
  h1_trk_endy_selected->SetTitle("End");
  h1_trk_endy_selected->Draw();
  c1_start_end->cd(6);
  h1_trk_endz_selected->SetTitle("End");
  h1_trk_endz_selected->Draw();
  
  c1_start_end->SaveAs("track_start_end.pdf");

  // track length
  TCanvas *c1_len = new TCanvas("c1_len", "c1_len", 800, 600);
  
  gPad->SetLogy();

  h1_trk_len->Draw();

  c1_len->SaveAs("track_length.pdf");

  // track anagle
  TCanvas *c1_theta = new TCanvas("c1_theta", "c1_theta", 1200, 600);
  
  c1_theta->Divide(2,1);

  c1_theta->cd(1);
  h1_trk_thetaxz->SetStats(0);
  h1_trk_thetaxz->Draw();

  c1_theta->cd(2);
  h1_trk_thetayz->SetStats(0);
  h1_trk_thetayz->Draw();

  c1_theta->SaveAs("track_angle_1d.pdf");

  
  // h2
  TCanvas *c1_theta_2d = new TCanvas("c1_theta_2d", "c1_theta_2d", 800, 600);

  h2_trk_angle->SetStats(0);
  h2_trk_angle->Draw("colz");
  
  c1_theta_2d->SaveAs("track_angle_2d.pdf");

  // hp2d
  TCanvas *c1_theta_hp2d = new TCanvas("c1_theta_hp2d", "c1_theta_hp2d", 1200, 400);
  
  c1_theta_hp2d->Divide(3,1);

  for (int p=0; p<3; p++) {
    c1_theta_hp2d->cd(p+1);
    hp2d_angle[p]->SetStats(0);
    hp2d_angle[p]->Draw("colz");
  }
  c1_theta_hp2d->SaveAs("track_angle_sn2d.pdf");



  // rms
  TCanvas *c1_rms = new TCanvas("c1_rms","c1_rms", 800, 600);
  
  auto lg_rms = new TLegend(0.35, 0.6, 0.9, 0.9);
  float mean_rms[3];
  for (int p=0; p<3; p++) {
    h1_rms[p]->SetStats(0);
    h1_rms[p]->GetXaxis()->SetTitle("Noise #sigma [ADC counts]");
    h1_rms[p]->GetYaxis()->SetTitle("");
    h1_rms[p]->SetLineColor(5-color_options[p]);
    h1_rms[p]->SetLineWidth(2);
    h1_rms[p]->Draw(draw_options[p]);

    TLegendEntry* leg = lg_rms->AddEntry(h1_rms[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    TF1 *f1 = fit_rms(h1_rms[p], 1);
    mean_rms[p] = f1->GetParameter(0);
    printf("rms[%d]: %.10f\n", p, mean_rms[p]);
    printf("sigma of rms[%d]: %.10f\n", p, f1->GetParameter(1));

  }
  
  TString rms_txt(Form("Noise #sigma (u/v/y): %.1f/%.1f/%.1f", mean_rms[0], mean_rms[1], mean_rms[2]));
  lg_rms->AddEntry(rms_txt, rms_txt.Data(), "");
  lg_rms->SetTextColor(color_options[0]);
  lg_rms->SetTextSize(0.05);
  lg_rms->SetFillStyle(0);
  lg_rms->SetLineWidth(0);

  lg_rms->Draw();
  c1_rms->SaveAs("noise_rms.pdf");


  // raw waveform
  TCanvas *c1_raw_wf = new TCanvas("c1_raw_wf", "c1_raw_wf", 800, 600);

  auto lg_wf = new TLegend(0.15, 0.6, 0.9, 0.9);
  float mean_s2n[3];
  float average_amax[3];
  for (int p=0; p<3; p++) {
    /*
    h1_s2n[p]->SetStats(0);
    h1_s2n[p]->GetXaxis()->SetTitle("Signal-to-Noise ratio [a.u.]");
    h1_s2n[p]->GetYaxis()->SetTitle("Counts [#]");
    h1_s2n[p]->GetYaxis()->SetRangeUser(0,40.);
    h1_s2n[p]->SetLineColor(5-color_options[p]);
    h1_s2n[p]->SetLineWidth(2);
    h1_s2n[p]->SetLineStyle(2);
    h1_s2n[p]->Draw(draw_options[p]);
    */

    //TLegendEntry* leg = lg_wf->AddEntry(h1_s2n[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    //leg->SetTextColor(5-color_options[p]);

    h1_s2n_aftercut[p]->SetStats(0);
    h1_s2n_aftercut[p]->GetXaxis()->SetTitle("Signal-to-Noise ratio [a.u.]");
    h1_s2n_aftercut[p]->GetYaxis()->SetTitle("Counts [#]");
    h1_s2n_aftercut[p]->GetYaxis()->SetRangeUser(0, 40.);
    h1_s2n_aftercut[p]->SetLineColor(5-color_options[p]);
    h1_s2n_aftercut[p]->SetLineWidth(2);
    h1_s2n_aftercut[p]->Draw(draw_options[p+3]);

    TLegendEntry* leg = lg_wf->AddEntry(h1_s2n_aftercut[p], TString::Format("plane %s (raw)", plane_options[p]), "");
    leg->SetTextColor(5-color_options[p]);

    average_amax[p] = h1_s2n_aftercut[p]->GetMean();
    
    cout << "plane: " << p << endl;
    TF1 *f1 = fit_s2n(h1_s2n_aftercut[p], 1);
    mean_s2n[p] = f1->GetParameter(0);
    printf("peak position[%d]: %.10f\n", p, mean_s2n[p]);
    printf("average [%d]: %.10f\n", p, average_amax[p]);

  }
  
  TString s2n_txt(Form("Signal-to-Noise (u/v/y): %.1f/%.1f/%.1f", mean_s2n[0], mean_s2n[1], mean_s2n[2]));
  lg_wf->AddEntry(s2n_txt, s2n_txt.Data(), "");
  lg_wf->SetTextColor(color_options[0]);
  lg_wf->SetTextSize(0.05);
  lg_wf->SetFillStyle(0);
  lg_wf->SetLineWidth(0);

  lg_wf->Draw();
  c1_raw_wf->SaveAs("rawwf_signal_to_noise.pdf");

  
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
  c1_dqdx->SaveAs("h1_dqdx.pdf");


  // dqdx vs time
  TCanvas *c1_dqdxtime[3];
  for (int p=0; p<3; p++) {
    
    c1_dqdxtime[p] = new TCanvas(TString::Format("c1_dqdxtime_plane_%d",p), TString::Format("c1_dqdxtime_plane_%d",p), 800, 600);
    h2_dqdxtime[p]->SetStats(0);
    h2_dqdxtime[p]->GetYaxis()->SetRangeUser(0, 600.);
    h2_dqdxtime[p]->Draw("colz");

    c1_dqdxtime[p]->SaveAs(TString::Format("h2_dqdxtime_plane_%d.pdf",p));
  }

  // angle dqdx
  TCanvas *c1_angle_dqdx[3];

  for (int p=0; p<3; p++) {
    
    c1_angle_dqdx[p] = new TCanvas(TString::Format("c1_angle_dqdx_plane_%d",p), TString::Format("c1_angle_dqdx_plane_%d",p), 800, 600);
    hp2d_angle_dqdx[p]->SetStats(0);
    hp2d_angle_dqdx[p]->Draw("colz");

    c1_angle_dqdx[p]->SaveAs(TString::Format("hp2d_angle_dqdx_plane_%d.pdf",p));
  }




}

