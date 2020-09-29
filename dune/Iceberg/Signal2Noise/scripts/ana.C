///
// To use this, run: 1) root; 2) .L ana.C 3) ana t 4) t.Loop();
//
#define ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TProfile2D.h>

using namespace std;

void ana::Loop()
{
  if (fChain == 0) return;

  // track select
  double trk_x_cut_1 = 1.;
  double trk_x_cut_2 = 40.;
  double trk_z_cut_1 = 3.;
  double trk_z_cut_2 = 114.;

  // track length
  double lbin_trk_len = 0.;
  double hbin_trk_len = 160; // [cm]
  int nbin_trk_len = 160;
  double trk_len_cut = 30.;

  // track angle
  double lbin_trk_thetaxz = -180.;
  double hbin_trk_thetaxz = 180.;
  int nbin_trk_thetaxz = 72.;

  double lbin_trk_thetayz = -180.;
  double hbin_trk_thetayz = 180.;
  int nbin_trk_thetayz = 72.;

  double radTodeg = 180.0/3.1415927;
  
  double trk_thetaxz_cut[3][4] = {
      {60., 120., -120., -60.},
      {60., 120., -120., -60.},
      {60., 120., -120., -60.}
    };
  double trk_thetayz_cut[3][4] = {
      {-84.29, -24.29, 95.71, 155.71}, // 35.71 deg to y+
      {24.29, 84.29, -155.71, -95.71}, // -35.71 deg to y+
      {-120.0, -60.0, 60.0, 120.}
    };


  // hit
  int max_hits = 1000;
  double threshold_negative = 8.;

  // signal-to-noise
  int nbin_s2n = 320;
  double lbin_s2n = 0.;
  double hbin_s2n = 160.;
  
  // rms
  int nbin_rms = 8000;
  double lbin_rms = 0.;
  double hbin_rms = 400.;

  // time for max pulse height
  int nbin_t_maxph = 2000;
  double lbin_t_maxph = 0.;
  double hbin_t_maxph = 2000.;

  // dqdx
  int nbin_dqdx = 200;
  double lbin_dqdx = 0.;
  double hbin_dqdx = 2000.;

  int nbin_time = 250;
  double lbin_time = 0.;
  double hbin_time = 250.;

  // output
  TFile output("output_ana.root","recreate");
  
  // histogram
  TH1F *h1_trk_startx = new TH1F("h1_trk_startx", "; x [cm]; Counts [#]", 120, -60., 60.);
  TH1F *h1_trk_starty = new TH1F("h1_trk_starty", "; y [cm]; Counts [#]", 200, -10, 190.);
  TH1F *h1_trk_startz = new TH1F("h1_trk_startz", "; z [cm]; Counts [#]", 150, -10, 140);

  TH1F *h1_trk_startx_selected = new TH1F("h1_trk_startx_selected", "; x [cm]; Counts [#]", 120, -60., 60.);
  TH1F *h1_trk_starty_selected = new TH1F("h1_trk_starty_selected", "; y [cm]; Counts [#]", 200, -10, 190);
  TH1F *h1_trk_startz_selected = new TH1F("h1_trk_startz_selected", "; z [cm]; Counts [#]", 150, -10, 140);

  TH1F *h1_trk_endx = new TH1F("h1_trk_endx", "; x [cm]; Counts [#]", 120, -60., 60.);
  TH1F *h1_trk_endy = new TH1F("h1_trk_endy", "; y [cm]; Counts [#]", 200, -10, 190);
  TH1F *h1_trk_endz = new TH1F("h1_trk_endz", "; z [cm]; Counts [#]", 150, -10, 140);

  TH1F *h1_trk_endx_selected = new TH1F("h1_trk_endx_selected", "; x [cm]; Counts [#]", 120, -60., 60.);
  TH1F *h1_trk_endy_selected = new TH1F("h1_trk_endy_selected", "; y [cm]; Counts [#]", 200, -10, 190);
  TH1F *h1_trk_endz_selected = new TH1F("h1_trk_endz_selected", "; z [cm]; Counts [#]", 150, -10, 140);

  TH1F *h1_trk_len = new TH1F("h1_trk_len", "; Track length [cm]; Counts [#]", nbin_trk_len, lbin_trk_len, hbin_trk_len);
  //TH1F *h1_trk_len_with_cut = new TH1F("h1_trk_len_with_cut", "; Track length [cm]; Counts [#]", nbin_trk_len, lbin_trk_len, hbin_trk_len);
  
  TH1F *h1_trk_thetaxz = new TH1F("h1_trk_thetaxz", "; #theta_{xz} [deg]; Counts [#]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz);
  TH1F *h1_trk_thetayz = new TH1F("h1_trk_thetayz", "; #theta_{yz} [deg]; Counts [#]", nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);

  TH1F *h1_trk_thetaxz_1 = new TH1F("h1_trk_thetaxz_1", "(x<0); #theta_{xz} [deg]; Counts [#]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz);
  TH1F *h1_trk_thetayz_1 = new TH1F("h1_trk_thetayz_1", "(x<0); #theta_{yz} [deg]; Counts [#]", nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);

  TH1F *h1_trk_thetaxz_2 = new TH1F("h1_trk_thetaxz_2", "(x>0); #theta_{xz} [deg]; Counts [#]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz);
  TH1F *h1_trk_thetayz_2 = new TH1F("h1_trk_thetayz_2", "(x>0); #theta_{yz} [deg]; Counts [#]", nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);

  TH1F *h1_trk_thetaxz_with_cut = new TH1F("h1_trk_thetaxz_with_cut", "; #theta_{xz} [deg]; Counts [#]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz);
  TH1F *h1_trk_thetayz_with_cut = new TH1F("h1_trk_thetayz_with_cut", "; #theta_{yz} [deg]; Counts [#]", nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);

  //TH1F *h1_cosgma = new TH1F("h1_cosgma", " #cos#gamma;", 20, -1, 1);

  TH2F *h2_trk_angle = new TH2F("h2_trk_angle", "; #theta_{xz}[deg]; #theta_{yz}[deg]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);
  TH2F *h2_trk_angle_1 = new TH2F("h2_trk_angle_1", "(x<0); #theta_{xz}[deg]; #theta_{yz}[deg]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);
  TH2F *h2_trk_angle_2 = new TH2F("h2_trk_angle_2", "(x>0); #theta_{xz}[deg]; #theta_{yz}[deg]", nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);
  
  TH1F* h1_cosgma[3];
  TH1F* h1_s2n[3];
  TH1F* h1_s2n_aftercut[3];
  TH1F* h1_rms[3];
  TH1F* h1_rmsfit[3];
  TH1F* h1_amp[3];
  TH1F* h1_t_maxpulseheight[3];

  TProfile2D *hp2d_angle[3];
  TProfile2D *hp2d_angle_1[3];
  TProfile2D *hp2d_angle_2[3];

  TH1F* h1_dqdx[3];
  TH2F* h2_dqdxtime[3];
  TProfile2D *hp2d_angle_dqdx[3];

  for (int p=0; p<3; p++) {
    h1_cosgma[p] = new TH1F(TString::Format("h1_cosgma_plane_%d",p), " #cos#gamma;", 20, -1, 1);

    h1_s2n[p] = new TH1F(TString::Format("h1_s2n_plane_%d", p), "; S/N ratio; Counts [#]", nbin_s2n, lbin_s2n, hbin_s2n);
    h1_rms[p] = new TH1F(TString::Format("h1_rms_plane_%d", p), "; rms; Counts [#]", nbin_rms, lbin_rms, hbin_rms);
    h1_rmsfit[p] = new TH1F(TString::Format("h1_rmsfit_plane_%d", p), "; rms; Counts [#]", nbin_rms, lbin_rms, hbin_rms);
    h1_amp[p] = new TH1F(TString::Format("h1_amp_plane_%d", p), "; amp; Counts [#]", 700, 0, 700);
    h1_t_maxpulseheight[p] = new TH1F(TString::Format("h1_t_maxpulseheight_plane_%d", p), "; Time [tick]; Counts [#]", nbin_t_maxph, lbin_t_maxph, hbin_t_maxph);
  
    h1_s2n_aftercut[p] = new TH1F(TString::Format("h1_s2n_aftercut_plane_%d", p), "; S/N ratio; Counts [#]", nbin_s2n, lbin_s2n, hbin_s2n);

    hp2d_angle[p] = new TProfile2D(TString::Format("hp2d_angle_plane_%d",p), TString::Format("Plane %d; #theta_{xz} [deg]; #theta_{yz} [deg]", p), nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);
    hp2d_angle_1[p] = new TProfile2D(TString::Format("hp2d_angle_1_plane_%d",p), TString::Format("Plane %d (x<0); #theta_{xz} [deg]; #theta_{yz} [deg]", p), nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);
    hp2d_angle_2[p] = new TProfile2D(TString::Format("hp2d_angle_2_plane_%d",p), TString::Format("Plane %d (x>0); #theta_{xz} [deg]; #theta_{yz} [deg]", p), nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);

    h1_dqdx[p] = new TH1F(TString::Format("h1_dqdx_plane_%d",p), "; dQ/dx [ADC/cm]; Counts [#]", nbin_dqdx, lbin_dqdx, hbin_dqdx);
    h2_dqdxtime[p] = new TH2F(TString::Format("h2_dqdxtime_plane_%d",p), "; Drift time [#mus]; dQ/dx [ADC/cm]", nbin_time, lbin_time, hbin_time, nbin_dqdx, lbin_dqdx, hbin_dqdx);

    hp2d_angle_dqdx[p] = new TProfile2D(TString::Format("hp2d_angle_dqdx_plane_%d",p), TString::Format("Plane %d; #theta_{xz} [deg]; #theta_{yz} [deg]", p), nbin_trk_thetaxz, lbin_trk_thetaxz, hbin_trk_thetaxz, nbin_trk_thetayz, lbin_trk_thetayz, hbin_trk_thetayz);
  }




  // entry
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries: " << nentries << endl;


  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    int n_trk_selected = 0;
    int n_trk_tpc_0 = 0;
    int n_trk_tpc_1 = 0;
    std::vector<int> trkid_selected;
    std::vector<int> trkid_tpc_0;
    std::vector<int> trkid_tpc_1;
    for (int i=0; i<ntrks; i++) {
      // track start
      h1_trk_startx->Fill(trkstart[i][0]);
      h1_trk_starty->Fill(trkstart[i][1]);
      h1_trk_startz->Fill(trkstart[i][2]);
      
      // track  end
      h1_trk_endx->Fill(trkend[i][0]);
      h1_trk_endy->Fill(trkend[i][1]);
      h1_trk_endz->Fill(trkend[i][2]);

      // track length
      h1_trk_len->Fill(trklen[i]);
      
      // track angle
      double thetaxz_deg = trkthetaxz[i]*radTodeg;
      double thetayz_deg = trkthetayz[i]*radTodeg;
      h1_trk_thetaxz->Fill(thetaxz_deg);
      h1_trk_thetayz->Fill(thetayz_deg);
      h2_trk_angle->Fill(trkthetaxz[i]*radTodeg, trkthetayz[i]*radTodeg);
      
      // x < 0
      if (trkstart[i][0] < 0 && trkend[i][0] < 0) {
        h1_trk_thetaxz_1->Fill(thetaxz_deg);
        h1_trk_thetayz_1->Fill(thetayz_deg);
        h2_trk_angle_1->Fill(trkthetaxz[i]*radTodeg, trkthetayz[i]*radTodeg);
      }
      // x>0
      if (trkstart[i][0] > 0 && trkend[i][0] > 0) {
        h1_trk_thetaxz_2->Fill(thetaxz_deg);
        h1_trk_thetayz_2->Fill(thetayz_deg);
        h2_trk_angle_2->Fill(trkthetaxz[i]*radTodeg, trkthetayz[i]*radTodeg);
      }
      
      // 2d profile anglular distribution
      for (int p=0; p<3; p++) {
        for (int j=0; j<max_hits; j++) {
          if (p != hit_plane[i][j]) continue;
          if (amp[i][j] == -9999.0) continue; // remove empty hits
          //if (amp[i][j] < threshold_negative || cosgma[i][j]*amp[i][j]<threshold_negative) continue; // remove negative hits, i.e., hits from anode gap
          double s2n = cosgma[i][j]*amp[i][j] / noiserms[i][j];
          hp2d_angle[p]->Fill(thetaxz_deg, thetayz_deg, s2n);
          
          if (trkstart[i][0] < 0 && trkend[i][0] < 0) {
            hp2d_angle_1[p]->Fill(thetaxz_deg, thetayz_deg, s2n);
          }

          if (trkstart[i][0] > 0 && trkend[i][0] > 0) {
            hp2d_angle_2[p]->Fill(thetaxz_deg, thetayz_deg, s2n);
          }
        }
      } // before track selection, check angular distribution of s2n
      

      // ---- select tracks ----
      // track length cut
      if (trklen[i] < trk_len_cut) continue;

      // track start/end cut
      if (! ( (trkstart[i][0]<trk_x_cut_2 && trkend[i][0]>trk_x_cut_1) ||
              (trkend[i][0]<trk_x_cut_2 && trkstart[i][0]>trk_x_cut_1) ||
              (trkstart[i][0]>-1.0*trk_x_cut_2 && trkend[i][0]<-1.0*trk_x_cut_1) ||
              (trkend[i][0]>-1.0*trk_x_cut_2 && trkstart[i][0]<-1.0*trk_x_cut_1) ) ) continue;
      if (! ( trkstart[i][2]<trk_z_cut_1 || trkstart[i][2] > trk_z_cut_2 ||
              trkend[i][2]<trk_z_cut_1 || trkend[i][2] > trk_z_cut_2 )) continue;
      
      if (trkstart[i][0]*trkend[i][0] < 0) continue;

      // track angle cut
      for (int p=0; p<3; p++) {
        
        // thetaxz
        if ( (thetaxz_deg > trk_thetaxz_cut[p][0] && thetaxz_deg < trk_thetaxz_cut[p][1]) || 
             (thetaxz_deg > trk_thetaxz_cut[p][2] && thetaxz_deg < trk_thetaxz_cut[p][3])) continue;

        // thetayz
        if (trkstart[i][0] < 0 && trkend[i][0] < 0) {
          if ( (thetayz_deg > trk_thetayz_cut[p][0] && thetayz_deg<trk_thetayz_cut[p][1]) || (thetayz_deg>trk_thetayz_cut[p][2] && thetayz_deg<trk_thetayz_cut[p][3]) )  continue;
        } // x<0

        if (trkstart[i][0] > 0 && trkend[i][0] > 0) {
          if ( (thetayz_deg < -1.0*trk_thetayz_cut[p][0] && thetayz_deg>-1.0*trk_thetayz_cut[p][1]) || (thetayz_deg<-1.0*trk_thetayz_cut[p][2] && thetayz_deg>-1.0*trk_thetayz_cut[p][3]) )  continue;
        } // x>0

        // hit
        for (int j=0; j<max_hits; j++) {
          if (hit_plane[i][j] == -1 || amp[i][j] == -9999.0 ) continue; // remove empty hits

          int hitplane = hit_plane[i][j];

          if (hitplane != p) continue;

          h1_cosgma[p]->Fill(cosgma[i][j]);

          h1_rms[p]->Fill(noiserms[i][j]);
          h1_rmsfit[p]->Fill(noisermsfit[i][j]);
          
          if (amp[i][j] == -9999.0 || tamp[i][j] == -1 || cosgma[i][j] == -99.) continue;
          if (amp[i][j] < threshold_negative || cosgma[i][j]*amp[i][j]<threshold_negative) continue; // remove negative hits, i.e., hits from anode gap
          double s2n = cosgma[i][j]*amp[i][j] / noiserms[i][j];
          h1_s2n[p]->Fill(s2n);
          h1_amp[p]->Fill((amp[i][j])*cosgma[i][j]);
          h1_t_maxpulseheight[p]->Fill(tamp[i][j]);
        
        } // loop over hit j
      
      } // loop over p

    } // loop over ntrks i
  } // loop over nentries jentry
  






  output.Write();

  cout << ".... done! ...." << endl;
}
