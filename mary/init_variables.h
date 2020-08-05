#ifndef INIT_VARIABLES_H
#define INIT_VARIABLES_H

using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;

TString fttl_main = "$$(pname)_$$(sys)_$$(vshort)";
TString fttly_ntne = "N_{tracks} / N_{events}";
TString fttly_n1 = "N_{prob > 0.5}";
TString fttly_n2 = "N_{prob > 0.2, eff > 0.05, abs(sd) < 5}";

auto fvar_prob     = variable("prob" , "prob" , "$$(cut0)", titles(fttl_main, "probability"                 , fttly_ntne), binning(200,0,1.01));
//auto fvar_eff      = variable("eff"  , "eff"  , "$$(cut0)", titles(fttl_main, "efficiency"                  , fttly_ntne), binning(200,0,1.01));
auto fvar_eff      = variable("eff"  , "eff"  , "$$(cut0)", titles(fttl_main, "efficiency"                  , fttly_ntne), binning(100,0,0.15));
auto fvar_sd       = variable("sd"   , "sd"   , "$$(cut0)", titles(fttl_main, "distance_{PID mean} (#sigma)", fttly_ntne), binning(200,-5,5));
auto fvar_pt_cm    = variable("pt_cm", "pt_cm", "$$(cut0)", titles(fttl_main, "p_{T,CM} (MeV/c^{2})"        , fttly_ntne), binning(400,0,2000));
auto fvar_fy_cm    = variable("fy_cm", "fy_cm", "$$(cut0)", titles(fttl_main, "y_{Frag,CM}"                 , fttly_ntne), binning(100,-2.,2.));
auto fvar_p_cm     = variable("p_cm" , "p_cm" , "$$(cut0)", titles(fttl_main, "p_{CM} (MeV/c^{2})"          , fttly_ntne), binning(400,0,3000));
auto fvar_ke_cm    = variable("ke_cm", "ke_cm", "$$(cut0)", titles(fttl_main, "K.E._{CM} (MeV)"             , fttly_ntne), binning(40,0,400));
auto fvar_p_lab    = variable("p_lab", "p_lab", "$$(cut0)", titles(fttl_main, "p_{Lab} (MeV/c^{2})"         , fttly_ntne), binning(400,0,2000));
auto fvar_dedx     = variable("dedx" , "dedx" , "$$(cut0)", titles(fttl_main, "dE/dx"                       , fttly_ntne), binning(400,0,1500));
auto fvar_dpoca    = variable("dpoca", "dpoca", "$$(cut0)", titles(fttl_main, "distance_{POCA-PV} (mm)"     , fttly_ntne), binning(200,0,10));
auto fvar_nr       = variable("nr"   , "nc"   , "$$(cut0)", titles(fttl_main, "N_{row-cluster}"             , fttly_ntne), binning(150,1,150));
auto fvar_nl       = variable("nl"   , "nl"   , "$$(cut0)", titles(fttl_main, "N_{layer-cluster}"           , fttly_ntne), binning(150,1,150));
auto fvar_tta_cm   = variable("tta_cm" ,  "theta_cm*57.295780", "$$(cut0)", titles(fttl_main, "#theta_{cm}" , fttly_ntne), binning(200,   0,180));
auto fvar_phi_cm   = variable("phi_cm" ,    "phi_cm*57.295780", "$$(cut0)", titles(fttl_main, "#phi_{cm}"   , fttly_ntne), binning(200,-180,180));
auto fvar_tta_lab  = variable("tta_lab", "theta_lab*57.295780", "$$(cut0)", titles(fttl_main, "#theta_{lab}", fttly_ntne), binning(200,   0,180));
auto fvar_phi_lab  = variable("phi_lab",   "phi_lab*57.295780", "$$(cut0)", titles(fttl_main, "#phi_{lab}"  , fttly_ntne), binning(200,-180,180));

auto fvar_asd      = variable("asd"    , "abs(sd)"      , "$$(cut0)", titles(fttl_main, "|distance_{PID mean}} (#sigma)", fttly_ntne), binning(200,0,5));
auto fvar_poz_lab  = variable("poz_lab", "p_lab/$$(z)"  , "$$(cut0)" ,titles(fttl_main, "p_{Lab}/Z (MeV/c^{2})", fttly_ntne), binning(400,0,2000));
auto fvar_keoa_cm  = variable("keoa_cm", "ke_cm/$$(a)"  , "$$(cut0)", titles(fttl_main, "KE_{CM}/A (MeV)"      , fttly_ntne), binning(40,0,400));
auto fvar_ptoa_cm  = variable("ptoa_cm", "pt_cm/$$(a)"  , "$$(cut0)", titles(fttl_main, "p_{T}/A (MeV/c)"      , fttly_ntne), binning(100,0,1000));
auto fvar_ny_cm    = variable("ny_cm"  , "fy_cm/$$(ynn)", "$$(cut0)", titles(fttl_main, "y_{CM}/y_{NN}"        , fttly_ntne), binning(100,-2.,2.));

auto fvar_n      = variable("n"        , "np+nd+nt+nhe3+nhe4", "$$(cutn)", titles("", "all_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_np     = variable("np"       , "np"                , "$$(cutn)", titles("",   "p_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nd     = variable("nd"       , "nd"                , "$$(cutn)", titles("",   "d_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nt     = variable("nt"       , "nt"                , "$$(cutn)", titles("",   "t_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nhe3   = variable("nhe3"     , "nhe3"              , "$$(cutn)", titles("", "he3_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nhe4   = variable("nhe4"     , "nhe4"              , "$$(cutn)", titles("", "he4_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_ng     = variable("n_good"   , "n_good"            , "$$(cutn)", titles("", "all_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngp    = variable("np_good"  , "np_good"           , "$$(cutn)", titles("",   "p_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngd    = variable("nd_good"  , "nd_good"           , "$$(cutn)", titles("",   "d_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngt    = variable("nt_good"  , "nt_good"           , "$$(cutn)", titles("",   "t_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_nghe3  = variable("nhe3_good", "nhe3_good"         , "$$(cutn)", titles("", "he3_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_nghe4  = variable("nhe4_good", "nhe4_good"         , "$$(cutn)", titles("", "he4_$$(sys)", fttly_n2), binning(100,0,100));

#endif
