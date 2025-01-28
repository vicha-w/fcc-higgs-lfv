
# Cross sections are in pb
MUTAU_SIGNAL = {
    "mutauM200": {
        "path": "../data/LFV_mutau/M200/output_powheg_M200.root",
        # "cross-section": 310.91,
        # Assume cross-section is 1 fb
        "cross-section": 1.0 / 1000,
    },
    "mutauM300": {
        "path": "../data/LFV_mutau/M300/output_powheg_M300.root",
        # "cross-section": 191.23,
        "cross-section": 1.0 / 1000,
    },
    "mutauM450": {
        "path": "../data/LFV_mutau/M450/output_powheg_M450.root",
        # "cross-section": 99.49,
        "cross-section": 1.0 / 1000,
    },
    "mutauM600": {
        "path": "../data/LFV_mutau/M600/output_powheg_M600.root",
        # "cross-section": 18.56,
        "cross-section": 1.0 / 1000,
    },
    "mutauM750": {
        "path": "../data/LFV_mutau/M750/output_powheg_M750.root",
        # "cross-section": 4.44,
        "cross-section": 1.0 / 1000,
    },
    "mutauM900": {
        "path": "../data/LFV_mutau/M900/output_powheg_M900.root",
        # "cross-section": 1.423,
        "cross-section": 1.0 / 1000,
    },
}


ETAU_SIGNAL = {
    "etauM200": {
        "path": "../data/output_pythia_etauM200.root",
        # "cross-section": 310.91,
        "cross-section": 1.0 / 1000,
    },
    "etauM300": {
        "path": "../data/LFV_etau/M300/output_powheg_M300.root",
        # "cross-section": 191.23,
        "cross-section": 1.0 / 1000,
    },
    "etauM450": {
        "path": "../data/LFV_etau/M450/output_powheg_M450.root",
        # "cross-section": 99.49,
        "cross-section": 1.0 / 1000,
    },
    "etauM600": {
        "path": "../data/LFV_etau/M600/output_powheg_M600.root",
        # "cross-section": 18.56,
        "cross-section": 1.0 / 1000,
    },
    "etauM750": {
        "path": "../data/LFV_etau/M750/output_powheg_M750.root",
        # "cross-section": 4.44,
        "cross-section": 1.0 / 1000,
    },
    "etauM900": {
        "path": "../data/LFV_etau/M900/output_powheg_M900.root",
        # "cross-section": 1.423,
        "cross-section": 1.0 / 1000,
    },
}


mutau_common_path = '/work/project/physics/psriling/FCC/DelphesPythia/FCChh/purple_LFV_mutau'
etau_common_path = '/work/project/physics/psriling/FCC/DelphesPythia/FCChh/purple_LFV_etau'
PURPLE_MUTAU_E = {}
PURPLE_ETAU_E = {}
for mass in [200, 300, 450, 600, 750, 900]:
    PURPLE_MUTAU_E[f"purple_mutauM{mass}"] = {
        "path": mutau_common_path + f"/M{mass}/output_powheg_M{mass}.root",
        "cross-section": 1.0 / 1000,
    }
    
    PURPLE_ETAU_E[f"purple_etauM{mass}"] = {
        "path": etau_common_path + f"/M{mass}/output_powheg_M{mass}.root",
        "cross-section": 1.0 / 1000,
    }
    


BACKGROUND = {
    'DYToLLJets': {
        'path': '/work/project/physics/psriling/FCC/core/output_backgrounds/DYToLLJets_S3_PFalse_N1000000',
        'cross-section': 160807.0
    },
    "DY1Jets": {
        "path": "/work/project/physics/psriling/FCC/core/output_backgrounds/DY1Jets_S3_PFalse_N1000000",
        "cross-section": 34722.71,
    },
    "DY2Jets": {
        "path": "/work/project/physics/psriling/FCC/core/output_backgrounds/DY2Jets_S3_PFalse_N1000000",
        "cross-section": 30898.09,
    },
    "DY3Jets": {
        "path": "/work/project/physics/psriling/FCC/core/output_backgrounds/DY3Jets_S3_PFalse_N1000000",
        "cross-section": 28478.54,
    },
    "DY4Jets": {
        "path": "/work/project/physics/psriling/FCC/core/output_backgrounds/DY4Jets_S3_PFalse_N1000000",
        "cross-section": 28366.36,
    },
    "WZTo1L1Nu2Q01j": {
        "path": "../data/bg_events/output_WZTo1L1Nu2Q01j_5f_NLO_FXFX.root",
        "cross-section": 349.4,
        
    },
    "WZTo2L2Q01j": {
        "path": "../data/bg_events/output_WZTo2L2Q01j_5f_NLO_FXFX.root",
        "cross-section": 171.7,
    },
    "ST_tch_top": {
        "path": "../data/bg_events/ST_tch/ST_tch_t.root",
        "cross-section": 1758.35
    },
    "ST_tch_antitop": {
        "path": "../data/bg_events/ST_tch/ST_tch_tb.root",
        "cross-section": 2316.52
    },
    "WAToLNuA0123j": {
        "path": "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/WAToLNuA0123j",
        "cross-section": 7995.77,
    },
    "ttbar_hvq": {
        "path": [
            "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti1M",
            "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti10M",
        ],
        "cross-section": 31429.12,
    },
    # "ttbj": {
    #     "path": "../data/ttbj/ttbj.root",
    #     "cross-section": 31429.12,
    # },
    "WJetsToLNu": {
        "path": "/work/project/physics/psriling/FCC/core/output_new/WJetsToLNu_S5_PFalse_N1000000/DelphesPythia8.root",
        "cross-section": 1218177.0,
        
    },
    "W1JetsToLNu": {
        "path": "/work/project/physics/psriling/FCC/core/output_new/W1JetsToLNu_S5_PFalse_N1000000/DelphesPythia8.root",
        "cross-section": 267200.0,
    },
    "W2JetsToLNu": {
        "path": "/work/project/physics/psriling/FCC/core/output_new/W2JetsToLNu_S5_PFalse_N1000000/DelphesPythia8.root",
        "cross-section": 232600.0,
    },
    "W3JetsToLNu": {
        "path": "/work/project/physics/psriling/FCC/core/output_new/W3JetsToLNu_S5_PFalse_N1000000/DelphesPythia8.root",
        "cross-section": 205600.0,
    },
    "W4JetsToLNu": {
        "path": "/work/project/physics/psriling/FCC/core/output_new/W4JetsToLNu_S5_PFalse_N1000000/DelphesPythia8.root",
        "cross-section": 197900.0,
    },
}
