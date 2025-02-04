
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

# TODO: add minimal version of file which only exclude needed branches -> 10x faster in cut framework, add as "path_minimal" in each dict

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
    "ttbar_hvq": {
        "path": [
            "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti1M",        # 1M
            "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti10M",       # 10M set1
            "/work/project/physics/psriling/FCC/DelphesPythia/FCChh/hvqMulti10M_set2",  # 10M set2
        ],
        "minimal": "/work/project/physics/psriling/FCC/extraction/ttbar_hvq",
        "cross-section": 31429.120490154033, # pb - POWHEG
    },
}
