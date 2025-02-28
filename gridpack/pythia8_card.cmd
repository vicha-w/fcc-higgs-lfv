! General parameters
! ==================
!
Main:numberOfEvents      = -1
!
! -------------------------------------------------------------------
!
HEPMCoutput:file         = $OUTFILE


JetMatching:qCut         = 1.9000000000e+01
JetMatching:doShowerKt   = off
JetMatching:nJetMax      = 4
SysCalc:fullCutVariation = off
!
Main:timesAllowErrors=10000
JetMatching:etaJetMax=5.0000000000e+00
JetMatching:merge=on
JetMatching:scheme=1
JetMatching:setMad=off
JetMatching:coneRadius=1.0000000000e+00
JetMatching:nQmatch=5
JetMatching:jetAlgorithm=2
JetMatching:slowJetPower=1
Check:epTolErr=1.0000000000e-02
PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118
SpaceShower:alphaSvalue=1.1800000000e-01
SpaceShower:alphaSorder=2
TimeShower:alphaSvalue=1.1800000000e-01
TimeShower:alphaSorder=2
Beams:setProductionScalesFromLHEF=off
Tune:preferLHAPDF=2
Tune:pp=14
Tune:ee=7
SLHA:minMassSM=1000.
ParticleDecays:limitTau0=on
ParticleDecays:tau0Max=10
ParticleDecays:allowPhotonRadiation=on
MultipartonInteractions:ecmPow=0.03344
MultipartonInteractions:bProfile=2
MultipartonInteractions:pT0Ref=1.41
MultipartonInteractions:coreRadius=0.7634
MultipartonInteractions:coreFraction=0.63
MultipartonInteractions:alphaSvalue=0.118
MultipartonInteractions:alphaSorder=2
ColourReconnection:range=5.176
SigmaTotal:zeroAXB=off
SigmaTotal:mode=0
SigmaTotal:sigmaEl=21.89
SigmaTotal:sigmaTot=100.309
SigmaProcess:alphaSvalue=0.118
SigmaProcess:alphaSorder=2
Init:showChangedSettings=on
Init:showChangedParticleData=off
Next:numberCount=100
Next:numberShowInfo=1
Next:numberShowProcess=1
Next:numberShowEvent=1
Random:setSeed=$SEED
!
! Additional technical parameters set by MG5_aMC.
!
! Tell Pythia8 that an LHEF input is used.
Beams:frameType=4
! 1.0 corresponds to HEPMC weight given in [mb]. We choose here the [pb] normalization.
HEPMCoutput:scaling=1.0000000000e+09
! Value of the merging scale below which one does not even write the HepMC event.
SysCalc:qWeed=1.0000000000e+01
!
!     ====================
!     Subrun definitions
!     ====================
!
LHEFInputs:nSubruns=1
Main:subrun=0
!
!  Definition of subrun 0
!
Beams:LHEF=$INFILE