
NPsPBPK.code <-'
$PROB
## Nanoparticles (NPs) PBPK model for female SD rats
- Author    : Wei-Chun Chou
- Advisor   : Zhoumeng Lin
- Date      : August, 2021
- Strucutre : GI tract, Plasma, Liver, Kidney, Spleen, rest of body

$PARAM @annotated
// #+ Body weight and fraction of blood flow to tissues
BW     : 0.273     : Kg,       body weight                       ; Measured data from Dr. Kreyling
QCC    : 15.0      : L/h/kg^0.75, Cardiac output; Value obtained ; Brown et al., 1997
QLuC   : 1         : % of QCC, Fraction of blood flow to lung    ; Brown et al., 1997
QGIC   : 0.101     : % of QCC, Fraction of blood flow to GI tract; Davis and Morris, 1993
QLC    : 0.021     : % of QCC, Fraction of blood flow to Liver   ; Brown et al., 1997
QKC    : 0.141     : % of QCC, Fraction of blood flow to kidney  ; Brown et al., 1997
QSC    : 0.0085    : % of QCC, Fraction of blood flow to spleen  ; Davis and Morris, 1993

// #+ Fraction of volume of tissues out of body weight
VLuC   : 0.0074    : % of BW,  Fraction of volume of lung        ; Measured data from Dr. Kreyling
VGIC   : 0.1014    : % of BW,  Fraction of volume of GI          ; Measured data from Dr. Kreyling
VLC    : 0.0346    : % of BW,  Fraction of volume of liver       ; Measured data from Dr. Kreyling
VKC    : 0.0091    : % of BW,  Fraction of volume of kidney      ; Measured data from Dr. Kreyling
VSC    : 0.0038    : % of BW,  Fraction of volume of spleen      ; Measured data from Dr. Kreyling
VBC    : 0.0717    : % of BW,  Fraction of volume of blood       ; Measured data from Dr. Kreyling

// #+ Fraction of blood volume in tissues
BVLu   : 0.36      : % of VLu, Fraction of blood volume in lung  ; Brown et al., 1997 
BVGI   : 0.04      : % of VGI, Fraction of blood volume in GI    ; Brown et al., 1997
BVL    : 0.21      : % of VL,  Fraction of blood volume in liver ; Brown et al., 1997
BVK    : 0.16      : % of VLK, Fraction of blood volume in kidney; Brown et al., 1997
BVS    : 0.22      : % of VLS, Fraction of blood volume in spleen; Brown et al., 1997
BVR    : 0.04      : % of VR,  Fraction of blood volume in rest of body; Brown et al., 1997

// #+ Partition coefficient
PLu    : 0.15      : Unitless, Partition coefficient of lung     ; Lin et al., 2016
PGI    : 0.15      : Unitless, Partition coefficient of GI       ; Lin et al., 2016
PL     : 0.08      : Unitless, Partition coefficient of Liver    ; Lin et al., 2016
PK     : 0.15      : Unitless, Partition coefficient of Kidney   ; Lin et al., 2016
PS     : 0.15      : Unitless, Partition coefficient of Spleen   ; Lin et al., 2016
PR     : 0.15      : Unitless, Partition coefficient of rest of body; Lin et al., 2016

// #+ Membrane-limited permeability
PALuC  : 0.001     : Unitless, Membrane-limited permeability coefficient of lung         ; Lin et al., 2016
PAGIC  : 0.001     : Unitless, Membrane-limited permeability coefficient of GI           ; Lin et al., 2016
PALC   : 0.001     : Unitless, Membrane-limited permeability coefficient of Liver        ; Lin et al., 2016
PAKC   : 0.001     : Unitless, Membrane-limited permeability coefficient of Kidney       ; Lin et al., 2016 
PASC   : 0.001     : Unitless, Membrane-limited permeability coefficient of Spleen       ; Lin et al., 2016
PARC   : 0.001     : Unitless, Membrane-limited permeability coefficient of Rest of body ; Lin et al., 2016

// #+ Endocytic parameters in liver
KLRESrelease       : 10       : 1/h,      Release rate constant of phagocytic cells 
KLRESmax           : 20       : 1/h,      Maximum uptake rate constant of phagocytic cells
KLRES50            : 24       : h,        Time reaching half maximum uptake rate
KLRESn             : 0.5      : Unitless, Hill coefficient, (unitless)
ALREScap           : 100      : ug/g,     tissue, Uptake capacity per tissue weight 

// #+ Endocytic parameters in spleen
KSRESrelease       : 1        : 1/h,      Release rate constant of phagocytic cells
KSRESmax           : 5        : 1/h,      Maximum uptake rate constant of phagocytic cells
KSRES50            : 24       : h,        Time reaching half maximum uptake rate
KSRESn             : 0.5      : Unitless, Hill coefficient, (unitless)
ASREScap           : 100      : ug/g,     tissue, Uptake capacity per tissue weight

// #+ Endocytic parameters in kidney
KKRESrelease       : 1        : 1/h,      Release rate constant of phagocytic cells
KKRESmax           : 5        : 1/h,      Maximum uptake rate constant of phagocytic cells
KKRES50            : 24       : h,        Time reaching half maximum uptake rate
KKRESn             : 0.5      : Unitless, Hill coefficient, (unitless)
AKREScap           : 100      : ug/g,     tissue, Uptake capacity per tissue weight

// #+ Endocytic parameters in pulmonary region of lung
KpulRESrelease     : 7        : 1/h,      Release rate constant of phagocytic cells
KpulRESmax         : 20       : 1/h,      Maximum uptake rate constant of phagocytic cells
KpulRES50          : 24       : h,        Time reaching half maximum uptake rate
KpulRESn           : 0.5      : Unitless, Hill coefficient, (unitless)
ApulREScap         : 100      : ug/g,     tissue, Uptake capacity per tissue weight

// #+ Endocytic parameters in interstitium of lung
KinRESrelease      : 7        : 1/h,      Release rate constant of phagocytic cells
KinRESmax          : 20       : 1/h,      Maximum uptake rate constant of phagocytic cells
KinRES50           : 24       : h,        Time reaching half maximum uptake rate
KinRESn            : 0.5      : Unitless, Hill coefficient, (unitless)
AinREScap          : 100      : ug/g,     tissue, Uptake capacity per tissue weight

// #+ Endocytic parameters in rest of body
KRRESrelease       : 0.2      : 1/h,      Release rate constant of phagocytic cells
KRRESmax           : 0.1      : 1/h,      Maximum uptake rate constant of phagocytic cells
KRRES50            : 24       : h,        Time reaching half maximum uptake rate
KRRESn             : 0.5      : Unitless, Hill coefficient, (unitless)
ARREScap           : 100      : ug/g,     tissue, Uptake capacity per tissue weight

// #+ Endocytic parameters in GI tract
KGIRESrelease      : 0.001    : 1/h,      Release rate constant of phagocytic cells
KGIRESmax          : 0.1      : 1/h,      Maximum uptake rate constant of phagocytic cells
AGIREScap          : 100      : ug/g,     tissue, Uptake capacity per tissue weight

// #+ Uptake and elimination parameters
KGIb               : 4e-5     : 1/h, absorption rate of GI tract
KbileC             : 0.05     : L/h/kg^0.75, Biliary clearance
KfecesC            : 1.32     : L/h/kg^0.75, Fecal clearance
KurineC            : 0.0012   : L/h/kg^0.75, Urinary clearance

// #+ Respiratory tract compartment--Membrane-limited model with phagocytic cells in pulmonary and interstitial regions 
// #+ Transport or clearance rate in lung (1/h)
Kuab               : 8.36e-5  : 1/h, Transport rate constant from upper airway to brain
Kpulrestra         : 8.65e-4  : 1/h, Transport rate constant from inactive pulmonary phagocytizing cells to trachealbronchial region
Kpulin             : 0.126    : 1/h, Transport rate constant from pulmonary region to interstitium of lung
Kinpul             : 1.12e-6  : 1/h, Transport rate constant from interstitium of lung to pulmonary region
KtragiC            : 5.52e-3  : 1/h, Transport rate constant from trachealbronchial region to the GI lumen
KuagiC             : 0.335    : 1/h, Transport rate constant from upper airways to the GI lumen

$MAIN
// #+ Blood flow to tissues (L/h)
double QC     = QCC*pow(BW,0.75);
double QGI    = QC*QGIC;
double QL     = QC*QLC;
double QK     = QC*QKC;
double QS     = QC*QSC;
double QR     = QC*(1 - QGIC - QLC - QKC - QSC);
double QBal   = QC - (QGI + QL + QK + QS + QR);

// #+ Tissue volumes (L)
double VLu    = BW*VLuC;
double VGI    = BW*VGIC;
double VL     = BW*VLC;
double VK     = BW*VKC;
double VS     = BW*VSC;
double VB     = BW*VBC;
double VR     = BW*(1 - VGIC - VLC - VKC - VSC - VBC);
double VBal   = BW - (VLu + VGI + VL + VK + VS + VR + VB);

// #+ Tissue volumes for different compartments (L)
double VLub   = VLu*BVLu;  
double VLut   = VLu-VLub;  
double VGIb   = VGI*BVGI;  
double VGIt   = VGI-VGIb;  
double VLb    = VL*BVL;    
double VLt    = VL-VLb;    
double VKb    = VK*BVK;    
double VKt    = VK-VKb;    
double VSb    = VS*BVS;    
double VSt    = VS-VSb;    
double VRb    = VR*BVR;    
double VRt    = VR-VRb;    

// #+ Permeability coefficient-surface area cross-product (L/h)
double PALu   = PALuC*QC;  
double PAGI   = PAGIC*QGI; 
double PAL    = PALC*QL; 
double PAK    = PAKC*QK; 
double PAS    = PASC*QS; 
double PAR    = PARC*QR; 

// #+ Endocytosis rate (1/h)
double KLRESUP     = KLRESmax*(1-(ALRES/(ALREScap*VL)));
double KSRESUP     = KSRESmax*(1-(ASRES/(ASREScap*VS)));
double KKRESUP     = KKRESmax*(1-(AKRES/(AKREScap*VK)));

double KpulRESUP   = KpulRESmax*(1-(ApulRES/(ApulREScap*VLu)));

double KinRESUP    = KinRESmax*(1-(AinRES/(AinREScap*VLu)));

double KGIRESUP    = KGIRESmax*(1-(AGIRES/(AGIREScap*VGI))); 

double KRRESUP     = KRRESmax*(1-(ARRES/(ARREScap*VR)));

// #+ Biliary excretion 
double Kbile       = KbileC*pow(BW, 0.75);
double Kfeces      = KfecesC*pow(BW, 0.75);
double Kurine      = KurineC*pow(BW, 0.75); 

// #+ Respiratory tract compartment
double Kuagi       = TIME<0.6 ? KuagiC + 0.8: KuagiC;
double Ktragi      = TIME<0.6 ? KtragiC + 0.8: KtragiC;

$INIT @annotated
AA                : 0  : mg, Amount of NPs in arterial blood compartment
AV                : 0  : mg, Amount of NPs in venous blood compartment
ARb               : 0  : mg, Amount of NPs in capillary blood of rest of tissues
ARt               : 0  : mg, Amount of NPs in rest of tissues compartment
ARRES             : 0  : mg, Amount of NPs in phagocytic cells of remaining tissues
Aua               : 0  : mg, Amount of NPs in upper airway region of lung compartment
Atra              : 0  : mg, Amount of NPs in trachealbronchial region of lung compartment
Apul              : 0  : mg, Amount of NPs in pulmonary region of lung compartment
ApulRES           : 0  : mg, Amount of NPs in phagocytic cells of pulmonary region
Ain               : 0  : mg, Amount of NPs in interstitium of rest of lung
AinRES            : 0  : mg, Amount of NPs in phagocytic cells of interstitium
ALub              : 0  : mg, Amount of NPs in capillary blood of lung
ASb               : 0  : mg, Amount of NPs in capillary blood of spleen
ASt               : 0  : mg, Amount of NPs in spleen compartment
ASRES             : 0  : mg, Amount of NPs in phagocytic cells of spleen
AGIb              : 0  : mg, Amount of NPs in capillary blood of GI
AGIt              : 0  : mg, Amount of NPs in GI tract
ALumen            : 0  : mg, Amount of NPs in GI tract lumen
AGIRES            : 0  : mg, Amount of NPs in phagocytic cells of GI
ALb               : 0  : mg, Amount of NPs in capillary blood of liver
ALt               : 0  : mg, Amount of NPs in liver compartment
ALRES             : 0  : mg, Amount of NPs in phagocytic cells of liver
AKb               : 0  : mg, Amount of NPs in capillary blood of Kidney
AKt               : 0  : mg, Amount of NPs in Kidney
AKRES             : 0  : mg, Amount of NPs in phagocytic cells of Kidney
Aurine            : 0  : mg, Amount of NPs in urinary excretion
Abile             : 0  : mg, Amount of NPs in biliary excretion
Afeces            : 0  : mg, Amount of NPs in feces excretion
AUCB              : 0  : mg/L*h, AUC in blood
AUCLu             : 0  : mg/L*h, AUC in lung
AUCS              : 0  : mg/L*h, AUC in spleen
AUCRt             : 0  : mg/L*h, AUC in rest of tissue
AUCGI             : 0  : mg/L*h, AUC in GI tract
AUCL              : 0  : mg/L*h, AUC in liver
AUCK              : 0  : mg/L*h, AUC in kidney

$ODE
// #+ Concentrations in the tissues (C) and in the venous plasma leaving each of the tissues (CV) (Unit: mg/L)
// #+ A:arterial blood; V: venous blood compartment; in:interstitium of rest of lung; Lu: lung
// #+ S: Spleen; GI: GI tract; R: rest of tissues

double CA        = AA/ (VB*0.2);
double CV        = AV/ (VB*0.8);
double Cin       = Ain/VLut;                      
double CVLu      = ALub/VLub;
double CVS       = ASb/VSb;
double CSt       = ASt/VSt;
double CVGI      = AGIb/VGIb;
double CGIt      = AGIt/VGIt;
double CVL       = ALb/VLb;
double CLt       = ALt/VLt;
double CVK       = AKb/VKb;
double CKt       = AKt/VKt;
double CVR       = ARb/VRb; 
double CRt       = ARt/VRt;

// #+ Equation for estimation of the rate of each compartment 
double RAA       = QC*CVLu - QC*CA;                                             
double RAV       = QL*CVL + QK*CVK + QR*CVR - QC*CV;                                
double RARb      = QR*(CA - CVR) - PAR*CVR + (PAR*CRt)/PR;
double RARt      = PAR*CVR - (PAR*CRt)/PR - (KRRESUP*ARt - KRRESrelease*ARRES) + Kuab*Aua;
double RARRES    = KRRESUP*ARt - KRRESrelease*ARRES;
double RAua      = -(Kuagi + Kuab)*Aua;
double Rtra      = Kpulrestra*ApulRES - Ktragi*Atra;
double RApul     = Kinpul*Ain - Kpulin*Apul - KpulRESUP*Apul + KpulRESrelease*ApulRES;
double RApulRES  = KpulRESUP*Apul - KpulRESrelease*ApulRES - Kpulrestra*ApulRES;
double RAin      = Kpulin*Apul - Kinpul*Ain + PALu*CVLu - (PALu*Cin)/PLu - KinRESUP*Ain + KinRESrelease*AinRES;
double RAinRES   = KinRESUP*Ain - KinRESrelease*AinRES;
double RALub     = QC*(CV-CVLu) - PALu*CVLu + (PALu*Cin)/PLu;
double RASb      = QS*(CA-CVS) - PAS*CVS + (PAS*CSt)/PS;
double RASt      = PAS*CVS - (PAS*CSt)/PS - KSRESUP*ASt + KSRESrelease*ASRES;
double RASRES    = KSRESUP*ASt-KSRESrelease*ASRES;
double RAGIb     = QGI*(CA-CVGI) - PAGI*CVGI + (PAGI*CGIt)/PGI + KGIb*ALumen;
double RAGIt     = PAGI*CVGI - (PAGI*CGIt)/PGI - (KGIRESUP*AGIt-KGIRESrelease*AGIRES);
double RALumen   = Kuagi*Aua + Ktragi*Atra + Kbile*CLt - (Kfeces + KGIb)*ALumen;
double RAGIRES   = KGIRESUP*AGIt - KGIRESrelease*AGIRES;
double RALb      = QL*(CA-CVL) + (QS*CVS) + (QGI*CVGI) - PAL*CVL + (PAL*CLt)/PL - KLRESUP*ALb + KLRESrelease*ALRES;
double RALt      = PAL*CVL - (PAL*CLt)/PL - Kbile*CLt;
double RALRES    = KLRESUP*ALb - KLRESrelease*ALRES;
double RKb       = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK - Kurine*CVK;
double RKt       = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES;
double RKRES     = KKRESUP*AKt-KKRESrelease*AKRES;
double Rurine    = Kurine*CVK;
double Rbile     = Kbile*CLt;
double Rfeces    = Kfeces*ALumen;

// #+ ODE equations for compartments in the female rat
dxdt_AA          = RAA;
dxdt_AV          = RAV;
dxdt_ARb         = RARb;
dxdt_ARt         = RARt;
dxdt_ARRES       = RARRES;
dxdt_Aua         = RAua;
dxdt_Atra        = Rtra;
dxdt_Apul        = RApul;
dxdt_ApulRES     = RApulRES;
dxdt_Ain         = RAin;
dxdt_AinRES      = RAinRES;
dxdt_ALub        = RALub;
dxdt_ASb         = RASb;
dxdt_ASt         = RASt;
dxdt_ASRES       = RASRES;
dxdt_AGIb        = RAGIb;
dxdt_AGIt        = RAGIt;
dxdt_ALumen      = RALumen;
dxdt_AGIRES      = RAGIRES;
dxdt_ALb         = RALb;
dxdt_ALt         = RALt;
dxdt_ALRES       = RALRES;
dxdt_AKb         = RKb;
dxdt_AKt         = RKt;
dxdt_AKRES       = RKRES;
dxdt_Aurine      = Rurine;
dxdt_Abile       = Rbile;
dxdt_Afeces      = Rfeces;

// #+ Total amount of NPs in tissues
double ABlood    = AA + AV;
double ALung     = Atra + Apul + Ain + ALub + ApulRES + AinRES;
double ASpleen   = ASb + ASt + ASRES;
double ARest     = ARb + ARt + ARRES;
double AGI       = AGIb + AGIt + ALumen + AGIRES;
double ALiver    = ALb + ALt + ALRES;
double AKidney   = AKb + AKt + AKRES;

// #+ Total concentrations of NPs in Tissues
double CBlood    = ABlood/VB;
double CLung     = ALung/VLu;
double CSPleen   = ASpleen/VS;
double CRest     = ARest/VR;
double CGI       = AGI/VGI;
double CLiver    = ALiver/VL;
double CKidney   = AKidney/VK;

// #+ AUC
dxdt_AUCB        = CBlood;
dxdt_AUCLu       = CLung;
dxdt_AUCS        = CSPleen;
dxdt_AUCRt       = CRest;
dxdt_AUCGI       = CGI;
dxdt_AUCL        = CLiver;
dxdt_AUCK        = CKidney;



// #+ Mass Balance
//double Tmass = ABlood + Aua + ALung + AGI + ALiver + ASpleen + AKidney + ARt + Aurine + Afeces;
//double BAL   = DOSE - Tmass

$TABLE
capture Blood     = ABlood;
capture Lung      = ALung;
capture GI        = AGI;
capture Spleen    = ASpleen;
capture Rest      = ARest;
capture Liver     = ALiver;
capture Kidney    = AKidney;
capture Urine     = Aurine;
capture Feces     = Afeces;
capture AUC_B     = AUCB;
capture AUC_Lu    = AUCLu;
capture AUC_S     = AUCS;
capture AUC_Rt    = AUCRt;
capture AUC_GI    = AUCGI;
capture AUC_L     = AUCL;
capture AUC_K     = AUCK;
'
