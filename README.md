# OPERA
OPERA is a free and open-source/open-data suite of QSAR models providing predictions on physicochemical properties, environmental fate and toxcicity endpoints as well as additional information including applicability domain and accuracy assessment. All models were built on curated data and standardized QSAR-ready chemical structures. OPERA is available in command line and user-friendly graphical interface for Windows and Linux operating systems. It can be installed as a standalone desktop application or embedded in a different tool/workflow. 


References:

[1] Mansouri K. et al. J Cheminform (2018) https://doi.org/10.1186/s13321-018-0263-1.

[2] Mansouri, K. et al. SAR and QSAR in Env. Res. (2016). https://doi.org/10.1080/1062936X.2016.1253611

[3] Williams A. J. et al. J Cheminform (2017) https://doi.org/10.1186/s13321-017-0247-6

[4] The CompTox Chemistry Dashboard https://comptox.epa.gov/dashboard

[5] JRC QSAR Model Database https://qsardb.jrc.ec.europa.eu/qmrf/endpoint

[6] Mansouri, K. et al. EHP (2016) https://doi.org/10.1289/ehp.1510267

[7] Mansouri, K. et al. J Cheminform (2019) https://doi.org/10.1186/s13321-019-0384-1

[8] Mansouri, K. et al. EHP (2020) https://doi.org/10.1289/EHP5580

[9] Mansouri, K. et al. EHP (2021) https://doi.org/10.1289/EHP8495

[10] Mansouri, K. et al. J Cheminform (2024) https://doi.org/10.1186/s13321-024-00814-3



Models:
----

  
          + Molecular descriptors:  
    - PaDEL (2.21) (https://doi.org/10.1002/jcc.21707 )
    - CDK (2.0) (https://doi.org/10.1186/s13321-017-0220-4)
 

          + List of models:
   
    - Caco-2 permeability (logPapp) model
    
    - FuB: Plasma fraction unbound (human, fraction)
    
    - Clint: hepatic intrinsic clearance (human, ul/min/10^6 cells)
    
    - pKa: acid dissociation constant

    - LogD: Octanol-water distribution constant. LogD is equivalent to logP for non-ionizable compounds.

    - CERAPP: Collaborative Estrogen Receptor Activity Prediction Project. Binding, Agonist, and Antagonist Estrogen Receptor activity (https://ehp.niehs.nih.gov/15-10267/)

    - CoMPARA: Collaborative Modeling Project for Androgen Receptor Activity. Binding, Agonist, and Antagonist Androgen Receptor activity (https://doi.org/10.13140/RG.2.2.19612.80009, https://doi.org/10.13140/RG.2.2.21850.03520)

    - CATMoS: Collaborative Acute Toxicity Modeling Suite. Very-Toxic, Non-Toxic, EPA categories, GHS categories, LD50 (Log mg/kg) (https://doi.org/10.1016/j.comtox.2018.08.002)
    
    - Structural Properties: MolWeight, nbAtoms, nbHeavyAtoms, nbC, nbO, nbN, nbAromAtom, nbRing, nbHeteroRing, Sp3Sp2HybRatio, nbRotBd, nbHBdAcc, ndHBdDon, nbLipinskiFailures, TopoPolSurfAir, MolarRefract, CombDipolPolarizability.
    
    - OH (LogOH) in cm3/molecule-sec: The OH rate constant for the atmospheric, gas-phase reaction between photochemically produced hydroxyl radicals and organic chemicals.

    - BCF (Log): Fish bioconcentration factor

    - Biodeg (LogHalfLife) in days: biodegradation half-life for compounds containing only carbon and hydrogen (i.e. hydrocarbons). 

    - Ready_biodeg (classification: 0/1): Ready biodegradability of organic chemicals. 

    - BP in deg C: Boiling Point at 760 mm Hg

    - HL (LogHL) in atm-m3/mole: The Henry’s Law constant (air/water partition coefficient) at 25C

    - Km (Log KmHL) half-lives in days: The whole body primary biotransformation rate constant for organic chemicals in fish. 

    - KOA (Log): The octanol/air partition coefficient.

    - LogP (Log): Octanol-water partition coefficient, log KOW, of chemicals.

    - MP in deg C: Melting Point

    - Koc (Log) in L/Kg: the soil adsorption coefficient of organic compounds.  The ratio of the amount of chemical adsorbed per unit weight of organic carbon in the soil or sediment to the concentration of the chemical in solution at equilibrium.

    - VP (Log) in mmHg: Vapor Pressure experimental values between 15 and 30 deg C (majority at 25-20C)

    - WS (Log) in Molar moles/L: Water solubility at 25C. 

    - RT in minutes: HPLC retention time.
