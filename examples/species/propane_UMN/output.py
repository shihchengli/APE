# Coordinates for propane in Input Orientation (angstroms):
#   C    1.2672    0.0000    0.2595
#   C    0.0000    0.0000   -0.5872
#   H    1.3053   -0.8807    0.9038
#   H    1.3053    0.8807    0.9038
#   H    2.1646    0.0000   -0.3600
#   C   -1.2672   -0.0000    0.2595
#   H    0.0000   -0.8739   -1.2430
#   H   -0.0000    0.8739   -1.2430
#   H   -1.3053   -0.8807    0.9038
#   H   -2.1646    0.0000   -0.3600
#   H   -1.3053    0.8807    0.9038

conformer(
    label = 'propane',
    E0 = (273289, 'J/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0626, 'amu')),
        NonlinearRotor(
            inertia = ([17.0849, 59.518, 67.2712], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(
            frequencies = ([247.72, 297.55, 381.05, 766.04, 887.89, 926.29, 946.74, 1074.84, 1189.44, 1224.37, 1328.88, 1380.86, 1418.85, 1431.46, 1499.66, 1502.64, 1504.91, 1517.68, 1521.16, 3039.01, 3040.88, 3047.33, 3074.94, 3106.69, 3116.67, 3117.22, 3119.35], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Thermodynamics for propane at 100 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0005672448
# Entropy (hartree/K): 0.0000005574
# Free energy (hartree): 0.0005115084
# Partition function: 0.1988481355

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0006747553
# Entropy (hartree/K): 0.0000002999
# Free energy (hartree): 0.0006447617
# Partition function: 0.1305506121

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0008747527
# Entropy (hartree/K): 0.0000000868
# Free energy (hartree): 0.0008660724
# Partition function: 0.0649045626

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0017511772
# Entropy (hartree/K): 0.0000000006
# Free energy (hartree): 0.0017511192
# Partition function: 0.0039676116

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0020216975
# Entropy (hartree/K): 0.0000000001
# Free energy (hartree): 0.0020216846
# Partition function: 0.0016884098

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0021090250
# Entropy (hartree/K): 0.0000000001
# Free energy (hartree): 0.0021090174
# Partition function: 0.0012814760

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0021563611
# Entropy (hartree/K): 0.0000000001
# Free energy (hartree): 0.0021563553
# Partition function: 0.0011035488

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0024545197
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0024545188
# Partition function: 0.0004304202

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0027080989
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0027080987
# Partition function: 0.0001932570

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0027920969
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0027920968
# Partition function: 0.0001482316

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0030215483
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0030215483
# Partition function: 0.0000718245

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0031393346
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0031393346
# Partition function: 0.0000495156

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0032386304
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0032386304
# Partition function: 0.0000361883

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0032681785
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0032681784
# Partition function: 0.0000329645

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0034159801
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034159801
# Partition function: 0.0000206706

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0034230533
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034230533
# Partition function: 0.0000202140

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0034287531
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034287531
# Partition function: 0.0000198534

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0034564266
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034564266
# Partition function: 0.0000181922

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0034648654
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034648654
# Partition function: 0.0000177138

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069455625
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455625
# Partition function: 0.0000000003

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069204253
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204253
# Partition function: 0.0000000003

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185284
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185284
# Partition function: 0.0000000003

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070564059
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070564059
# Partition function: 0.0000000002

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071174979
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174979
# Partition function: 0.0000000002

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071492938
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071492938
# Partition function: 0.0000000002

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071282071
# Entropy (hartree/K): -0.0000000000
# Free energy (hartree): 0.0071282071
# Partition function: 0.0000000002

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071238545
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071238545
# Partition function: 0.0000000002

# 	********** Final results **********


# Temperature (K): 100.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 0.2980809751
# Translational entropy (cal/mol/K): 31.8740751550
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.2980809751
# Rotational entropy (cal/mol/K): 18.0253188155
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 65.4657312262
# Internal (tor+vib) entropy (cal/mol/K): 0.5929699872
# Internal (tor+vib) Cv (cal/mol/K): 1.7775170353


# Total energy (kcal/mol): 66.0618931765
# Total enthalpy (kcal/mol): 66.2606138379
# Enthalpy H(100.000000 K)-H(0 K) (kcal/mol):  0.8414521993
# Total entropy (cal/mol/K): 50.4923639577
# Total Cv (cal/mol/K): 7.7391365381



# 	********** HOhf results **********


# Translational energy (kcal/mol): 0.2980809751
# Rotational energy (kcal/mol): 0.2980809751
# Vibrational energy (kcal/mol): 65.3831742376
# gas constant (RT): 0.1987206501
# Translational entropy (cal/mol/K): 31.8740751550
# Rotational entropy (cal/mol/K): 18.0253188155
# Vibrational entropy (cal/mol/K): 0.4648323619


# Total energy (kcal/mol): 65.9793361879
# Total enthalpy (kcal/mol): 66.1780568380
# Enthalpy H(100.000000 K)-H(0 K) (kcal/mol): 0.8320568380
# Total entropy (cal/mol/K): 50.3642263324
# Total Cv (cal/mol/K): 9.4804127129
# Thermodynamics for propane at 200 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0007838437
# Entropy (hartree/K): 0.0000020116
# Free energy (hartree): 0.0003815333
# Partition function: 0.5475000375

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0008545522
# Entropy (hartree/K): 0.0000014898
# Free energy (hartree): 0.0005565938
# Partition function: 0.4152846599

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0009877477
# Entropy (hartree/K): 0.0000008144
# Free energy (hartree): 0.0008248586
# Partition function: 0.2718928296

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0017648481
# Entropy (hartree/K): 0.0000000810
# Free energy (hartree): 0.0017486552
# Partition function: 0.0632345030

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0020286262
# Entropy (hartree/K): 0.0000000401
# Free energy (hartree): 0.0020205964
# Partition function: 0.0411609276

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0021144820
# Entropy (hartree/K): 0.0000000314
# Free energy (hartree): 0.0021081966
# Partition function: 0.0358441291

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0021611511
# Entropy (hartree/K): 0.0000000275
# Free energy (hartree): 0.0021556511
# Partition function: 0.0332566628

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0024565976
# Entropy (hartree/K): 0.0000000117
# Free energy (hartree): 0.0024542513
# Partition function: 0.0207553321

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0027091604
# Entropy (hartree/K): 0.0000000059
# Free energy (hartree): 0.0027079743
# Partition function: 0.0139044192

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0027929280
# Entropy (hartree/K): 0.0000000046
# Free energy (hartree): 0.0027920025
# Partition function: 0.0121768529

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0030219872
# Entropy (hartree/K): 0.0000000024
# Free energy (hartree): 0.0030215022
# Partition function: 0.0084755528

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0031396561
# Entropy (hartree/K): 0.0000000018
# Free energy (hartree): 0.0031393021
# Partition function: 0.0070370904

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0032388619
# Entropy (hartree/K): 0.0000000013
# Free energy (hartree): 0.0032386078
# Partition function: 0.0060158861

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0032683915
# Entropy (hartree/K): 0.0000000012
# Free energy (hartree): 0.0032681578
# Partition function: 0.0057416572

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0034161212
# Entropy (hartree/K): 0.0000000008
# Free energy (hartree): 0.0034159670
# Partition function: 0.0045465826

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0034231918
# Entropy (hartree/K): 0.0000000008
# Free energy (hartree): 0.0034230405
# Partition function: 0.0044960884

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0034288890
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0034287405
# Partition function: 0.0044558067

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0034565523
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0034564151
# Partition function: 0.0042653041

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0034649884
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0034648542
# Partition function: 0.0042088491

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069455625
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455625
# Partition function: 0.0000172765

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069204253
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204253
# Partition function: 0.0000179760

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185284
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185284
# Partition function: 0.0000180299

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070564059
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070564059
# Partition function: 0.0000145028

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071174979
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174979
# Partition function: 0.0000131692

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071492938
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071492938
# Partition function: 0.0000125244

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071282071
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282071
# Partition function: 0.0000129484

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071238545
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071238545
# Partition function: 0.0000130377

# 	********** Final results **********


# Temperature (K): 200.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 0.5961619503
# Translational entropy (cal/mol/K): 35.3176416133
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.5961619503
# Rotational entropy (cal/mol/K): 20.0914586905
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 65.8083999876
# Internal (tor+vib) entropy (cal/mol/K): 2.8416029337
# Internal (tor+vib) Cv (cal/mol/K): 5.0055986032


# Total energy (kcal/mol): 67.0007238882
# Total enthalpy (kcal/mol): 67.3981652109
# Enthalpy H(200.000000 K)-H(0 K) (kcal/mol):  1.9790035724
# Total entropy (cal/mol/K): 58.2507032376
# Total Cv (cal/mol/K): 10.9672181061



# 	********** HOhf results **********


# Translational energy (kcal/mol): 0.5961619503
# Rotational energy (kcal/mol): 0.5961619503
# Vibrational energy (kcal/mol): 65.7010182789
# gas constant (RT): 0.3974413002
# Translational entropy (cal/mol/K): 35.3176416133
# Rotational entropy (cal/mol/K): 20.0914586905
# Vibrational entropy (cal/mol/K): 2.5405990110


# Total energy (kcal/mol): 66.8933421795
# Total enthalpy (kcal/mol): 67.2907834797
# Enthalpy H(200.000000 K)-H(0 K) (kcal/mol): 1.9447834797
# Total entropy (cal/mol/K): 57.9496993148
# Total Cv (cal/mol/K): 12.7200305784
# Thermodynamics for propane at 298.15 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000003 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0010619421
# Entropy (hartree/K): 0.0000031382
# Free energy (hartree): 0.0001262911
# Partition function: 0.8748029866

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0011168697
# Entropy (hartree/K): 0.0000025502
# Free energy (hartree): 0.0003565283
# Partition function: 0.6855019735

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0011967811
# Entropy (hartree/K): 0.0000016562
# Free energy (hartree): 0.0007029906
# Partition function: 0.4749494972

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0018380836
# Entropy (hartree/K): 0.0000003691
# Free energy (hartree): 0.0017280337
# Partition function: 0.1603841017

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0020788453
# Entropy (hartree/K): 0.0000002363
# Free energy (hartree): 0.0020083973
# Partition function: 0.1191799694

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0021582933
# Entropy (hartree/K): 0.0000002021
# Free energy (hartree): 0.0020980346
# Partition function: 0.1083859671

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0022017903
# Entropy (hartree/K): 0.0000001856
# Free energy (hartree): 0.0021464506
# Partition function: 0.1029682401

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0024814809
# Entropy (hartree/K): 0.0000001077
# Free energy (hartree): 0.0024493593
# Partition function: 0.0747093941

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0027257880
# Entropy (hartree/K): 0.0000000697
# Free energy (hartree): 0.0027050145
# Partition function: 0.0569878494

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0028072569
# Entropy (hartree/K): 0.0000000594
# Free energy (hartree): 0.0027895357
# Partition function: 0.0521081031

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0030316747
# Entropy (hartree/K): 0.0000000393
# Free energy (hartree): 0.0030199654
# Partition function: 0.0408239664

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0031476465
# Entropy (hartree/K): 0.0000000321
# Free energy (hartree): 0.0031380809
# Partition function: 0.0360235030

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0032453770
# Entropy (hartree/K): 0.0000000259
# Free energy (hartree): 0.0032376488
# Partition function: 0.0324181327

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0032745782
# Entropy (hartree/K): 0.0000000246
# Free energy (hartree): 0.0032672555
# Partition function: 0.0314173720

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0034209029
# Entropy (hartree/K): 0.0000000188
# Free energy (hartree): 0.0034153001
# Partition function: 0.0268580352

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0034279165
# Entropy (hartree/K): 0.0000000186
# Free energy (hartree): 0.0034223828
# Partition function: 0.0266573158

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0034335593
# Entropy (hartree/K): 0.0000000183
# Free energy (hartree): 0.0034280917
# Partition function: 0.0264966226

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0034609992
# Entropy (hartree/K): 0.0000000174
# Free energy (hartree): 0.0034558023
# Partition function: 0.0257302811

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0034693748
# Entropy (hartree/K): 0.0000000172
# Free energy (hartree): 0.0034642511
# Partition function: 0.0255010680

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069455679
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455621
# Partition function: 0.0006386600

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069204318
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204249
# Partition function: 0.0006558916

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185358
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185279
# Partition function: 0.0006572107

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070564101
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070564056
# Partition function: 0.0005679176

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071175017
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174977
# Partition function: 0.0005323350

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071492974
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071492936
# Partition function: 0.0005147068

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071282109
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282069
# Partition function: 0.0005263312

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071238586
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071238542
# Partition function: 0.0005287632

# 	********** Final results **********


# Temperature (K): 298.15
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 0.8887284274
# Translational entropy (cal/mol/K): 37.3012679085
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.8887284274
# Rotational entropy (cal/mol/K): 21.2816344676
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 66.4777157124
# Internal (tor+vib) entropy (cal/mol/K): 5.5138076238
# Internal (tor+vib) Cv (cal/mol/K): 8.8373562788


# Total energy (kcal/mol): 68.2551725672
# Total enthalpy (kcal/mol): 68.8476582191
# Enthalpy H(298.150000 K)-H(0 K) (kcal/mol):  3.4284965805
# Total entropy (cal/mol/K): 64.0967099999
# Total Cv (cal/mol/K): 14.7989757817



# 	********** HOhf results **********


# Translational energy (kcal/mol): 0.8887284274
# Rotational energy (kcal/mol): 0.8887284274
# Vibrational energy (kcal/mol): 66.3474139282
# gas constant (RT): 0.5924856183
# Translational entropy (cal/mol/K): 37.3012679085
# Rotational entropy (cal/mol/K): 21.2816344676
# Vibrational entropy (cal/mol/K): 5.1196683348


# Total energy (kcal/mol): 68.1248707829
# Total enthalpy (kcal/mol): 68.7173564012
# Enthalpy H(298.150000 K)-H(0 K) (kcal/mol): 3.3713564012
# Total entropy (cal/mol/K): 63.7025707109
# Total Cv (cal/mol/K): 16.5430222039
# Thermodynamics for propane at 300 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000004 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0010674705
# Entropy (hartree/K): 0.0000031567
# Free energy (hartree): 0.0001204683
# Partition function: 0.8809073631

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0011222262
# Entropy (hartree/K): 0.0000025681
# Free energy (hartree): 0.0003517938
# Partition function: 0.6905326985

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0012012486
# Entropy (hartree/K): 0.0000016711
# Free energy (hartree): 0.0006999128
# Partition function: 0.4786834123

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0018401519
# Entropy (hartree/K): 0.0000003760
# Free energy (hartree): 0.0017273444
# Partition function: 0.1623221908

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0020803932
# Entropy (hartree/K): 0.0000002415
# Free energy (hartree): 0.0020079554
# Partition function: 0.1208097609

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0021596841
# Entropy (hartree/K): 0.0000002068
# Free energy (hartree): 0.0020976564
# Partition function: 0.1099251182

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0022031011
# Entropy (hartree/K): 0.0000001900
# Free energy (hartree): 0.0021461031
# Partition function: 0.1044601004

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0024823664
# Entropy (hartree/K): 0.0000001107
# Free energy (hartree): 0.0024491573
# Partition function: 0.0759302963

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0027264259
# Entropy (hartree/K): 0.0000000718
# Free energy (hartree): 0.0027048836
# Partition function: 0.0580115897

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0028078213
# Entropy (hartree/K): 0.0000000613
# Free energy (hartree): 0.0027894240
# Partition function: 0.0530724019

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0030320826
# Entropy (hartree/K): 0.0000000406
# Free energy (hartree): 0.0030198915
# Partition function: 0.0416404111

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0031479936
# Entropy (hartree/K): 0.0000000332
# Free energy (hartree): 0.0031380205
# Partition function: 0.0367717769

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0032456693
# Entropy (hartree/K): 0.0000000269
# Free energy (hartree): 0.0032375999
# Partition function: 0.0331126395

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0032748580
# Entropy (hartree/K): 0.0000000255
# Free energy (hartree): 0.0032672092
# Partition function: 0.0320965588

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0034211277
# Entropy (hartree/K): 0.0000000195
# Free energy (hartree): 0.0034152646
# Partition function: 0.0274648874

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0034281390
# Entropy (hartree/K): 0.0000000193
# Free energy (hartree): 0.0034223478
# Partition function: 0.0272608813

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0034337797
# Entropy (hartree/K): 0.0000000191
# Free energy (hartree): 0.0034280571
# Partition function: 0.0270975482

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0034612105
# Entropy (hartree/K): 0.0000000181
# Free energy (hartree): 0.0034557694
# Partition function: 0.0263185421

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0034695837
# Entropy (hartree/K): 0.0000000179
# Free energy (hartree): 0.0034642187
# Partition function: 0.0260855153

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069455684
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455621
# Partition function: 0.0006682987

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069204324
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204248
# Partition function: 0.0006862172

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185365
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185278
# Partition function: 0.0006875888

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070564105
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070564056
# Partition function: 0.0005947037

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071175021
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174976
# Partition function: 0.0005576652

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071492977
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071492936
# Partition function: 0.0005393102

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071282113
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282068
# Partition function: 0.0005514143

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071238590
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071238542
# Partition function: 0.0005539464

# 	********** Final results **********


# Temperature (K): 300.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 0.8942429254
# Translational entropy (cal/mol/K): 37.3319988602
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.8942429254
# Rotational entropy (cal/mol/K): 21.3000730386
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 66.4941417227
# Internal (tor+vib) entropy (cal/mol/K): 5.5687302461
# Internal (tor+vib) Cv (cal/mol/K): 8.9205451235


# Total energy (kcal/mol): 68.2826275735
# Total enthalpy (kcal/mol): 68.8787895576
# Enthalpy H(300.000000 K)-H(0 K) (kcal/mol):  3.4596279191
# Total entropy (cal/mol/K): 64.2008021449
# Total Cv (cal/mol/K): 14.8821646264



# 	********** HOhf results **********


# Translational energy (kcal/mol): 0.8942429254
# Rotational energy (kcal/mol): 0.8942429254
# Vibrational energy (kcal/mol): 66.3633893816
# gas constant (RT): 0.5961619503
# Translational entropy (cal/mol/K): 37.3319988602
# Rotational entropy (cal/mol/K): 21.3000730386
# Vibrational entropy (cal/mol/K): 5.1730844538


# Total energy (kcal/mol): 68.1518752325
# Total enthalpy (kcal/mol): 68.7480371827
# Enthalpy H(300.000000 K)-H(0 K) (kcal/mol): 3.4020371827
# Total entropy (cal/mol/K): 63.8051563526
# Total Cv (cal/mol/K): 16.6254390390
# Thermodynamics for propane at 400 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000407 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0013769920
# Entropy (hartree/K): 0.0000040457
# Free energy (hartree): -0.0002413039
# Partition function: 1.2098472896

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000043 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0014264563
# Entropy (hartree/K): 0.0000034415
# Free energy (hartree): 0.0000498634
# Partition function: 0.9614006693

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0014597521
# Entropy (hartree/K): 0.0000024128
# Free energy (hartree): 0.0004946509
# Partition function: 0.6767214531

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0019843040
# Entropy (hartree/K): 0.0000007866
# Free energy (hartree): 0.0016696656
# Partition function: 0.2676446120

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0021961917
# Entropy (hartree/K): 0.0000005704
# Free energy (hartree): 0.0019680352
# Partition function: 0.2114769649

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0022663507
# Entropy (hartree/K): 0.0000005095
# Free energy (hartree): 0.0020625625
# Partition function: 0.1962702901

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0023050201
# Entropy (hartree/K): 0.0000004791
# Free energy (hartree): 0.0021133865
# Partition function: 0.1885513648

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0025576899
# Entropy (hartree/K): 0.0000003237
# Free energy (hartree): 0.0024282287
# Partition function: 0.1470573198

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0027849579
# Entropy (hartree/K): 0.0000002368
# Free energy (hartree): 0.0026902194
# Partition function: 0.1195814105

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0028610518
# Entropy (hartree/K): 0.0000002113
# Free energy (hartree): 0.0027765435
# Partition function: 0.1117037154

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0030734684
# Entropy (hartree/K): 0.0000001569
# Free energy (hartree): 0.0030106981
# Partition function: 0.0928513396

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0031845196
# Entropy (hartree/K): 0.0000001358
# Free energy (hartree): 0.0031302190
# Partition function: 0.0844910443

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0032776085
# Entropy (hartree/K): 0.0000001164
# Free energy (hartree): 0.0032310402
# Partition function: 0.0780268984

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0033057296
# Entropy (hartree/K): 0.0000001120
# Free energy (hartree): 0.0032609301
# Partition function: 0.0762073102

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0034471743
# Entropy (hartree/K): 0.0000000924
# Free energy (hartree): 0.0034102098
# Partition function: 0.0677355072

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0034539801
# Entropy (hartree/K): 0.0000000916
# Free energy (hartree): 0.0034173435
# Partition function: 0.0673551198

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0034594228
# Entropy (hartree/K): 0.0000000908
# Free energy (hartree): 0.0034231012
# Partition function: 0.0670496612

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0034860327
# Entropy (hartree/K): 0.0000000875
# Free energy (hartree): 0.0034510133
# Partition function: 0.0655883938

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0034941819
# Entropy (hartree/K): 0.0000000867
# Free energy (hartree): 0.0034595167
# Partition function: 0.0651495777

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069457962
# Entropy (hartree/K): 0.0000000006
# Free energy (hartree): 0.0069455412
# Partition function: 0.0041565700

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069206899
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0069204009
# Partition function: 0.0042398880

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069188203
# Entropy (hartree/K): 0.0000000008
# Free energy (hartree): 0.0069185012
# Partition function: 0.0042462514

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070565990
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0070563887
# Partition function: 0.0038083009

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071176771
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071174820
# Partition function: 0.0036289879

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071494638
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071492788
# Partition function: 0.0035390284

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071283860
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071281913
# Partition function: 0.0035984368

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071240436
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071238377
# Partition function: 0.0036108255

# 	********** Final results **********


# Temperature (K): 400.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 1.1923239006
# Translational entropy (cal/mol/K): 38.7612080717
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 1.1923239006
# Rotational entropy (cal/mol/K): 22.1575985655
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 67.6219077018
# Internal (tor+vib) entropy (cal/mol/K): 8.7801364331
# Internal (tor+vib) Cv (cal/mol/K): 13.6945851259


# Total energy (kcal/mol): 70.0065555030
# Total enthalpy (kcal/mol): 70.8014381484
# Enthalpy H(400.000000 K)-H(0 K) (kcal/mol):  5.3822765099
# Total entropy (cal/mol/K): 69.6989430703
# Total Cv (cal/mol/K): 19.6562046288



# 	********** HOhf results **********


# Translational energy (kcal/mol): 1.1923239006
# Rotational energy (kcal/mol): 1.1923239006
# Vibrational energy (kcal/mol): 67.4635420211
# gas constant (RT): 0.7948826004
# Translational entropy (cal/mol/K): 38.7612080717
# Rotational entropy (cal/mol/K): 22.1575985655
# Vibrational entropy (cal/mol/K): 8.3055671947


# Total energy (kcal/mol): 69.8481898223
# Total enthalpy (kcal/mol): 70.6430724227
# Enthalpy H(400.000000 K)-H(0 K) (kcal/mol): 5.2970724227
# Total entropy (cal/mol/K): 69.2243738319
# Total Cv (cal/mol/K): 21.3240893399
# Thermodynamics for propane at 500 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000007540 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0017052192
# Entropy (hartree/K): 0.0000047774
# Free energy (hartree): -0.0006834816
# Partition function: 1.5397997064

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000001151 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0017538592
# Entropy (hartree/K): 0.0000041712
# Free energy (hartree): -0.0003317386
# Partition function: 1.2330728160

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0017395780
# Entropy (hartree/K): 0.0000030365
# Free energy (hartree): 0.0002213055
# Partition function: 0.8695622499

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0021777796
# Entropy (hartree/K): 0.0000012167
# Free energy (hartree): 0.0015694415
# Partition function: 0.3711385731

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0023648827
# Entropy (hartree/K): 0.0000009450
# Free energy (hartree): 0.0018924064
# Partition function: 0.3026591368

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0024261856
# Entropy (hartree/K): 0.0000008642
# Free energy (hartree): 0.0019940707
# Partition function: 0.2838373040

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0024601346
# Entropy (hartree/K): 0.0000008233
# Free energy (hartree): 0.0020484886
# Partition function: 0.2742482298

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0026844495
# Entropy (hartree/K): 0.0000006045
# Free energy (hartree): 0.0023821767
# Partition function: 0.2221367303

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0028921509
# Entropy (hartree/K): 0.0000004741
# Free energy (hartree): 0.0026551045
# Partition function: 0.1869658034

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0029615182
# Entropy (hartree/K): 0.0000004335
# Free energy (hartree): 0.0027447513
# Partition function: 0.1766745574

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0031580806
# Entropy (hartree/K): 0.0000003439
# Free energy (hartree): 0.0029861315
# Partition function: 0.1516940764

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0032622063
# Entropy (hartree/K): 0.0000003073
# Free energy (hartree): 0.0031085458
# Partition function: 0.1404083837

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0033483215
# Entropy (hartree/K): 0.0000002725
# Free energy (hartree): 0.0032120721
# Partition function: 0.1315218774

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0033747909
# Entropy (hartree/K): 0.0000002644
# Free energy (hartree): 0.0032425870
# Partition function: 0.1290115042

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0035085413
# Entropy (hartree/K): 0.0000002277
# Free energy (hartree): 0.0033946728
# Partition function: 0.1171964813

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0035150097
# Entropy (hartree/K): 0.0000002262
# Free energy (hartree): 0.0034019255
# Partition function: 0.1166608983

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0035201230
# Entropy (hartree/K): 0.0000002247
# Free energy (hartree): 0.0034077980
# Partition function: 0.1162290309

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0035453672
# Entropy (hartree/K): 0.0000002184
# Free energy (hartree): 0.0034361850
# Partition function: 0.1141638675

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0035531440
# Entropy (hartree/K): 0.0000002167
# Free energy (hartree): 0.0034448178
# Partition function: 0.1135431351

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069476705
# Entropy (hartree/K): 0.0000000047
# Free energy (hartree): 0.0069453228
# Partition function: 0.0124465209

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069227470
# Entropy (hartree/K): 0.0000000052
# Free energy (hartree): 0.0069201581
# Partition function: 0.0126459111

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069210353
# Entropy (hartree/K): 0.0000000056
# Free energy (hartree): 0.0069182369
# Partition function: 0.0126612637

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070582216
# Entropy (hartree/K): 0.0000000040
# Free energy (hartree): 0.0070562034
# Partition function: 0.0116047512

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071192108
# Entropy (hartree/K): 0.0000000038
# Free energy (hartree): 0.0071173082
# Partition function: 0.0111654464

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071509379
# Entropy (hartree/K): 0.0000000037
# Free energy (hartree): 0.0071491127
# Partition function: 0.0109434136

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071299180
# Entropy (hartree/K): 0.0000000038
# Free energy (hartree): 0.0071280177
# Partition function: 0.0110901831

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071256410
# Entropy (hartree/K): 0.0000000040
# Free energy (hartree): 0.0071236556
# Partition function: 0.0111207770

# 	********** Final results **********


# Temperature (K): 500.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 1.4904048757
# Translational entropy (cal/mol/K): 39.8697888612
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 1.4904048757
# Rotational entropy (cal/mol/K): 22.8227470392
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 69.2310704922
# Internal (tor+vib) entropy (cal/mol/K): 12.3512048791
# Internal (tor+vib) Cv (cal/mol/K): 18.4310670188


# Total energy (kcal/mol): 72.2118802436
# Total enthalpy (kcal/mol): 73.2054835504
# Enthalpy H(500.000000 K)-H(0 K) (kcal/mol):  7.7863219119
# Total entropy (cal/mol/K): 75.0437407795
# Total Cv (cal/mol/K): 24.3926865217



# 	********** HOhf results **********


# Translational energy (kcal/mol): 1.4904048757
# Rotational energy (kcal/mol): 1.4904048757
# Vibrational energy (kcal/mol): 69.0344260958
# gas constant (RT): 0.9936032505
# Translational entropy (cal/mol/K): 39.8697888612
# Rotational entropy (cal/mol/K): 22.8227470392
# Vibrational entropy (cal/mol/K): 11.7917806347


# Total energy (kcal/mol): 72.0152358472
# Total enthalpy (kcal/mol): 73.0088390977
# Enthalpy H(500.000000 K)-H(0 K) (kcal/mol): 7.6628390977
# Total entropy (cal/mol/K): 74.4843165352
# Total Cv (cal/mol/K): 25.9252492424
# Thermodynamics for propane at 600 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000055602 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0020517166
# Entropy (hartree/K): 0.0000054086
# Free energy (hartree): -0.0011934664
# Partition function: 1.8740667505

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000011070 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0021012364
# Entropy (hartree/K): 0.0000048040
# Free energy (hartree): -0.0007811677
# Partition function: 1.5085088139

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0020307541
# Entropy (hartree/K): 0.0000035672
# Free energy (hartree): -0.0001095590
# Partition function: 1.0593546519

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0024035839
# Entropy (hartree/K): 0.0000016276
# Free energy (hartree): 0.0014269950
# Partition function: 0.4718874397

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0025712624
# Entropy (hartree/K): 0.0000013204
# Free energy (hartree): 0.0017790451
# Partition function: 0.3920775314

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0026248770
# Entropy (hartree/K): 0.0000012256
# Free energy (hartree): 0.0018895248
# Partition function: 0.3699304985

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0026546851
# Entropy (hartree/K): 0.0000011771
# Free energy (hartree): 0.0019484363
# Partition function: 0.3586369242

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0028529204
# Entropy (hartree/K): 0.0000009107
# Free energy (hartree): 0.0023065000
# Partition function: 0.2970394011

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0030420684
# Entropy (hartree/K): 0.0000007464
# Free energy (hartree): 0.0025942497
# Partition function: 0.2552962620

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0031045748
# Entropy (hartree/K): 0.0000006933
# Free energy (hartree): 0.0026886011
# Partition function: 0.2429287992

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0032845842
# Entropy (hartree/K): 0.0000005735
# Free energy (hartree): 0.0029405068
# Partition function: 0.2127659305

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0033812567
# Entropy (hartree/K): 0.0000005233
# Free energy (hartree): 0.0030672798
# Partition function: 0.1990334978

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0034593457
# Entropy (hartree/K): 0.0000004739
# Free energy (hartree): 0.0031750336
# Partition function: 0.1880604104

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0034839220
# Entropy (hartree/K): 0.0000004623
# Free energy (hartree): 0.0032065338
# Partition function: 0.1849684004

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0036087278
# Entropy (hartree/K): 0.0000004093
# Free energy (hartree): 0.0033631182
# Partition function: 0.1703365219

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0036147987
# Entropy (hartree/K): 0.0000004071
# Free energy (hartree): 0.0033705646
# Partition function: 0.1696702794

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0036195139
# Entropy (hartree/K): 0.0000004048
# Free energy (hartree): 0.0033766257
# Partition function: 0.1691299122

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0036431218
# Entropy (hartree/K): 0.0000003955
# Free energy (hartree): 0.0034057927
# Partition function: 0.1665535401

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0036504584
# Entropy (hartree/K): 0.0000003930
# Free energy (hartree): 0.0034146375
# Partition function: 0.1657800390

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069546984
# Entropy (hartree/K): 0.0000000173
# Free energy (hartree): 0.0069443165
# Partition function: 0.0258684256

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069303076
# Entropy (hartree/K): 0.0000000187
# Free energy (hartree): 0.0069190607
# Partition function: 0.0262145608

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069290472
# Entropy (hartree/K): 0.0000000200
# Free energy (hartree): 0.0069170614
# Partition function: 0.0262421594

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070644977
# Entropy (hartree/K): 0.0000000153
# Free energy (hartree): 0.0070553229
# Partition function: 0.0244004460

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071252154
# Entropy (hartree/K): 0.0000000146
# Free energy (hartree): 0.0071164725
# Partition function: 0.0236276802

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071567582
# Entropy (hartree/K): 0.0000000141
# Free energy (hartree): 0.0071483071
# Partition function: 0.0232351139

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071359175
# Entropy (hartree/K): 0.0000000146
# Free energy (hartree): 0.0071271829
# Partition function: 0.0234948717

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071318411
# Entropy (hartree/K): 0.0000000151
# Free energy (hartree): 0.0071227879
# Partition function: 0.0235492794

# 	********** Final results **********


# Temperature (K): 600.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 1.7884858509
# Translational entropy (cal/mol/K): 40.7755653185
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 1.7884858509
# Rotational entropy (cal/mol/K): 23.3662129136
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 71.2924185561
# Internal (tor+vib) entropy (cal/mol/K): 16.0976355620
# Internal (tor+vib) Cv (cal/mol/K): 22.7078064808


# Total energy (kcal/mol): 74.8693902578
# Total enthalpy (kcal/mol): 76.0617142259
# Enthalpy H(600.000000 K)-H(0 K) (kcal/mol):  10.6425525874
# Total entropy (cal/mol/K): 80.2394137941
# Total Cv (cal/mol/K): 28.6694259837



# 	********** HOhf results **********


# Translational energy (kcal/mol): 1.7884858509
# Rotational energy (kcal/mol): 1.7884858509
# Vibrational energy (kcal/mol): 71.0420133855
# gas constant (RT): 1.1923239006
# Translational entropy (cal/mol/K): 40.7755653185
# Rotational entropy (cal/mol/K): 23.3662129136
# Vibrational entropy (cal/mol/K): 15.4406613413


# Total energy (kcal/mol): 74.6189850873
# Total enthalpy (kcal/mol): 75.8113089878
# Enthalpy H(600.000000 K)-H(0 K) (kcal/mol): 10.4653089878
# Total entropy (cal/mol/K): 79.5824395734
# Total Cv (cal/mol/K): 30.0337003190
# Thermodynamics for propane at 700 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000238669 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0024154931
# Entropy (hartree/K): 0.0000059691
# Free energy (hartree): -0.0017628652
# Partition function: 2.2149731648

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000057875 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0024668758
# Entropy (hartree/K): 0.0000053673
# Free energy (hartree): -0.0012902363
# Partition function: 1.7896753094

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0023285397
# Entropy (hartree/K): 0.0000040261
# Free energy (hartree): -0.0004897484
# Partition function: 1.2472346098

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0026506826
# Entropy (hartree/K): 0.0000020082
# Free energy (hartree): 0.0012449356
# Partition function: 0.5702965703

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0028038773
# Entropy (hartree/K): 0.0000016785
# Free energy (hartree): 0.0016289222
# Partition function: 0.4795928893

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0028510011
# Entropy (hartree/K): 0.0000015737
# Free energy (hartree): 0.0017494079
# Partition function: 0.4542218355

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0028773193
# Entropy (hartree/K): 0.0000015198
# Free energy (hartree): 0.0018134545
# Partition function: 0.4412862783

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0030526981
# Entropy (hartree/K): 0.0000012181
# Free energy (hartree): 0.0022000079
# Partition function: 0.3706717922

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0032257515
# Entropy (hartree/K): 0.0000010289
# Free energy (hartree): 0.0025055012
# Partition function: 0.3229530417

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0032818405
# Entropy (hartree/K): 0.0000009660
# Free energy (hartree): 0.0026056754
# Partition function: 0.3086838801

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0034463312
# Entropy (hartree/K): 0.0000008222
# Free energy (hartree): 0.0028708112
# Partition function: 0.2738863321

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0035359647
# Entropy (hartree/K): 0.0000007611
# Free energy (hartree): 0.0030031670
# Partition function: 0.2580121299

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0036058414
# Entropy (hartree/K): 0.0000006990
# Free energy (hartree): 0.0031165152
# Partition function: 0.2451510237

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0036285226
# Entropy (hartree/K): 0.0000006846
# Free energy (hartree): 0.0031493202
# Partition function: 0.2415498609

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0037443613
# Entropy (hartree/K): 0.0000006178
# Free energy (hartree): 0.0033119152
# Partition function: 0.2244668862

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0037500327
# Entropy (hartree/K): 0.0000006149
# Free energy (hartree): 0.0033196228
# Partition function: 0.2236877858

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0037543312
# Entropy (hartree/K): 0.0000006120
# Free energy (hartree): 0.0033259411
# Partition function: 0.2230511320

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0037762582
# Entropy (hartree/K): 0.0000006001
# Free energy (hartree): 0.0033561679
# Partition function: 0.2200303560

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0037831551
# Entropy (hartree/K): 0.0000005969
# Free energy (hartree): 0.0033652992
# Partition function: 0.2191258680

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0069716245
# Entropy (hartree/K): 0.0000000432
# Free energy (hartree): 0.0069414179
# Partition function: 0.0436601276

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069482599
# Entropy (hartree/K): 0.0000000462
# Free energy (hartree): 0.0069159442
# Partition function: 0.0441647352

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069478586
# Entropy (hartree/K): 0.0000000487
# Free energy (hartree): 0.0069137598
# Partition function: 0.0442082765

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0070799511
# Entropy (hartree/K): 0.0000000389
# Free energy (hartree): 0.0070527312
# Partition function: 0.0415219044

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071401300
# Entropy (hartree/K): 0.0000000373
# Free energy (hartree): 0.0071139914
# Partition function: 0.0403901616

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071713030
# Entropy (hartree/K): 0.0000000363
# Free energy (hartree): 0.0071459010
# Partition function: 0.0398129244

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071508227
# Entropy (hartree/K): 0.0000000373
# Free energy (hartree): 0.0071247038
# Partition function: 0.0401954502

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071471480
# Entropy (hartree/K): 0.0000000385
# Free energy (hartree): 0.0071202270
# Partition function: 0.0402767075

# 	********** Final results **********


# Temperature (K): 700.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 2.0865668260
# Translational entropy (cal/mol/K): 41.5413884012
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.0865668260
# Rotational entropy (cal/mol/K): 23.8257067632
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 73.7549440447
# Internal (tor+vib) entropy (cal/mol/K): 19.8862110310
# Internal (tor+vib) Cv (cal/mol/K): 26.4566203415


# Total energy (kcal/mol): 77.9280776967
# Total enthalpy (kcal/mol): 79.3191223262
# Enthalpy H(700.000000 K)-H(0 K) (kcal/mol):  13.8999606877
# Total entropy (cal/mol/K): 85.2533061954
# Total Cv (cal/mol/K): 32.4182398444



# 	********** HOhf results **********


# Translational energy (kcal/mol): 2.0865668260
# Rotational energy (kcal/mol): 2.0865668260
# Vibrational energy (kcal/mol): 73.4338723279
# gas constant (RT): 1.3910445507
# Translational entropy (cal/mol/K): 41.5413884012
# Rotational entropy (cal/mol/K): 23.8257067632
# Vibrational entropy (cal/mol/K): 19.1206293568


# Total energy (kcal/mol): 77.6070059799
# Total enthalpy (kcal/mol): 78.9980505305
# Enthalpy H(700.000000 K)-H(0 K) (kcal/mol): 13.6520505305
# Total entropy (cal/mol/K): 84.4877245212
# Total Cv (cal/mol/K): 33.6186876451
# Thermodynamics for propane at 800 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000725569 		 0.0000000035
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705523
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0027934652
# Entropy (hartree/K): 0.0000064736
# Free energy (hartree): -0.0023854266
# Partition function: 2.5640080020

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000204807 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0028480181
# Entropy (hartree/K): 0.0000058760
# Free energy (hartree): -0.0018528169
# Partition function: 2.0778647592

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0026304322
# Entropy (hartree/K): 0.0000044292
# Free energy (hartree): -0.0009129268
# Partition function: 1.4338296520

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0029121930
# Entropy (hartree/K): 0.0000023572
# Free energy (hartree): 0.0010264076
# Partition function: 0.6668821683

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0030550090
# Entropy (hartree/K): 0.0000020136
# Free energy (hartree): 0.0014441168
# Partition function: 0.5655144738

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0030966353
# Entropy (hartree/K): 0.0000019015
# Free energy (hartree): 0.0015754684
# Partition function: 0.5369414248

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0031200409
# Entropy (hartree/K): 0.0000018437
# Free energy (hartree): 0.0016451099
# Partition function: 0.5223825760

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0032755913
# Entropy (hartree/K): 0.0000015155
# Free energy (hartree): 0.0020632176
# Partition function: 0.4429094832

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0034352730
# Entropy (hartree/K): 0.0000013084
# Free energy (hartree): 0.0023885754
# Partition function: 0.3895299519

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0034855494
# Entropy (hartree/K): 0.0000012376
# Free energy (hartree): 0.0024954514
# Partition function: 0.3734390583

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0036361753
# Entropy (hartree/K): 0.0000010753
# Free energy (hartree): 0.0027759322
# Partition function: 0.3343017059

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0037195920
# Entropy (hartree/K): 0.0000010060
# Free energy (hartree): 0.0029148264
# Partition function: 0.3164672750

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0037814615
# Entropy (hartree/K): 0.0000009332
# Free energy (hartree): 0.0030349349
# Partition function: 0.3018139604

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0038023613
# Entropy (hartree/K): 0.0000009163
# Free energy (hartree): 0.0030693094
# Partition function: 0.2977465342

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0039098484
# Entropy (hartree/K): 0.0000008384
# Free energy (hartree): 0.0032391633
# Partition function: 0.2784387481

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0039151511
# Entropy (hartree/K): 0.0000008350
# Free energy (hartree): 0.0032471874
# Partition function: 0.2775582569

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0039190413
# Entropy (hartree/K): 0.0000008315
# Free energy (hartree): 0.0032538223
# Partition function: 0.2768323059

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0039393667
# Entropy (hartree/K): 0.0000008175
# Free energy (hartree): 0.0032853455
# Partition function: 0.2734090839

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0039458643
# Entropy (hartree/K): 0.0000008138
# Free energy (hartree): 0.0032948239
# Partition function: 0.2723880865

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0070028410
# Entropy (hartree/K): 0.0000000846
# Free energy (hartree): 0.0069351634
# Partition function: 0.0647361032

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0069810303
# Entropy (hartree/K): 0.0000000897
# Free energy (hartree): 0.0069092885
# Partition function: 0.0654006609

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069819204
# Entropy (hartree/K): 0.0000000939
# Free energy (hartree): 0.0069067662
# Partition function: 0.0654658057

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0071089180
# Entropy (hartree/K): 0.0000000773
# Free energy (hartree): 0.0070470479
# Partition function: 0.0619393814

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0071682699
# Entropy (hartree/K): 0.0000000747
# Free energy (hartree): 0.0071085157
# Partition function: 0.0604546639

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0071988683
# Entropy (hartree/K): 0.0000000729
# Free energy (hartree): 0.0071405673
# Partition function: 0.0596946469

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0071789514
# Entropy (hartree/K): 0.0000000746
# Free energy (hartree): 0.0071192316
# Partition function: 0.0601994940

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0071759033
# Entropy (hartree/K): 0.0000000766
# Free energy (hartree): 0.0071146001
# Partition function: 0.0603096473

# 	********** Final results **********


# Temperature (K): 800.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 2.3846478011
# Translational entropy (cal/mol/K): 42.2047745300
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.3846478011
# Rotational entropy (cal/mol/K): 24.2237384405
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 76.5673145055
# Internal (tor+vib) entropy (cal/mol/K): 23.6367697228
# Internal (tor+vib) Cv (cal/mol/K): 29.7128215281


# Total energy (kcal/mol): 81.3366101078
# Total enthalpy (kcal/mol): 82.9263753987
# Enthalpy H(800.000000 K)-H(0 K) (kcal/mol):  17.5072137601
# Total entropy (cal/mol/K): 90.0652826932
# Total Cv (cal/mol/K): 35.6744410310



# 	********** HOhf results **********


# Translational energy (kcal/mol): 2.3846478011
# Rotational energy (kcal/mol): 2.3846478011
# Vibrational energy (kcal/mol): 76.1607181951
# gas constant (RT): 1.5897652008
# Translational entropy (cal/mol/K): 42.2047745300
# Rotational entropy (cal/mol/K): 24.2237384405
# Vibrational entropy (cal/mol/K): 22.7571789325


# Total energy (kcal/mol): 80.9300137974
# Total enthalpy (kcal/mol): 82.5197789982
# Enthalpy H(800.000000 K)-H(0 K) (kcal/mol): 17.1737789982
# Total entropy (cal/mol/K): 89.1856919029
# Total Cv (cal/mol/K): 36.7453727869
# Thermodynamics for propane at 900 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0001746086 		 0.0000000035
# 
# 	 53 		 0.0001499986 		 0.0000000028
# 
# 	 54 		 0.0001285582 		 0.0000000037
# 
# 	 55 		 0.0001099337 		 0.0000000018
# 
# 	 56 		 0.0000938011 		 0.0000000034
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705497
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0031840555
# Entropy (hartree/K): 0.0000069339
# Free energy (hartree): -0.0030564599
# Partition function: 2.9223608648

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000556048 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480633
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0032407551
# Entropy (hartree/K): 0.0000063385
# Free energy (hartree): -0.0024639059
# Partition function: 2.3737901410

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000001 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0029349980
# Entropy (hartree/K): 0.0000047879
# Free energy (hartree): -0.0013741125
# Partition function: 1.6195075746

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0031837516
# Entropy (hartree/K): 0.0000026770
# Free energy (hartree): 0.0007744657
# Partition function: 0.7620602923

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0033195288
# Entropy (hartree/K): 0.0000023250
# Free energy (hartree): 0.0012269905
# Partition function: 0.6501821022

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0033564244
# Entropy (hartree/K): 0.0000022073
# Free energy (hartree): 0.0013698487
# Partition function: 0.6183960712

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0033773958
# Entropy (hartree/K): 0.0000021466
# Free energy (hartree): 0.0014454203
# Partition function: 0.6022146899

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0035156543
# Entropy (hartree/K): 0.0000017980
# Free energy (hartree): 0.0018974098
# Partition function: 0.5138999756

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0036645032
# Entropy (hartree/K): 0.0000015782
# Free energy (hartree): 0.0022441524
# Partition function: 0.4550331540

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0037095562
# Entropy (hartree/K): 0.0000015013
# Free energy (hartree): 0.0023584228
# Partition function: 0.4371503872

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0038480537
# Entropy (hartree/K): 0.0000013246
# Free energy (hartree): 0.0026558814
# Partition function: 0.3938266068

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0039261889
# Entropy (hartree/K): 0.0000012491
# Free energy (hartree): 0.0028020371
# Partition function: 0.3741401387

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0039803786
# Entropy (hartree/K): 0.0000011672
# Free energy (hartree): 0.0029298902
# Partition function: 0.3577275970

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0039996505
# Entropy (hartree/K): 0.0000011484
# Free energy (hartree): 0.0029660484
# Partition function: 0.3532179535

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0040996512
# Entropy (hartree/K): 0.0000010617
# Free energy (hartree): 0.0031441573
# Partition function: 0.3318204134

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0041046308
# Entropy (hartree/K): 0.0000010579
# Free energy (hartree): 0.0031525414
# Partition function: 0.3308457390

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0041081319
# Entropy (hartree/K): 0.0000010540
# Free energy (hartree): 0.0031595430
# Partition function: 0.3300339883

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0041269875
# Entropy (hartree/K): 0.0000010383
# Free energy (hartree): 0.0031925557
# Partition function: 0.3262333092

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0041331447
# Entropy (hartree/K): 0.0000010341
# Free energy (hartree): 0.0032024278
# Partition function: 0.3251052748

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0070514138
# Entropy (hartree/K): 0.0000001416
# Free energy (hartree): 0.0069239783
# Partition function: 0.0880938770

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0070316357
# Entropy (hartree/K): 0.0000001491
# Free energy (hartree): 0.0068974779
# Partition function: 0.0889167902

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0070342094
# Entropy (hartree/K): 0.0000001553
# Free energy (hartree): 0.0068944317
# Partition function: 0.0890118737

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0071545433
# Entropy (hartree/K): 0.0000001309
# Free energy (hartree): 0.0070367592
# Partition function: 0.0846760377

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0072128161
# Entropy (hartree/K): 0.0000001270
# Free energy (hartree): 0.0070985540
# Partition function: 0.0828599058

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0072426523
# Entropy (hartree/K): 0.0000001242
# Free energy (hartree): 0.0071308313
# Partition function: 0.0819268244

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0072234900
# Entropy (hartree/K): 0.0000001269
# Free energy (hartree): 0.0071092747
# Partition function: 0.0825488134

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0072212829
# Entropy (hartree/K): 0.0000001299
# Free energy (hartree): 0.0071043969
# Partition function: 0.0826902111

# 	********** Final results **********


# Temperature (K): 900.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 2.6827287763
# Translational entropy (cal/mol/K): 42.7899225654
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.6827287763
# Rotational entropy (cal/mol/K): 24.5748272617
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 79.6846009897
# Internal (tor+vib) entropy (cal/mol/K): 27.3053559830
# Internal (tor+vib) Cv (cal/mol/K): 32.5388731216


# Total energy (kcal/mol): 85.0500585423
# Total enthalpy (kcal/mol): 86.8385444945
# Enthalpy H(900.000000 K)-H(0 K) (kcal/mol):  21.4193828560
# Total entropy (cal/mol/K): 94.6701058101
# Total Cv (cal/mol/K): 38.5004926244



# 	********** HOhf results **********


# Translational energy (kcal/mol): 2.6827287763
# Rotational energy (kcal/mol): 2.6827287763
# Vibrational energy (kcal/mol): 79.1800954627
# gas constant (RT): 1.7884858509
# Translational entropy (cal/mol/K): 42.7899225654
# Rotational entropy (cal/mol/K): 24.5748272617
# Vibrational entropy (cal/mol/K): 26.3103329630


# Total energy (kcal/mol): 84.5455530153
# Total enthalpy (kcal/mol): 86.3340388662
# Enthalpy H(900.000000 K)-H(0 K) (kcal/mol): 20.9880388662
# Total entropy (cal/mol/K): 93.6750827901
# Total Cv (cal/mol/K): 39.4789631898
# Thermodynamics for propane at 1000 K:


# 	********** Mode 1 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0003559905 		 0.0000000035
# 
# 	 53 		 0.0003101955 		 0.0000000028
# 
# 	 54 		 0.0002697308 		 0.0000000037
# 
# 	 55 		 0.0002340715 		 0.0000000018
# 
# 	 56 		 0.0002027280 		 0.0000000034
# 
# 	 57 		 0.0001752466 		 0.0000000011
# 
# 	 58 		 0.0001512092 		 0.0000000027
# 
# 	 59 		 0.0001302328 		 0.0000000006
# 
# 	 60 		 0.0001119687 		 0.0000000019
# 
# 	 61 		 0.0000961007 		 0.0000000004
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 223.7365705472
# Zero point vibrational energy (hartree): 0.0005244858
# Energy (hartree): 0.0035829317
# Entropy (hartree/K): 0.0000073546
# Free energy (hartree): -0.0037716277
# Partition function: 3.2903190205

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0001250498 		 0.0000000017
# 
# 	 53 		 0.0001057174 		 0.0000000018
# 
# 	 54 		 0.0000891523 		 0.0000000017
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 277.0315480634
# Zero point vibrational energy (hartree): 0.0006507202
# Energy (hartree): 0.0036421443
# Entropy (hartree/K): 0.0000067615
# Free energy (hartree): -0.0031193789
# Partition function: 2.6778677361

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000006 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.0156819052
# Zero point vibrational energy (hartree): 0.0008674121
# Energy (hartree): 0.0032413610
# Entropy (hartree/K): 0.0000051107
# Free energy (hartree): -0.0018693115
# Partition function: 1.8044958456

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 771.3642446129
# Zero point vibrational energy (hartree): 0.0017511240
# Energy (hartree): 0.0034625075
# Entropy (hartree/K): 0.0000029706
# Free energy (hartree): 0.0004918810
# Partition function: 0.8561381224

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 884.9873684896
# Zero point vibrational energy (hartree): 0.0020216856
# Energy (hartree): 0.0035939794
# Entropy (hartree/K): 0.0000026141
# Free energy (hartree): 0.0009798537
# Partition function: 0.7338776252

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 924.1996576628
# Zero point vibrational energy (hartree): 0.0021090179
# Energy (hartree): 0.0036267177
# Entropy (hartree/K): 0.0000024920
# Free energy (hartree): 0.0011347129
# Partition function: 0.6988538777

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 945.4852312923
# Zero point vibrational energy (hartree): 0.0021563557
# Energy (hartree): 0.0036456445
# Entropy (hartree/K): 0.0000024292
# Free energy (hartree): 0.0012164641
# Partition function: 0.6810438621

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.0720224070
# Zero point vibrational energy (hartree): 0.0024545188
# Energy (hartree): 0.0037686484
# Entropy (hartree/K): 0.0000020645
# Free energy (hartree): 0.0017041471
# Partition function: 0.5838410619

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1186.5113691648
# Zero point vibrational energy (hartree): 0.0027080987
# Energy (hartree): 0.0039088934
# Entropy (hartree/K): 0.0000018355
# Free energy (hartree): 0.0020733582
# Partition function: 0.5195906401

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1224.9600356474
# Zero point vibrational energy (hartree): 0.0027920968
# Energy (hartree): 0.0039492371
# Entropy (hartree/K): 0.0000017537
# Free energy (hartree): 0.0021955762
# Partition function: 0.4999199090

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.5670345652
# Zero point vibrational energy (hartree): 0.0030215483
# Energy (hartree): 0.0040772108
# Entropy (hartree/K): 0.0000015659
# Free energy (hartree): 0.0025112753
# Partition function: 0.4524865962

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1372.8492394814
# Zero point vibrational energy (hartree): 0.0031393346
# Energy (hartree): 0.0041509780
# Entropy (hartree/K): 0.0000014857
# Free energy (hartree): 0.0026652311
# Partition function: 0.4310149550

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1423.5021379763
# Zero point vibrational energy (hartree): 0.0032386304
# Energy (hartree): 0.0041978149
# Entropy (hartree/K): 0.0000013961
# Free energy (hartree): 0.0028016660
# Partition function: 0.4128399948

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.2766360607
# Zero point vibrational energy (hartree): 0.0032681784
# Energy (hartree): 0.0042156150
# Entropy (hartree/K): 0.0000013758
# Free energy (hartree): 0.0028397808
# Partition function: 0.4079009688

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1499.5371239225
# Zero point vibrational energy (hartree): 0.0034159801
# Energy (hartree): 0.0043090408
# Entropy (hartree/K): 0.0000012821
# Free energy (hartree): 0.0030269291
# Partition function: 0.3844938267

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4630848767
# Zero point vibrational energy (hartree): 0.0034230533
# Energy (hartree): 0.0043137477
# Entropy (hartree/K): 0.0000012780
# Free energy (hartree): 0.0030357066
# Partition function: 0.3834295902

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.2898594367
# Zero point vibrational energy (hartree): 0.0034287531
# Energy (hartree): 0.0043168815
# Entropy (hartree/K): 0.0000012738
# Free energy (hartree): 0.0030431170
# Partition function: 0.3825334072

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1517.2514196842
# Zero point vibrational energy (hartree): 0.0034564266
# Energy (hartree): 0.0043344120
# Entropy (hartree/K): 0.0000012566
# Free energy (hartree): 0.0030777750
# Partition function: 0.3783697378

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1520.5875415127
# Zero point vibrational energy (hartree): 0.0034648654
# Energy (hartree): 0.0043402948
# Entropy (hartree/K): 0.0000012522
# Free energy (hartree): 0.0030880753
# Partition function: 0.3771410693

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000039
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7700682024
# Zero point vibrational energy (hartree): 0.0069455625
# Energy (hartree): 0.0071189056
# Entropy (hartree/K): 0.0000002125
# Free energy (hartree): 0.0069063792
# Partition function: 0.1129443959

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8963255471
# Zero point vibrational energy (hartree): 0.0069204253
# Energy (hartree): 0.0071015547
# Entropy (hartree/K): 0.0000002225
# Free energy (hartree): 0.0068790050
# Partition function: 0.1139249286

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.7536676315
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0071061389
# Entropy (hartree/K): 0.0000002309
# Free energy (hartree): 0.0068752293
# Partition function: 0.1140608375

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000196
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.2199275173
# Zero point vibrational energy (hartree): 0.0070564059
# Energy (hartree): 0.0072185328
# Entropy (hartree/K): 0.0000001981
# Free energy (hartree): 0.0070204157
# Partition function: 0.1089496411

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.0124723804
# Zero point vibrational energy (hartree): 0.0071174979
# Energy (hartree): 0.0072755414
# Entropy (hartree/K): 0.0000001929
# Free energy (hartree): 0.0070826682
# Partition function: 0.1068288489

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.0289415956
# Zero point vibrational energy (hartree): 0.0071492938
# Energy (hartree): 0.0073044643
# Entropy (hartree/K): 0.0000001892
# Free energy (hartree): 0.0071152645
# Partition function: 0.1057348908

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000056
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4685873933
# Zero point vibrational energy (hartree): 0.0071282071
# Energy (hartree): 0.0072862200
# Entropy (hartree/K): 0.0000001928
# Free energy (hartree): 0.0070933941
# Partition function: 0.1064676350

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3121.5716242493
# Zero point vibrational energy (hartree): 0.0071238545
# Energy (hartree): 0.0072850403
# Entropy (hartree/K): 0.0000001969
# Free energy (hartree): 0.0070881661
# Partition function: 0.1066435459

# 	********** Final results **********


# Temperature (K): 1000.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.4191616385
# Translational energy (kcal/mol): 2.9808097514
# Translational entropy (cal/mol/K): 43.3133553195
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.9808097514
# Rotational entropy (cal/mol/K): 24.8888869142
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 83.0662336267
# Internal (tor+vib) entropy (cal/mol/K): 30.8663809581
# Internal (tor+vib) Cv (cal/mol/K): 34.9847307031


# Total energy (kcal/mol): 89.0278531296
# Total enthalpy (kcal/mol): 91.0150597432
# Enthalpy H(1000.000000 K)-H(0 K) (kcal/mol):  25.5958981047
# Total entropy (cal/mol/K): 99.0686231918
# Total Cv (cal/mol/K): 40.9463502060



# 	********** HOhf results **********


# Translational energy (kcal/mol): 2.9808097514
# Rotational energy (kcal/mol): 2.9808097514
# Vibrational energy (kcal/mol): 82.4553637596
# gas constant (RT): 1.9872065010
# Translational entropy (cal/mol/K): 43.3133553195
# Rotational entropy (cal/mol/K): 24.8888869142
# Vibrational entropy (cal/mol/K): 29.7589599326


# Total energy (kcal/mol): 88.4169832625
# Total enthalpy (kcal/mol): 90.4041897634
# Enthalpy H(1000.000000 K)-H(0 K) (kcal/mol): 25.0581897634
# Total entropy (cal/mol/K): 97.9612021663
# Total Cv (cal/mol/K): 41.8707892271

