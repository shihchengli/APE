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
# 	 52 		 0.0000018967 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0006402857
# Entropy (hartree/K): 0.0000002675
# Free energy (hartree): 0.0006135403
# Partition function: 0.1440773752

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000018967 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0006402848
# Entropy (hartree/K): 0.0000002675
# Free energy (hartree): 0.0006135392
# Partition function: 0.1440778619

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0008764635
# Entropy (hartree/K): 0.0000000860
# Free energy (hartree): 0.0008678614
# Partition function: 0.0645389376

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0017406082
# Entropy (hartree/K): 0.0000000006
# Free energy (hartree): 0.0017405464
# Partition function: 0.0041023116

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0020246095
# Entropy (hartree/K): 0.0000000001
# Free energy (hartree): 0.0020245968
# Partition function: 0.0016729543

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0021179203
# Entropy (hartree/K): 0.0000000001
# Free energy (hartree): 0.0021179132
# Partition function: 0.0012459790

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0021841320
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0021841272
# Partition function: 0.0010108935

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0024551908
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0024551898
# Partition function: 0.0004295090

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0027283383
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0027283382
# Partition function: 0.0001812921

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0028419784
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0028419783
# Partition function: 0.0001266291

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0030180882
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0030180882
# Partition function: 0.0000726136

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0031566485
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0031566485
# Partition function: 0.0000468811

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0032406334
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0032406334
# Partition function: 0.0000359601

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0032683420
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0032683420
# Partition function: 0.0000329475

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0034203834
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034203834
# Partition function: 0.0000203851

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0034315030
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034315030
# Partition function: 0.0000196818

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0034606954
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034606954
# Partition function: 0.0000179486

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0034964439
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0034964439
# Partition function: 0.0000160326

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0035487164
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0035487164
# Partition function: 0.0000135931

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069455467
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455467
# Partition function: 0.0000000003

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069204270
# Entropy (hartree/K): -0.0000000000
# Free energy (hartree): 0.0069204270
# Partition function: 0.0000000003

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185284
# Entropy (hartree/K): -0.0000000000
# Free energy (hartree): 0.0069185284
# Partition function: 0.0000000003

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070574601
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070574601
# Partition function: 0.0000000002

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071174023
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174023
# Partition function: 0.0000000002

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071495835
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071495835
# Partition function: 0.0000000002

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071282627
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282627
# Partition function: 0.0000000002

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071245410
# Entropy (hartree/K): -0.0000000000
# Free energy (hartree): 0.0071245410
# Partition function: 0.0000000002

# 	********** Final results **********


# Temperature (K): 100.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 0.2980809751
# Translational entropy (cal/mol/K): 31.8740751550
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.2980809751
# Rotational entropy (cal/mol/K): 18.0253188155
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 65.6707653949
# Internal (tor+vib) entropy (cal/mol/K): 0.3901899623
# Internal (tor+vib) Cv (cal/mol/K): 2.1613476950


# Total energy (kcal/mol): 66.2669273452
# Total enthalpy (kcal/mol): 66.4656480065
# Enthalpy H(100.000000 K)-H(0 K) (kcal/mol):  0.9651836513
# Total entropy (cal/mol/K): 50.2895839328
# Total Cv (cal/mol/K): 8.1229671979



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
# 	 52 		 0.0000021230 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0008676567
# Entropy (hartree/K): 0.0000017892
# Free energy (hartree): 0.0005098199
# Partition function: 0.4471142905

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000021230 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0008676553
# Entropy (hartree/K): 0.0000017892
# Free energy (hartree): 0.0005098189
# Partition function: 0.4471150083

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0009890720
# Entropy (hartree/K): 0.0000008110
# Free energy (hartree): 0.0008268672
# Partition function: 0.2710319551

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0017546793
# Entropy (hartree/K): 0.0000000834
# Free energy (hartree): 0.0017379942
# Partition function: 0.0643079039

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0020314936
# Entropy (hartree/K): 0.0000000399
# Free energy (hartree): 0.0020235169
# Partition function: 0.0409715663

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0021231859
# Entropy (hartree/K): 0.0000000303
# Free energy (hartree): 0.0021171262
# Partition function: 0.0353423180

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0021885293
# Entropy (hartree/K): 0.0000000252
# Free energy (hartree): 0.0021834901
# Partition function: 0.0318265524

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0024572643
# Entropy (hartree/K): 0.0000000117
# Free energy (hartree): 0.0024549231
# Partition function: 0.0207333311

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0027293227
# Entropy (hartree/K): 0.0000000055
# Free energy (hartree): 0.0027282240
# Partition function: 0.0134669032

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0028427062
# Entropy (hartree/K): 0.0000000040
# Free energy (hartree): 0.0028418971
# Partition function: 0.0112544045

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0030185274
# Entropy (hartree/K): 0.0000000024
# Free energy (hartree): 0.0030180421
# Partition function: 0.0085219825

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0031569410
# Entropy (hartree/K): 0.0000000016
# Free energy (hartree): 0.0031566192
# Partition function: 0.0068472922

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0032408632
# Entropy (hartree/K): 0.0000000013
# Free energy (hartree): 0.0032406110
# Partition function: 0.0059968887

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0032685549
# Entropy (hartree/K): 0.0000000012
# Free energy (hartree): 0.0032683214
# Partition function: 0.0057401748

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0034205219
# Entropy (hartree/K): 0.0000000008
# Free energy (hartree): 0.0034203706
# Partition function: 0.0045150812

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0034316388
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0034314904
# Partition function: 0.0044365026

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0034608195
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0034606841
# Partition function: 0.0042366520

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0034965568
# Entropy (hartree/K): 0.0000000006
# Free energy (hartree): 0.0034964337
# Partition function: 0.0040041413

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0035488135
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0035487078
# Partition function: 0.0036869334

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069455467
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455467
# Partition function: 0.0000172769

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069204270
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204270
# Partition function: 0.0000179759

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185284
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185284
# Partition function: 0.0000180299

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070574601
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070574601
# Partition function: 0.0000144786

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071174023
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174023
# Partition function: 0.0000131712

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071495835
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071495835
# Partition function: 0.0000125187

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071282627
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282627
# Partition function: 0.0000129473

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071245410
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071245410
# Partition function: 0.0000130236

# 	********** Final results **********


# Temperature (K): 200.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 0.5961619503
# Translational entropy (cal/mol/K): 35.3176416133
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.5961619503
# Rotational entropy (cal/mol/K): 20.0914586905
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 66.0494902770
# Internal (tor+vib) entropy (cal/mol/K): 2.8860566122
# Internal (tor+vib) Cv (cal/mol/K): 5.4384671179


# Total energy (kcal/mol): 67.2418141776
# Total enthalpy (kcal/mol): 67.6392555003
# Enthalpy H(200.000000 K)-H(0 K) (kcal/mol):  2.1387911451
# Total entropy (cal/mol/K): 58.2951569161
# Total Cv (cal/mol/K): 11.4000866207



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
# 	 52 		 0.0000018812 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0011800354
# Entropy (hartree/K): 0.0000030527
# Free energy (hartree): 0.0002698712
# Partition function: 0.7513946914

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000018812 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0011800340
# Entropy (hartree/K): 0.0000030527
# Free energy (hartree): 0.0002698704
# Partition function: 0.7513953496

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0011977935
# Entropy (hartree/K): 0.0000016515
# Free energy (hartree): 0.0007054064
# Partition function: 0.4737358580

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0018290923
# Entropy (hartree/K): 0.0000003763
# Free energy (hartree): 0.0017169005
# Partition function: 0.1622864240

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0020815303
# Entropy (hartree/K): 0.0000002353
# Free energy (hartree): 0.0020113779
# Partition function: 0.1188043318

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0021660896
# Entropy (hartree/K): 0.0000001974
# Free energy (hartree): 0.0021072392
# Partition function: 0.1073344766

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0022271886
# Entropy (hartree/K): 0.0000001755
# Free energy (hartree): 0.0021748710
# Partition function: 0.0999150352

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0024821163
# Entropy (hartree/K): 0.0000001076
# Free energy (hartree): 0.0024500388
# Partition function: 0.0746556485

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0027452038
# Entropy (hartree/K): 0.0000000663
# Free energy (hartree): 0.0027254262
# Partition function: 0.0557690875

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0028559200
# Entropy (hartree/K): 0.0000000545
# Free energy (hartree): 0.0028396621
# Partition function: 0.0494138599

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0030282185
# Entropy (hartree/K): 0.0000000393
# Free energy (hartree): 0.0030165045
# Partition function: 0.0409738790

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0031644763
# Entropy (hartree/K): 0.0000000302
# Free energy (hartree): 0.0031554800
# Partition function: 0.0353657532

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0032473468
# Entropy (hartree/K): 0.0000000258
# Free energy (hartree): 0.0032396574
# Partition function: 0.0323492387

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0032747393
# Entropy (hartree/K): 0.0000000246
# Free energy (hartree): 0.0032674195
# Partition function: 0.0314119171

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0034252465
# Entropy (hartree/K): 0.0000000186
# Free energy (hartree): 0.0034197129
# Partition function: 0.0267328013

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0034363075
# Entropy (hartree/K): 0.0000000183
# Free energy (hartree): 0.0034308419
# Partition function: 0.0264195564

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0034652311
# Entropy (hartree/K): 0.0000000173
# Free energy (hartree): 0.0034600769
# Partition function: 0.0256140556

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0035007122
# Entropy (hartree/K): 0.0000000162
# Free energy (hartree): 0.0034958674
# Partition function: 0.0246612956

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0035525930
# Entropy (hartree/K): 0.0000000147
# Free energy (hartree): 0.0035482005
# Partition function: 0.0233315924

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069455522
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455463
# Partition function: 0.0006386707

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069204334
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204265
# Partition function: 0.0006558904

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185358
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185279
# Partition function: 0.0006572107

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070574643
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070574598
# Partition function: 0.0005672839

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071174061
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174021
# Partition function: 0.0005323889

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071495870
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071495833
# Partition function: 0.0005145490

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071282665
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282624
# Partition function: 0.0005263003

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071245451
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071245407
# Partition function: 0.0005283789

# 	********** Final results **********


# Temperature (K): 298.15
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 0.8887284274
# Translational entropy (cal/mol/K): 37.3012679085
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.8887284274
# Rotational entropy (cal/mol/K): 21.2816344676
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 66.7680585979
# Internal (tor+vib) entropy (cal/mol/K): 5.7573268912
# Internal (tor+vib) Cv (cal/mol/K): 9.3593527273


# Total energy (kcal/mol): 68.5455154527
# Total enthalpy (kcal/mol): 69.1380011045
# Enthalpy H(298.150000 K)-H(0 K) (kcal/mol):  3.6375367493
# Total entropy (cal/mol/K): 64.3402292673
# Total Cv (cal/mol/K): 15.3209722302



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
# 	 52 		 0.0000018767 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0011863406
# Entropy (hartree/K): 0.0000030738
# Free energy (hartree): 0.0002642042
# Partition function: 0.7572236542

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000018767 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0011863392
# Entropy (hartree/K): 0.0000030738
# Free energy (hartree): 0.0002642033
# Partition function: 0.7572243104

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0012022562
# Entropy (hartree/K): 0.0000016664
# Free energy (hartree): 0.0007023373
# Partition function: 0.4774633814

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0018311863
# Entropy (hartree/K): 0.0000003833
# Free energy (hartree): 0.0017161979
# Partition function: 0.1642378696

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0020830739
# Entropy (hartree/K): 0.0000002405
# Free energy (hartree): 0.0020109379
# Partition function: 0.1204310967

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0021674572
# Entropy (hartree/K): 0.0000002020
# Free energy (hartree): 0.0021068698
# Partition function: 0.1088642308

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0022284481
# Entropy (hartree/K): 0.0000001797
# Free energy (hartree): 0.0021745424
# Partition function: 0.1013794576

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0024830009
# Entropy (hartree/K): 0.0000001105
# Free energy (hartree): 0.0024498371
# Partition function: 0.0758759871

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0027458180
# Entropy (hartree/K): 0.0000000684
# Free energy (hartree): 0.0027253016
# Partition function: 0.0567781289

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0028564479
# Entropy (hartree/K): 0.0000000563
# Free energy (hartree): 0.0028395596
# Partition function: 0.0503442920

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0030286264
# Entropy (hartree/K): 0.0000000407
# Free energy (hartree): 0.0030164306
# Partition function: 0.0417923784

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0031648067
# Entropy (hartree/K): 0.0000000313
# Free energy (hartree): 0.0031554232
# Partition function: 0.0361043308

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0032476379
# Entropy (hartree/K): 0.0000000268
# Free energy (hartree): 0.0032396088
# Partition function: 0.0330426946

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0032750189
# Entropy (hartree/K): 0.0000000255
# Free energy (hartree): 0.0032673732
# Partition function: 0.0320910196

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0034254690
# Entropy (hartree/K): 0.0000000193
# Free energy (hartree): 0.0034196779
# Partition function: 0.0273375992

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0034365278
# Entropy (hartree/K): 0.0000000191
# Free energy (hartree): 0.0034308073
# Partition function: 0.0270192192

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0034654410
# Entropy (hartree/K): 0.0000000180
# Free energy (hartree): 0.0034600443
# Partition function: 0.0262003834

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0035009116
# Entropy (hartree/K): 0.0000000169
# Free energy (hartree): 0.0034958368
# Partition function: 0.0252316594

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0035527769
# Entropy (hartree/K): 0.0000000153
# Free energy (hartree): 0.0035481727
# Partition function: 0.0238792918

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069455527
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069455463
# Partition function: 0.0006683097

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069204340
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069204265
# Partition function: 0.0006862160

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069185365
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0069185278
# Partition function: 0.0006875888

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070574647
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0070574598
# Partition function: 0.0005940441

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071174065
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071174020
# Partition function: 0.0005577213

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071495874
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071495832
# Partition function: 0.0005391458

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071282668
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071282624
# Partition function: 0.0005513821

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071245455
# Entropy (hartree/K): 0.0000000000
# Free energy (hartree): 0.0071245407
# Partition function: 0.0005535463

# 	********** Final results **********


# Temperature (K): 300.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 0.8942429254
# Translational entropy (cal/mol/K): 37.3319988602
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 0.8942429254
# Rotational entropy (cal/mol/K): 21.3000730386
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 66.7854489564
# Internal (tor+vib) entropy (cal/mol/K): 5.8154739648
# Internal (tor+vib) Cv (cal/mol/K): 9.4410658002


# Total energy (kcal/mol): 68.5739348073
# Total enthalpy (kcal/mol): 69.1700967914
# Enthalpy H(300.000000 K)-H(0 K) (kcal/mol):  3.6696324361
# Total entropy (cal/mol/K): 64.4475458636
# Total Cv (cal/mol/K): 15.4026853031



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
# 	 52 		 0.0000017406 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0015339190
# Entropy (hartree/K): 0.0000040732
# Free energy (hartree): -0.0000953674
# Partition function: 1.0781930196

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000017406 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0015339176
# Entropy (hartree/K): 0.0000040732
# Free energy (hartree): -0.0000953680
# Partition function: 1.0781935562

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0014605418
# Entropy (hartree/K): 0.0000024074
# Free energy (hartree): 0.0004975823
# Partition function: 0.6751571789

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0019767392
# Entropy (hartree/K): 0.0000007979
# Free energy (hartree): 0.0016575803
# Partition function: 0.2702103166

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0021986295
# Entropy (hartree/K): 0.0000005687
# Free energy (hartree): 0.0019711549
# Partition function: 0.2109567840

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0022726872
# Entropy (hartree/K): 0.0000005006
# Free energy (hartree): 0.0020724644
# Partition function: 0.1947420402

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0023272199
# Entropy (hartree/K): 0.0000004598
# Free energy (hartree): 0.0021433149
# Partition function: 0.1841487461

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0025582645
# Entropy (hartree/K): 0.0000003233
# Free energy (hartree): 0.0024289321
# Partition function: 0.1469756861

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0028026378
# Entropy (hartree/K): 0.0000002286
# Free energy (hartree): 0.0027112164
# Partition function: 0.1176155857

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0029069928
# Entropy (hartree/K): 0.0000001986
# Free energy (hartree): 0.0028275523
# Partition function: 0.1072949633

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0030700154
# Entropy (hartree/K): 0.0000001570
# Free energy (hartree): 0.0030072349
# Partition function: 0.0931055402

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0031999385
# Entropy (hartree/K): 0.0000001298
# Free energy (hartree): 0.0031480029
# Partition function: 0.0833131429

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0032794753
# Entropy (hartree/K): 0.0000001160
# Free energy (hartree): 0.0032330759
# Partition function: 0.0779016020

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0033058827
# Entropy (hartree/K): 0.0000001120
# Free energy (hartree): 0.0032610961
# Partition function: 0.0761973238

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0034513084
# Entropy (hartree/K): 0.0000000916
# Free energy (hartree): 0.0034146738
# Partition function: 0.0674972237

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0034621665
# Entropy (hartree/K): 0.0000000908
# Free energy (hartree): 0.0034258528
# Partition function: 0.0669041751

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0034901321
# Entropy (hartree/K): 0.0000000870
# Free energy (hartree): 0.0034553199
# Partition function: 0.0653657880

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0035246404
# Entropy (hartree/K): 0.0000000832
# Free energy (hartree): 0.0034913432
# Partition function: 0.0635330875

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0035750536
# Entropy (hartree/K): 0.0000000776
# Free energy (hartree): 0.0035440217
# Partition function: 0.0609451619

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069457805
# Entropy (hartree/K): 0.0000000006
# Free energy (hartree): 0.0069455254
# Partition function: 0.0041566217

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069206915
# Entropy (hartree/K): 0.0000000007
# Free energy (hartree): 0.0069204026
# Partition function: 0.0042398825

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069188204
# Entropy (hartree/K): 0.0000000008
# Free energy (hartree): 0.0069185012
# Partition function: 0.0042462514

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070576528
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0070574429
# Partition function: 0.0038051329

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071175815
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071173864
# Partition function: 0.0036292618

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071497535
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071495685
# Partition function: 0.0035382192

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071284416
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071282468
# Partition function: 0.0035982790

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071247298
# Entropy (hartree/K): 0.0000000005
# Free energy (hartree): 0.0071245242
# Partition function: 0.0036088691

# 	********** Final results **********


# Temperature (K): 400.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 1.1923239006
# Translational entropy (cal/mol/K): 38.7612080717
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 1.1923239006
# Rotational entropy (cal/mol/K): 22.1575985655
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 67.9552743285
# Internal (tor+vib) entropy (cal/mol/K): 9.1496099827
# Internal (tor+vib) Cv (cal/mol/K): 13.9634091457


# Total energy (kcal/mol): 70.3399221296
# Total enthalpy (kcal/mol): 71.1348047751
# Enthalpy H(400.000000 K)-H(0 K) (kcal/mol):  5.6343404198
# Total entropy (cal/mol/K): 70.0684166199
# Total Cv (cal/mol/K): 19.9250286485



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
# 	 52 		 0.0000024556 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0018754025
# Entropy (hartree/K): 0.0000048359
# Free energy (hartree): -0.0005425373
# Partition function: 1.4086602283

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000024556 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0018754013
# Entropy (hartree/K): 0.0000048359
# Free energy (hartree): -0.0005425377
# Partition function: 1.4086606424

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0017401911
# Entropy (hartree/K): 0.0000030308
# Free energy (hartree): 0.0002247933
# Partition function: 0.8676489715

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0021715183
# Entropy (hartree/K): 0.0000012309
# Free energy (hartree): 0.0015560724
# Partition function: 0.3742854508

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0023671072
# Entropy (hartree/K): 0.0000009428
# Free energy (hartree): 0.0018957223
# Partition function: 0.3020259935

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0024309009
# Entropy (hartree/K): 0.0000008517
# Free energy (hartree): 0.0020050493
# Partition function: 0.2818761170

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0024789305
# Entropy (hartree/K): 0.0000007964
# Free energy (hartree): 0.0020807424
# Partition function: 0.2687183304

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0026849547
# Entropy (hartree/K): 0.0000006041
# Free energy (hartree): 0.0023829203
# Partition function: 0.2220324443

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0029076034
# Entropy (hartree/K): 0.0000004608
# Free energy (hartree): 0.0026771823
# Partition function: 0.1843769894

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0030039872
# Entropy (hartree/K): 0.0000004131
# Free energy (hartree): 0.0027974213
# Partition function: 0.1708943731

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0031546054
# Entropy (hartree/K): 0.0000003439
# Free energy (hartree): 0.0029826675
# Partition function: 0.1520263016

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0032754703
# Entropy (hartree/K): 0.0000002966
# Free energy (hartree): 0.0031271580
# Partition function: 0.1387676095

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0033500288
# Entropy (hartree/K): 0.0000002717
# Free energy (hartree): 0.0032141677
# Partition function: 0.1313479286

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0033749317
# Entropy (hartree/K): 0.0000002643
# Free energy (hartree): 0.0032427575
# Partition function: 0.1289976062

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0035123272
# Entropy (hartree/K): 0.0000002261
# Free energy (hartree): 0.0033992573
# Partition function: 0.1168576524

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0035228687
# Entropy (hartree/K): 0.0000002246
# Free energy (hartree): 0.0034105516
# Partition function: 0.1160270809

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0035492473
# Entropy (hartree/K): 0.0000002174
# Free energy (hartree): 0.0034405673
# Partition function: 0.1138483393

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0035821320
# Entropy (hartree/K): 0.0000002100
# Free energy (hartree): 0.0034771467
# Partition function: 0.1112483967

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0036300535
# Entropy (hartree/K): 0.0000001988
# Free energy (hartree): 0.0035306629
# Partition function: 0.1075512410

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069476548
# Entropy (hartree/K): 0.0000000047
# Free energy (hartree): 0.0069453071
# Partition function: 0.0124466448

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069227487
# Entropy (hartree/K): 0.0000000052
# Free energy (hartree): 0.0069201597
# Partition function: 0.0126458978

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069210358
# Entropy (hartree/K): 0.0000000056
# Free energy (hartree): 0.0069182369
# Partition function: 0.0126612641

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070592732
# Entropy (hartree/K): 0.0000000040
# Free energy (hartree): 0.0070572579
# Partition function: 0.0115970254

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071191154
# Entropy (hartree/K): 0.0000000038
# Free energy (hartree): 0.0071172126
# Partition function: 0.0111661207

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071512272
# Entropy (hartree/K): 0.0000000036
# Free energy (hartree): 0.0071494024
# Partition function: 0.0109414115

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071299737
# Entropy (hartree/K): 0.0000000038
# Free energy (hartree): 0.0071280732
# Partition function: 0.0110897941

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071263252
# Entropy (hartree/K): 0.0000000040
# Free energy (hartree): 0.0071243424
# Partition function: 0.0111159545

# 	********** Final results **********


# Temperature (K): 500.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 1.4904048757
# Translational entropy (cal/mol/K): 39.8697888612
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 1.4904048757
# Rotational entropy (cal/mol/K): 22.8227470392
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 69.5688536209
# Internal (tor+vib) entropy (cal/mol/K): 12.7324926479
# Internal (tor+vib) Cv (cal/mol/K): 18.2291620721


# Total energy (kcal/mol): 72.5496633724
# Total enthalpy (kcal/mol): 73.5432666791
# Enthalpy H(500.000000 K)-H(0 K) (kcal/mol):  8.0428023239
# Total entropy (cal/mol/K): 75.4250285484
# Total Cv (cal/mol/K): 24.1907815749



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
# 	 52 		 0.0000069393 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0021965880
# Entropy (hartree/K): 0.0000054221
# Free energy (hartree): -0.0010566760
# Partition function: 1.7438916716

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000069393 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0021965869
# Entropy (hartree/K): 0.0000054221
# Free energy (hartree): -0.0010566763
# Partition function: 1.7438919755

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0020311947
# Entropy (hartree/K): 0.0000035611
# Free energy (hartree): -0.0001054802
# Partition function: 1.0570830584

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0023984939
# Entropy (hartree/K): 0.0000016440
# Free energy (hartree): 0.0014120923
# Partition function: 0.4756030839

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0025733292
# Entropy (hartree/K): 0.0000013179
# Free energy (hartree): 0.0017825951
# Partition function: 0.3913457024

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0026278811
# Entropy (hartree/K): 0.0000012099
# Free energy (hartree): 0.0019019154
# Partition function: 0.3675259993

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0026701162
# Entropy (hartree/K): 0.0000011440
# Free energy (hartree): 0.0019836991
# Partition function: 0.3520425500

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0028533566
# Entropy (hartree/K): 0.0000009101
# Free energy (hartree): 0.0023072977
# Partition function: 0.2969147189

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0030551186
# Entropy (hartree/K): 0.0000007287
# Free energy (hartree): 0.0026178769
# Partition function: 0.2521413662

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0031434942
# Entropy (hartree/K): 0.0000006664
# Free energy (hartree): 0.0027436471
# Partition function: 0.2359920520

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0032810300
# Entropy (hartree/K): 0.0000005733
# Free energy (hartree): 0.0029370512
# Partition function: 0.2131532238

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0033919162
# Entropy (hartree/K): 0.0000005079
# Free energy (hartree): 0.0030872007
# Partition function: 0.1969576872

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0034608645
# Entropy (hartree/K): 0.0000004727
# Free energy (hartree): 0.0031772243
# Partition function: 0.1878437170

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0034840482
# Entropy (hartree/K): 0.0000004622
# Free energy (hartree): 0.0032067117
# Partition function: 0.1849510847

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0036120835
# Entropy (hartree/K): 0.0000004070
# Free energy (hartree): 0.0033679018
# Partition function: 0.1699082240

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0036222849
# Entropy (hartree/K): 0.0000004048
# Free energy (hartree): 0.0033793790
# Partition function: 0.1688850110

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0036467383
# Entropy (hartree/K): 0.0000003941
# Free energy (hartree): 0.0034102998
# Partition function: 0.1661589329

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0036776602
# Entropy (hartree/K): 0.0000003831
# Free energy (hartree): 0.0034477995
# Partition function: 0.1629118106

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0037225369
# Entropy (hartree/K): 0.0000003664
# Free energy (hartree): 0.0035027155
# Partition function: 0.1582707544

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069546829
# Entropy (hartree/K): 0.0000000173
# Free energy (hartree): 0.0069443007
# Partition function: 0.0258686407

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069303092
# Entropy (hartree/K): 0.0000000187
# Free energy (hartree): 0.0069190624
# Partition function: 0.0262145378

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069290486
# Entropy (hartree/K): 0.0000000200
# Free energy (hartree): 0.0069170611
# Partition function: 0.0262421623

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070655424
# Entropy (hartree/K): 0.0000000153
# Free energy (hartree): 0.0070563786
# Partition function: 0.0243868935

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071251205
# Entropy (hartree/K): 0.0000000146
# Free energy (hartree): 0.0071163768
# Partition function: 0.0236288703

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071570465
# Entropy (hartree/K): 0.0000000141
# Free energy (hartree): 0.0071485970
# Partition function: 0.0232315696

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071359736
# Entropy (hartree/K): 0.0000000146
# Free energy (hartree): 0.0071272383
# Partition function: 0.0234941858

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071325191
# Entropy (hartree/K): 0.0000000151
# Free energy (hartree): 0.0071234757
# Partition function: 0.0235407564

# 	********** Final results **********


# Temperature (K): 600.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 1.7884858509
# Translational entropy (cal/mol/K): 40.7755653185
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 1.7884858509
# Rotational entropy (cal/mol/K): 23.3662129136
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 71.5835039004
# Internal (tor+vib) entropy (cal/mol/K): 16.3952552904
# Internal (tor+vib) Cv (cal/mol/K): 21.9757349522


# Total energy (kcal/mol): 75.1604756022
# Total enthalpy (kcal/mol): 76.3527995703
# Enthalpy H(600.000000 K)-H(0 K) (kcal/mol):  10.8523352151
# Total entropy (cal/mol/K): 80.5370335225
# Total Cv (cal/mol/K): 27.9373544550



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
# 	 52 		 0.0000226781 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0024945936
# Entropy (hartree/K): 0.0000058819
# Free energy (hartree): -0.0016227618
# Partition function: 2.0793153997

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000226781 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0024945927
# Entropy (hartree/K): 0.0000058819
# Free energy (hartree): -0.0016227621
# Partition function: 2.0793156087

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0023287913
# Entropy (hartree/K): 0.0000040198
# Free energy (hartree): -0.0004850487
# Partition function: 1.2445931702

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0026466546
# Entropy (hartree/K): 0.0000020262
# Free energy (hartree): 0.0012283121
# Partition function: 0.5745892735

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0028058429
# Entropy (hartree/K): 0.0000016759
# Free energy (hartree): 0.0016327280
# Partition function: 0.4787702115

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0028522053
# Entropy (hartree/K): 0.0000015553
# Free energy (hartree): 0.0017635039
# Partition function: 0.4513426909

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0028894441
# Entropy (hartree/K): 0.0000014817
# Free energy (hartree): 0.0018522844
# Partition function: 0.4336238266

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0030530688
# Entropy (hartree/K): 0.0000012174
# Free energy (hartree): 0.0022008711
# Partition function: 0.3705274801

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0032363610
# Entropy (hartree/K): 0.0000010075
# Free energy (hartree): 0.0025310840
# Partition function: 0.3192474027

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0033174699
# Entropy (hartree/K): 0.0000009340
# Free energy (hartree): 0.0026636740
# Partition function: 0.3007123428

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0034426175
# Entropy (hartree/K): 0.0000008218
# Free energy (hartree): 0.0028673835
# Partition function: 0.2743101620

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0035437487
# Entropy (hartree/K): 0.0000007413
# Free energy (hartree): 0.0030248557
# Partition function: 0.2555000641

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0036071618
# Entropy (hartree/K): 0.0000006976
# Free energy (hartree): 0.0031188335
# Partition function: 0.2448947847

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0036286338
# Entropy (hartree/K): 0.0000006845
# Free energy (hartree): 0.0031495078
# Partition function: 0.2415294117

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0037472501
# Entropy (hartree/K): 0.0000006147
# Free energy (hartree): 0.0033169735
# Partition function: 0.2239552802

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0037571676
# Entropy (hartree/K): 0.0000006121
# Free energy (hartree): 0.0033286869
# Partition function: 0.2227750143

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0037796036
# Entropy (hartree/K): 0.0000005982
# Free energy (hartree): 0.0033608450
# Partition function: 0.2195666119

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0038084947
# Entropy (hartree/K): 0.0000005841
# Free energy (hartree): 0.0033996018
# Partition function: 0.2157611864

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0038502035
# Entropy (hartree/K): 0.0000005625
# Free energy (hartree): 0.0034564419
# Partition function: 0.2102991884

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0069716094
# Entropy (hartree/K): 0.0000000432
# Free energy (hartree): 0.0069414020
# Partition function: 0.0436604405

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069482615
# Entropy (hartree/K): 0.0000000462
# Free energy (hartree): 0.0069159459
# Partition function: 0.0441647020

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069478620
# Entropy (hartree/K): 0.0000000487
# Free energy (hartree): 0.0069137592
# Partition function: 0.0442082891

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0070809822
# Entropy (hartree/K): 0.0000000388
# Free energy (hartree): 0.0070537897
# Partition function: 0.0415020837

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071400361
# Entropy (hartree/K): 0.0000000373
# Free energy (hartree): 0.0071138955
# Partition function: 0.0403919092

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071715894
# Entropy (hartree/K): 0.0000000363
# Free energy (hartree): 0.0071461913
# Partition function: 0.0398077119

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071508797
# Entropy (hartree/K): 0.0000000373
# Free energy (hartree): 0.0071247591
# Partition function: 0.0401944474

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071478136
# Entropy (hartree/K): 0.0000000384
# Free energy (hartree): 0.0071209173
# Partition function: 0.0402641668

# 	********** Final results **********


# Temperature (K): 700.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 2.0865668260
# Translational entropy (cal/mol/K): 41.5413884012
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.0865668260
# Rotational entropy (cal/mol/K): 23.8257067632
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 73.9475671244
# Internal (tor+vib) entropy (cal/mol/K): 20.0330308873
# Internal (tor+vib) Cv (cal/mol/K): 25.2302811314


# Total energy (kcal/mol): 78.1207007764
# Total enthalpy (kcal/mol): 79.5117454059
# Enthalpy H(700.000000 K)-H(0 K) (kcal/mol):  14.0112810507
# Total entropy (cal/mol/K): 85.4001260518
# Total Cv (cal/mol/K): 31.1919006343



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
# 	 52 		 0.0000633745 		 0.0019683777
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8224049998
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0027709136
# Entropy (hartree/K): 0.0000062512
# Free energy (hartree): -0.0022300576
# Partition function: 2.4114895565

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000633746 		 0.0019683791
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8214966544
# Zero point vibrational energy (hartree): 0.0005082679
# Energy (hartree): 0.0027709128
# Entropy (hartree/K): 0.0000062512
# Free energy (hartree): -0.0022300577
# Partition function: 2.4114896853

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0026304673
# Entropy (hartree/K): 0.0000044226
# Free energy (hartree): -0.0009075772
# Partition function: 1.4308052006

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0029091419
# Entropy (hartree/K): 0.0000023765
# Free energy (hartree): 0.0010079172
# Partition function: 0.6717672070

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0030569235
# Entropy (hartree/K): 0.0000020109
# Free energy (hartree): 0.0014481896
# Partition function: 0.5646060829

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0030959326
# Entropy (hartree/K): 0.0000018805
# Free energy (hartree): 0.0015915348
# Partition function: 0.5335470636

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0031288616
# Entropy (hartree/K): 0.0000018011
# Free energy (hartree): 0.0016879799
# Partition function: 0.5136174012

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0032759001
# Entropy (hartree/K): 0.0000015147
# Free energy (hartree): 0.0020641555
# Partition function: 0.4427455573

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0034434332
# Entropy (hartree/K): 0.0000012837
# Free energy (hartree): 0.0024164643
# Partition function: 0.3852654215

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0035182494
# Entropy (hartree/K): 0.0000012017
# Free energy (hartree): 0.0025568500
# Partition function: 0.3644974965

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0036322088
# Entropy (hartree/K): 0.0000010746
# Free energy (hartree): 0.0027725615
# Partition function: 0.3347467928

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0037242958
# Entropy (hartree/K): 0.0000009820
# Free energy (hartree): 0.0029387094
# Partition function: 0.3134979277

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0037825816
# Entropy (hartree/K): 0.0000009315
# Free energy (hartree): 0.0030374094
# Partition function: 0.3015193170

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0038024576
# Entropy (hartree/K): 0.0000009162
# Free energy (hartree): 0.0030695090
# Partition function: 0.2977230738

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0039122562
# Entropy (hartree/K): 0.0000008346
# Free energy (hartree): 0.0032445642
# Partition function: 0.2778457951

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0039219947
# Entropy (hartree/K): 0.0000008318
# Free energy (hartree): 0.0032565479
# Partition function: 0.2765346368

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0039424546
# Entropy (hartree/K): 0.0000008153
# Free energy (hartree): 0.0032902306
# Partition function: 0.2728823911

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0039693966
# Entropy (hartree/K): 0.0000007986
# Free energy (hartree): 0.0033305314
# Partition function: 0.2685758720

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0040080640
# Entropy (hartree/K): 0.0000007729
# Free energy (hartree): 0.0033897431
# Partition function: 0.2623715173

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0070028265
# Entropy (hartree/K): 0.0000000846
# Free energy (hartree): 0.0069351474
# Partition function: 0.0647365130

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0069810319
# Entropy (hartree/K): 0.0000000897
# Free energy (hartree): 0.0069092901
# Partition function: 0.0654006179

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0069819267
# Entropy (hartree/K): 0.0000000940
# Free energy (hartree): 0.0069067648
# Partition function: 0.0654658415

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0071099282
# Entropy (hartree/K): 0.0000000773
# Free energy (hartree): 0.0070481116
# Partition function: 0.0619133811

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0071681776
# Entropy (hartree/K): 0.0000000747
# Free energy (hartree): 0.0071084194
# Partition function: 0.0604569618

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0071991518
# Entropy (hartree/K): 0.0000000729
# Free energy (hartree): 0.0071408583
# Partition function: 0.0596877908

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0071790096
# Entropy (hartree/K): 0.0000000747
# Free energy (hartree): 0.0071192866
# Partition function: 0.0601981873

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0071765499
# Entropy (hartree/K): 0.0000000766
# Free energy (hartree): 0.0071152952
# Partition function: 0.0602931029

# 	********** Final results **********


# Temperature (K): 800.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004643552
# Translational energy (kcal/mol): 2.3846478011
# Translational entropy (cal/mol/K): 42.2047745300
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.3846478011
# Rotational entropy (cal/mol/K): 24.2237384405
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 76.6158056847
# Internal (tor+vib) entropy (cal/mol/K): 23.5917456698
# Internal (tor+vib) Cv (cal/mol/K): 28.0714005513


# Total energy (kcal/mol): 81.3851012870
# Total enthalpy (kcal/mol): 82.9748665779
# Enthalpy H(800.000000 K)-H(0 K) (kcal/mol):  17.4744022227
# Total entropy (cal/mol/K): 90.0202586402
# Total Cv (cal/mol/K): 34.0330200541



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
# 	 52 		 0.0001486174 		 0.0019683777
# 
# 	 53 		 0.0001490170 		 0.0001291716
# 
# 	 54 		 0.0000861283 		 0.0021405194
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8246746909
# Zero point vibrational energy (hartree): 0.0005082576
# Energy (hartree): 0.0030296222
# Entropy (hartree/K): 0.0000065563
# Free energy (hartree): -0.0028710544
# Partition function: 2.7383082168

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0001486174 		 0.0019683791
# 
# 	 53 		 0.0001490171 		 0.0001291724
# 
# 	 54 		 0.0000861284 		 0.0021405201
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8237663468
# Zero point vibrational energy (hartree): 0.0005082575
# Energy (hartree): 0.0030296216
# Entropy (hartree/K): 0.0000065563
# Free energy (hartree): -0.0028710544
# Partition function: 2.7383082780

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000001 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0029347834
# Entropy (hartree/K): 0.0000047810
# Free energy (hartree): -0.0013680839
# Partition function: 1.6160856542

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0031816100
# Entropy (hartree/K): 0.0000026974
# Free energy (hartree): 0.0007539901
# Partition function: 0.7675546885

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0033214357
# Entropy (hartree/K): 0.0000023223
# Free energy (hartree): 0.0012313339
# Partition function: 0.6491920169

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0033536897
# Entropy (hartree/K): 0.0000021840
# Free energy (hartree): 0.0013881319
# Partition function: 0.6144418346

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0033828518
# Entropy (hartree/K): 0.0000021001
# Free energy (hartree): 0.0014927477
# Partition function: 0.5922973038

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0035159040
# Entropy (hartree/K): 0.0000017972
# Free energy (hartree): 0.0018984298
# Partition function: 0.5137160889

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0036701946
# Entropy (hartree/K): 0.0000015506
# Free energy (hartree): 0.0022746553
# Partition function: 0.4501892298

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0037396822
# Entropy (hartree/K): 0.0000014624
# Free energy (hartree): 0.0024235668
# Partition function: 0.4272719980

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0038437361
# Entropy (hartree/K): 0.0000013235
# Free energy (hartree): 0.0026526052
# Partition function: 0.3942795754

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0039276206
# Entropy (hartree/K): 0.0000012212
# Free energy (hartree): 0.0028285121
# Partition function: 0.3706808335

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0039812984
# Entropy (hartree/K): 0.0000011653
# Free energy (hartree): 0.0029325460
# Partition function: 0.3573944163

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0039997323
# Entropy (hartree/K): 0.0000011483
# Free energy (hartree): 0.0029662618
# Partition function: 0.3531915051

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0041015719
# Entropy (hartree/K): 0.0000010573
# Free energy (hartree): 0.0031499615
# Partition function: 0.3311453561

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0041112598
# Entropy (hartree/K): 0.0000010545
# Free energy (hartree): 0.0031622303
# Partition function: 0.3297229595

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0041298405
# Entropy (hartree/K): 0.0000010357
# Free energy (hartree): 0.0031976798
# Partition function: 0.3256473137

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0041549767
# Entropy (hartree/K): 0.0000010169
# Free energy (hartree): 0.0032397605
# Partition function: 0.3208746298

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0041908434
# Entropy (hartree/K): 0.0000009879
# Free energy (hartree): 0.0033017112
# Partition function: 0.3139753284

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0070514002
# Entropy (hartree/K): 0.0000001416
# Free energy (hartree): 0.0069239621
# Partition function: 0.0880943803

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0070316373
# Entropy (hartree/K): 0.0000001491
# Free energy (hartree): 0.0068974796
# Partition function: 0.0889167380

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0070342194
# Entropy (hartree/K): 0.0000001553
# Free energy (hartree): 0.0068944292
# Partition function: 0.0890119535

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0071555259
# Entropy (hartree/K): 0.0000001308
# Free energy (hartree): 0.0070378312
# Partition function: 0.0846441960

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0072127259
# Entropy (hartree/K): 0.0000001270
# Free energy (hartree): 0.0070984570
# Partition function: 0.0828627233

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0072429319
# Entropy (hartree/K): 0.0000001242
# Free energy (hartree): 0.0071311234
# Partition function: 0.0819184267

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0072235500
# Entropy (hartree/K): 0.0000001269
# Free energy (hartree): 0.0071093292
# Partition function: 0.0825472354

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0072219041
# Entropy (hartree/K): 0.0000001298
# Free energy (hartree): 0.0071050995
# Partition function: 0.0826698288

# 	********** Final results **********


# Temperature (K): 900.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004513761
# Translational energy (kcal/mol): 2.6827287763
# Translational entropy (cal/mol/K): 42.7899225654
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.6827287763
# Rotational entropy (cal/mol/K): 24.5748272617
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 79.5519990394
# Internal (tor+vib) entropy (cal/mol/K): 27.0474225526
# Internal (tor+vib) Cv (cal/mol/K): 30.5736105812


# Total energy (kcal/mol): 84.9174565920
# Total enthalpy (kcal/mol): 86.7059425442
# Enthalpy H(900.000000 K)-H(0 K) (kcal/mol):  21.2054911680
# Total entropy (cal/mol/K): 94.4121723797
# Total Cv (cal/mol/K): 36.5352300840



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
# 	 52 		 0.0003016568 		 0.0019683777
# 
# 	 53 		 0.0003023964 		 0.0001291716
# 
# 	 54 		 0.0001835243 		 0.0021405194
# 
# 	 55 		 0.0001836232 		 0.0001064918
# 
# 	 56 		 0.0001097204 		 0.0021769421
# 
# 	 57 		 0.0001094344 		 0.0000881719
# 
# 	 58 		 0.0000645165 		 0.0021095520
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8291558486
# Zero point vibrational energy (hartree): 0.0005082371
# Energy (hartree): 0.0032745528
# Entropy (hartree/K): 0.0000068148
# Free energy (hartree): -0.0035402752
# Partition function: 3.0585142631

# 	********** Mode 2 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0003016569 		 0.0019683791
# 
# 	 53 		 0.0003023964 		 0.0001291724
# 
# 	 54 		 0.0001835244 		 0.0021405201
# 
# 	 55 		 0.0001836232 		 0.0001064922
# 
# 	 56 		 0.0001097205 		 0.0021769423
# 
# 	 57 		 0.0001094344 		 0.0000881721
# 
# 	 58 		 0.0000645166 		 0.0021095520
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 40.8282475054
# Zero point vibrational energy (hartree): 0.0005082371
# Energy (hartree): 0.0032745522
# Entropy (hartree/K): 0.0000068148
# Free energy (hartree): -0.0035402752
# Partition function: 3.0585142675

# 	********** Mode 3 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000005 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 380.7559734379
# Zero point vibrational energy (hartree): 0.0008691869
# Energy (hartree): 0.0032408615
# Entropy (hartree/K): 0.0000051034
# Free energy (hartree): -0.0018625744
# Partition function: 1.8006610046

# 	********** Mode 4 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 766.4771170053
# Zero point vibrational energy (hartree): 0.0017405515
# Energy (hartree): 0.0034612218
# Entropy (hartree/K): 0.0000029919
# Free energy (hartree): 0.0004693220
# Partition function: 0.8622586330

# 	********** Mode 5 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 886.0538241084
# Zero point vibrational energy (hartree): 0.0020245978
# Energy (hartree): 0.0035959166
# Entropy (hartree/K): 0.0000026115
# Free energy (hartree): 0.0009844666
# Partition function: 0.7328094178

# 	********** Mode 6 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 930.0283117176
# Zero point vibrational energy (hartree): 0.0021179137
# Energy (hartree): 0.0036218123
# Entropy (hartree/K): 0.0000024664
# Free energy (hartree): 0.0011554465
# Partition function: 0.6942933228

# 	********** Mode 7 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 959.3969027048
# Zero point vibrational energy (hartree): 0.0021841275
# Energy (hartree): 0.0036476235
# Entropy (hartree/K): 0.0000023790
# Free energy (hartree): 0.0012686290
# Partition function: 0.6699173710

# 	********** Mode 8 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1080.4084647657
# Zero point vibrational energy (hartree): 0.0024551899
# Energy (hartree): 0.0037688408
# Entropy (hartree/K): 0.0000020636
# Free energy (hartree): 0.0017052558
# Partition function: 0.5836366932

# 	********** Mode 9 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1198.3829620108
# Zero point vibrational energy (hartree): 0.0027283382
# Energy (hartree): 0.0039120752
# Entropy (hartree/K): 0.0000018053
# Free energy (hartree): 0.0021067520
# Partition function: 0.5141403870

# 	********** Mode 10 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1245.7584869351
# Zero point vibrational energy (hartree): 0.0028419783
# Energy (hartree): 0.0039771013
# Entropy (hartree/K): 0.0000017124
# Free energy (hartree): 0.0022647350
# Partition function: 0.4891207151

# 	********** Mode 11 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1324.4630260928
# Zero point vibrational energy (hartree): 0.0030180882
# Energy (hartree): 0.0040724432
# Entropy (hartree/K): 0.0000015643
# Free energy (hartree): 0.0025081379
# Partition function: 0.4529350941

# 	********** Mode 12 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1387.4183361429
# Zero point vibrational energy (hartree): 0.0031566485
# Energy (hartree): 0.0041489388
# Entropy (hartree/K): 0.0000014543
# Free energy (hartree): 0.0026946728
# Partition function: 0.4270264017

# 	********** Mode 13 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1424.6965035293
# Zero point vibrational energy (hartree): 0.0032406334
# Energy (hartree): 0.0041985333
# Entropy (hartree/K): 0.0000013940
# Free energy (hartree): 0.0028045255
# Partition function: 0.4124673909

# 	********** Mode 14 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1436.3711703578
# Zero point vibrational energy (hartree): 0.0032683420
# Energy (hartree): 0.0042156829
# Entropy (hartree/K): 0.0000013757
# Free energy (hartree): 0.0028400096
# Partition function: 0.4078714994

# 	********** Mode 15 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000001
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1502.4634886093
# Zero point vibrational energy (hartree): 0.0034203834
# Energy (hartree): 0.0043104688
# Entropy (hartree/K): 0.0000012773
# Free energy (hartree): 0.0030331911
# Partition function: 0.3837342773

# 	********** Mode 16 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000002
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1505.3814982738
# Zero point vibrational energy (hartree): 0.0034315030
# Energy (hartree): 0.0043202433
# Entropy (hartree/K): 0.0000012745
# Free energy (hartree): 0.0030457433
# Partition function: 0.3822162936

# 	********** Mode 17 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1519.1937641984
# Zero point vibrational energy (hartree): 0.0034606954
# Energy (hartree): 0.0043370559
# Entropy (hartree/K): 0.0000012539
# Free energy (hartree): 0.0030831630
# Partition function: 0.3777265358

# 	********** Mode 18 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1533.7548162262
# Zero point vibrational energy (hartree): 0.0034964439
# Energy (hartree): 0.0043605422
# Entropy (hartree/K): 0.0000012333
# Free energy (hartree): 0.0031272162
# Partition function: 0.3725084000

# 	********** Mode 19 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000000
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 1556.7729586366
# Zero point vibrational energy (hartree): 0.0035487164
# Energy (hartree): 0.0043938857
# Entropy (hartree/K): 0.0000012017
# Free energy (hartree): 0.0031922030
# Partition function: 0.3649419872

# 	********** Mode 20 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000038
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3056.7549842723
# Zero point vibrational energy (hartree): 0.0069455467
# Energy (hartree): 0.0071188929
# Entropy (hartree/K): 0.0000002125
# Free energy (hartree): 0.0069063626
# Partition function: 0.1129449888

# 	********** Mode 21 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000165
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3018.8966083269
# Zero point vibrational energy (hartree): 0.0069204270
# Energy (hartree): 0.0071015563
# Entropy (hartree/K): 0.0000002225
# Free energy (hartree): 0.0068790066
# Partition function: 0.1139248682

# 	********** Mode 22 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000459
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 2988.6880859248
# Zero point vibrational energy (hartree): 0.0069185284
# Energy (hartree): 0.0071061534
# Entropy (hartree/K): 0.0000002309
# Free energy (hartree): 0.0068752251
# Partition function: 0.1140609884

# 	********** Mode 23 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000197
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3115.7842095654
# Zero point vibrational energy (hartree): 0.0070574601
# Energy (hartree): 0.0072194826
# Entropy (hartree/K): 0.0000001980
# Free energy (hartree): 0.0070214994
# Partition function: 0.1089123673

# 	********** Mode 24 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000112
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3137.9693648417
# Zero point vibrational energy (hartree): 0.0071174023
# Energy (hartree): 0.0072754536
# Entropy (hartree/K): 0.0000001929
# Free energy (hartree): 0.0070825704
# Partition function: 0.1068321474

# 	********** Mode 25 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000153
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3154.1143118738
# Zero point vibrational energy (hartree): 0.0071495835
# Energy (hartree): 0.0073047395
# Entropy (hartree/K): 0.0000001892
# Free energy (hartree): 0.0071155583
# Partition function: 0.1057250822

# 	********** Mode 26 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000055
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3138.4338596295
# Zero point vibrational energy (hartree): 0.0071282627
# Energy (hartree): 0.0072862822
# Entropy (hartree/K): 0.0000001928
# Free energy (hartree): 0.0070934479
# Partition function: 0.1064658277

# 	********** Mode 27 **********
# 
# 	 51 		-		-
# 
# 	 52 		 0.0000000000 		 0.0000000022
# Convergence criterion met
# ------------------------------------
# Frequency (cm-1): 3122.0859558210
# Zero point vibrational energy (hartree): 0.0071245410
# Energy (hartree): 0.0072856310
# Entropy (hartree/K): 0.0000001968
# Free energy (hartree): 0.0070888793
# Partition function: 0.1066195296

# 	********** Final results **********


# Temperature (K): 1000.00
# Pressure (Pa): 100000
# Zero point vibrational energy (kcal/mol): 65.5004257514
# Translational energy (kcal/mol): 2.9808097514
# Translational entropy (cal/mol/K): 43.3133553195
# Translational Cv (cal/mol/K): 2.9808097514
# Rotational energy (kcal/mol): 2.9808097514
# Rotational entropy (cal/mol/K): 24.8888869142
# Rotational Cv (cal/mol/K): 2.9808097514
# Internal (rot+vib) energy (kcal/mol): 82.7249220109
# Internal (tor+vib) entropy (cal/mol/K): 30.3888114412
# Internal (tor+vib) Cv (cal/mol/K): 32.7852224652


# Total energy (kcal/mol): 88.6865415138
# Total enthalpy (kcal/mol): 90.6737481274
# Enthalpy H(1000.000000 K)-H(0 K) (kcal/mol):  25.1733223759
# Total entropy (cal/mol/K): 98.5910536749
# Total Cv (cal/mol/K): 38.7468419681



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

