# -*- coding: utf-8 -*-
"""
Created on Thu May 12 15:05:48 2022

@author: John Klumpp

Description: This code calculates tables of inhalation dose coefficients for a wide variety of radionuclides to workers
and members of the public, including children and infants. It uses KDEP to calculate the fractional deposition of the
inhaled material, nd thus accepts a wide range of parameters describing the inhaled aerosol, the individual doing the breathing,
and the activity being performed at the time of inhalation (e.g., sleep, excercise, etc.). 

This code uses the ICRP 66 respiratory tract model, and the ICRP 60 and 70 series dosimetry system for workers and 
members of the public, respectively. Hoever, the code is easily modified to use the OIR/EIR dosimetry systems, as KDEP
already accomodates the ICRP 130 respiratory tract deposition model. The option to use the OIR/EIR model will be added
as soon as EIR is complete. 

Here is a brief description of how it works. First it reads the input files, Nuclides.csv, 
SizesInp.csv, and CaseParams.csv. It then writes input files for KDEP, calls KDEP, and reads 
the KDEP output files. This contains the fractional depisition of inhaled particles in the respiratory 
tract. Next, the script reads "HDB" tables, which contain organ doses per unit of deposition in each region 
of the respiratory tract. This code adds up the total organ doses per unit of intake using these tables and 
the KDEP output. Finally, the tissue-weighted sum of the organ doses is calculated to give the committed 
effective dose. Committed Effective doses and organ doses are written to tables in the OutputFiles directory.

This code was designed to work in concert with QUIC to estimate doses to individuals downind of a radoactive plume.
Many of the design decisions were informed by the need to facilitate integration with QUIC.
"""

# -*- coding: utf-8 -*-


import csv
import pandas as pd
import numpy as np
import os, subprocess

#os.chdir('C:\Dropbox\Workspaces\Python\Python3\DCAL-KDEP')
BaseDir = os.getcwd()
DCPAKDir = os.path.join(BaseDir, 'DC_PAK_Files')
KDEPDir = os.path.join(BaseDir, 'KDEP')
InputDir = os.path.join(KDEPDir, 'InputFiles')

activities = { 1: "Sleep", 2: "Sitting", 3:"LightExercise",4:"HeavyExercise",5:"StandardWorker",6:"HeavyWorker",7:"MemberOfPublicFullDay"}
subjects = { 1: "AdultMale", 2: "AdultFemale", 3:"15YearOldMale",4:"15YearOldFemale",5:"10YearOld",6:"5YearOld",7:"1YearOld",8:"3MonthOld"}
ages = { '1': 9125, '2': 9125, '3':5475, '4':5475, '5':3650, '6':1825 ,'7':365, '8':100}

#Remainder tissue masses from ICRP Publication 71 Table 11, Page 28. Dictionary key is subject.
remainder_tissues_dict = {1:[28000.,1400.,640.,310.,100.,180.,20.,80.,14.,15.5],2:[28000.,1400.,640.,310.,100.,180.,20.,80.,14.,15.5],\
                          3:[22000.,1410.,516.,248.,64.9,123.,28.4,80.,10.5,12.1], 4:[22000.,1410.,516.,248.,64.9,123.,28.4,80.,10.5,12.1],\
                          5:[11000.,1360.,286.,173.,30.,77.4,31.4,4.2,7.2,7.1],6:[5000.,1260.,169.,116.,23.6,48.3,29.6,2.7,5.3,4.3], \
                          7:[2500.,884.,84.9,62.9,10.3,25.5,22.9,1.5,3.5,2.2], 8:[760.,352.,32.6,22.9,2.8,9.1,11.3,3.9,5.8,1.3]}
    
main_tissues_wt = np.array([0.2,0.12,0.12,0.12,0.12,0.05,0.05,0.05,0.05,0.05,0.01,0.01])

#Read Input Files

print("Reading Nuclides.csv")
nuclides, solubilities,sizes,distributions =[],[],[],[]
NuclideLines,SizesLines = [],[]

with open("Nuclides.csv") as csvfile:
    lines = csv.reader(csvfile,delimiter=',')
    for line in lines:
        if (line[0][0]=='!'): #Skip comment lines
            pass
        else:
            NuclideLines= NuclideLines + [line]
            nuclides = nuclides+[line[0].strip()]
            solubilities = solubilities + [line[1].strip()]

#Remove column headers
nuclides = nuclides[1:]
solubilities = solubilities[1:]

print("Reading SizesInp.csv")
with open("SizesInp.csv") as csvfile:
    lines = csv.reader(csvfile,delimiter=',')
    for line in lines:
        SizesLines= SizesLines + [line]
        if (line[0][0]=='!'): #Skip comment lines
            pass
        else:
            distributions = distributions + [line[0]]
            sizes = sizes + [line[1]]

#Remove column headers
sizes = sizes[1:]
distributions = distributions[1:]

monodispursed, nose_breather,subject,activity,rho,shape_factor,wind_speed,Atm_Pressure,Chronic,ICRP130,age = ([] for i in range(11))
CaseParamsLines,CaseParamsHeader = [],[]

print("Reading CaseParams.csv")
with open("CaseParams.csv") as csvfile:
    lines = csv.reader(csvfile,delimiter=',')
    for line in lines:
        if (line[0][0]=='!'): #Skip comment lines
            CaseParamsHeader = CaseParamsHeader+[line]
            pass
        else:
            CaseParamsLines = CaseParamsLines + [line]

for line in CaseParamsLines[1:]:
    monodispursed = monodispursed + [line[0]]
    nose_breather = nose_breather + [line[1]]
    subject = subject + [int(line[2])]
    activity = activity + [int(line[3])]
    rho = rho + [line[4]]
    shape_factor = shape_factor + [line[5]]
    wind_speed = wind_speed + [line[6]]
    Atm_Pressure = Atm_Pressure + [line[7]]
    Chronic = Chronic + [line[8]]
    ICRP130 = ICRP130 + [line[9]]
    age = age + [ages[line[2]]]

#ReadHDB
print("Reading HDB Tables")
os.chdir(DCPAKDir)
dfET1 = pd.read_table("ET1.HDB", delim_whitespace=True,skiprows=1)
dfET2 = pd.read_table("ET2.HDB", delim_whitespace=True,skiprows=1)

dfBBEGEL = pd.read_table("BBE-GEL.HDB", delim_whitespace=True,skiprows=1)
dfBBESEQ = pd.read_table("BBE-SEQ.HDB", delim_whitespace=True,skiprows=1)
dfBBESOL = pd.read_table("BBE-SOL.HDB", delim_whitespace=True,skiprows=1)

dfBBISEQ = pd.read_table("BBI-SEQ.HDB", delim_whitespace=True,skiprows=1)
dfBBISOL = pd.read_table("BBI-SOL.HDB", delim_whitespace=True,skiprows=1)
dfBBIGEL = pd.read_table("BBI-GEL.HDB", delim_whitespace=True,skiprows=1)

dfAI = pd.read_table("AI.HDB", delim_whitespace=True,skiprows=1)

os.chdir(BaseDir)


#Call KDEP
os.chdir(InputDir)
print(os.getcwd())

print("Writing Sizes.csv")
with open('Sizes.csv', 'wb') as f: #Python 2 version
#with open('Sizes.csv', 'w+',newline='') as f: #Python3 version
    writer = csv.writer(f)
    writer.writerows(SizesLines)

for subject, activity,age,CaseParamsLine in zip(subject,activity,age,CaseParamsLines):
    print("Writing input.csv")
    print(subject,activity,age)
    os.chdir(InputDir)
    Lines =[]
    Lines = CaseParamsHeader+[CaseParamsLine]
    print(CaseParamsLine)
    with open('input.csv', 'wb') as f: #Python 2 Version
#    with open('input.csv', 'w',newline='') as f: #Python 3 Version
    #    f.write(CaseParamsLines)
        writer = csv.writer(f)
        writer.writerows(Lines)


    os.chdir(KDEPDir)
    print(os.getcwd())

    print("Running KDEP")
    subprocess.call('kdep.exe')
    print("KDEP Finished")


    AMADs, AMTDs,ET1s,ET2s,Bs,bbs,AIs,Totals,FsBs,Fsbbs,BBi_GELs, BBi_SOLs, BBi_SEQs, bbe_GELs, bbe_SOLs, bbe_SEQs =([] for i in range(16))
    KDEPLines = []
    print("Reading kdep.csv")
    with open("kdep.csv") as csvfile:
        lines = csv.reader(csvfile,delimiter=',')
        for line in lines:
            KDEPLines = KDEPLines + [line]

    n = len(KDEPLines)-4
    for line in KDEPLines[1:n]:
            AMADs = AMADs+[float(line[0])]
            AMTDs = AMTDs + [float(line[1])]
            ET1s = ET1s+[float(line[2])]
            ET2s = ET2s+[float(line[3])]
            Bs = Bs+[float(line[4])]
            bbs = bbs+[float(line[5])]
            AIs = AIs+[float(line[6])]
            Totals = Totals+[float(line[7])]
            FsBs = FsBs+[float(line[8])]
            Fsbbs = Fsbbs+[float(line[9])]
            BBi_GELs = BBi_GELs + [float(line[4])*(0.993-float(line[8]))]
            BBi_SOLs = BBi_SOLs + [float(line[4])*(float(line[8]))]
            BBi_SEQs = BBi_SEQs + [float(line[4])*(0.007)]
            bbe_GELs = bbe_GELs + [float(line[5])*(0.993-float(line[9]))]
            bbe_SOLs = bbe_SOLs + [float(line[5])*(float(line[9]))]
            bbe_SEQs = bbe_SEQs + [float(line[5])*(0.007)]


    for nuclide,solubility in zip(nuclides,solubilities):
        #CalcEffDose
        #find the right radionuclide
        print("Calculating Effective Dose")
        print(nuclide,solubility)
        queryline = "Nuclide==\""+nuclide+"\" and Age=="+str(age)+ " and Class==\""+solubility+"\""

        dfET1Eq = dfET1.query(queryline).iloc[:,3:]

        if (age==9125 and len(dfET1Eq)==0):
            print("Age not found, trying 7300")
            queryline = "Nuclide==\""+nuclide+"\" and Age=="+str(7300)+ " and Class==\""+solubility+"\""
            dfET1Eq = dfET1.query(queryline).iloc[:,3:]
        if (len(dfET1Eq)==0):
            print("Error: entry not found. Queryline = "+queryline)
            exit()

        f = float(dfET1Eq.loc[:,"f"])

        dfET2Eq = dfET2.query(queryline).iloc[:,4:]
        dfBBEGELEq = dfBBEGEL.query(queryline).iloc[:,4:]
        dfBBESEQEq = dfBBESEQ.query(queryline).iloc[:,4:]
        dfBBESOLEq = dfBBESOL.query(queryline).iloc[:,4:]
        dfBBISEQEq = dfBBISEQ.query(queryline).iloc[:,4:]
        dfBBISOLEq = dfBBISOL.query(queryline).iloc[:,4:]
        dfBBIGELEq = dfBBIGEL.query(queryline).iloc[:,4:]
        dfAIEq = dfAI.query(queryline).iloc[:,4:]

        CED = []
        data = []
        Columns = ["AMAD","AMTD","Nuclide","Subject","Activity","Solubility","f1","CED","Gonads","R-Marrow","Colon","Lungs","Stomach","Bladder", \
                   "Breasts","Liver","Oesophagus","Thyroid","Skin","Bone Surface","Muscle","Brain","SI-Wall","Kidneys","Pancreas","Spleen","Thymus", \
                   "Uterus","Adrenals","Extrathoracic Tissues","Remainder"]


        #Take sum weighted by lung depostion
        for AMAD,AMTD,ET1,ET2,bbe_GEL,bbe_SEQ,bbe_SOL,BBi_SEQ,BBi_SOL,BBi_GEL,AI in zip(AMADs,AMTDs,ET1s,ET2s,bbe_GELs,bbe_SEQs,bbe_SOLs,BBi_SEQs,BBi_SOLs,BBi_GELs,AIs):
             
            dfSumEq = dfET1Eq*ET1+dfET2Eq*ET2+dfBBEGELEq*bbe_GEL+dfBBESEQEq*bbe_SEQ+dfBBESOLEq*bbe_SOL+dfBBISEQEq*BBi_SEQ+dfBBISOLEq*BBi_SOL+dfBBIGELEq*BBi_GEL+dfAIEq*AI
            print("AMAD",AMAD)

            #Add up tissues to determine CED
            colon = float(dfSumEq.loc[:,"ULI-Wall"]*0.57+dfSumEq.loc[:,"LLI-Wall"]*0.43)
            Ovaries = float(dfSumEq.loc[:,"Ovaries"])
            Testes = float(dfSumEq.loc[:,"Testes"])
            gonads = max(Ovaries,Testes)

            #Main Tissues are Gonads, Bone Marrow, Colon, Lung, Stomach, Bladder, Breast, Liver, Oesophagus, Thyroid, Skin, and Bone Surface
            main_tissues = np.array([gonads,float(dfSumEq.loc[:,"R-Marrow"]),colon,float(dfSumEq.loc[:,"Lungs"]),float(dfSumEq.loc[:,"St-Wall"]), \
                                     float(dfSumEq.loc[:,"UB-Wall"]),float(dfSumEq.loc[:,"Breasts"]),float(dfSumEq.loc[:,"Liver"]), \
                                     float(dfSumEq.loc[:,"Thymus"]),float(dfSumEq.loc[:,"Thyroid"]),float(dfSumEq.loc[:,"Skin"]),\
                                     float(dfSumEq.loc[:,"B-Surface"])])
            main_tissues_eff = np.dot(main_tissues,main_tissues_wt)

            #Remainder Tissues are Muscle, Brain, Small Intestine, Kidneys, Pancreas, Spleen, Thymus, Uterus, Adrenals, and Extrathoracic Airways
            remainder_tissues = np.array([float(dfSumEq.loc[:,"Muscle"]),float(dfSumEq.loc[:,"Brain"]),float(dfSumEq.loc[:,"SI-Wall"]),float(dfSumEq.loc[:,"Kidneys"]), \
                                          float(dfSumEq.loc[:,"Pancreas"]),float(dfSumEq.loc[:,"Spleen"]),float(dfSumEq.loc[:,"Thymus"]),float(dfSumEq.loc[:,"Uterus"]),\
                                          float(dfSumEq.loc[:,"Adrenals"]),float(dfSumEq.loc[:,"ET-Region"])])
            remainder_masses = np.array(remainder_tissues_dict[subject])
            remainder_mass = sum(remainder_masses)

            if max(remainder_tissues)<max(main_tissues):
                print("Remainder formulation: standard")
                remainder_tissues = remainder_tissues
                remainder_mass_fractions = remainder_masses/remainder_mass
                remainder_eff = np.dot(remainder_mass_fractions,remainder_tissues)*0.05

            else:
                print("Remainder formulation: split")
                i_max_remainder = np.argmax(remainder_tissues)
                max_remainder_eff = max(remainder_tissues)*0.025
                remainder_masses[i_max_remainder] = 0
                remainder_mass_fractions = remainder_masses/remainder_mass
                remainder_eff = np.dot(remainder_mass_fractions,remainder_tissues)*0.025
                remainder_eff = remainder_eff+max_remainder_eff
            person = subjects[subject]
            act = activities[activity]
            remainder_eq = remainder_eff/0.05
            e50 = [main_tissues_eff+remainder_eff]
            output = np.concatenate(([AMAD],[AMTD],[nuclide],[person],[act],[solubility],[f],e50,main_tissues,remainder_tissues,[remainder_eq]))
            data = data + [list(output)]


        Filename = str(nuclide)+'-'+person+'-'+act+'-'+str(solubility)+'.csv'
        results = pd.DataFrame(data,columns=Columns)
        filepath = os.path.join(BaseDir, 'OutputFiles',Filename)
        results.to_csv(filepath)