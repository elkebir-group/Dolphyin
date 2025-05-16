# Make imports
import os
import sys
import random

# For all perfect datasets, get options ovr matrix dimension, loss chance, and replicate
# Copy over the data, but with random errors

mNs = [
        "m25_n25", "m25_n50", "m25_n100", 
        "m50_n25", "m50_n50", "m50_n100", 
        "m100_n25", "m100_n50", "m100_n100"
    ]
lossVals = ["0.1", "0.2", "0.4"]
replicates = [("s" + (str)(k)) for k in range(1, 21)]

# mNs = [
#         "m100_n100"
#     ]
# lossVals = ["0.4"]
# replicates = ["s20"]

for mn in mNs:
    for replicate in replicates:
        for lossVal in lossVals:

            # For each perfect dataset combination, get the input and output location
            dataFileName = "../../../SphyRCodebase/SPhyR/data/k_dollo/" + \
                mn + "_" + replicate + "_k1_loss" + lossVal + ".B"
            modDataFileName = "errors/" + mn + "_" + replicate + "_k1_loss" + lossVal + ".txt"
            fileIn = open(dataFileName, "r") 
            fileOut = open(modDataFileName, "w")

            # Set the random seed
            random.seed((int)(replicate[1:]))

            # Read all lines and write 1s if original has a 1 and no FN
            lineNum = 0
            numFalse, numTotal = 0, 0

            for line in fileIn:
                lineNum += 1
                if lineNum <= 2:
                    fileOut.write(line)
                    continue
                entries = line[:-1].split(" ")
                newLine = ""
                for entry in entries:
                    if entry == "1":
                        if random.random() < 0.05:
                            newLine += "0 "
                            numFalse += 1
                        else:
                            newLine += "1 "
                        numTotal += 1
                    else:
                        newLine += "0 "
                fileOut.write(newLine[:-1] + "\n")
            fileOut.close()
            print("FNs induced:", numFalse / numTotal)
