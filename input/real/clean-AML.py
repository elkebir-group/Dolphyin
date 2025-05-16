# Make imports
import os
import os.path
import sys

# Define executables
dolphyinLocation = "./Dolphyin/dolphyin.o"

def getAMLID(i):
    if i < 10:
        return "0" + (str)(i)
    return (str)(i)

for i in range(1, 120):

    # Read the data in
    AMLID = getAMLID(i)
    dataFileName = "AML_input/AML-" + AMLID + "-001_input.csv"
    if not os.path.isfile(dataFileName):
        print(dataFileName, "is not a file!")
        continue
    infile = open(dataFileName, "r")
    data = [line[:-1].split(",") for line in infile]

    # Clean the data
    cleanData = [row for row in data if "-1" not in row]
    numCells = len(cleanData)
    numChars = len(cleanData[0])

    # Write the data, cleaned
    cleanedDataFileName = "AML_input/AML-" + AMLID + "-001-cleaned.csv"
    outfile = open(cleanedDataFileName, "w")
    outfile.write("# " + (str)(numCells) + " CELLS \n")
    outfile.write("# " + (str)(numChars) + " CHARACTERS \n")

    for line in cleanData:
        for char in line:
            if char == "-1":
                outfile.write("0")
            else:
                outfile.write(char)
            outfile.write(" ")
        outfile.write("\n")
    outfile.close()
    print("Wrote", dataFileName, "without -1s!")