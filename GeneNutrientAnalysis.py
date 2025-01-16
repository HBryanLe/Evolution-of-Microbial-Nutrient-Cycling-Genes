# Script by Bryan Le, Carleton College, July 2022

# This Python script looks into how marine microbiomes respond to nutrient-poor conditions
# by identifying genes correlated with specific environmental stressors. Using data from 
# TARA Ocean Samples, the script performs linear regression analyses between selected 
# nutrients (nitrate, nitrogen dioxide, phosphate, and silicon) and gene abundances. 

# To call script:
# 	>> python GeneNutrientAnalysis.py...
# 		Arguments:
#           OM.CompanionTables.csv (required)
#           OM-RGC_Kegg_Database.tsv (required)
# 			OM-RGC_v2_gene_profile_metaT.tsv (required)



##########################################################################################
### Functions:

def linear_regression(nutrient, abundance):
    import numpy as np

    r = np.corrcoef(nutrient, abundance)
    r = r[0, 1]
    r_squared = r**2
    r = round(r, 2)
    r_squared = round(r_squared, 2)

    return r, r_squared



##########################################################################################
### Main Script:

def main():
    import sys
    import csv
    import math

    inputFile_1 = sys.argv[1]
    inputFile_2 = sys.argv[2]
    inputFile_3 = sys.argv[3]
    base = 2

    # Step 1: Creating Output File
    outputFile = open(
        "Coverage_LogTransformed_GeneNutrientAnalysis.csv", "w", newline=""
    )
    outputDictWriter = csv.DictWriter(
        outputFile,
        [
            "COG",
            "r_mean_nitrate",
            "r2_mean_nitrate",
            "r_nitrogen_dioxide",
            "r2_nitrogen_dioxide",
            "r_phosphate",
            "r2_phosphate",
            "r_nitrate_nitrite",
            "r2_nitrate_nitrite",
            "r_silicon",
            "r2_silicon",
        ],
    )
    outputDictWriter.writeheader()

    # Step 2: Associating TARA Sample IDs with their Respective Nutrient Concentrations
    dict_mn = {}
    dict_no2 = {}
    dict_po4 = {}
    dict_nn = {}
    dict_si = {}

    companion_file = open(inputFile_1, "r")
    companion_reader = companion_file.readlines()

    for row in companion_reader[1:]:
        col = row.split(",")

        if col[1].strip("\n") != "NA":
            dict_mn[col[0]] = float(col[1].strip("\n"))

        if col[2].strip("\n") != "NA":
            dict_no2[col[0]] = float(col[2].strip("\n"))

        if col[3].strip("\n") != "NA":
            dict_po4[col[0]] = float(col[3].strip("\n"))

        if col[4].strip("\n") != "NA":
            dict_nn[col[0]] = float(col[4].strip("\n"))

        if col[5].strip("\n") != "NA":
            dict_si[col[0]] = float(col[5].strip("\n"))

    # Step 3: Associating OM-RGC Gene IDs to their Respective KO Numbers
    dict_ID = {}

    ID_file = open(inputFile_2, "r")
    ID_reader = ID_file.readlines()

    for row in ID_reader[1:]:
        col = row.split("\t")
        dict_ID[col[0]] = col[1].strip("\n")

    # Step 4: Creating List of TARA Sample IDs
    TARA_ID = []

    abundance_file = open(inputFile_3, "r")
    abundance_reader = abundance_file.readlines()

    for row in abundance_reader:
        col = row.split("\t")
        for item in col:
            TARA_ID.append(item.strip("\n"))
        del TARA_ID[0]
        break

    # Step 5: Using Arrays to Allign Gene Abundance and Nutrient Concentration Respective to each TARA Sample
    count = 0

    for row in abundance_reader[1:]:
        list_mn = []
        list_no2 = []
        list_po4 = []
        list_nn = []
        list_si = []

        list_mn_abundance = []
        list_no2_abundance = []
        list_po4_abundance = []
        list_nn_abundance = []
        list_si_abundance = []

        col = row.split("\t")
        gene_ID = col[0].strip()
        if col[0] in dict_ID:
            COG = dict_ID.get(gene_ID)

            for i in TARA_ID:
                if i in dict_mn and col[TARA_ID.index(i) + 1].strip() > "0":
                    list_mn.append(float(dict_mn.get(i)))
                    list_mn_abundance.append(
                        math.log((float(col[TARA_ID.index(i) + 1].strip())), base)
                    )
                    print((float(col[TARA_ID.index(i) + 1].strip())))

                if i in dict_no2 and col[TARA_ID.index(i) + 1].strip() > "0":
                    list_no2.append(float(dict_no2.get(i)))
                    list_no2_abundance.append(
                        math.log((float(col[TARA_ID.index(i) + 1].strip())), base)
                    )

                if i in dict_po4 and col[TARA_ID.index(i) + 1].strip() > "0":
                    list_po4.append(float(dict_po4.get(i)))
                    list_po4_abundance.append(
                        math.log((float(col[TARA_ID.index(i) + 1].strip())), base)
                    )

                if i in dict_nn and col[TARA_ID.index(i) + 1].strip() > "0":
                    list_nn.append(float(dict_nn.get(i)))
                    list_nn_abundance.append(
                        math.log((float(col[TARA_ID.index(i) + 1].strip())), base)
                    )

                if i in dict_si and col[TARA_ID.index(i) + 1].strip() > "0":
                    list_si.append(float(dict_si.get(i)))
                    list_si_abundance.append(
                        math.log((float(col[TARA_ID.index(i) + 1].strip())), base)
                    )

            # Step 6: Performing Linear Regression & Calculating Statistical Values
            r_mean_nitrate, r2_mean_nitrate = linear_regression(
                list_mn, list_mn_abundance
            )
            r_nitrogen_dioxide, r2_nitrogen_dioxide = linear_regression(
                list_no2, list_no2_abundance
            )
            r_phosphate, r2_phosphate = linear_regression(list_po4, list_po4_abundance)
            r_nitrate_nitrite, r2_nitrate_nitrite = linear_regression(
                list_nn, list_nn_abundance
            )
            r_silicon, r2_silicon = linear_regression(list_si, list_si_abundance)

            # Step 7: Writing to Output File
            outputDictWriter.writerow(
                {
                    "COG": COG,
                    "r_mean_nitrate": r_mean_nitrate,
                    "r2_mean_nitrate": r2_mean_nitrate,
                    "r_nitrogen_dioxide": r_nitrogen_dioxide,
                    "r2_nitrogen_dioxide": r2_nitrogen_dioxide,
                    "r_phosphate": r_phosphate,
                    "r2_phosphate": r2_phosphate,
                    "r_nitrate_nitrite": r_nitrate_nitrite,
                    "r2_nitrate_nitrite": r2_nitrate_nitrite,
                    "r_silicon": r_silicon,
                    "r2_silicon": r2_silicon,
                }
            )

        count += 1
        percent = (count / 46775154) * 100
        percent_rounded = str(round(percent, 5))
        percent_rounded = percent_rounded + "% Complete"
        print(percent_rounded)

    outputFile.close()


if __name__ == "__main__":
    main()
