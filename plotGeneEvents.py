# Script by Bryan Le, Carleton College, July 2022

# To call script:
#       >> python3 plotGeneEvents.py

# This Python script visualizes frequency of gene events (e.g. gene loss, gain, and 
# horizontal gene transfer) within the prokaryotic lineage across Earth's history. We used 
# the Gene Events file for each gene to construct histograms showing the proportion of 
# event frequencies over time. This analysis provides a temporal perspective on the birth, 
# loss, and transfer of genes between branches



##########################################################################################
def main():
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import os
    
    # Step 1: Looping through Gene Event File & Collecting Data on Dates of each Event
    for i in os.listdir("INSERT PATHWAY TO DIRECTORY CONTAINING GENE EVENT FILES"):
        gene_name = i[:-37]
        geneFilePath = "INSERT PATHWAY TO DIRECTORY CONTAINING GENE EVENT FILES" + i

        pathway = i[-36:-33]
        if pathway == "tst" and pathway != ("dsm" or "org"):
            border_color = "#85c242"
        if pathway == "org" and pathway != ("dsm" or "tst"):
            border_color = "#de6467"
        if pathway == "dsm" and pathway != ("tst" or "org"):
            border_color = "#8360a9"

        file = open(geneFilePath, "r")
        data = file.read().splitlines()

        hgt_events = []
        spe_events = []
        los_events = []
        dup_events = []
        total_events = []
        num_los_events = 0
        num_hgt_events = 0
        num_dup_events = 0
        num_spe_events = 0
        num_total_events = -1

        for i in data:
            num_total_events += 1
            line = i.split("\t")
            if line[0] == "hgt" and line[5] != "?":
                hgt_events.append(float(line[5]))
                total_events.append(float(line[5]))
                num_hgt_events += 1
            elif line[0] == "spe" and line[5] != "?":
                spe_events.append(float(line[5]))
            elif line[0] == "dup" and line[5] != "?":
                dup_events.append(float(line[5]))
            elif line[0] == "los" and line[5] != "?":
                los_events.append(float(line[5]))

        # Step 2: Plotting Arrays for HGT and Total Event Dates on Histogram
        plt.style.use("ggplot")
        bins = [
            0,
            250,
            500,
            750,
            1000,
            1250,
            1500,
            1750,
            2000,
            2250,
            2500,
            2750,
            3000,
            3250,
            3500,
            3750,
            4000,
        ]
        hgt_color = "#00cecb"
        hgt_alpha = float(0.80)
        spe_color = "#CB21D4"
        spe_alpha = float(0.60)
        los_color = "#FF4A47"
        los_alpha = float(0.60)
        dup_color = "#FFB60A"
        dup_alpha = float(0.60)

        fig, ax = plt.subplots(linewidth=10, edgecolor=border_color)
        plt.gca().invert_xaxis()
        ax.set_title(gene_name, fontsize=18, y=1.025)
        ax.set_xlabel("Million Years Ago", fontsize=12)
        ax.set_ylabel("Proportion of Total Events", fontsize=12)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=num_total_events))
        ax.set_facecolor("#e7e7e7")
        ax.tick_params(axis="x", labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
        # ax.grid(which = "major", linewidth = 0.945)
        # ax.grid(which = "minor", linewidth = 0.0945)

        ax.hist(
            x=(hgt_events, spe_events, los_events, dup_events),
            bins=bins,
            alpha=0.80,
            color=(hgt_color, spe_color, los_color, dup_color),
            width=250,
            label=(
                "HGT events",
                "Speciation Events",
                "Loss Events",
                "Duplication events",
            ),
            stacked=True,
        )
        # ax.hist(hgt_events, bins = bins, alpha = hgt_alpha, color = hgt_color, width = 94.5, label = "hgt")

        # ax.legend(loc = 2, facecolor = "#ffffff", framealpha = 0.85)
        plt.axis("off")
        outputFile = gene_name + "_cirGeneEvents_Histogram.pdf"
        outputPath = (
            "INSERT PATH TO OUTPUT DIRECTORY HERE" + outputFile
        )
        plt.savefig(outputPath, edgecolor=fig.get_edgecolor())

if __name__ == "__main__":
    main()
