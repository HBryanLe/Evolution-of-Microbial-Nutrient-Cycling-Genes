# Evolution-of-Microbial-Nutrient-Cycling-Genes

This repository contains some of Python scripts developed during my undergraduate research in Dr. Anderson’s lab at Carleton College. My research focused on microbial evolution and the adaptation of marine microbiomes to changing environmental conditions on early Earth. The work aimed to uncover the mechanisms that enabled microbial survival during historical atmospheric and nutrient transitions. 

*** Note: This is repository is not a complete representation of what I completed while in Dr. Anderson's lab.  Rather, these are some select projects that I worked on.

Through this research, I contributed to advancing our understanding of:
- Gene abundance patterns in marine metagenomes under nutrient stress.
- The evolutionary dynamics of gene loss, gain, and horizontal transfer in prokaryotic lineages.

The repository contains scripts that address these research questions and provide tools for analyzing gene events and nutrient-gene correlations in microbial communities.

## Overview of Research

1. **Studying Marine Microbiome Adaptations**  
   One of my primary projects involved analyzing the relationship between gene abundance and nutrient availability in marine environments. Using data from TARA Ocean samples, I identified genes linked to low nitrogen and phosphorus conditions, providing insights into microbial survival strategies during nutrient-poor periods in Earth's history.

2. **Tracing Gene Events in Prokaryotes**  
   I also worked on comparative phylogenetics to explore the evolutionary history of sulfur-cycling metabolic genes. By tracking gene loss, gain, and horizontal gene transfer events, I helped uncover how microbial lineages responded to shifts in nutrient availability across geological time.

## Scripts

### 1. **Plot Gene Events (plot_Gene_Events)**
   - **Description**: This script visualizes the frequency of gene events (e.g., gene loss, gain, and horizontal gene transfer) within prokaryotic lineages across Earth's history.
   - **Functionality**:
     - Reads Gene Events files for each gene.
     - Constructs histograms showing the proportion of event frequencies over time.
   - **Significance**: Provides a temporal perspective on the birth, loss, and transfer of genes between branches, enhancing our understanding of microbial evolution.

### 2. **Marine Microbiome Nutrient Cycling Gene Response (GeneNutrientAnalysis.py)**
   - **Description**: This script investigates how marine microbiomes respond to nutrient-poor conditions by identifying genes correlated with specific environmental stressors.
   - **Functionality**:
     - Analyzes data from TARA Ocean samples.
     - Performs linear regression analyses between selected nutrients (nitrate, nitrogen dioxide, phosphate, and silicon) and gene abundances.
   - **Significance**: Identifies genes linked to nutrient stress, shedding light on microbial adaptation mechanisms during historical atmospheric transitions.
