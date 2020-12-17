Hidden biodiversity: search for uncultured protists in metagenomes

Background: metagenomics allows to reveal the hidden diversity of uncultivated organisms, including unicellular creatures (protists). There are groups of protists that are quite difficult to isolate from the environment and even more difficult to cultivate. The discovery of hidden groups of protists and previously undescribed representatives of known groups will help us learn more about their biodiversity, evolution, etc.

Project objective: search for a hidden diversity of mushroom-like protists and amoebas, new supergroups of protists in metagenomes.

Project tasks:
  1. Get metagenomes from NCBI database
  2. Predict and select marker genes from raw contigs
  3. Make taxonomic annotation
  4. Perform phylogenetic analysis
  
Workflow:
  1. Get a total stats table about the number and volume of metagenomes from different habitats using BioPython;
  2. Download metagenomes of different habitats using BioPython;
  3. Predict genes in rRNA assemblies using Barrnap (https://github.com/tseemann/barrnap);
  4. Conduct taxonomic annotation of the obtained sequences using SINA (https://github.com/epruesse/SINA), get predicted and annotated .fasta;
  5. Selection and filtering of 18S and 28S rRNA genes using command-line tools;
  6. Perform phylogenetic analysis using FastTree (http://www.microbesonline.org/fasttree/);
  7. Make beautiful and easy-to-read tree plots using ETE Toolkit for Python (http://etetoolkit.org/);
  8. Perform BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi) to verify our results.
  You can find more info about tools using given links.
  You can find all intermediate results in Github repository.
