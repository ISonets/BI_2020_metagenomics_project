# Hidden biodiversity: search for uncultured protists in metagenomes

Background: metagenomics allows to reveal the hidden diversity of uncultivated organisms, including unicellular creatures (protists). There are groups of protists that are quite difficult to isolate from the environment and even more difficult to cultivate. The discovery of hidden groups of protists and previously undescribed representatives of known groups will help us learn more about their biodiversity, evolution, etc.

Project objective: search for a hidden diversity of mushroom-like protists and amoebas, new supergroups of protists in metagenomes.

## Project tasks:
  1. Get metagenomes from NCBI database
  2. Predict and select marker genes from raw contigs
  3. Make taxonomic annotation
  4. Perform phylogenetic analysis
  
## Workflow:
  1. Get a total stats table about the number and volume of metagenomes from different habitats using BioPython;
  2. Download metagenomes of different habitats using BioPython;
  3. Predict genes in rRNA assemblies using Barrnap (https://github.com/tseemann/barrnap);
  4. Selection and filtering of 18S rRNA genes using command-line tools;
  5. Conduct taxonomic annotation of the obtained sequences using SINA (https://github.com/epruesse/SINA) and SILVA SSU Ref 99 database(https://www.arb-silva.de/download/arb-files/), get predicted and annotated .fasta;
  6. Filtering SINA results from undesired fungi using command-line tools;
  7. Removing alignment gaps after SINA using BioPython;
  8. Realignment of degapped sequences using MAFFT online server (https://mafft.cbrc.jp/alignment/server/);
  9. Perform phylogenetic analysis using FastTree (http://www.microbesonline.org/fasttree/);
  10. Make beautiful and easy-to-read tree plots using ETE Toolkit for Python (http://etetoolkit.org/);
  11. Perform BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi) to verify our results.
  
  To sum up, we developed new pipeline that allows to find new species. This pipeline is unique and first to be developed to seacrh new species.
  
  Pipeline was tested on previously found microsporidia *Enterospora canceri* before lauhching.
  
  You can find more info about tools using given links.
  You can find all intermediate results in Github repository, as also a thesis.
  If you have any questions, please contact ignatsonets@gmail.com and/or vanypyankov@gmail.com
  
  Thank you for your attention!
  
  
# Pipeline with commands
## 1.Script for search metagenomes.py
  ```
  python3 Script for search metagenomes.py -o List_with_metagenomes
  ```
  
  You can see the instructions for use in the "python scripts" folder

## 2. Download_assembly.py
  ```
  python3 Download_assembly.py -t search_metagenome -i input_file.csv -d path_to_dir
  ```
  
  You can see the instructions for use in the "python scripts" folder
  
## 3. Barrnap
  Barrnap was installed accroding to the manual in Anaconda.
  ```
  for i in raw_data/*.fna;  do barrnap --threads 4 --kingdom arc --kingdom euk --outseq ${i%.*}.fasta < "$i"; done
  ```

  ### 3.1. Extracting only 18S sequences
    18S:
    ```
    for i in *.fasta; do samtools faidx ${i}; done
    for i in *.fasta.fai; do grep 18S ${i} | cut -f 1 > ${i%.**}_18S_names; done
    for i in *.fasta; do more ${i}_18S_names | xargs samtools faidx ${i} > ${i%.*}_18S_only.fasta; done
    ```
    
## 4. SINA
  SINA was lauched from pre-compiled tarball archive.
```
./sina -i ../18S_all_add.fasta -o 18S_all_aligned.fasta --db ../../../SILVA_138.1_SSURef_NR99_12_06_20_opt.arb  --search --meta-fmt CSV --lca-fields tax_slv
```
  
## 5. Results filtreing

  ### 5.1.Removing Dikarya/Mucoromycota from aligned 18S sequences
  ``` samtools faidx 18S_all_aligned.fasta
  grep Fungi 18S_all_aligned.csv | cut -f 1 > 18S_fungi.csv
  sed -i '/\b\(Dikarya\|Mucoromycota\)\b/d' 18S_fungi.csv 
  cat 18S_fungi.csv | cut -d \, -f 1 > 18S_fungi_namesonly
  more 18S_fungi_namesonly | xargs samtools faidx 18S_all_aligned.fasta > 18S_desired_fungi.fasta 
  ```
    
  ### 5.2. Adding additional annotated fungi SSU set
  ```
  cat 18S_desired_fungi.fasta ssu.fst > 18S_allweneed_v2.fasta
  ``` 
  
  ### 5.3. Extracting 'Unclassified' sequences
  ```
  grep Unclassified 18S_all_aligned.csv | cut -f 1 > 18S_unclassified.csv
  cat 18S_unclassified.csv | cut -f 1 > 18S_unclassified_namesonly
  more 18S_unclassified_namesonly | xargs samtools faidx 18S_all_aligned.fasta > 18S_unclassified.fasta
  ```

## 6. Removing alignment gaps using Python script(available in GitHub repository)
  ```
  python3 removegap.py 18S_desired_fungi.fasta 18S_desired_fungi_nogaps.fasta
  python3 removegap.py 18S_allweneed_v2.fasta 18S_allweneed_v2_nogaps.fasta
  python3 removegap.py 18S_unclassified.fasta 18S_unclassified_nogaps.fasta
  ```
  You can see the instructions for use in the "python scripts" folder
  
## 7. Realignment using MAFFT online version

  Parameters:
  ```
  -direction of nucleotide sequences: same as input
  -output order: aligned
  -strategy: auto (depends on data size)
  -align unrelated sequences: try to align gappy regions anyway
  -scoring matrix for nuclotide sequences: 200PAM / k = 2
  -gap opening penalty: 1.53
  -offset value: 0
  -guide tree: default
  -mafft-homologs: use UniRef50
   ```
   
   MAFFT was performed for 3 sequences:
   
   ```
  -18S_desired_fungi_nogaps.fasta
  -18S_allweneed_v2_nogaps.fasta
  -18S_unclassified_nogaps.fasta
   ```
  
   Results were renamed:
   
   ```
  -18S_allweneed_v2_nogaps_aligned.fasta
  -18S_desired_fungi_nogaps_aligned.fasta
  -18S_unclassified_nogaps_aligned.fasta
   ```
  
## 8. FastTree
  FastTree was installed according to the manual in Anaconda.
  
  ```
  fasttree -nt 18S_allweneed_v2_nogaps_aligned.fasta > 18S_allweneed_v2.tre
  fasttree -nt 18S_desired_fungi_nogaps_aligned.fasta > 18S_desired_fungi.tre
  fasttree -nt 18S_unclassified_nogaps_aligned.fasta > 18S_unclassified.tre
  ```

## 9. ete3
  Trees have been built using package ete3 in Python.
  
  You can see the instructions for use in the "python scripts" folder

