## Materials

This repository contains additional matrials for the following preprint:

*Micom: metagenome-scale modeling to infer metabolic interactions in the microbiota.*<br>
Christian Diener and Osbaldo Resendis-Antonio<br>
https://doi.org/10.1101/361907

All scripts are provided as Python 3 code. libraries required to reproduce
figures and data processing csan be installed using the provided requirements
file:

```bash
pip install --user -r requirements.txt
```

### Subfolders

See each subfolder for detailed desription of contained files.

[*figures*](/figures)<br>
Scripts to reproduce the figures in the publication.

[*results*](/results)<br>
Intermediate results required for figures.

[*workflows*](/workflows)<br>
Long running parallel analyses that produce the results provided in the
previous folder.

### Scripts

**genera.py**<br>
Reads the genus-level abundances and combines them with the AGORA models
for each genus.<br>
Output: `genera.csv`

**taxa_stats.py**<br>
Counts the fraction of reads assigned to each taxonomic rank and the fraction
of reads that has at least one represenative in the AGORA models.<br>
Output: a text table

### Files

**recent.csv**
Provides a list of samples and the corresponding [SRA](https://www.ncbi.nlm.nih.gov/sra) abundances.




