# <center> **BinderSpace Tutorial**

Data information:
The input data is peptide sequence data that bind to Bcl-2 family of proteins, the dataset id is: GSE118147. For more details, refer to this link [GSE118147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118147).


 

*   positive sequence: GSM3319530_m1r.csv' (peptide sequences that bind to Mcl-1)
*    negative sequence: GSM3319530_x100r.csv and GSM3319530_x1.csv (peptide sequences that bind to Bcl-xL)



>




## 1. Installing

> Downloading the motifSearch.py, VisualizeMotif.ipynb and put it in the same folder with your
data.



## 2. Run motifSearch.py

1.   Prepare input file:
#####The input file is a CSV file contains DNA or peptide sequences. The sequences column header should be “Sequence”. Here we use 22_nt peptide sequences (3480 rows) from Z-keating dataset, and we also get 1470 negative sequences from peptide that bind to other proteins (pos_seqs.csv, neg_seqs.csv). The sequences should have same length and the sequences should not contain letters beyond 22 amino acids or 4 DNA bases.
<center>
<table>
<tr><th> positive sequences </th><th> negative sequences</th></tr>

<tr><td>

Sequence              |            
----------------------|       
GRTEHQTGQELIRIGDEIQAYY      
GRTEHQIGQELIRIGDEDDAYY
GRTEHQIGQELARIGDEDNAYY
GRTEHQVSQELVRIGDEIDAYY
GRTEHHVSQELVRIGDETQAYY

</td><td>

Sequence |
---------|
GRTERQAVQELARIGDEVDAYY
GRTERHVVQELTRIGDEHHAYY
GRTELQIGQELKRIGDEFNAYY
GRTEDDIGQELVRIGDESNAYY
GRTEAQVVQELVRIGDEYQAYY

</td></tr> </table>





















2.   Run motif search:
#####Parameter setting: no gap, motif minimal length is 8, minmal frequency is 0.5% * total number of positive sequences, use 20 processors 
#####command line:
python MotifSearch.py -i pos_seqs.csv -n neg_seqs.csv -o keating -m 8 -f 0.005 -g 2 -a 2 






3.   Output:
#####The output keating_length*.csv files contains 8 columns: motifs, frequency of this motif occurs in positive sequences (positive_occurrence %), frequency of this motif occurs in negative sequences, percentage of the positive occurrence, percentage of the negative occurrence, index of positive sequences, index of negative sequences


4.   Other command options: run “python MotifSearch.py“ to see the parameter options



> *  -i posfile -- the file with positive sequences (REQUIRED)
*   -n negfile -- the file with negative sequences
*  -o motif_file  --the name of the output file (default: 'motifs')
*   -f minfreq -- the fraction of the minimal motif occurence in positive sequences (default: f = 0.001,
minfreq = f * (the total number of positive sequences))
*   -m minlength  --the minimal length of the motifs (default: 3)
*   -l maxlength  --the maximal length of the motifs (default: the length of the positive sequences)
*   -c category --choose from 'dna' or 'amino_acids',(default: 'amino_acids')
*   -g maxgaps --maximal number of gaps (default: 0)
*  -a maxgaplength --maximal gap length (default: 1)
*   -p processor --number of processors used for running the program (default: 20) italicized text

##3. PCA/tSNE visualization for one of the motifs (this part only works for dataset that has negative control)

Open **VisualizeMotif.ipynb** by jupyter notenook. Choose one motif file (‘keating_length 11.csv’) and from that motif file choose one motif ('L*RIGDE*DAY'), and we will check the locations of those sequences that contain that motif in PCA and t-SNE space. 
<br>
**Input:**

> pos_seq = ‘pos_seqs.csv’
<br>
neg_seq = ‘neg_seqs.csv’
<br>
motif_file = ‘keating_length 11.csv’
<br>
motif = 'L*RIGDE*DAY'
<br>
bases_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

**Output:**

> ![picture](https://drive.google.com/uc?id=1NTJWKAbM4UbWfiKNFomPj37a6Bv_ROZY)





```python

```
