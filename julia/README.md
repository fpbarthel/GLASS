# SubClonalSelection Pipeline

## Files used
|File|Description|
|---|---|
|GLASS_genotypes.csv| Contains information for every mutation called in a sample.|
|GLASS\_genotype\_comparison| Contains shared/private status for every mutation called in a tumour pair, including whether it is shared or private and also the TITAN copy number estimate for the position.|
|GLASS\_genotype\_comparison\_extracted.tsv| The extracted GLASS\_genotype\_comparison file.|
|silver_set.csv| Contains a list of paired samples to include in the analysis.|
|titanparams_synapse.tsv| Contains TITAN ploidy and purity estimates for each sample.|

## Pipeline

### 1. Extract necessary columns from the genotypes comparisons file

(The position column counts as 2 due to the comma).

```
cat GLASS_genotype_comparison.tsv | cut -d, -f 1,9,10,12,21,22 | tr "," "\\t" > GLASS_genotype_comparison_extracted.tsv
```

### 2. Run extrac_vafs.py 

```
python extrac_vafs.py -c GLASS_genotype_comparison_extracted.tsv -g GLASS_genotypes.csv -s silver_set.csv -t titanparams_synapse.tsv -o ./
```

### 3. Add minimum VAF column

Manualy add a column for minimum VAF to the metadata file output from extrac_vafs.py by inspecting the histogram outputs for each sample and choosing the VAF for the highest point of the left most peak.

### 4. Run the analysis through qsubsec
See https://www.ncbi.nlm.nih.gov/pubmed/26635140 for details on qsubsec.

```
qsubsec subclonalselection.qsubsec subclonalselection.tff -s
```

### 5. Subsample inputs 

For runs that don't finish within 48h due to large numbers of mutations, subsample their VAF inputs and rerun.

```
shuf -n 20000 XXX.txt > XXX.txt
``` 

### 6. Remove runs with high error in the model results

Remove any runs with "New ϵ is within 7.0% of previous population, stop ABC SMC" warning in logs.