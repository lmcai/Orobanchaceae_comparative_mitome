# MTPT annotation
1. Calculate mtpt proportion for each species: `mtpt_characterization.sh`

This script outputs the number of site in the mito assembly with pt origin (similarity based)

2. Examine bias in pt position for mtpt

Use `mtpt_pt_bias.sh` to blast mtpt to Rehmannia plastome, then `cat` all bed files into one `all_mtpt.bed`. Use `bedtools merge` to count site coverage

```
bedtools coverage -d -a pt.bed -b all_mtpt.bed
```
3. Examine retention pattens in holoparasites

`cat` all holoparasite `*.comparative_mtpt.bed` into `holoparasite_mtpt.bed `. Then use similar method to count coverage:
```
bedtools coverage -d -a pt.bed -b holoparasite_mtpt.bed >holoparasite_pt_mtpt_coverage.bed
```
Identify pt regions with high transferability
```
awk '{if ($5>13) print}' holoparasite_pt_mtpt_coverage.bed | awk 'OFS="\t"{print $1, $4-1, $4}' >temp.bed 
bedtools merge -i temp.bed >holoparasite_mtpt_high_coverage.bed 
rm temp.bed
NC_034308.1	30353	30475
NC_034308.1	42130	42294
NC_034308.1	45432	45571
NC_034308.1	67036	67231
NC_034308.1	100873	101736
NC_034308.1	136491	137354

```

Then use `bedtools intersect` to identify overlap with Rehmannia gene regions

```
bedtools intersect -a Rehmannia_glutinosa.gene.bed -b holoparasite_mtpt_high_coverage.bed
NC_034308.1	101143	101220	trnI-CAU

```
The majority of these universal MTPT is from non-coding regions in plastomes.
