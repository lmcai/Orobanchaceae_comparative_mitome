#I. MTPT
1. Calculate mtpt proportion for each species: `mtpt_characterization.sh`

This script outputs the number of site in the mito assembly with pt origin (similarity based)

2. Examine bias in pt position for mtpt

Use `mtpt_pt_bias.sh` to blast mtpt to Rehmannia plastome, then `cat` all bed files into one `all_mtpt.bed`. Use `bedtools merge` to count site coverage

```
bedtools coverage -d -a pt.bed -b all_mtpt.bed
```

#II. Repeat
