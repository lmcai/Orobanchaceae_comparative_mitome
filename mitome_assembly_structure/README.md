# I. MTPT
1. Calculate mtpt proportion for each species: `mtpt_characterization.sh`

This script outputs the number of site in the mito assembly with pt origin (similarity based)

2. Examine bias in pt position for mtpt

Use `mtpt_pt_bias.sh` to blast mtpt to Rehmannia plastome, then `cat` all bed files into one `all_mtpt.bed`. Use `bedtools merge` to count site coverage

```
bedtools coverage -d -a pt.bed -b all_mtpt.bed
```

# II. Repeat
Use `ROUSFinder2.0.py` to identify repeat with size >15 bp. The batch file to execute it across all species is `ROUSFinder_batch.sh`

ROUSFinder will incorrectly identify long single-copy sequences as repeats somehow. Use `summarize_ROUSfinder_result.sh` to remove these and calculate the total size of repeats.