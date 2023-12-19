# ResinGlycoside_Pipeline
 Scripts for processing LCMS data and annotating resin glycoside like peaks

## Requirements
1. Unix workspace
2. Python
3. Pandas, NumPy

Before using these scripts, obtain exported data from MS-DIAL as described in the Methods paper. You need two types of files from your side:
1. Peak area files
2. MGF files.
These can be individual samples or aligned files.

## How to run the scripts

**Step 1:** Download this entire folder to your Unix workspace. This folder has all the demo files and standard input files needed for the scripts to work. 

**Step 2:** Run the following command to receive instructions for running the script. Once you run the first line, it will output the instructions block.

```
python parse_MSDIAL_wrapper.py

---------- PRE-SCRIPT INSTRUCTIONS -----------
* Make four sub-folders called 0_RgFrags, 1_PeakArea, 2_MGF, 3_predictMotifs
       mkdir 0_RgFrags 1_PeakArea 2_MGF 3_predictMotifs
* Put mz132, mz146, mz162 frag pairs in 0_RgFrags
* Put all Area files in 1_PeakArea
* Put all GnpsMgf files in 2_MGF
* Put motifs.groups in 3_predictMotifs
* Ensure all scripts are in this working directory
---------- REQUIRED ARGUMENTS -----------
-filelist: File containing a list of files to process
---------- INSTRUCTIONS FOR HOW TO FORMAT -filelist -----------
      Tab-delimited file should have following columns
#PeakAreaFile-name      Species-name    Samples-index   Blank-indexPost-curation-index     Fillper-index   S/N-index       MSMS-index  GnpsMGF-file
      Columns should start with index0 and should be csv
---------------------------------------------------------------

```

**Step 3:** Run a demo version with pre-populated filelist.tab using the following command

```
python parse_MSDIAL_wrapper.py -filelist filelist.tab
```

**Step 4:** Ensure that the script ran properly. If everything is correct, the last file that the wrapper makes is **filelist.tab.nid1.id2.tab**. Open that file and ensure it is populated. Also check the **filelist.tab.*.stats** files to ensure they are populated.

**Step 5:** If everything looks good, populate the filelist.tab file as per instructed in the help menu of the above script, as per your sample organization. If you want, you can rename filelist.tab to anything else.

**Step 6:** Run the wrapper once again with the correct file name

**Step 7:** Sit back and relax

## OUTPUTS
This script produces the following output files (if your input filelist name was filelist.tab)

| Filename | Location   | Description |
| -------- | ------------ | ------------ |
| filelist.tab.mz132.combined.stats | Main folder | # of RG-like peaks containing the 132 fragment (pentose) |
| filelist.tab.mz146.combined.stats | Main folder | # of RG-like peaks containing the 146 fragment (deoxyhexose) |
| filelist.tab.mz162.combined.stats | Main folder | # of RG-like peaks containing the 162 fragment (hexose) |
| filelist.tab.nid1.id2.tab | Main folder | A well-formatted file containing annotations of individual MS/MS peaks in each RG-like peaks in the alignment file |
| <fn>.mz132* and mz146* and mz162* | 1_PeakArea | Intermediate files created while parsing the peak area files |
| <fn>.mz132* and mz146* and mz162* | 2_MGF | Intermediate files created while parsing the mgf files |

* The stats files can be combined to obtain a comprehensive view of RG-like peaks in each sample/series/species
* The *.nid1.id2.tab file can be separated into different files of each sample/series/species. For now, all alignment files listed in filelist.tab are included in this file.

For further questions, email Gaurav Moghe

