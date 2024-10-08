# NXF - TGS Oxford Nanopore pipeline
Process raw Oxford Nanopore data, generate html report, and do de novo assembly (plasmid, amplicon, bacterial genome)

## About
Nextflow pipeline for processing raw Oxford Nanopore run data - merge/rename files, generate read reports, run ONT workflows 
(wf-clone-validation, wf-bacterial-genomes or wf-amplicon). Raw reads are mapped back to the assemblies and the alignments are visualised in IGV html reports.


## Running

Only `nextflow` and `docker` are required. The required inputs for the pipeline are:
- path to samplesheet - path to a csv or excel file with columns `sample`, `barcode` and `user`. Other columns can also be there.
- path to raw ONT data (`path/to/fastq_pass`)
  
The pipelines runs the following steps:
- `merge_reads` - merge all reads belonging to a barcode and rename according to the provided samplesheet
- `report` - generate a `csv` and `html` reports about read counts, quality etc. One report per user is generated.
- `assembly` - perform assembly using one of the epi2me pipelines wf-clone-validation, wf-bacterial-genomes or wf-amplicon (default wf-clone-validation)
- `mapping` and `IGV` - map raw reads to the assemblies and generate html IGV reports (one report per sample).

```bash
nextflow run angelovangel/nxf-tgs --samplesheet path/to/samplesheet.csv --fastq path/to/fastq_pass
```

It is possible to run only `merge_reads` or `merge_reads` + `report`. For this, use the `-entry` parameter (note the single dash) in nextflow like so:

```bash
nextflow run angelovangel/nxf-tgs --samplesheet path/to/samplesheet.csv --fastq path/to/fastq_pass -entry report
```
This will run `merge_reads` + `report` and no asembly and mapping...

## Run in EPI2ME

You can import and run the pipeline using GUI in [EPI2ME](https://labs.epi2me.io/downloads/). There is also a test dataset to try if all works - after importing select 'Options' --> 'Run demo analysis'. 