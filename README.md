This repository describes how version 2 of the *Heliconius melpomene* genome was produced, and includes bespoke code written to generate markers and linkage maps, incorporate haplotype scaffolds and PacBio sequence into the genome and order scaffolds. The bespoke code is provided in the `code` folder for transparency only. All config files mentioned can be found in the `config` folder. Full path names have been omitted.

A preprint about this genome assembly is now [available on bioRxiv](http://www.biorxiv.org/content/early/2015/10/15/029199). The genome itself is available in a full distribution with annotation, maps and other information from  [butterflygenome.org](http://butterflygenome.org/node/4) and as genome and gene sequences from [LepBase v1.0](http://ensembl.lepbase.org/index.html) (Hmel2).

# Generate Markers

VCF files were generated containing SNPs for the Heliconius melpomene mapping cross, containing F0 grandmother, F1 parents and 69 offspring against *Heliconius melpomene* genome version 1.1 (Hmel1-1) (see Methods for details of alignment and SNP calling). Separate VCF files were created for the primary and haplotype scaffolds (`Hmel_cross.Hmel1-1_primaryScaffolds.vcf`, `Hmel_cross.Hmel_haplotype_scaffolds.vcf`). SNPs were converted to markers using the `scaffoldgenome.pl` script, parallelised using `vcf-parallel.pl`:

```
vcf-parallel.pl -s scaffoldgenome.pl -v Hmel_cross.Hmel1-1_primaryScaffolds.vcf -e "-g Hmel_cross_marker_types.txt -r Hmel_cross_marker_types.rms_pval.txt" -b -o Hmel_cross.linkage_map -t 30
```

```
vcf-parallel.pl -s scaffoldgenome.pl -v Hmel_cross.Hmel_haplotype_scaffolds.vcf -e "-g Hmel_cross_marker_types.txt -r Hmel_cross_marker_types.rms_pval.txt" -b -o Hmel_cross.linkage_map -t 15
```

These commands produce an SQLite3 database, `Hmel_cross.linkage_map.db`, containing SNP markers in the `markers` table and scaffold regions for all scaffolds in the `blocks` table.

`Hmel_cross_marker_types.txt` (included in the repository) is a config file describing the cross (parents and sexes, as well as some poor-quality individuals to ignore) and marker types to process (see Table S1). For a particular marker type, the expected parent calls and offspring calls are given, along with a name for the type (Type) and whether recombination can be detected with this marker or not (Recombination, Y/N). Calls are defined using alleles A or B for homozygotes and H for heterozygotes. Parent calls are defined by the order given in the Parents line. For offspring, where multiple calls are possible, they are separated by commas. Multiples of an allele can be given by adding a number before the allele. For example, Intercross calls segregating by Mendelian inheritance in a 1:2:1 ratio are defined as A,2H,B.

`Hmel_cross_marker_types.rms_pval.txt` (included in the repository) is a cache of results for the root mean square (RMS) test, generated using the `generate_rms_distributions.pl` script. This script is run by `scaffoldgenome.pl` if no cached values are provided, but they can be generated separately as follows:

```
generate_rms_distributions.pl -g Hmel_cross_marker_types.txt -s 1000000
```

The `-s` option defines how many samples to draw from the distribution.

Scaffold regions (blocks) defined by `scaffoldgenome.pl` in the `blocks` table in `Hmel_cross.linkage_map.db` are then cleaned using `clean_blocks.pl`, writing new blocks to the table `cleanblocks` (see Methods section "Identification of maternal chromosome prints and paternal markers" for details);

```
clean_blocks.pl -i Hmel_cross.linkage_map.db -o Hmel_cross.linkage_map
```

# Build Linkage Maps


Linkage maps were constructed by the `build_linkage_maps.pl` script, which takes the `Hmel_cross.linkage_map.db` database as input, cleans the markers, separates them into linkage groups, and passes them to `MSTMap` to build maps for each linkage group. `MSTMap.exe` must be in the path. It takes two additional optional files, `known_patterns.txt`, which contains known patterns inferred manually from inspection of SNPs but which were rejected by `scaffoldgenome.pl` due to segregation distortion or rejected by `build_linkage_maps` as they appeared at the ends of linkage groups, and `correct_patterns.txt`, which contains patterns inferred incorrectly by `scaffoldgenome.pl`, to be corrected by `build_linkage_maps.pl`. This script outputs new tables `mapblocks`, `chromosome_map` and `scaffold_map` to the input database.

```
build_linkage_maps.pl -i Hmel_cross.linkage_map.db -k known_patterns.txt -c correct_patterns.txt
```

Many individual scaffold regions were called incorrectly; either they were called with the incorrect marker, or they were not called at all. Many of these regions were identified manually during fixing of misassemblies and ordering of scaffolds. These regions were corrected using `clean_errors.py`, providing the file `block_errors_additions.txt` as input. The maps for several chromosomes were oriented differently to `Hmel1.1` and so were reversed at this stage, which required updating markers for all scaffold regions for these chromosomes. The file `chromosomes_to_reverse.txt` lists the reversed chromosomes. The new scaffold blocks were written to a new SQLite3 database, `Hmel_cross.linkage_map.clean.db`.

```
clean_errors.py -d Hmel_cross.linkage_map.db -e block_errors_additions.txt -r chromosomes_to_reverse.txt  -o Hmel_cross.linkage_map.clean.db
```

# Correct and Merge Draft Genome

The `Hmel1.1` primary and haplotype scaffolds were concatenated into one FASTA file and repeat masked using `RepeatMasker` (see Methods). This FASTA file was then used to create a new version of the genome, `Hmel1-2.fasta`, fixing the misassemblies described in `draft_misassemblies.txt`. The 'revise_genome.py' script fixes these misassemblies, including revisiting scaffolds broken due to misassemblies during assembly of Hmel1.1, re-merging those misassembled scaffolds and re-breaking them using more accurate breakpoints inferred from the new whole genome markers. Scaffolds which were correctly broken and did not need to be remerged were listed in `real_broken_scaffolds.tsv`. The `-p` and `-n` options define the scaffold prefix and scaffold number to begin naming of new merged scaffolds.

```
revise_genome.py -f Hmel1-1_primaryScaffolds.haplotypes.fasta.masked -a Hmel1-1_primary_haplotype_scaffolds.agp -c draft_misassemblies.txt -r real_broken_scaffolds.tsv -o Hmel1-2 -p HE -n 679999
```

This script produces a revised genome `Hmel1-2.fasta` and a TSV file describing the transfer of `Hmel1-1` to `Hmel1-2`, `Hmel1-2.tsv`. The corrected `Hmel1-2` genome, containing primary and haplotype scaffolds, was then collapsed to a haploid version using `batchhm.pl` (see HaploManager repository for details).

```
batchhm.pl -i Hmel1-2.fasta -c heliconius_template -o Hmel1-2_merge -p RV -b bad_merges_Hmel1-2.txt -t Hmel1-2.tsv -a /whale-data/jd626/hmel_genome/Hmel1-1_primary_haplotype.gff
```

# Assemble and Correct PacBio Genome


The PacBio reads `heliconius_melpomene_melpomene.subreads.fastq` were corrected using `PBcR` with the original genome strain data as follows (see Methods for further details). Original genome strain data was downloaded from NCBI (see Experiment and Run accessions in filenames) and saved in folders `SRX124669_Illumina`, `SRX124544_454_shotgun` and `SRX124545_454_3kb`.
```
PBcR -partitions 334 -genomeSize 292000000 -libraryname helmel_genome_strain -s PBcR_genome_strain.spec -fastq heliconius_melpomene_melpomene.subreads.fastq SRX124669_Illumina/SRR424576.frg SRX124544_454_shotgun/SRR424191.frg SRX124544_454_shotgun/SRR424192.frg SRX124544_454_shotgun/SRR424193.frg SRX124544_454_shotgun/SRR424194.frg SRX124544_454_shotgun/SRR424195.frg SRX124544_454_shotgun/SRR424196.frg SRX124544_454_shotgun/SRR424197.frg SRX124544_454_shotgun/SRR424198.frg SRX124544_454_shotgun/SRR424199.frg SRX124544_454_shotgun/SRR424200.frg SRX124544_454_shotgun/SRR424201.frg SRX124544_454_shotgun/SRR424202.frg SRX124544_454_shotgun/SRR424203.frg /disk2/jd626/Hmel_genome_DNAseq/SRX124544_454_shotgun/SRR424207.frg SRX124544_454_shotgun/SRR424208.frg SRX124545_454_3kb/SRR424204.frg SRX124545_454_3kb/SRR424205.frg SRX124545_454_3kb/SRR424206.frg SRX124545_454_3kb/SRR424209.frg assemble=1
```

The PacBio reads were also self-corrected:

```
PBcR -length 200 -threads 50 -genomeSize 292000000 -s PBcR_selfcorrected.spec -l PBcR_selfcorrected fastqFile=heliconius_melpomene_melpomene.subreads.fastq
```
The two read sets were assembled using `FALCON`:

```
$ cat input.fofn
PBcR_selfcorrected.fasta
helmel_genome_strain.fasta
```

```
fc_run.py FALCON.cfg &> FALCON.log &
```

The resulting genome `p_ctg.fa` was masked with `RepeatMasker` and misassemblies identified iteratively through merging with the draft genome and comparing to the linkage map. The final merge used a corrected genome, fixed with `revise_genome.py`:

```
revise_genome.py -f p_ctg.fa.masked -c pacbio_misassemblies.txt -o pacbio_revised -p F -n 100000
```

The PacBio genome was then collapsed to a haploid version with `batchhm.pl`:
```
batchhm.pl -i pacbio_revised.fasta -c heliconius_template -o hm_pacbio_revised -p PF
```


# Merge Draft and PacBio Genomes and Reassemble

The haploid draft genome and haploid PacBio genome were concatenated and then merged with HaploMerger:

```
cat Hmel1-2_merge.fa hm_pacbio_revised_merge.fa > Hmel1-2_pacbio.fa
runhm.pl -i Hmel1-2_pacbio.fa -c heliconius_template -o Hmel1-2_pacbio -p HM -s PB -a Hmel1-2.gff -g -x
```

The `runhm.pl` script generates `Hmel1-2_pacbio_refined.fa`, `Hmel1-2_pacbio_map.db`, `Hmel1-2_pacbio.gff`, `Hmel1-2_pacbio_orig.tsv` and `HMSc_pacbio.nodes` which were then used to order the merged genome using `reassemble_genome.py`:

```
reassemble_genome.py -d Hmel1-2_pacbio_map.db -f Hmel1-2_pacbio_refined.fa -g Hmel1-2_pacbio.gff -o Hmel1-2_pacbio_orig.tsv -p PB -n HMSc_pacbio.nodes -r Hmel2
```

# Transfer Annotation

The `Hmel1` annotation was transferred to `Hmel2` by merging `Hmel1-2` and `Hmel2` with HaploMerger to generate chain files, filtering the chain to retain only between-genome matches, using `CrossMap` to transfer features with the chains, and then using `transfer_features.py` to transfer as many features directly as possible, falling back on the `CrossMap` transfers where necessary.

```
runhm.pl -i Hmel1-2_Hmel2.fa -c heliconius_template -o Hmel1-2_Hmel2 -p H12
filter_chain.py -c Hmel1-2_Hmel2/genome.genomex.result/all.tbest.chain.gz -p Hmel -o Hmel1-2_Hmel2.chain.filtered.gz
CrossMap.py gff Hmel1-2_Hmel2.chain.filtered.gz Hmel1-1_primary_haplotype_mtDNA.gff >Hmel1-2_Hmel2.crossmap.out
transfer_features.py -m Hmel2_transfer_merge.tsv -c Hmel1-2_Hmel2.crossmap.out -f Hmel1-1_primary_haplotype_mtDNA.fa -n Hmel2.fa -o Hmel2
```
