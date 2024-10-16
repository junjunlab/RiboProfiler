# RiboProfiler update news

---

## RiboProfiler version 0.1.0 update

- Add **get_track_df** for gene track plot.
- Add **metagene_plot** for metagene analysis.
- Add **triAmino_motif_plot** for enriched motif logo plot.

## RiboProfiler version 0.1.1 update

- Add **get_pausing_site** to detect pausing site across transcipt. The method
is using PausePred described with a sliding window. More details please refer to
[PausePred and Rfeet: webtools for inferring ribosome pauses and visualizing footprint density from ribosome profiling data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140459/)
- Add 5'end or 3'end assginment of reads for **pre_qc_data** function.
- Add in-frame information for aggregaton plot.

## RiboProfiler version 0.1.2 update

- Add **peptide_motif_score2** to calculate average pausec score according to 
**The ubiquitin conjugase Rad6 mediates ribosome pausing during oxidative stress** reference.
- Add **show_cds_region_only** parameter for **track_plot**.
- Add **calculatePolarity2** to calculate polarity score for nomalized data.
- Add **frame_col** parameter for **metagene_plot**.

## RiboProfiler version 0.2.0 update
- Constructed into a ribosomeObj Object, methods can be applied on ribosomeObj Object.
- Implementary with RiboSeqTools package functions to deal with selective ribosome profilling data.
- Add **ribo_bam_to_bw** to transform ribo bamfiles to bigwig.
- Add **XAM_version** for **ribo_bam_to_bw** and **pre_qc_data* to choose suitable XAM version to analysis.

## RiboProfiler version 0.2.2 update
- Add rpm and average normalization for **get_nomalized_counts**.
- Add **periodicity_check** function which applied Discrete Fourier transform method to check periodicity.
- Add **single_gene_plot** function to visualize single gene profile of ribosome data.

## RiboProfiler version 0.2.3 update
- Add **codonRelativeDist.py** to calculate relative codon occupancy along transcript.
- Add **bam_to_bw* to convert bam files to bigwigs.
- Add **get_counts_from_bam** to extract counts for ribo or rna bam files.
