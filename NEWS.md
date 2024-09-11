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
