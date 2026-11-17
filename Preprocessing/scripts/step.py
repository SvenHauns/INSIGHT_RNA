"""Snakemake script adapter (runs inside the rule's Conda environment)."""
from steps import prepare_reference, index_reference, trim, align, finalize

s = snakemake
settings = dict(s.params.settings)
action = s.params.action
if action == "reference":
    prepare_reference(s.input.fasta, s.input.gtf, s.output.fasta, s.output.report)
elif action == "index":
    index_reference(s.input.fasta, s.input.gtf, s.output.index, s.threads, settings, s.log[0])
elif action == "trim":
    trim(s.input.r1, s.input.r2, s.output.r1, s.output.r2, s.output.json, s.output.html, s.threads, settings, s.log[0])
elif action == "align":
    align(list(s.input.r1), list(s.input.r2), s.input.index, s.input.fasta, s.wildcards.sample,
          s.output.bam, s.output.star_log, s.threads, settings, s.log[0])
elif action == "finalize":
    finalize(s.input.bam, s.input.fasta, s.wildcards.sample, s.output.bam, s.output.flagstat,
             s.output.idxstats, s.output.provenance, s.threads, settings, list(s.params.runs), s.log[0])
else:
    raise ValueError(f"Unknown preprocessing action: {action}")
