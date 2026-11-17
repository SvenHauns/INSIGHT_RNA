"""Quality-filtered fragment observations shared by selection and profiling."""
from dms_common import DEFAULT_MIN_BASE_QUALITY, mutation_counters


def observations(read, start, end, min_base_quality):
    if any(getattr(read, flag, False) for flag in
           ('is_unmapped', 'is_secondary', 'is_supplementary', 'is_qcfail')):
        return {}
    sequence = read.query_sequence
    qualities = read.query_qualities
    if sequence is None or qualities is None:
        return {}
    result = {}
    for query, position, reference in read.get_aligned_pairs(matches_only=True, with_seq=True):
        if not start <= position < end:
            continue
        observed = sequence[query].upper()
        reference = reference.upper()
        quality = qualities[query]
        if quality >= min_base_quality and observed in 'ACGT' and reference in 'ACGT':
            result[position] = (reference, observed, quality)
    return result


def merge_mates(first, second):
    """Use higher-quality overlapping calls; discard equal-quality disagreements."""
    merged = dict(first)
    for position, call in second.items():
        previous = merged.get(position)
        if previous is None or call[2] > previous[2]:
            merged[position] = call
        elif call[2] == previous[2] and call[:2] != previous[:2]:
            del merged[position]
    return merged


def fragment_calls(bam, chrom, start, end, min_base_quality):
    pending = {}
    for read in bam.fetch(chrom, start=start, stop=end):
        if any(getattr(read, flag, False) for flag in
               ('is_unmapped', 'is_secondary', 'is_supplementary', 'is_qcfail')):
            continue
        calls = observations(read, start, end, min_base_quality)
        if not read.is_paired or not (read.is_read1 or read.is_read2):
            yield calls
            continue
        group = read.get_tag('RG') if read.has_tag('RG') else ''
        key = (group, read.query_name)
        mate = 1 if read.is_read1 else 2
        previous = pending.pop(key, None)
        if previous is None:
            pending[key] = (mate, calls)
        elif previous[0] != mate:
            yield merge_mates(previous[1], calls)
        else:
            # Repeated names for the same mate are not independent templates.
            raise ValueError(f'Duplicate primary mate for read group/name {key}')
    for _, calls in pending.values():
        yield calls


def count_interval(bam, chrom, start, end, strand, min_base_quality=DEFAULT_MIN_BASE_QUALITY):
    if start < 0 or end <= start or strand not in {'+', '-'}:
        raise ValueError('Expected a positive BED0 interval and explicit strand')
    coverage, mutations = [0] * (end - start), [0] * (end - start)
    stable, changed = mutation_counters()
    complement = str.maketrans('ACGT', 'TGCA')
    for calls in fragment_calls(bam, chrom, start, end, min_base_quality):
        for position, (reference, observed, _) in calls.items():
            index = position - start
            coverage[index] += 1
            if strand == '-':
                reference, observed = reference.translate(complement), observed.translate(complement)
            if reference == observed:
                stable[reference] += 1
            else:
                mutations[index] += 1
                changed[reference + '->' + observed] += 1
    return coverage, mutations, stable, changed
