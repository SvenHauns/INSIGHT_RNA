"""One ordered target stream for coverage lists and sequence-extraction BEDs."""


def target_records(records):
    # Preserve biological transcript identities even when coordinates coincide.
    for group in records.values():
        for record in group:
            starts, ends = record[0][1]
            if len(starts) != len(ends) or any(s < 0 or e <= s for s, e in zip(starts, ends)):
                raise ValueError('Invalid BED0 target intervals')
            if sum(e - s for s, e in zip(starts, ends)) > 30:
                yield record
