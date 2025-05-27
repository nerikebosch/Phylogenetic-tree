def merge_alignment(prev_center, new_center, new_seq, aligned_others):
    """
        Merge a new sequence into an existing MSA using a new center sequence.

        Args:
            prev_center (str): Previously aligned central sequence.
            new_center (str): Newly aligned central sequence (with gaps).
            new_seq (str): Newly aligned sequence corresponding to new_center.
            aligned_others (list): List of other sequences already aligned to prev_center.

        Returns:
            tuple: (new_center, realigned_new_seq, updated_others)
    """

    # Project both the previous center and new_seq onto new_center
    updated_prev = project_onto_master(new_center, prev_center)
    updated_new = project_onto_master(new_center, new_seq)

    # Project all previously aligned sequences (aligned to prev_center) onto new_center
    updated_others = [project_onto_master(new_center, s) for s in aligned_others]

    # Insert the realigned previous center at the beginning of others (it’s part of the alignment)
    updated_others.append(updated_new)


    return new_center, updated_new, updated_others

def project_onto_master(master: str, seq):
    """
        Project an unaligned sequence onto a gapped master sequence.

        Inserts gaps into `seq` where `master` has them to maintain column consistency.

        Args:
            master (str): Gapped reference sequence.
            seq (str): Sequence to align with the master.

        Returns:
            str: Gapped version of `seq` aligned to `master`.
    """

    out = []
    p = 0
    for c in master:
        if c == "-":

            if len(master) > len(seq):
                out.append("-")
            else:
                continue
        else:
            if p < len(seq):
                out.append(seq[p])
                p += 1

    while p < len(seq):
        out.append(seq[p])
        p += 1

    return "".join(out)