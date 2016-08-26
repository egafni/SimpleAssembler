#!/usr/bin/env python
import argparse


def iter_fasta_sequences(fasta_path):
    """
    :yields: the sequences in the fasta file
    """
    with open(fasta_path) as fp:
        for i, line in enumerate(fp):
            if i % 2 != 0:
                yield line.strip()


def assemble_contig(fragments):
    """
    Assembles fragments into a contig.  Assumes each fragment overlaps one other fragment by more than half of it's sequence.

    :param list[str] fragments: A list of fragments of type str.
    :return: The assembled contig.
    """
    fragments = list(fragments)  # avoid in place mutation
    while len(fragments) > 1:
        fragment = fragments.pop()

        for candidate_fragment in fragments:
            merged_fragment = assemble_two_fragments(fragment, candidate_fragment)
            if merged_fragment is not None:
                fragments.append(merged_fragment)
                fragments.remove(candidate_fragment)
                break

        assert merged_fragment, 'No fragments overlap with %s' % fragment
    return fragments[0]


def assemble_two_fragments(contig_a, contig_b):
    """
    Assembles two contigs if they overlap by more than half of the bases of at least one of the contigs.

    note:: An alternative is to do a smith waterman and then look for left/right matches on right_seq.
    :param str super_contig: The super contig.
    :param str contig: The Contig.
    :return: The assembled contig if assembly was successful, else None.


    >>> assemble_two_fragments('CT', 'CTA')
    'CTA'
    >>> assemble_two_fragments('CTCT', 'CTA')
    'CTCTA'
    >>> assemble_two_fragments('CT', 'CTAAAA')
    'CTAAAA'
    >>> assemble_two_fragments('CTTT', 'CTTA') is None
    True
    >>> assemble_two_fragments('CLMFWFZV', 'LMFWFZVDHWXUZKT')
    'CLMFWFZVDHWXUZKT'
    """
    contig, super_contig = sorted((contig_a, contig_b), key=len)
    assert len(super_contig) > 1 and len(contig) > 1
    # try to assemble `contig` onto left side of `super_contig`
    i = 0
    while i < len(contig) / 2 + 1:
        if super_contig[:len(contig) - i or None] == contig[i:]:
            return contig[:i] + super_contig
        i += 1

    # try to assemble `contig` onto right side of `super_contig`
    i = 0
    while i < len(contig) / 2 + 1:
        if super_contig[-len(contig) + i:] == contig[:-i or None]:
            return super_contig + contig[-i:]
        i += 1

    return None


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('in_fasta', help="input fasta file")
    args = p.parse_args()

    fasta_fragments = iter_fasta_sequences(args.in_fasta)

    print assemble_contig(fasta_fragments)

