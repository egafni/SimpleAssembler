"""
Simulate fragmenting a contig

.. note::
    This simulation does not *guarantee* unique overlaps, which would break the assembly,
    so there is a chance this test will fail incorrectly, although it is pretty unlikely.

"""
from simple_assembler import assemble_contig
import random


def random_sequence(size, chars='ACTG'):
    return ''.join(random.choice(chars) for _ in range(size))


def create_contig_and_fragments(contig, overlap_size, fragment_size):
    """
    Creates a contig and overlapping fragments

    :param str contig: original sequence to create test data from
    :param int overlap_size: number of bases fragments should overlap
    :param int fragment_size: length of bases

    note:: returned contig is probably going to be smaller than the input contig so that the
    last fragment isn't too short
    """
    assert overlap_size < fragment_size
    assert fragment_size < len(contig)
    step_size = fragment_size - overlap_size

    fragments = []
    i = 0
    while i + fragment_size <= len(contig):
        fragments.append(contig[i: i + fragment_size])
        i += step_size

    random.shuffle(fragments)
    return contig[:i - step_size + fragment_size], fragments


def test_assemble_contig():
    contig_scaffold = random_sequence(1000)
    for overlap_size, fragment_size in [(9, 16), (10, 15), (13, 20), (13, 25)]:
        # Note you need enough separation between overlap_size and fragment_size, or
        # it will be too likely for there to be multiple ways to scaffold
        true_contig, fragments = create_contig_and_fragments(contig_scaffold, overlap_size, fragment_size)
        assert assemble_contig(fragments) == true_contig

test_assemble_contig()

print 'All tests passed!'
