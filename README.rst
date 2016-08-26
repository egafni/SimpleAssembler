Introduction
===============
To sequence the DNA for a given individual we typically fragment each chromosome to many small pieces that can be sequenced in parallel and then re-assemble the sequenced fragments into one long DNA sequence. In this task we ask that you take on a specific subtask of this process.

Challenge
===========

The input to the problem is at most 50 DNA sequences (i.e, the character set is limited to T/C/G/A) whose length does not exceed 1000 characters. The sequences are given in FASTA format (https://en.wikipedia.org/wiki/FASTA_format). These sequences are all different fragments of one chromosome.

The specific set of sequences you will get satisfy a very unique property:  there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length. An example set of input strings is attached.

The output of your program should be this unique sequence that contains each of the given input strings as a substring.

In addition to the code you wrote, we also ask for a README describing your general approach as well as any additional code you wrote to evaluate your solution. We would prefer your code to be written in Python, Go, Scala, Javascript, or Java.



Example input
=============

.. code-block:: bash

    >Frag_56
    ATTAGACCTG
    >Frag_57
    CCTGCCGGAA
    >Frag_58
    AGACCTGCCG
    >Frag_59
    GCCGGAATAC

Example output
===============

.. code-block:: bash

    ATTAGACCTGCCGGAATAC

Solution
=========

The algorithm works by reading all of the input fragments into memory and storing into a List.
A fragment is popped off the List, and the rest of the fragments are iterated over as candidate fragments
for assembly.  If the fragment can be assembled with the candidate fragment (using more than half it's sequence),
the candidate fragment is removed from the List, and the new assembled fragment is added back
in.  This reduces the total size of the List by one (since two fragments have been merged into one), and the
process is repeated until there is only one fragment left.  The final fragment is the original contig.

This is not a very memory efficient algorithm since all fragments are stored in memory.  It's run time is also
not great, and is O(|Fragments|^2) since the worst case scenario is each fragment is compared to all remaining fragments.



Testing
========

.. code-block:: bash

    $ python -m doctest simple_assembler.py
    $ python test_simple_assembler.py
