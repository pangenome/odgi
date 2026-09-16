.. _privacy:

##################################
Privacy-preserving variation graphs
##################################

Pangenome graphs summarize the genomes of many individuals, so sharing them openly risks leaking private information about the people whose sequences built them.
The NLnet-funded `privvg <https://privvg.github.io/>`_ project ("privacy-preserving variation graphs") addressed this by developing practical differential privacy models for variation graphs, work that was first prototyped in `vg <https://github.com/vgteam/vg>`_ and has now landed in `odgi` as :ref:`odgi priv`.
:ref:`odgi priv` applies the exponential mechanism to sample shared sub-haplotypes from the graph, and the strongest ε-differential privacy guarantees are obtained by publishing only the FASTA sequences of these sampled paths via :ref:`odgi paths` **-f**.

Acknowledgments
===============

Development of the differential privacy model and :ref:`odgi priv` was funded through the `NGI0 Discovery Fund <https://nlnet.nl/discovery>`_, a fund established by the `NLnet Foundation <https://nlnet.nl/project/VariationGraph/>`_ with financial support from the European Commission's `Next Generation Internet <https://ngi.eu/>`_ programme under grant agreement No 825322.
