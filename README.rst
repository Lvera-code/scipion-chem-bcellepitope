================================
B-Cell Epitope Prediction scipion plugin
================================

Scipion framework plugin wrapping the B-Cell-Epitope-Prediction pipeline
(https://github.com/Lvera-code/B-Cell-Epitope-Prediction): local BepiPred-3.0
and EpiDope antigenicity prediction, an annotated BepiPred/EpiDope consensus,
a BLASTp immune tolerance filter, and NetMHCIIpan-4.3 T-helper promiscuity
prediction.

All 5 phases of the original pipeline are implemented as independent,
chainable protocols: ``ProtBCellEpitopeBepiPredPredict``,
``ProtBCellEpitopeEpiDopePredict``, ``ProtBCellEpitopeConsensus``,
``ProtBCellEpitopeBlastFilter`` and ``ProtBCellEpitopeNetMHCIIpan``.

None of BepiPred-3.0, EpiDope, NetMHCIIpan-4.3 or the BLAST human proteome
database are bundled with this plugin: they must each be downloaded/built
separately (see below) and pointed to via ``scipion.conf``.

================================
Download BepiPred-3.0
================================

BepiPred-3.0 is **academic-use only** software (DTU Health Tech). Request it
from:
https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=bepipred&version=3.0&packageversion=3.0b&platform=src

The request form requires a valid **institutional/academic email address**
(a personal Gmail/Outlook address will be rejected) - use your university or
research center email. The download link is emailed to that address.

Unzip it and create a dedicated Python environment for it (BepiPred-3.0 pins
old dependency versions - torch==1.12.0, numpy==1.20.2 - that clash with
Scipion's own environment). Then, in ``scipion.conf``, set:

.. code-block::

      BEPIPRED_HOME = <path to the BepiPred3_src folder>
      BEPIPRED_PYTHON_BIN = <path to the python binary inside that dedicated env>

================================
Download EpiDope
================================

EpiDope is **open-source (MIT license)** - unlike BepiPred-3.0 and
NetMHCIIpan-4.3, it does **not** require an academic request or institutional
email. Just clone it:
https://github.com/rnajena/EpiDope

Install it in a dedicated conda environment using the official
``epidope.yml`` (do NOT install its dependencies by hand: the resolution of
that old stack - Python 3.6, TensorFlow 1.13, ELMo/AllenNLP - is fragile):

.. code-block::

      git clone https://github.com/rnajena/EpiDope.git /tmp/EpiDope
      conda env create -f /tmp/EpiDope/epidope.yml -p <path to your EpiDope conda prefix>

Then, in ``scipion.conf``, set:

.. code-block::

      EPIDOPE_HOME = <path to that EpiDope conda prefix>

================================
Install NCBI BLAST+ and the human proteome database
================================

NCBI BLAST+ itself has no license restriction, but this plugin needs a local,
pre-indexed BLAST database of the human proteome (used to filter out
epitope candidates with high homology to human proteins). Download the human
proteome FASTA (e.g. from UniProt, reference proteome UP000005640) and index
it locally:

.. code-block::

      makeblastdb -in human_proteome.fasta -dbtype prot -out human_proteome_db

Then, in ``scipion.conf``, set:

.. code-block::

      BLASTP_BIN = <path to the blastp binary, or just "blastp" if it is on PATH>
      BLAST_HUMAN_DB = <path prefix to the indexed database, e.g. .../human_proteome_db>

================================
Download NetMHCIIpan-4.3
================================

NetMHCIIpan-4.3 is **academic-use only** software (DTU Health Tech). Request
it from:
https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/

As with BepiPred-3.0, the request form requires a valid **institutional/
academic email address**. Unzip it, then edit the ``NMHOME`` line inside the
``netMHCIIpan`` wrapper script with the absolute installation path (a manual
step required by DTU's own install instructions). Then, in ``scipion.conf``,
set:

.. code-block::

      NETMHCIIPAN_HOME = <path to the netMHCIIpan-4.3 folder>

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-bcellepitope.git
            cd scipion-chem-bcellepitope
            scipion3 installp -p . --devel
