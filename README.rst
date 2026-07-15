================================
B-Cell Epitope Prediction scipion plugin
================================

Scipion framework plugin wrapping the B-Cell-Epitope-Prediction pipeline
(https://github.com/Lvera-code/B-Cell-Epitope-Prediction): local BepiPred-3.0
and EpiDope antigenicity prediction, BLASTp immune tolerance filtering and
NetMHCIIpan-4.3 T-helper promiscuity prediction.

**Work in progress.** Currently implemented: BepiPred-3.0
(``ProtBCellEpitopeBepiPredPredict``) and EpiDope
(``ProtBCellEpitopeEpiDopePredict``) antigenicity prediction.

================================
Download BepiPred-3.0
================================

BepiPred-3.0 is academic-use only software. Download it from
https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=bepipred&version=3.0&packageversion=3.0b&platform=src

Unzip it and create a dedicated Python environment for it (BepiPred-3.0 pins
old dependency versions - torch==1.12.0, numpy==1.20.2 - that clash with
Scipion's own environment). Then, in ``scipion.conf``, set:

.. code-block::

      BEPIPRED_HOME = <path to the BepiPred3_src folder>
      BEPIPRED_PYTHON_BIN = <path to the python binary inside that dedicated env>

================================
Download EpiDope
================================

EpiDope is open-source (MIT license, no academic request needed):
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

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-bcellepitope.git
            cd scipion-chem-bcellepitope
            scipion3 installp -p . --devel
