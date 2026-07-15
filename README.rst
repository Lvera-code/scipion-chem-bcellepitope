================================
B-Cell Epitope Prediction scipion plugin
================================

Scipion framework plugin wrapping the B-Cell-Epitope-Prediction pipeline
(https://github.com/Lvera-code/B-Cell-Epitope-Prediction): local BepiPred-3.0
and EpiDope antigenicity prediction, BLASTp immune tolerance filtering and
NetMHCIIpan-4.3 T-helper promiscuity prediction.

**Work in progress.** Currently only the BepiPred-3.0 antigenicity protocol
is implemented (``ProtBCellEpitopeBepiPredPredict``).

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

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-bcellepitope.git
            cd scipion-chem-bcellepitope
            scipion3 installp -p . --devel
