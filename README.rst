================================
DEPRECATED: B-Cell Epitope Prediction scipion plugin
================================

This plugin is **deprecated** and pending archival. Following feedback from
the CNB-CSIC Scipion-chem maintainers, the 5-phase pipeline it used to wrap
(https://github.com/Lvera-code/BCell-Epitope-Prediction) has been split
across the existing Scipion-chem ecosystem instead of living in its own
plugin, to avoid duplicating protocols that already exist upstream.

Where each phase moved to
================================

**Phase 2/3 (BepiPred-3.0 antigenicity)**
    Use the existing ``scipion-chem-bepipred`` plugin directly:
    https://github.com/scipion-chem/scipion-chem-bepipred . A generic
    "gap-tolerant sliding window" epitope-extraction mode (an alternative to
    its existing threshold + soft-extension algorithm) was added on a branch
    pending PR: https://github.com/Lvera-code/scipion-chem-bepipred/tree/feat/gap-tolerant-window-mode

**Phase 2/3 (EpiDope antigenicity)**
    Moved into ``scipion-chem`` itself (``pwchem/protocols/Sequences/protocol_epidope.py``),
    since EpiDope is open source (MIT) and installs automatically via
    ``defineBinaries``, unlike the two academic-license engines below. On a
    branch pending PR: https://github.com/Lvera-code/scipion-chem/tree/feat/epidope-protocol

**Phase 3 (BepiPred + EpiDope consensus)**
    No new protocol: use the existing generic ``Operate sequence ROIs``
    protocol (``ProtOperateSeqROI``) in ``scipion-chem``
    (``pwchem/protocols/Sequences/protocol_operate_sequence_rois.py``) with
    the BepiPred and EpiDope ``outputROIs`` as the two input sets:
    operation=Union, minOverlap=0, keepNonOverlaping=True, keepAttributes=True.
    Note this keeps each output ROI's ``_meanScore`` (from whichever input
    ROI it originated from) but, unlike the retired ``ProtBCellEpitopeConsensus``,
    it does not explicitly label a ROI as 'Consenso'/'BepiPred'/'EpiDope'.

**Phase 4 (BLASTp immune tolerance filter)**
    Use the existing ``scipion-chem-blast`` plugin directly:
    https://github.com/scipion-chem/scipion-chem-blast . It did not support
    batching multiple ROI queries with length-tiered E-values, so a generic
    "batch mode" (accepts a whole SetOfSequences/SetOfSequenceROIs, with
    optional auto length-tiered E-value) was added on a branch pending PR:
    https://github.com/Lvera-code/scipion-chem-blast/tree/feat/batch-roi-input

**Phase 5 (NetMHCIIpan-4.3 T-helper promiscuity)**
    Split into its own standalone plugin, since it requires a separate
    academic-license installation (not redistributable, unlike EpiDope):
    https://github.com/Lvera-code/scipion-chem-netmhciipan . This is the only
    piece that stays as an independently maintained plugin; once validated it
    is meant to be transferred to the ``scipion-chem`` GitHub organization.

Status
================================

All of the above are on branches, not yet merged/PR'd upstream: they need to
be installed and tested end-to-end before opening pull requests to
``scipion-chem``. Until then, this repository is kept around for reference
only and should not be installed.
