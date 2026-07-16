"""DEPRECATED: this plugin has been split up, see README.rst.

Kept only so the entry point still resolves for anyone with it already
installed; it no longer defines any protocol. This package is pending
archival once the migration below is validated.
"""

from pwchem import Plugin as pwchemPlugin


class Plugin(pwchemPlugin):

    @classmethod
    def _defineVariables(cls):
        pass

    @classmethod
    def defineBinaries(cls, env):
        pass
