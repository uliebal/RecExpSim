
from dataclasses import dataclass

from ...events import Event
from .records.gene import Gene


@dataclass
class InsertGeneEvent ( Event ) :
    gene: Gene


@dataclass
class RemoveGeneEvent ( Event ) :
    gene: Gene
