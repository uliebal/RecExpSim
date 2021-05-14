"""
Events that serve as communication for and between organism modules.
"""

from __future__ import annotations
from abc import ABC
from typing import Literal, NamedTuple
from dataclasses import dataclass

from ..record.gene.base import Gene


EventType = str # Literal[ 'initialize', 'insert_gene' ]


class Event :
    typ: EventType


@dataclass
class InitializeEvent ( Event ) :
    typ: Literal['initialize'] = 'initialize'


@dataclass
class InsertGeneEvent ( Event ) :
    gene: Gene
    typ: Literal['insert_gene'] = 'insert_gene'


@dataclass
class RemoveGeneEvent ( Event ) :
    gene: Gene
    typ: Literal['remove_gene'] = 'remove_gene'


@dataclass
class PromStrCalcEvent ( Event ) :
    gene: Gene
    result: float
    typ: Literal['prom_str_calc'] = 'prom_str_calc'
