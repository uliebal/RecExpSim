"""
This file contains helper methods to visualize data structures and aid in debugging.
"""

from math import floor, ceil

from ..common import Scaffold

def print_scaffold ( scaffold:Scaffold ) -> None :
    """
    Visually print a scaffold.
    """
    r2_exists = scaffold.r2_sequence != None

    # Calculate the gap size.
    gap = 0
    if scaffold.read_method == 'single-read' :
        gap = scaffold.expected_len - len(scaffold.r1_sequence)
    elif scaffold.read_method == 'paired-end' :
        gap = scaffold.expected_len - len(scaffold.r1_sequence) - len(scaffold.r2_sequence)

    print( "| {r1s}~{lspad}(ca.{gap}){rspad}~{r2s}\n| {r1q}~{lspad}~~~~{cpad}{rspad}~~{r2q}".format(
        r1s="".join(scaffold.r1_sequence),
        r1q="".join(map(str,[ min(9,floor(q*10)) for q in scaffold.r1_quality ])),
        r2s="".join(reversed(scaffold.r2_sequence)) if r2_exists else "",
        r2q="".join(map(str,[ min(9,floor(q*10)) for q in reversed(scaffold.r2_quality) ])) if r2_exists else "",
        lspad="~" * ceil(gap / 2),
        rspad="~" * floor(gap / 2),
        gap=gap,
        cpad="~" * len(str(gap))
    ))


