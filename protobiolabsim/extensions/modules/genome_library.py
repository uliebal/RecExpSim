"""
WIP - Work in Progress
This Module tries to fuse both the gene list and sequence together.
Changes to one will affect the other part.
"""


from __future__ import annotations
from copy import copy
from typing import Optional, Set, List, NamedTuple

from numpy.random import Generator
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Blast.Record import Alignment
from Bio.Align import PairwiseAligner

from ...organism import Organism
from ...module import Module
from ..records.gene.gene import Gene
from ..records.gene.localized_gene import LocalizedGene
from ..events import InsertGeneEvent, RemoveGeneEvent



class PrimerMatch (NamedTuple) :
    loc_start: int # Match goes from start (inclusive) to end (exclusive).
    loc_end: int
    success: float # (0,1)



class GenomeLibrary (Module) :
    """
    GenomeLibrary is a module that models the genome in a simplified way to ease interpretation.
    It stores a sequence of base pairs and provides a Gene interface that reads from the sequence
    at the gene's locus.

    The Sequence is stored as a base-pair list, whereas genes are stored as locations inside
    that sequence.
    """

    # bg_size: int # Amount of background gene bases.

    # bg_gc_content: float # Ratio of GC-content in the background genome. Range: [0,1]

    locgenes: List[LocalizedGene]

    sequence: Seq



    def __init__ (
        self, org:Organism,
        sequence:Optional[Seq] = None, # (option 1: pre-defined sequence)
        bg_size:Optional[int] = None, # (option 2: new sequence)
        bg_gc_content: Optional[float] = None, # (option 2: new sequence)
        locgenes:List[LocalizedGene] = []
    ) :
        super().__init__(org)

        # Two options: sequence can come either pre-defined or is generated on the spot.
        if sequence is not None :
            self.sequence = sequence
        elif bg_size is not None and bg_gc_content is not None :
            self.sequence = make_background_seq( rnd=org.make_generator(), size=bg_size, gc_content=bg_gc_content )
        else :
            raise Exception ("Creating GenomeLibrary module without either pre-defined or generated sequence.")

        # Add the locgenes and check if they are matching.
        self.locgenes = []
        for lg in locgenes :
            new_loc_gene = LocalizedGene(
                name=lg.name, seq=self.sequence,
                start_loc=lg._start_loc, prom_len=lg.prom_len, orf_len=lg.orf_len
            )
            self.locgenes.append( new_loc_gene )

        # self.org.observe( InsertGeneEvent, self.listen_insert_gene )
        #self.org.observe( RemoveGeneEvent, self.listen_remove_gene )


    def clone ( self, org:Organism ) -> GenomeLibrary :
        return GenomeLibrary(
            org=org,
            sequence=self.sequence,
            locgenes=self.locgenes
        )


    @property
    def genes ( self ) -> List[Gene] :
        return self.locgenes


    # def listen_insert_gene ( self, event:InsertGeneEvent ) -> None :
    #     self.genes[event.gene] = LocGene( gene=event.gene, loc=event.loc )
    #     print( "Added gene={} loc={} to the GenomeLibrary.".format(
    #         event.gene.get_name(), event.loc
    #     ))



    # def listen_remove_gene ( self, event:RemoveGeneEvent ) -> None :
    #     """ Remove a gene from the library and sequence. """
    #     del self.genes[event.gene]
    #     print( "Removed gene={} from the GenomeLibrary.".format(event.gene.get_name()) )



    def insert_gene ( self, gene:Gene, loc:int ) :
        """
        Insert a gene and its sequence at a specific location in the sequence.
        """
        # Insert the sequence at the insertion point.
        self.sequence = self.sequence[0:loc] + gene.prom + gene.orf + self.sequence[loc:None]

        # Remove genes that start before the insertion and end afterwards (knockout)
        is_knocked_out = lambda lg : lg.start_loc < loc and lg.end_loc > loc
        self.locgenes[:] = [ lg for lg in self.locgenes if not is_knocked_out(lg) ]
            # Replace the list of locgenes with the filtered locgenes, while keeping the list reference.

        # Shift the location of all genes after the insertion point.
        new_gene_len = len(gene.prom) + len(gene.orf)
        for lg in self.locgenes :
            lg.seq = self.sequence # Update to use this new sequence (since they are immutable)
            if lg.start_loc >= loc :
                lg.start_loc = lg.start_loc + new_gene_len

        # Annotate the new sequence as a localized gene.
        new_loc_gene = LocalizedGene(
            name=gene.name, seq=self.sequence,
            start_loc=loc, prom_len=len(gene.prom), orf_len=len(gene.orf)
        )
        self.locgenes.append( new_loc_gene )



    def calc_primer_matches ( self, primer:Seq ) -> List[PrimerMatch] :
        """
        Return the insertion sites a primer can have, alongside with success rate.
        The insertion site is located right after a primer match.
        """
        matches: List[PrimerMatch] = []
        targets: List[Seq] = [ self.sequence ] # List of sequences to check.
        for target in targets :
            if len(target) > 0 :
                match = find_best_primer_match( target, primer )
                if match is not None :
                    matches.append( match )
                    targets.append( target[None:match.loc_start] ) # Also search left of match.
                    targets.append( target[match.loc_end:None] ) # Also search right of match.
        return matches






def make_background_seq ( rnd:Generator, size:int, gc_content:float ) -> Seq :
    gc = gc_content
    at = 1 - gc_content
    seq_array = rnd.choice( ["A","C","G","T"], size=size, p=[at/2,gc/2,gc/2,at/2] ) # 'ATCG' is ArrayLike
    return Seq( "".join(seq_array) )



def check_primer_integrity ( primer:Seq ) -> bool :
    """ Non-deterministic check if the primer itself is well built. """
    return True



def find_best_primer_match ( sequence:Seq, primer:Seq ) -> Optional[PrimerMatch] :
    """ Find best match for a template sequence. """

    # The primer matches with the complement sequence.
    comp_seq = sequence.complement()

    # Build the aligner for matching.
    aligner = PairwiseAligner()
    aligner.mode = 'global' # target_end_gap_score, query_end_gap_score
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.gap_score = -1 # open_gap_score, extend_gap_score
    aligner.query_end_gap_score = 0

    # Get the alignment with the best score, anywhere on the sequence.
    all_alignments = aligner.align( comp_seq, primer )
    best_alignment = next(iter(sorted(all_alignments)))
    #print("score:" + str(best_alignment.score) + "\n" + str(best_alignment))
    if best_alignment is not None and best_alignment.score > 0 :
        return PrimerMatch(
            loc_start= best_alignment.aligned[0][0][0], # start of first chunk in target
            loc_end= best_alignment.aligned[0][-1][1], # end of last chunk in target
            success= best_alignment.score / len(primer)
        )
    return None





# def build_sequence ( self ) -> Seq :
#     rnd = self.org.make_generator()
#     gc = self.bg_gc_content
#     at = 1 - self.bg_gc_content

#     # Write the sequence from begin to end. Genes have their sequence copied and random bases
#     # are used to fill non-coding positions. Overlapping genes have undefined behaviour.
#     # TODO: When gene positions overlap it actually writes both gene sequences in sequence but
#     # in future this behaviour should be improved. This deals with validity of genes in
#     # overlapping regions.
#     locgenes = self.genes.values()
#     locgenes.sort( key=lambda el: el.loc )

#     seq = Seq()
#     cur = 0
#     bg_size = 0 # Current amount of background bases.
#     for locgene in locgenes :

#         # If the gene starts later, then fill the space with random (background) bases.
#         filler = locgene.loc - cur
#         if filler > 0 :
#             new_bases = rnd.choice( 'ATCG', size=filler, p=[at,at,gc,gc] )
#             bg_size += filler
#             seq += "".join(new_bases)

#         # Copy the gene sequence.
#         seq = seq + locgene.gene.get_prom() + locgene.gene.get_orf()

#     # At the end, if the minimum of background bases was not filled, to that now.
#     filler = min( 0, self.bg_min_size - bg_size )
#     if filler > 0 :
#         new_bases = rnd.choice( 'ATCG', size=filler, p=[at,at,gc,gc] )
#         seq += "".join(new_bases)

#     return seq