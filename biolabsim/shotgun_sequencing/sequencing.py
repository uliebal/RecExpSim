
from typing import List, Optional, NamedTuple, Literal
from recordclass import RecordClass

from .fake_genome import SimplifiedGenome
from ..common import Sequence, ReadMethod, Scaffold
from ..random import pick_normal, pick_float, pick_integer



class RatedSequence (RecordClass) :
    """
    Internal structure to store sequences and their quality score.
    """
    seq : Optional[Sequence]
    qual : Optional[List[float]]



class Sequencer :
    """
    Function that, after setup, can take in a genome and perform whole-genome sequencing on it.
    The result will simulate a sequencing machine and provide a set of fragmented genome libraries
    that can later be assembled together to reconstruct the real genome of the host.

    Sequencing machines can be highly customizable and allow a multitude of parameters. Some
    considerations about the type of parameters and their default values are given below:

      - Usually a library size selection is done and this simulator will perform similar effects to
        those described under "Double-sided size selection and bead clean-up".
        A simplification is done and the blue graph in the link is converted to a normal distribution.
        https://emea.support.illumina.com/bulletins/2020/07/library-size-selection-using-sample-purification-beads.html

      - The read length of each fragment can also be configured, but a limit at roughly 150bp is set
        for each library, with some versions achieving up to 300bp. This limit is for each read.
        For example, with a value of 150bp, running the simulation on single-read method will yield
        libraries of 150bp each, whereas a paired-end method will yield two linked sequences
        of 150bp each.

      - Coverage of the extracted libraries can be controlled by reagents and number of cycles in
        the cloning step. The distribution of covered regions, however, depends on how the fragments
        attach to the flow cell prior to cluster forming. Attachment to the flow cell will follow
        random sampling from the available fragments. Given the amount of fragments this can be
        sampling to "sampling with reposition". The average coverage and read length; the genome
        size is known before-hand; and the number of reads will be calculated from all these values
        by using the Lander/Waterman equation.
        https://www.illumina.com/documents/products/technotes/technote_coverage_calculation.pdf

    The sequences that are output will only contain the target region of each library.
    Adapter, index and tag regions are discarded.
    """

    # Average number of base pairs each resulting library will have.
    library_size_mean : float

    # Standard deviation in the length of resulting OpenReadFrames.
    library_size_sd : float

    # Method used when reading fragments.
    read_method : ReadMethod

    # Read length of each fragment.
    read_length : int

    # Average coverage of each base. Influences amount of frames that are returned.
    average_coverage : float

    # Probability that a base read is wrong.
    # If this occurs, the base that was read will be randomly and equally selected (could even be
    # the correct one).
    # TODO: See how this behaves in reality.
    read_error_prob : float



    def __init__ (
        self,
        library_size_mean:float = 400,
        library_size_sd:float = 75,
        read_method:ReadMethod = 'single-read',
        read_length:int = 150,
        average_coverage:float = 10,
        read_error_prob:float = 0.01  # Notes say its typically around 0.5% - 2.0%
    ) :
        self.library_size_mean = library_size_mean
        self.library_size_sd = library_size_sd
        self.read_method = read_method
        self.read_length = read_length
        self.average_coverage = average_coverage
        self.read_error_prob = read_error_prob



    def apply ( self, genome:SimplifiedGenome ) -> List[Scaffold] :
        """
        TODO: Decide between Host and Strain.
        """

        # Go for the approach of selecting randomly a start point and length, then check if it is
        # possible. Repeat if not.
        out_scaffolds:List[Scaffold] = []
        num_scaffolds = self.calc_num_scaffolds_obtained(genome)
        while len(out_scaffolds) < num_scaffolds :
            start = pick_integer( low=0, high=len(genome.template_strand) ) # [low,high)
            library_size = int(pick_normal( self.library_size_mean, self.library_size_sd ))
            end = start + library_size

            # Repeat this sampling if it surpassed the possible range.
            if not ( end <= len(genome.template_strand) ) :
                continue

            # Define the extracted the slices. R2 is only extracted if paired-end method is used.
            r1 = RatedSequence( seq=None, qual=None )
            r2 = RatedSequence( seq=None, qual=None )

            # Extract the slices.
            read_size = min( library_size, self.read_length ) # Read is restricted by max read length.
            r1.seq = genome.template_strand[ start : start + read_size ]
            if self.read_method == 'paired-end' :
                r2.seq = genome.template_strand[ end - read_size : end ]

            # Apply quality estimation and substitution errors for both
            for r in [ r1, r2 ] : # Same code for both possible reads.
                if r.seq != None : # Catch non-existent R2 reads.
                    r.qual = [ 1 - self.read_error_prob for _ in range(len(r.seq)) ] # TODO: Uniform quality.
                    for i in range(len(r.seq)) :
                        if pick_float() < self.read_error_prob :
                            err_bases = [ b for b in ['A','T','C','G'] if b != r.seq[i] ]
                                # Always substitute for a base that is not the original.
                            sel_base = err_bases[ pick_integer(0,len(err_bases)) ]
                            r.seq[i] = sel_base
                            r.qual[i] = 0.25 # TODO: Remove this. Debugging.

            # Create the scaffold and add it to the final results.
            cur_scaffold = Scaffold(
                read_method=self.read_method, expected_len=library_size,
                r1_sequence=r1.seq, r1_quality=r1.qual,
                r2_sequence=r2.seq, r2_quality=r2.qual
            )
            out_scaffolds.append(cur_scaffold)

        return out_scaffolds



    def calc_num_scaffolds_obtained ( self, genome:SimplifiedGenome ) -> int :
        """
        Number of frames that will be obtained. Calculated based on Lander/Waterman equation with
        regards to the current parameters and the host genome.
        This will account for the read method used. Paired-end methods will output half of it.
        TODO: Change between Host and Strain.
        """
        method_factor = 2 if self.read_method == 'paired-end' else 1
        return round(
            ( len(genome) * self.average_coverage )
            /
            ( self.read_length * method_factor )
        )
