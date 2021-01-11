
from typing import List

import numpy as np

from biolabsim.random import set_seed, pick_integer
from biolabsim.common import Sequence
from biolabsim.shotgun_sequencing.sequencing import Sequencer
from biolabsim.shotgun_sequencing.visualization import print_scaffold, print_assembly_evaluation
from biolabsim.shotgun_sequencing.assembly import LocalizedSequence, to_consensus, GreedyContigAssembler
from biolabsim.shotgun_sequencing.storage import write_scaffolds_to_file
from biolabsim.organism.ecol import Ecol
from biolabsim.organism.fabricated import FabricatedHost


# Set random seed
set_seed(2020)

# Get a good Host to be sequenced.
#wt_host = Ecol()
wt_host = FabricatedHost( genome_size=80, gc_content=0.6 )
#mut_host = wt_host.clone_with_mutation(name="Mut1")
genome = wt_host.get_genome()

# Create the Sequencer. It will fragment the genome.
shotgun = Sequencer(
    library_size_mean=30, library_size_sd=5, read_method='single-read', read_length=20,
    average_coverage=4
)
scafs = shotgun.apply( genome )

# Print the obtained scaffolds.
#for i in range(len(scafs)) :
#    print_scaffold(scafs[i])

# Store them in files.
write_scaffolds_to_file( scafs, "output/tryseq_R1.fastq" )

# Method that places multiple sequences together.
def reconstruct ( loci:List[int], seqs:List[Sequence] ) -> List[LocalizedSequence] :
    return [
        LocalizedSequence( seqs[k], v )
        for k,v in enumerate(loci)
    ]

# Place the obtained R1 sequences randomly together and test the consensus sequence.
rnd_loci = [ pick_integer( -5, 80-30+5 ) for _ in range(len(scafs)) ]
rnd_locseqs = reconstruct( rnd_loci, [ scaf.r1_seqrecord.seq for scaf in scafs ] )

print("\nEvaluation using Random Locations <=> Real Genome:")
print_assembly_evaluation( rnd_locseqs, genome )

# Try to assemble the scaffolds using the greedy assembler.
assembler = GreedyContigAssembler()
asb_locseqs = assembler.apply_internal(scafs)

print("\nEvaluation using Greedy Contig Assembler <=> Real Genome:")
print_assembly_evaluation( asb_locseqs, genome )



# # START OF EDITABLE SECTION (will now show multiple examples)

# # The student exercise is generally to create a well-formed SequenceAssembly.
# # This example does it manually but its best to use an algorithm for it.
# print("--- DEFAULT RECONSTRUCTION ---")
# default_reconst = reconstruct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ])
# print_assembly_evaluation( default_reconst, genome )
# print()



# print("--- RANDOM RECONSTRUCTION ---")
# randgen = np.random.default_rng(2020)
# random_loci = [ randgen.integers(0,len(genome.template_strand)) for _ in range(12) ]
# random_reconst = reconstruct(random_loci)
# print_assembly_evaluation( random_reconst, genome )
# print()



# print("--- MANUAL RECONSTRUCTION ---")
# manual_reconst = reconstruct([ 1, 35, 19, 15, 19, 57, 4, 28, 58, 15, 28, 52 ])
# print_assembly_evaluation( manual_reconst, genome )
# print()




# # SERVER                                    CLIENT
# #
# # genome = Genome()
# # fragments = Sequencer.apply(genome)
# #          ------   fragments.txt    ------>
# #                                            client will reconstruct
# #          <----- [loci,loci,...].txt ------
# # compare( fragments, loci, genome )
# #          -------     score         ----->


# # SERVER                                     CLIENT
# #
# #         <----   assembler-code   ------
# # genome = Genome()
# # fragments = Sequencer.apply(genome)
# # reconstruction = "assember-code".apply(fragments)
# # compare( reconstruction, genome )
# #         -------     score       ----->


# # biolabsim.vs42.rwht-aachen.de


# # [Beispiel Anfrage]  --> downloads fragments.txt
# #
# # [     ] [Antwort Abgeben]   --> schicken