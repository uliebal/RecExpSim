
from typing import List

import numpy as np

from biolabsim.random import set_seed, pick_integer
from biolabsim.common import Sequence
from biolabsim.shotgun_sequencing.fake_genome import RandomGenome
from biolabsim.shotgun_sequencing.sequencing import Sequencer
from biolabsim.shotgun_sequencing.visualization import print_scaffold
from biolabsim.shotgun_sequencing.assembly import LocalizedSequence, to_consensus, GreedyContigAssembler
from biolabsim.shotgun_sequencing.evaluation import print_assembly_evaluation


# Set random seed
set_seed(2020)

# Create the genome.
genome = RandomGenome( gc_content=0.6, num_bp=80 )

# Creat the Sequencer. Its fragments the genome
shotgun = Sequencer(
    library_size_mean=30, library_size_sd=5, read_method='single-read', read_length=20,
    average_coverage=4, read_error_prob=0.01
)
scafs = shotgun.apply( genome )

# Print the obtained scaffolds.
for i in range(len(scafs)) :
    print("|\n| [{}]".format(i+1))
    print_scaffold(scafs[i])

# Method that places multiple sequences together.
def reconstruct ( loci:List[int], seqs:List[Sequence] ) -> List[LocalizedSequence] :
    return [
        LocalizedSequence( seqs[k], v )
        for k,v in enumerate(loci)
    ]

# Place the obtained R1 sequences randomly together and test the consensus sequence.
rnd_loci = [ pick_integer( -5, 80-30+5 ) for _ in range(len(scafs)) ]
rnd_locseqs = reconstruct( rnd_loci, [ scaf.r1_sequence for scaf in scafs ] )
rnd_consens = to_consensus(rnd_locseqs)

print("\nEvaluation using Random Locations <=> Random Consensus:")
print_assembly_evaluation( rnd_locseqs, rnd_consens )

print("\nEvaluation using Random Locations <=> Real Genome:")
print_assembly_evaluation( rnd_locseqs, genome.template_strand )

# Try to assemble the scaffolds using the greedy assembler.
assembler = GreedyContigAssembler()
asb_locseqs = assembler.apply_internal(scafs)
asb_consens = to_consensus(asb_locseqs)

print("\nEvaluation using Greedy Contig Assembler <=> Greedy Consensus:")
print_assembly_evaluation( asb_locseqs, asb_consens )

print("\nEvaluation using Greedy Contig Assembler <=> Real Genome:")
print_assembly_evaluation( asb_locseqs, genome.template_strand )



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