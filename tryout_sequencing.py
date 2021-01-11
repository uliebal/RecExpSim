
from typing import List

import numpy as np

from biolabsim.random import set_seed, pick_integer
from biolabsim.common import Sequence
from biolabsim.shotgun.datatype import LocalizedSequence, estimate_from_overlap
from biolabsim.shotgun.sequencing import Sequencer
from biolabsim.shotgun.visualization import print_scaffold, pretty_print_scaffold, \
    print_assembly_evaluation, print_estimation_evaluation
from biolabsim.shotgun.assembly import GreedyContigAssembler
from biolabsim.shotgun.storage import write_scaffolds_to_file
from biolabsim.shotgun.evaluation import evaluate_sequence
from biolabsim.organism.ecol import Ecol
from biolabsim.organism.fabricated import FabricatedHost


# Set random seed
set_seed(2021)

# Get a good Host to be sequenced.
#wt_host = Ecol()
wt_host = FabricatedHost( genome_size=500, gc_content=0.6 )
#mut_host = wt_host.clone_with_mutation(name="Mut1")
genome = wt_host.get_genome()

# Create the Sequencer. It will fragment the genome.
shotgun = Sequencer(
    library_size_mean=150, library_size_sd=10, read_method='paired-end', read_length=50,
    average_coverage=30, call_error_beta=2.85
)
scafs = shotgun.apply( genome )

# Print the obtained scaffolds.
for i in range(len(scafs)) :
    pretty_print_scaffold(scafs[i])

# Store them in files.
write_scaffolds_to_file( scafs, "output/tryseq_R1.fastq", "output/tryseq_R2.fastq" )

# Method that places multiple sequences together.
def reconstruct ( loci:List[int], seqs:List[Sequence] ) -> List[LocalizedSequence] :
    return [
        LocalizedSequence( seqs[k], v )
        for k,v in enumerate(loci)
    ]

# Place the obtained R1 sequences randomly together and test the consensus sequence.
rnd_loci = [ pick_integer( -5, 80-25+5 ) for _ in range(len(scafs)) ]
rnd_locseqs = reconstruct( rnd_loci, [ scaf.r1_seqrecord.seq for scaf in scafs ] )
rnd_estseq = estimate_from_overlap(rnd_locseqs)

print("\nRandom Positioning - Shannon Entropy: " + str(rnd_estseq.calc_shannon_entropy()))
print("\nRandom Positioning - Positioning <=> Real Genome:")
print_assembly_evaluation( rnd_locseqs, genome )
print("\nRandom Positioning - Estimation <=> Real Genome:")
print_estimation_evaluation( rnd_estseq, genome )
print("\nRandom Positioning - Evaluation: " + str(evaluate_sequence(rnd_estseq,genome)))

# Try to assemble the scaffolds using the greedy assembler.
gca_assembler = GreedyContigAssembler()
gca_locseqs = gca_assembler.apply_internal(scafs)
gca_estseq = gca_assembler.apply(scafs)

print("\nGreedy Contig Assembler - Shannon Entropy: " + str(gca_estseq.calc_shannon_entropy()))
print("\nGreedy Contig Assembler - Positioning <=> Real Genome:")
print_assembly_evaluation( gca_locseqs, genome )
print("\nGreedy Contig Assembler - Estimation <=> Real Genome:")
print_estimation_evaluation( gca_estseq, genome )
print("\nGreedy Contig Assembler - Evaluation: " + str(evaluate_sequence(gca_estseq,genome)))

# Try to assemble using the external SPAdes assembler.






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