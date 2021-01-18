
from typing import List

import numpy as np

from biolabsim import set_seed, FabricatedHost, Sequencer, print_scaffold, write_scaffolds_to_file, \
    RandomAssembler, estimate_from_overlap, print_assembly_evaluation, print_estimation_evaluation, \
    GreedyContigAssembler, evaluate_sequence


# Set random seed
set_seed(2021)

# Get a good Host to be sequenced.
#wt_host = Ecol()
#target_host = wt_host.clone_with_mutation(name="Mut1")
target_host = FabricatedHost( genome_size=120, gc_content=0.6 )
genome = target_host.get_genome()

# Create the Sequencer. It will fragment the genome.
shotgun = Sequencer(
    library_size_mean=80, library_size_sd=5, read_method='paired-end', read_length=35,
    average_coverage=5, call_error_beta=2.85
)
scafs = shotgun.apply( genome )

# Print the obtained scaffolds.
print("\n")
for i in range(len(scafs)) :
    print_scaffold(scafs[i])
print("\n")

# Store them in files.
write_scaffolds_to_file( scafs, "output/tryseq_R1.fastq", "output/tryseq_R2.fastq" )

# Place the obtained R1 sequences randomly together and test the consensus sequence.
rnd_assembler = RandomAssembler(expected_genome_size=120)
rnd_locseqs = rnd_assembler.apply_internal(scafs)
rnd_estseq = estimate_from_overlap(rnd_locseqs) # rnd_assembler.apply(scafs)

print("\nRandom Positioning - Shannon Entropy: {:.4f}".format(rnd_estseq.calc_shannon_entropy()))
print("\nRandom Positioning - Positioning <=> Real Genome:")
print_assembly_evaluation( rnd_locseqs, genome )
print("\nRandom Positioning - Estimation <=> Real Genome:")
print_estimation_evaluation( rnd_estseq, genome )
print("\nRandom Positioning - Evaluation: {:.4f}".format(evaluate_sequence(rnd_estseq,genome)))

# Try to assemble the scaffolds using the greedy assembler.
gca_assembler = GreedyContigAssembler()
gca_locseqs = gca_assembler.apply_internal(scafs)
gca_estseq = estimate_from_overlap(gca_locseqs) # gca_assembler.apply(scafs)

print("\nGreedy Contig Assembler - Shannon Entropy: {:.4f}".format(gca_estseq.calc_shannon_entropy()))
print("\nGreedy Contig Assembler - Positioning <=> Real Genome:")
print_assembly_evaluation( gca_locseqs, genome )
print("\nGreedy Contig Assembler - Estimation <=> Real Genome:")
print_estimation_evaluation( gca_estseq, genome )
print("\nGreedy Contig Assembler - Evaluation: {:.4f}".format(evaluate_sequence(gca_estseq,genome)))


