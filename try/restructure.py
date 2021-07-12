
from __future__ import annotations
from Bio.Seq import Seq

# Use biolabsim from parent folder.
import os
import sys
sys.path.append( os.path.abspath(os.path.join('.')) )
sys.path.append( os.path.abspath(os.path.join('..')) )

from protobiolabsim.catalog.recexpsim2 import RecExperiment2, Ecol
from protobiolabsim.extensions.records.gene.crafted_gene import CraftedGene

exp = RecExperiment2()

# Org 1 has a single gene.

org1 = Ecol( exp=exp )
org1.print_status()

gen1 = CraftedGene( name="Cx76", orf=Seq("ATGAACTGGTTACCGGATCGCTTATCGCTACG"), prom=Seq("TTCGTG") )
org1.genome.insert_gene( gene=gen1, loc=20 ) # Model import would call multiple insert gene.
org1.print_status()

gen2 = CraftedGene( name="Bh44", orf=Seq("ATGTGGTCTATGCATTTATG"), prom=Seq("CCCGGGCCC") )
org1.genome.insert_gene( gene=gen2, loc=70 ) # Model import would call multiple insert gene.
org1.print_status()

gen3 = CraftedGene( name="Pre6", orf=Seq("ATGTTATTAG"), prom=Seq("GCG") )
org1.genome.insert_gene( gene=gen3, loc=10 ) # Model import would call multiple insert gene.
org1.print_status()

gen4 = CraftedGene( name="X7", orf=Seq("ATGTATAGT"), prom=Seq("TAATAG") )
org1.genome.insert_gene( gene=gen4, loc=50 ) # Model import would call multiple insert gene.
org1.print_status()

pri = Seq("CGGGT")
matches = org1.genome.calc_primer_matches(pri)
print("Primer Matches:\n" + "\n".join([ "  > {} | {}:{}".format( m.success, m.loc_start, m.loc_end ) for m in matches ]))

pri = pri
gen5 = CraftedGene( name="Cx76", orf=Seq("TTTTTTTTTTT"), prom=Seq("GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG") )
( org2, outcome ) = org1.clone_with_recombination( primer=pri, gene=gen5, tm=25 )
print("Org2 Cloning Outcome: " + str(outcome))
org2.print_status()





# # Org 2 has an additional gene.

# gen2 = CraftedGene( name="Hut3", orf=Seq("GGGTTG"), prom=Seq("CCT") )

# org2 = org1.clone()

# org2.insert_gene( gen2 )

# print( "Organism 2 DNA: {}".format( org2.genseq.print_sequence() ) )

# # Org 3 has a mutation on the second gene.

# gen3 = gen2.make_variant( orf_loc=3, orf_sub="A" )

# org3 = org2.clone()

# org3.replace_gene( gen2, gen3 )

# print( "Organism 3 DNA: {}".format( org3.genseq.print_sequence() ) )

# # Show the gene registry.

# # print("Gene Registry:")
# # for gene in exp.gene_reg.records :
# #     print("  {}  {}  {}".format(gene.get_name(),gene.get_prom(),gene.get_orf()))



print("Finished")




# for org in exp.orgs :
#     org.genlib.check_included(gene_rec)

#     # if generec in genlib_mod.genes

# org.genlib.genes['ATP'].substitute_promoter( "ATTACG" )
# org.genseq.substitute_bases( ch=0, pos=1360, base="C" )
# org.genseq.substitute_seq( ch=0, pos=[1360:1363], sub="CAA" )


# org.genseq.match_and_replace( match="ATTCA", sub="ATTGA" )
# org.genseq.match_and_replace( match="ATTCA", offset=4, sub="G" )


