
from typing import List
from recordclass import RecordClass



# REFACTOR: The `Host.var_Library` variable is a tricky one. To maintain the single-responsibility
# principle for the Host it would be best to keep `var_Library` outside of the Host, since a single
# Host should not be responsible to keep track of all mutations that have been cloned from it.
# Such a task would be the responsibility of the scientist and maybe they can use a helper
# structure for that. This example here shows such a helper structure. The scientist can
# perform multiple experiments on the host and register the output. That output can then be kept
# track of in the Library.
#
# Maybe the Host itself will need to change its structure to correctly encapsulate what we want
# to simulate. Let's assume that each host has a `target_promoter` attribute that symbolizes
# the "target promoter we want to test for any given gene".
#
# Example Code:
#
#     # Initial preparation of the host and library.
#
#     my_library = Library()
#     wt_host = Ecol()
#
#     # Start an experiment on the cloned host and iteratively fill out the library entry.
#
#     mut_host = wt_host.clone_mutation_with_promoter( name='Mut1', prom="ACGT", prim="A", tm=38 )
#         # if experiment fails, the `mut_host.target_promoter` is None
#     my_entry = Entry( clone_id= mut_host.name )
#
#     ( ok, prom_seq, prom_gc_cont, prom_str ) = measure_target_promoter( mut_host )
#     if ok == True :
#         my_entry.promoter_sequence = prom_seq
#         my_entry.promoter_gc_content = prom_gc_cont
#         my_entry.promoter_strength = prom_str
#
#     ( ok, expr_temp, expr_biomass, expr_rate ) = measure_target_expression( mut_host, 38, 1, 30 )
#     if ok == True :
#         my_entry.expression_temperature = expr_temp
#         my_entry.expression_biomass = expr_biomass
#         my_entry.expression_rate = expr_rate
#
#     # Record the results in the library. Maybe later do some historical analysis.
#
#     my_library.add_entry( my_entry )



class Entry (RecordClass) :
    """
    A Library Entry holds historical information about Host properties that has been analyzed.
    """
    clone_id : str
    promoter_sequence : str # or biolabsim.common.Sequence
    promoter_gc_content
    promoter_strength
    expression_temperature
    expression_biomass
    expression_rate



class Library :
    """
    Library is helpful data structure to keep track of information from multiple strains.
    """

    entries : List[Entry]


    def add_entry ( new_entry:Entry ) :
        self.entries.append(new_entry)
