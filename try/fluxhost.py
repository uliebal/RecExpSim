# This example will show interplay of a cobra model and gene expression on the same host.
# By editing the promoter of a gene, we change the expression of said gene. That changed expression
# will reflect on the metabolic model.

import cobra
from cobra.core.solution import Solution as CobraSolution
from Bio.Seq import Seq

from silvio import DATADIR, Experiment, Host, GenomeLibrary, GenomeExpression, MetabolicFlux, \
    AlterGenePromoterEvent, coalesce



#
# START OF THE DEFINITION - CATALOG
#



class EcolFluxHost (Host) :

    genome: GenomeLibrary
    genexpr: GenomeExpression
    metflux: MetabolicFlux



    def make ( self ) -> None :

        self.genome = GenomeLibrary()
        self.genome.make( bg_size=300, bg_gc_content=0.45, bg_rnd=self.make_generator() )
        self.genome.bind( host=self )

        self.genexpr = GenomeExpression()
        self.genexpr.make(
            opt_primer_len=20, infl_prom_str=40, species_prom_str=0.057,
            regressor_file = DATADIR / 'expression_predictor' / 'Ecol-Promoter-predictor.pkl',
            addparams_file = DATADIR / 'expression_predictor' / 'Ecol-Promoter-AddParams.pkl',
        )
        self.genexpr.bind( host=self, genome=self.genome )

        path = DATADIR / 'metabolic_model' / 'e_coli_core.xml'
        model = cobra.io.read_sbml_model(str(path))
        self.metflux = MetabolicFlux()
        self.metflux.make( model=model )
        self.metflux.bind( host=self )



    def copy ( self, ref:'EcolFluxHost' ) -> None :

        self.genome = GenomeLibrary()
        self.genome.copy( ref=ref.genome )
        self.genome.bind( host=self )

        self.genexpr = GenomeExpression()
        self.genexpr.copy( ref=ref.genexpr )
        self.genexpr.bind( host=self, genome=self.genome )

        self.metflux = MetabolicFlux()
        self.metflux.copy( ref=ref.metflux )
        self.metflux.bind( host=self )



    def sync ( self ) -> None :
        self.sync_modules([ self.genome, self.genexpr, self.metflux ])



    def edit_gene ( self, gene_name:str, new_promoter:Seq ) -> None :
        found_genes = [ gene for gene in self.genome.locgenes if gene_name == gene.name ]
        if len(found_genes) > 0 :
            self.emit( AlterGenePromoterEvent(found_genes[0],new_promoter) )



    def optimize ( self ) -> CobraSolution :
        return self.metflux.optimize()



class FluxExp (Experiment) :

    def create_host ( self ) -> EcolFluxHost :
        new_host = EcolFluxHost()
        new_host.make()
        new_host.sync()
        self.bind_host(new_host)
        return new_host

    def clone_host ( self, host:Host ) -> EcolFluxHost :
        new_host = EcolFluxHost( ref=host )
        new_host.copy( ref=host )
        self.bind_host(new_host)
        return new_host



#
# START OF THE SCRIPT - NOTEBOOK
#

exp = FluxExp()

# Create first ECol host from scratch.
print("###\n### First Ecol\n###")
ecol1 = exp.create_host()
ecol1.print_event_log()
sol1 = ecol1.optimize()

# Create second Ecol host as a clone from the first. Slightly lower the expressions.
print("###\n### Second Ecol\n###")
ecol2 = exp.clone_host(ecol1)
for gn in ["glcB", "aceB", "sucD", "frmA","puuA"] :
    ecol2.edit_gene( gn, "GCCCAAAAAAAAAGCAAACACGTAAAGGAAAAAATGCACG" )
ecol2.print_event_log()
sol2 = ecol2.optimize()

# Create third Ecol host as a clone from the second. Reduce the expression by a lot.
print("###\n### Third Ecol\n###")
ecol3 = exp.clone_host(ecol1)
for gn in ["lpd","sdhD","cbdB","fumC","nuoK","gltD","gltP"] :
    ecol3.edit_gene( gn, "ATGGGGGGGGGGGGGGGGGGGGG" )
ecol3.print_event_log()
sol3 = ecol3.optimize()

print("Solution 1 = ", sol1)
print("Solution 2 = ", sol2)
print("Solution 3 = ", sol3)

pass
