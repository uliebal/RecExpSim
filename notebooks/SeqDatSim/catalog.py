"""
An experiment that uses hosts which implements all currently existing modules.
"""

from __future__ import annotations
from typing import Optional

from Bio.Seq import Seq

from silvio import (
    DATADIR, # config
    Experiment, ExperimentException, Host, HostException, # base classes
    alldef, coalesce, first, Generator, # utilities
    DataOutcome, DataWithPlotOutcome, combine_data, # outcome
    InsertGeneEvent, # events
    ShotgunSequencer, RandomContigAssembler, GreedyContigAssembler, ReadMethod, # tools
    GenomeLibrary, # modules
    Gene, # records
)



class SeqExperiment (Experiment) :
    """
    Experiment to perform shotgun sequencing and assembly.
    """

    # Keep track of how many hosts were created
    host_counter: int



    def __init__ ( self, seed:Optional[int] = None ) :
        super().__init__(seed=seed)
        self.host_counter = 0



    def create_host ( self, name:Optional[str], bg_size:int, bg_gc_content:float ) -> SeqHost:
        self.host_counter += 1
        seed = self.rnd_gen.pick_seed() # The experiment provides stable seed generation for hosts.
        chosen_name = coalesce( name, "host" + str(self.host_counter) )
        new_host = SeqHost( name=chosen_name, seed=seed )
        new_host.make( bg_size=bg_size, bg_gc_content=bg_gc_content )
        new_host.sync()
        self.bind_host(new_host)
        return new_host


    def create_sequencer (
        self,
        library_size_mean:float, library_size_sd:float,
        read_method:ReadMethod, read_length:int,
        average_coverage:float, call_error_beta:float
    ) -> ShotgunSequencer :
        seed = self.rnd_gen.pick_seed() # The experiment provides stable seed generation for tools.
        new_tool = ShotgunSequencer(
            library_size_mean=library_size_mean, library_size_sd=library_size_sd,
            read_method=read_method, read_length=read_length,
            average_coverage=average_coverage, call_error_beta=call_error_beta,
            name=None, seed=seed
        )
        return new_tool



    def create_random_assembler ( self, expected_genome_size:int ) -> RandomContigAssembler :
        seed = self.rnd_gen.pick_seed() # The experiment provides stable seed generation for tools.
        new_tool = RandomContigAssembler(
            expected_genome_size=expected_genome_size,
            name=None, seed=seed
        )
        return new_tool



    def create_greedy_assembler ( self ) -> GreedyContigAssembler :
        seed = self.rnd_gen.pick_seed() # The experiment provides stable seed generation for tools.
        new_tool = GreedyContigAssembler( name=None, seed=seed )
        return new_tool



    def find_host_or_abort ( self, host_name:str ) -> SeqHost :
        """
        Find a host by a given name or abort with an exception.
        Throws: ExperimentException
        """
        host:Optional[SeqHost] = first( self.hosts, lambda h: h.name == host_name )
            # Find the first host matching the name.
        if host is None :
            raise ExperimentException("Experiment has not found the host '{}'".format(host_name))
        return host



    def clone_host ( self, host_name:str ) -> SeqHost :
        """
        Make an exact clone of the host.
        """
        host = self.find_host_or_abort(host_name)
        new_host = SeqHost( ref=host ) # Clone the host but with different seed.
        new_host.copy( ref=host )
        self.bind_host(new_host)
        return new_host



    def print_status ( self ) -> None :
        print("Experiment:")
        print("  hosts = [ {} ]".format( " , ".join([h.name for h in self.hosts]) ))



class SeqHost (Host) :
    """
    Host with a genome that contains both genes and background genes are important.
    """

    genome: GenomeLibrary



    def make ( self, bg_size:int = None, bg_gc_content:float = None ) -> None :
        self.genome = GenomeLibrary()
        self.genome.make( bg_size=bg_size, bg_gc_content=bg_gc_content, bg_rnd=self.make_generator() )
        self.genome.bind( host=self )



    def copy ( self, ref:SeqHost ) -> None :
        self.genome = GenomeLibrary()
        self.genome.copy( ref=ref.genome )
        self.genome.bind( host=self )



    def sync ( self ) -> None :
        self.sync_modules([ self.genome ])



    def insert_gene ( self, gene:Gene ) -> None :
        self.emit( InsertGeneEvent( gene=gene, locus=0 ) )



    def get_genome_sequence ( self ) -> Seq :
        return self.genome.sequence



    def print_status ( self ) -> None :
        print("Host [{}]:".format( self.name ))
        print("  seed plus counter = {} + {}".format( self.rnd_seed, self.rnd_counter ))
        print("  Gene List: {} genes".format(len(self.genome.genes)))
        for gene in self.genome.genes :
            print("  - {} = {} * {}".format(gene.name, gene.prom, gene.orf))
        print("  Event History: {} events".format(len(self.event_log)))
        for el in self.event_log :
            print("  - {}".format(el))
