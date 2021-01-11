
from typing import Optional, List
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..common import Sequence, Scaffold



def write_scaffolds_to_file ( scaffolds:List[Scaffold], r1_path:str, r2_path:Optional[str] = None ) :
    """
    Write a list of scaffolds into the filesystem as a FASTQ file.
    This method will detect paired-end sequences and store them in two different files in their
    reversed form.
    """

    # Resolve to the absolute paths from what was given from the user.
    r1_abspath:Path = Path(r1_path).resolve()
    r2_abspath:Optional[Path] = None
    if r2_path is not None :
        r2_path = Path(r2_path).resolve()

    # Fill both lists of SeqRecords that will be extracted from the scaffolds.
    r1_seqrecs:List[SeqRecord] = []
    r2_seqrecs:List[SeqRecord] = []
    for scaf in scaffolds :
        if scaf.r1_seqrecord is not None :
            r1_seqrecs.append( scaf.r1_seqrecord )
        if scaf.r2_seqrecord is not None :
            r2_seqrecs.append( scaf.r2_seqrecord )

    # Save the files.
    SeqIO.write( r1_seqrecs, r1_abspath, "fastq" )
    print("R1 file stored in: " + str(r1_abspath))
    if len(r2_seqrecs) > 0 :
        if r2_abspath is None :
            raise Exception(
                "Tried to write_scaffolds_to_file. There is at least one paired-end sequence"
                " but no file to store it has been provided."
            )
        SeqIO.write( r2_seqrecs, r2_abspath, "fastq" )
        print("R2 file stored in: " + str(r2_abspath))

