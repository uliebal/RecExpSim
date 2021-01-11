
from __future__ import annotations

import numpy as np
from cobra.io import read_sbml_model
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Avoid circular dependencies because of typing.
    from .host import Host
    from .strain import Strain



def load_local_sbml_model ( path ) :
    """
    Load a SBML Model from the local filesystem.
    """
    return read_sbml_model(str(path))



def Help_LoadCobra(Path=False):
    '''
    This function loads cobra models from input path, or by default the Ecoli core model. Either from a local copy or dowloading it from the Bigg database.

    Input:
        Path (optional): string, relative path address to local copy

    Output:
        Model:           cobraby model, E. coli core reactions
    '''

    import os

    if Path:
        Model = read_sbml_model(Path)
    else:
        os.system('wget http://bigg.ucsd.edu/static/models/e_coli_core.xml')
        # generating cobra variable from SBML/xml/mat file
        Path = 'e_coli_core.xml'
        Model = read_sbml_model(Path)

    return Model

def Help_GeneAnnotator(HostName, Model):
    '''
    Function to characterize genes.
    Input
        HostName:     name of the host
        Model:    cobra model
    Output
        Genes_df: dataframe, gene name from enzyme id in model, expression strength, promoter sequence, ORF, flux
    '''
    import pandas as pd
    import multiprocessing
    from joblib import Parallel, delayed
    from .genome import make_GeneJoiner

    num_cores = multiprocessing.cpu_count()
#     use_core = min([num_cores, n])

    Result = Parallel(n_jobs=num_cores)(delayed(make_GeneJoiner)(HostName, Model, myRct.id) for myRct in Model.reactions)
    Genes_df = pd.DataFrame(Result)

    return Genes_df

def Help_StrainCharacterizer(HostName, GenesDF, Genome_WT, Genome_MT, Model):
    '''
    Function to scan a manipulated genome for changes in gene expression of enzymes.
    Input:
        Mutant:         parent class
        StrainID
    Output:
        RctNew_df:      dataframe, containing reactions with updated expression
    '''
    import pandas as pd
    import multiprocessing
    from joblib import Parallel, delayed
    from .expression import make_UpdateExpression

    num_cores = multiprocessing.cpu_count()
#     use_core = min([num_cores, n])

    Result = Parallel(n_jobs=num_cores)(delayed(make_UpdateExpression)(HostName, GenesDF, Genome_WT, Genome_MT, myRct.id) for myRct in Model.reactions)
    RctNew_df = pd.DataFrame(Result)

    return RctNew_df

def Help_FluxCalculator ( HostName:str, StrainWT:Strain, StrainMut:Optional[Strain] = None ) :
    '''
    Calculation of flux values.

    This method can work in 2 modes:
      [StrainWT only] : The metabolic model of the WT strain is used for calculation.
      [StrainWT + StrainMut] : An additional step of resetting boundaries is done on the model
        before the fluxes are calculated.

    TODO: The "reset boundary" step is mutating the `StrainMut.model`, mutation might not be intended.
    '''

    # adding flux values
    # setup of flux boundaries. For the reference boundary changes are set to 'False',
    # for mutant strains, ractions with altered promoter sequence will change enzyme levels and boundaries must be changed accordingly, their variable is 'True'

    if StrainMut is not None :
        print('resetting boundaries')
        # finding reactions for which the expression has changed
        RctNewDF, Set_Boundary, _ = measure_EnzymeLevel1(HostName, StrainWT, StrainMut)
        # Defining the model with the two combinations of either
        # increasing lower bound (increased forward, decreased reverse reaction)
        # decreasing upper bound (decreased forward, increased reverse reaction)
        with StrainWT.model as myModel:
            # Comb.1: positive flux with increased expression -> increasing lower bound
            for Indx in Set_Boundary['lower']:
                myModel.reactions[Indx].lower_bound = RctNewDF.loc[Indx, 'NewExpr'] * RctNewDF.loc[Indx, 'Expr2Flux']
            # Comb.2: positive flux with decreased expression -> decreasing upper bound
            for Indx in Set_Boundary['upper']:
                myModel.reactions[Indx].upper_bound = RctNewDF.loc[Indx, 'NewExpr'] * RctNewDF.loc[Indx, 'Expr2Flux']

            Fluxes = myModel.optimize()

    else:
        Fluxes = StrainWT.model.optimize()


    return Fluxes.fluxes.values, Fluxes.objective_value

def Help_Expr2Flux(GenesDF):
    '''
    Function to correlate expression strength with flux.
    Input:
        GenesDF:        DataFrame, details of enzymes
    Output:
        Corr_ExprFlux:  array, correlation factor
    '''
    Corr_ExprFlux = GenesDF['Fluxes'].values/GenesDF['Expression'].values

    return Corr_ExprFlux


def make_AdaptModel(StrainID):
    '''
    Function to adapt the Cobra-model boundaries according to the expression strength of the genome promoters.
    Input:
        :         parent class, subfield _Mutant_GenesDF and genome used.
        StrainID:       string, identifier for expressed genome in subclass var_strains
    Output:
        Model:          cobraby model, E. coli core reactions with updated boundary conditions
    '''


    return Model


def measure_EnzymeLevel1(HostName:str, StrainWT:Strain, StrainMut:Strain):
    '''
    The function differences of expression levels of enzymes between two strains.
    '''
    import numpy as np
    from .metabolism import Help_StrainCharacterizer

    RefGenDF = StrainWT.genes_df
    RefGenome = str(StrainWT.genome)
    MutGenome = str(StrainMut.genome)
    RefModel = StrainWT.model

    RctNewDF = Help_StrainCharacterizer(HostName, RefGenDF, RefGenome, MutGenome, RefModel)

    # finding reactions for which the expression has changed
    RctNew = RctNewDF[RctNewDF['RctFlag']==True].index.values
    # For reactions with reduced expression and positive flux: reduce the upper limit,
    # For reactions with increased expression and positive flux: increase the lower limit
    # for reactions with negative flux the limits are exchanged.
    FluxPos = RctNew[tuple([RctNewDF.loc[RctNew, 'RefFlux']>0])]
    FluxNeg = RctNew[tuple([RctNewDF.loc[RctNew, 'RefFlux']<0])]
    # Finding increased and decreased fluxes
    FluxInc = RctNew[RctNewDF.loc[RctNew, 'NewExpr'].values / RctNewDF.loc[RctNew, 'RefExpr'].values>1]
    FluxDec = RctNew[RctNewDF.loc[RctNew, 'NewExpr'].values / RctNewDF.loc[RctNew, 'RefExpr'].values<1]

    # Comb.1: positive flux with increased expression -> increasing lower bound
    PosIncInd = np.intersect1d(FluxPos,FluxInc)
    # Comb.2: positive flux with decreased expression -> decreasing upper bound
    PosDecInd = np.intersect1d(FluxPos,FluxDec)
    # Comb.3: negative flux with increased expression -> decreasing lower bound
    NegIncInd = np.intersect1d(FluxNeg,FluxInc)
    # Comb.4: positive flux with increased expression -> increasing lower bound
    NegDecInd = np.intersect1d(FluxNeg,FluxDec)

    Expr_Change = {'increase': np.hstack([PosIncInd,NegIncInd]),'decrease': np.hstack([PosDecInd,NegDecInd])}
    Set_Boundary = {'lower': np.hstack([PosIncInd,NegDecInd]),'upper': np.hstack([PosDecInd,NegIncInd])}

    return RctNewDF, Set_Boundary, Expr_Change
