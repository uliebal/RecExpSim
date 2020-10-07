def Help_LoadCobra(Path=False):
    '''
    This function loads cobra models from input path, or by default the Ecoli core model. Either from a local copy or dowloading it from the Bigg database.
    
    Input:
        Path (optional): string, relative path address to local copy
        
    Output:
        Model:           cobraby model, E. coli core reactions
    '''
    
    import os
    from cobra.io import read_sbml_model    
    
    if Path:
        Model = read_sbml_model(Path)
    else:
        os.system('wget http://bigg.ucsd.edu/static/models/e_coli_core.xml')
        # generating cobra variable from SBML/xml/mat file
        Path = 'e_coli_core.xml'
        Model = read_sbml_model(Path)
        
    return Model

def Help_GeneAnnotator(Host, Model):
    '''
    Function to characterize genes.
    Input
        Model:    cobra model
    Output
        Genes_df: dataframe, gene name from enzyme id in model, expression strength, promoter sequence, ORF, flux
    '''
    import pandas as pd
    import multiprocessing
    from joblib import Parallel, delayed
    from BioLabSim.ModuleVirtualOrganism.Genome import make_GeneJoiner
    
    num_cores = multiprocessing.cpu_count()
#     use_core = min([num_cores, n])
  
    Result = Parallel(n_jobs=num_cores)(delayed(make_GeneJoiner)(Host, Model, myRct.id) for myRct in Model.reactions)
    Genes_df = pd.DataFrame(Result)
    
    return Genes_df

def Help_StrainCharacterizer(Host, GenesDF, Genome_WT, Genome_MT, Model):
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
    from BioLabSim.ModuleVirtualOrganism.Expression import make_UpdateExpression
    
    num_cores = multiprocessing.cpu_count()
#     use_core = min([num_cores, n])
  
    Result = Parallel(n_jobs=num_cores)(delayed(make_UpdateExpression)(Host, GenesDF, Genome_WT, Genome_MT, myRct.id) for myRct in Model.reactions)
    RctNew_df = pd.DataFrame(Result)

    return RctNew_df

def Help_FluxCalculator(Strain, Host, Test_Strain=False):
    '''
    Calculation of flux values.
    '''
    import numpy as np
    from BioLabSim.ModuleMeasureOrganism.Fluxes import measure_EnzymeLevel1
    
    # adding flux values
    # setup of flux boundaries. For the reference boundary changes are set to 'False', 
    # for mutant strains, ractions with altered promoter sequence will change enzyme levels and boundaries must be changed accordingly, their variable is 'True'

    Model = Strain.var_Model
    
    if Test_Strain:
        print('resetting boundaries')
        # finding reactions for which the expression has changed
        RctNewDF, Set_Boundary, _ = measure_EnzymeLevel1(Host, Strain)
        # Defining the model with the two combinations of either 
        # increasing lower bound (increased forward, decreased reverse reaction)
        # decreasing upper bound (decreased forward, increased reverse reaction)
        with Model as myModel:
            # Comb.1: positive flux with increased expression -> increasing lower bound
            for Indx in Set_Boundary['lower']:
                myModel.reactions[Indx].lower_bound = RctNewDF.loc[Indx, 'NewExpr'] * RctNewDF.loc[Indx, 'Expr2Flux']
            # Comb.2: positive flux with decreased expression -> decreasing upper bound
            for Indx in Set_Boundary['upper']:
                myModel.reactions[Indx].upper_bound = RctNewDF.loc[Indx, 'NewExpr'] * RctNewDF.loc[Indx, 'Expr2Flux']           
                
            Fluxes = myModel.optimize()
        
    else:
        Fluxes = Model.optimize()
        

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
        :         parent class, subfield _Mutant_GenesDF and info_Genome used.
        StrainID:       string, identifier for expressed genome in subclass var_strains
    Output:
        Model:          cobraby model, E. coli core reactions with updated boundary conditions
    '''
    
    
    return Model