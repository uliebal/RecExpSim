# Functions to report flux associated information

def measure_EnzymeLevel1(Host, Strain):
    '''
    The function differences of expression levels of enzymes between two strains.
    '''
    import numpy as np
    from BioLabSim.ModuleVirtualOrganism.Metabolism import Help_StrainCharacterizer 
    
    Host_Strain = Host.var_Host
    RefGenDF = Host.WT.info_GenesDF
    RefGenome = Host.WT.info_Genome
    MutGenome = Strain.info_Genome
    RefModel = Host.WT.var_Model
    
    RctNewDF = Help_StrainCharacterizer(Host_Strain, RefGenDF, RefGenome, MutGenome, RefModel)
    
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