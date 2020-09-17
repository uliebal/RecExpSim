def Help_PromoterStrength(Host, Sequence, Scaler=100, Similarity_Thresh=.4, Predict_File=None):
    '''Expression of the recombinant protein.
        Arguments:
            Host:       class, contains optimal growth temperature, production phase
            Sequence:     string, Sequence for which to determine promoter strength
            Scaler:       int, multiplied to the regression result for higher values
            Predict_File: string, address of regression file
        Output: 
            Expression: float, expression rate
    '''
    import os
    import numpy as np
    import joblib
    import pickle
    from BioLabSim.ModuleVirtualOrganism.Genome import list_onehot, list_integer
    
    if Sequence_ReferenceDistance(Sequence) > Similarity_Thresh:
        Expression = 0
    else:
        if Predict_File!=None:
            Regressor_File = Predict_File
        else:    
            Data_Folder = 'ExpressionPredictor'
            if Host == 'Ecol':
                Regressor_File = os.path.join(Data_Folder,'Ecol-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ecol-Promoter-AddParams.pkl')
#                 Scaler_DictName = 'Ecol Promoter Activity_Scaler'
            elif Host == 'Pput':
                Regressor_File = os.path.join(Data_Folder,'Ptai-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ptai-Promoter-AddParams.pkl')
#                 Scaler_DictName = 'Ptai Promoter Activity_Scaler'
            else:
                print('Non-recognized host name. Rename host to either "Ecol" or "Pput."')

        Predictor = joblib.load(Regressor_File)
        Params = pickle.load(open(Add_Params, 'rb'))
        Positions_removed = Params['Positions_removed']
#         Expr_Scaler = Params[Scaler_DictName]

        X = np.array(list_onehot(np.delete(list_integer(Sequence),Positions_removed, axis=0))).reshape(1,-1)  
        GC_cont = (Sequence.count('G') + Sequence.count('C'))/len(Sequence)
        X = np.array([np.append(X,GC_cont)])
        Y = Predictor.predict(X)
        Expression = round(float(Y)*Scaler,3)

    return Expression

def make_UpdateExpression(Host, GenesDF, Genome_WT, Genome_Mut, RctID):
    '''
    Function to evaluate if expression for the input gene is different compared to the reference strain. 
    Input:
        Mutant:         parent class, subfield _Mutant_GenesDF and info_Genome used.
        StrainID:       string, identifier for expressed genome in subclass var_strains
        RctID:          string, identifier for the reaction to be tested
    Output:
        OutDict:        dictionary containing 'RctFlag' boolean, if True, the expression relative to the reference has changed; 'Gene_Activity' float, expression strength; 'RctID' string, identifier for the reaction tested
    '''
    
    import numpy as np
    
    # Extraction of ORF and reference promoter sequence of input reaction from reference annotation database
    Rct_Idx = GenesDF[GenesDF['RctID']==RctID].index.values
    Rct_ORF = ''.join(str(n) for n in GenesDF['ORF'].iloc[Rct_Idx].to_list())
    Rct_RefProm = GenesDF['Promoter'].iloc[Rct_Idx].values

    # Extraction of reaction association promoter sequence from input genome
    Gene_ORF = Genome_WT.find(Rct_ORF)
    Gene_Promoter = Genome_Mut[Gene_ORF-40:Gene_ORF:1]
    RefExpr = float(GenesDF['Expression'].iloc[Rct_Idx].values)

    # Comparison of current promoter with reference promoter
    # If promoters differ, then new expression strength is calculated
    if Gene_Promoter != Rct_RefProm:
#         print('New promoter.')
        RctFlag = True
        NewExpr = float(Help_PromoterStrength(Host, Gene_Promoter, Similarity_Thresh=.8))
    else:
        RctFlag = False
        NewExpr = RefExpr
#     print(Gene_Activity)
    OutDict = {'RctFlag':RctFlag, 'RctID':RctID, 'RefExpr':RefExpr, 'NewExpr':NewExpr, 'RefFlux':float(GenesDF.loc[Rct_Idx, 'Fluxes']), 'Expr2Flux':float(GenesDF.loc[Rct_Idx, 'Expr2Flux']) }
    
    return OutDict


def make_Promoter(RefFile=False, WeightFile=False):
    '''
    Generating a single promoter. The promoter sequences are selected from the exploration space used to train the regressor for promoter strength. The algorithm needs to know which nucleotides at each position were sampled strongly enough for predictions. Positions not sufficiently well sampled are filled with entries from a reference sequence.
    Input
        RefFile:               string, path to reference sequence in txt-file, used to fill positions with weak sampling
        WeightFile:            string, path to pickle with dataframe, columns represent bases (A-C-G-T), rows represent positions, index is position relative to start codon (e.g. -40, -1)
    Output
        Gene_Promoter:        string, sequence of the promoter, same length as reference sequence
    '''
    import numpy as np
    import pickle
    import os
    import random
    # generating promoter sequence

    # the reference sequence contains the most common tested nucleotides at each position
    if not RefFile:   
        RefFile = os.path.join('Models','RefSeq.txt')
    with open(RefFile) as f:
        RefSeq = f.read()
    
    # loading information that determines the exploratory space of the regressor
    # weight file has boolean representations for each base on each position whether it was part of the training set or not
    if not WeightFile:
        WeightFile = os.path.join('Models','NucleotideWeightTable.pkl')
    with open(WeightFile, 'rb') as handle:
        Nucleotides_Weight = pickle.load(handle)
    Nucleotides_Weight.index = Nucleotides_Weight.index+40

    # setup of default sequence with reference
    SeqDef = np.asarray([Letter for Letter in RefSeq])

    # changing positions within the exploratory space of the regressor
    Bases = ['A','C','G','T']
    for idx in range(len(Nucleotides_Weight)):
        np.put(SeqDef,Nucleotides_Weight.index[idx],random.choices(Bases,Nucleotides_Weight.iloc[idx].values))
    Gene_Promoter = ''.join(SeqDef)
    
    return Gene_Promoter

def Sequence_ReferenceDistance(SeqObj, RefSeq=None):
    '''Returns the genetic sequence distance to a reference sequence.
    Input:
           SeqDF: list, the sequence in conventional letter format
    Output:
           SequenceDistance: float, genetic distances as determined from the sum of difference in bases divided by total base number, i.e. max difference is 1, identical sequence =0
    '''
    import numpy as np
    
    if RefSeq != None:
        RefSeq = SeqObj[0]
    else:
        RefSeq = 'GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG'
        
    Num_Samp = len(SeqObj)
    SequenceDistance = np.sum([int(seq1!=seq2) for seq1,seq2 in zip(RefSeq, SeqObj)], dtype='float')/len(SeqObj)
    
    return SequenceDistance
