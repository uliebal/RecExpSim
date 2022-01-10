
def measure_BaseCompare(Seq1, Seq2):
    '''
    Comparison of two sequences from start to end. returns positions of base differences.
    '''
    SeqDiff = [[count, Pos] for count, Pos in enumerate(zip(Seq1, Seq2)) if Pos[0] != Pos[1]]

    return SeqDiff


# TODO: This file also occurs in `biolabsim.simulation.expression` with a very similar code.
def Help_PromoterStrength(Host, Sequence, Scaler=1, Similarity_Thresh=.4, Predict_File=None):
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

    from ..auxfun import Sequence_ReferenceDistance, list_onehot, list_integer

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
