print("\n--> Initializing SAPPHIRE...\n")
print("\n--> Importing packages...\n")
print("")
print("---------------------")

import sys
import re
import numpy as np
import csv
import operator

from keras.models import Sequential, load_model
from keras import backend as K

from bisect import bisect

from Bio import SeqIO
print("")
print("---------------------")

prob_rand_file_sig70 = "sig70NN_pos_probabilities.txt"

PROB_POS_SIG70 = 0.0056

def one_hot_encode_bases(nparr):
    result = []
    for base in nparr:
        if base == "A":
            result.append(np.array([1,0,0,0]))
        elif base == "C":
            result.append(np.array([0,1,0,0]))
        elif base == "G":
            result.append(np.array([0,0,1,0]))
        elif base == "T":
            result.append(np.array([0,0,0,1]))
        elif base == "a":
            result.append(np.array([1,0,0,0]))
        elif base == "c":
            result.append(np.array([0,1,0,0]))
        elif base == "g":
            result.append(np.array([0,0,1,0]))
        elif base == "t":
            result.append(np.array([0,0,0,1]))
        else:
            raise Exception("NoACGT")

    return np.array(result)

def predict_sig(model, sequence):
    # SPACER 18
    input = []
    format1 = [0, 1, 2, 3, 4, 5, 24, 25, 26, 27, 28, 29]
    for pos in format1:
        input.append(sequence[pos])

    bases_encoded = one_hot_encode_bases(np.array(input))
    prediction1 = model.predict(np.array([bases_encoded]))[0][0]

    # SPACER 17
    input = []
    format2 = [0, 1, 2, 3, 4, 5, 23, 24, 25, 26, 27, 28]
    for pos in format2:
        input.append(sequence[pos])

    bases_encoded = one_hot_encode_bases(np.array(input))
    prediction2 = model.predict(np.array([bases_encoded]))[0][0]

    # SPACER 16
    input = []
    format3 = [0, 1, 2, 3, 4, 5, 22, 23, 24, 25, 26, 27]
    for pos in format3:
        input.append(sequence[pos])

    bases_encoded = one_hot_encode_bases(np.array(input))
    prediction3 = model.predict(np.array([bases_encoded]))[0][0]

    return [prediction1,prediction2,prediction3]

def calc_p_value(prediction, prob_rand_list):
    """
    Calculates the probability of this predicted score in case the Null hypothesis (the analysed sequence is random)
    is true
    """
    return round(PROB_POS_SIG70*float(10000-bisect(prob_rand_list, prediction))/10000,8)

def evaluate_sequence(model, sequence, prob_rand_list, seq_name):
    predictions_list = []

    k = 0
    if seq_name[19:] == "_reverse_complement":
        k = len(sequence)

    for i in range(0, len(sequence)-30):
        test_sequence = sequence[i:i+30]

        local_predictions = predict_sig(model, test_sequence)
        if local_predictions[0] > 0.5:
            predictions_list.append({"location": abs(k-(i+35)), "p-value":calc_p_value(local_predictions[0], prob_rand_list), "sequence": str(sequence), "spacer": 18, "spacer_seq": str(sequence[i+6:i+24]), '-35box': str(sequence[i+0:i+6]), '-10box': str(sequence[i+24:i+30])})
        if local_predictions[1] > 0.5:
            predictions_list.append({"location": abs(k-(i+34)), "p-value":calc_p_value(local_predictions[1], prob_rand_list), "sequence": str(sequence), "spacer": 17, "spacer_seq": str(sequence[i+6:i+23]), '-35box': str(sequence[i+0:i+6]), '-10box': str(sequence[i+23:i+29])})
        if local_predictions[2] > 0.5:
            predictions_list.append({"location": abs(k-(i+33)), "p-value":calc_p_value(local_predictions[2], prob_rand_list), "sequence": str(sequence), "spacer": 16, "spacer_seq": str(sequence[i+6:i+22]), '-35box': str(sequence[i+0:i+6]), '-10box': str(sequence[i+22:i+28])})

    return predictions_list

def check_valid_sequence(sequence):
    valid = "acgtACGT"
    if not all(c in valid for c in sequence):
        raise Exception("NoACGT")
    if len(sequence) < 30:
        raise Exception("No30")

def specificity(y_pred, y_true):
    neg_y_true = 1 - y_true
    neg_y_pred = 1 - y_pred
    fp = K.sum(neg_y_true * y_pred)
    tn = K.sum(neg_y_true * neg_y_pred)
    specificity = tn / (tn + fp + K.epsilon())
    return specificity

def select_unique_highest(sequences_results):
    sequences_results_out = []
    for j in range(0, len(sequences_results["sequences_results"])):
        seq_preds = sequences_results["sequences_results"][j]["seq_preds"]
        seq_preds_out = []
        i=0
        while i < len(seq_preds):
            if i < len(seq_preds)-1:
                if seq_preds[i+1]["location"]-seq_preds[i]["location"]<2:
                    if i < len(seq_preds) - 2:
                        if seq_preds[i+2]["location"]-seq_preds[i+1]["location"]<2:
                            if seq_preds[i+2]["p-value"] < min(seq_preds[i + 1]["p-value"],seq_preds[i]["p-value"]):
                                seq_preds_out.append(seq_preds[i+2])
                            elif seq_preds[i+1]["p-value"] < min(seq_preds[i + 2]["p-value"],seq_preds[i]["p-value"]):
                                seq_preds_out.append(seq_preds[i+1])
                            else:
                                seq_preds_out.append(seq_preds[i])
                            i+=3
                        else:
                            if seq_preds[i]["p-value"] < seq_preds[i+1]["p-value"]:
                                seq_preds_out.append(seq_preds[i])
                            else:
                                seq_preds_out.append(seq_preds[i+1])

                            i+=2
                    else:
                        if seq_preds[i]["p-value"] < seq_preds[i+1]["p-value"]:
                            seq_preds_out.append(seq_preds[i])
                        else:
                            seq_preds_out.append(seq_preds[i+1])

                        i+=2
                else:
                    seq_preds_out.append(seq_preds[i])
                    i += 1
            else:
                seq_preds_out.append(seq_preds[i])
                i += 1
        sequences_results_out.append({"seq_name":sequences_results["sequences_results"][j]["seq_name"], "seq_len":sequences_results["sequences_results"][j]["seq_len"], "seq_preds":seq_preds_out})
    return {"sequences_results": sequences_results_out}


def main():
    try:
        prob_rand_list = []
        with open(prob_rand_file_sig70, 'r') as f:
            for line in f:
                prob_rand_list.append(float(line))

        print("\n--> Loading NN model...\n")

        model = load_model("sig70_model.h5", custom_objects={"specificity": specificity})

        include_reverse = False
        if len(sys.argv) == 3:
            if sys.argv[1] == "-r":
                include_reverse = True
                print("\n--> Loading input fasta file...\n")
                target_file = sys.argv[2]
            else:
                raise Exception('Commandformat')
        elif len(sys.argv) == 2:
            print("\n--> Loading input fasta file...\n")
            target_file = sys.argv[1]
        else:
            raise Exception('Commandformat')

        try:
            fasta_sequences = SeqIO.parse(open(target_file), 'fasta')
        except:
            raise Exception('Fileformat')

        fasta_sequences_list = list(fasta_sequences)

        if len(fasta_sequences_list)==0:
            raise Exception('Fileformat')


        if include_reverse:
            new_list = []
            for sequence in fasta_sequences_list:
                new_list.append(sequence)
                new_list.append(sequence.reverse_complement(name=sequence.name+"_reverse_complement"))

            fasta_sequences_list = new_list

        sequences_results = {"sequences_results": []}
        print("\n--> Evaluating sequences...\n")
        for sequence in fasta_sequences_list:
            check_valid_sequence(sequence)
            results = evaluate_sequence(model, sequence.seq, prob_rand_list, sequence.name)
            if not results == []:
                sequences_results["sequences_results"].append({"seq_name":sequence.name, "seq_len":len(sequence), "seq_preds":results})

        sequences_results = select_unique_highest(sequences_results)

        sequences_results = select_unique_highest(sequences_results)

        tar_out_file="SAPPHIRE_out.tsv"

        with open(tar_out_file, 'w') as outfile:
            outfile.write("Sequence_name\tstrand\tTSS\tp-value\t-35box\tspacer\t-10box\n")
            for result in sequences_results["sequences_results"]:
                strand = "+1"
                name = result['seq_name']
                if name[-19:] == "_reverse_complement":
                    strand = "-1"
                    name = name[:-19]
                for seq_preds in result['seq_preds']:
                    outfile.write(name + "\t" + strand + "\t" + str(seq_preds['location']) + "\t" + str(
                        format(seq_preds['p-value'], 'f')) + "\t" + str(seq_preds['-35box']) + "\t" + str(
                        seq_preds["spacer_seq"]) + "\t" + str(seq_preds['-10box']) + "\n")

        data = csv.reader(open(tar_out_file), delimiter='\t')
        headers = next(data)
        sortedlist = sorted(data, key=operator.itemgetter(3))

        # Opening and sorting of output file
        with open(tar_out_file, 'w') as outfile:
            fileWriter = csv.writer(outfile, delimiter='\t')
            fileWriter.writerow(headers)
            for row in sortedlist:
                fileWriter.writerow(row)
        print("\n--> DONE!\n")
        print("Found promoters in " + str(len(sequences_results["sequences_results"]))+" of the submitted sequences")
        print("Output written to SAPPHIRE_out.tsv")


    except Exception as e:
        print("\n--> An error occured: "+str(e)+"\n")
        if str(e)[:10]=='Fileformat':
            print("The provided input is not valid. Make sure following properties hold for the input file:")
            print("    * Your input fasta file is correctly formatted")
            print("    * Your input DNA sequences only contain A, C, G or T")
            print("    * Your input DNA sequences are longer than 30 basepairs")
        elif str(e)[:13]=='Commandformat':
            print("Invalid command. Correct usage:")
            print("  python3 SAPPHIRE.py input.fasta")
            print("  INCLUDE REVERSE COMPLEMENT:")
            print("  python3 SAPPHIRE.py -r input.fasta")
        else:
            raise e


main()

exit()
