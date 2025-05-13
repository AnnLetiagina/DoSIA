
import pandas as pd

r2_17_24 = 0.4558 
r2_21_28 = 0.5967
r2_25_32 = 0.511
r2_29_36 = 0.6641
r2_33_40 = 0.5455
r2_37_44 = 0.3797
r2_41_48 = 0.2674
r2_45_52 = 0.3145
r2_sum = r2_17_24 + r2_21_28 + r2_25_32 + r2_29_36 + r2_33_40 + r2_37_44 + r2_41_48 + r2_45_52
    
pred_17_24 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\17-24 log2 bc predicted.xlsx')
pred_21_28 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\21-28 log2 bc predicted.xlsx')
pred_25_32 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\25-32 log2 barc corr predicted.xlsx')
pred_29_36 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\29-36 log2 barc corr predicted.xlsx')
pred_33_40 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\33-40 log2 bc predicted.xlsx')
pred_37_44 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\37-44 log2 bc predicted.xlsx')
pred_41_48 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\41-48 log2 bc predicted.xlsx')
pred_45_52 = pd.read_excel('C:\\Users\\Alexey\\Documents\\Аспирантура\\ML models\\45-52 log2 bc predicted.xlsx')

mut_to_pred = 'catgcgtcaattttacgcatgattatctttaacgtacgtc'
mut_to_pred = mut_to_pred.upper()

if len(mut_to_pred) == 36:
    
    mut_17_24 = mut_to_pred[:8]
    mut_21_28 = mut_to_pred[4:12]
    mut_25_32 = mut_to_pred[8:16]
    mut_29_36 = mut_to_pred[12:20]
    mut_33_40 = mut_to_pred[16:24]
    mut_37_44 = mut_to_pred[20:28]
    mut_41_48 = mut_to_pred[24:32]
    mut_45_52 = mut_to_pred[28:]
    
    expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
    expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
    expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
    expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
    expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
    expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
    expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
    expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
    
    pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
    print(pred_expr)
else:
    print('Can predict only 36 bp length sequences.')
    

pred_17_24.prediction_label[0]
pred_17_24.Mutation[0:100]

ascen_pred_17_24 = pred_17_24.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_21_28 = pred_21_28.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_25_32 = pred_25_32.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_29_36 = pred_29_36.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_33_40 = pred_33_40.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_37_44 = pred_37_44.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_41_48 = pred_41_48.sort_values(by = 'prediction_label', 
                                          ignore_index = True)
ascen_pred_45_52 = pred_45_52.sort_values(by = 'prediction_label', 
                                          ignore_index = True)


def find_max(mut_start, lib, mode = 'forward', end = ''):
    mut_list = []
    for mutation in lib['Mutation']:
        if len(mut_list) >= 21:
            break
        if mode == 'forward':
            if end:
                if mut_start == mutation[0:4] and mutation[-len(end):] == end:
                    mut_list.append(mutation)
            else:
                if mut_start == mutation[0:4]:
                    mut_list.append(mutation)
        elif mode == 'reverse':
            if mut_start == mutation[4:]:
                mut_list.append(mutation)
    return mut_list

checked_mutations = {'21-28': [], '25-32': [], '29-36': [], '33-40': [], '37-44': [], '41-48': [], '45-52': []}
best_mutations = []
short_best_mutations = []

#Set len(mut_list) >= 21
for mutation in pred_17_24.Mutation[0:146]:
    best_mutation = mutation
    mutation_end = mutation[4:]
    if mutation_end not in checked_mutations['21-28']:
        checked_mutations['21-28'].append(mutation_end) 
        list_21_28 = find_max(mutation_end, pred_21_28)
        for mut_21_28 in list_21_28:
            mut_25_32_start = mut_21_28[4:]
            best_mutation += mut_25_32_start
            if mut_25_32_start not in checked_mutations['25-32']:
                checked_mutations['25-32'].append(mut_25_32_start)
                list_25_32 = find_max(mut_25_32_start, pred_25_32)
                for mut_25_32 in list_25_32:
                    mut_29_36_start = mut_25_32[4:]
                    best_mutation += mut_29_36_start
                    if mut_29_36_start not in checked_mutations['29-36']:
                        checked_mutations['29-36'].append(mut_29_36_start)
                        list_29_36 = find_max(mut_29_36_start, pred_29_36)
                        for mut_29_36 in list_29_36:
                            mut_33_40_start = mut_29_36[4:]
                            best_mutation += mut_33_40_start
                            if mut_33_40_start not in checked_mutations['33-40']:
                                checked_mutations['33-40'].append(mut_33_40_start)
                                list_33_40 = find_max(mut_33_40_start, pred_33_40)
                                for mut_33_40 in list_33_40:
                                    mut_37_44_start = mut_33_40[4:]
                                    best_mutation += mut_37_44_start
                                    if mut_37_44_start not in checked_mutations['37-44']:
                                        checked_mutations['37-44'].append(mut_37_44_start)
                                        list_37_44 = find_max(mut_37_44_start, pred_37_44)
                                        for mut_37_44 in list_37_44:
                                            mut_41_48_start = mut_37_44[4:]
                                            best_mutation += mut_41_48_start
                                            if mut_41_48_start not in checked_mutations['41-48']:
                                                checked_mutations['41-48'].append(mut_41_48_start)
                                                list_41_48 = find_max(mut_41_48_start, pred_41_48)
                                                for mut_41_48 in list_41_48:
                                                    mut_45_52_start = mut_41_48[4:]
                                                    best_mutation += mut_45_52_start
                                                    if mut_45_52_start not in checked_mutations['45-52']:
                                                        checked_mutations['45-52'].append(mut_45_52_start)
                                                        list_45_52 = find_max(mut_45_52_start, pred_45_52, end = 'GTA')
                                                        for mut_45_52 in list_45_52:
                                                            best_mutation += mut_45_52[4:]
                                                            best_mutations.append(best_mutation)
                                                            best_mutation = best_mutation[:-4]
                                                        else:
                                                            best_mutation = best_mutation[:-4]
                                                    else:
                                                        short_best_mutations.append(best_mutation)
                                                        best_mutation = best_mutation[:-4]
                                                else:
                                                    best_mutation = best_mutation[:-4]
                                            else:
                                                short_best_mutations.append(best_mutation)
                                                best_mutation = best_mutation[:-4]
                                        else:
                                            best_mutation = best_mutation[:-4]
                                    else:
                                        short_best_mutations.append(best_mutation)
                                        best_mutation = best_mutation[:-4]
                                else:
                                    best_mutation = best_mutation[:-4]
                            else:
                                short_best_mutations.append(best_mutation)
                                best_mutation = best_mutation[:-4]
                        else:
                            best_mutation = best_mutation[:-4]
                    else:
                        short_best_mutations.append(best_mutation)
                        best_mutation = best_mutation[:-4]
                else:
                    best_mutation = best_mutation[:-4]
            else:
                short_best_mutations.append(best_mutation)
                best_mutation = best_mutation[:-4]
        else:
            best_mutation = best_mutation[:-4]


elongated_mutations = []
for short_mutation in short_best_mutations:
    checked_motive = short_mutation[-4:]
    short_mot_len = len(short_mutation)
    for best_mutation in best_mutations:
        if best_mutation[short_mot_len - 4: short_mot_len] == checked_motive:
            elongated_mutations.append(short_mutation + best_mutation[short_mot_len:])
            
best_mutations += elongated_mutations


best_expr = 0        
for best_mut in best_mutations:
    mut_to_pred = best_mut
    if len(mut_to_pred) == 36:
        
        mut_17_24 = mut_to_pred[:8]
        mut_21_28 = mut_to_pred[4:12]
        mut_25_32 = mut_to_pred[8:16]
        mut_29_36 = mut_to_pred[12:20]
        mut_33_40 = mut_to_pred[16:24]
        mut_37_44 = mut_to_pred[20:28]
        mut_41_48 = mut_to_pred[24:32]
        mut_45_52 = mut_to_pred[28:]
        
        expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
        expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
        expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
        expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
        expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
        expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
        expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
        expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
        
        pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
        if pred_expr >= best_expr:
            best_expr = pred_expr
            print(pred_expr, mut_to_pred)



checked_low_mutations = {'21-28': [], '25-32': [], '29-36': [], '33-40': [], '37-44': [], '41-48': [], '45-52': []}
worst_mutations = []
short_worst_mutations = []

for mutation in ascen_pred_17_24.Mutation[0:100]:
    worst_mutation = mutation
    mutation_end = mutation[4:]
    if mutation_end not in checked_low_mutations['21-28']:
        checked_low_mutations['21-28'].append(mutation_end) 
        list_21_28 = find_max(mutation_end, ascen_pred_21_28)
        for mut_21_28 in list_21_28:
            mut_25_32_start = mut_21_28[4:]
            worst_mutation += mut_25_32_start
            if mut_25_32_start not in checked_low_mutations['25-32']:
                checked_low_mutations['25-32'].append(mut_25_32_start)
                list_25_32 = find_max(mut_25_32_start, ascen_pred_25_32)
                for mut_25_32 in list_25_32:
                    mut_29_36_start = mut_25_32[4:]
                    worst_mutation += mut_29_36_start
                    if mut_29_36_start not in checked_low_mutations['29-36']:
                        checked_low_mutations['29-36'].append(mut_29_36_start)
                        list_29_36 = find_max(mut_29_36_start, ascen_pred_29_36)
                        for mut_29_36 in list_29_36:
                            mut_33_40_start = mut_29_36[4:]
                            worst_mutation += mut_33_40_start
                            if mut_33_40_start not in checked_low_mutations['33-40']:
                                checked_low_mutations['33-40'].append(mut_33_40_start)
                                list_33_40 = find_max(mut_33_40_start, ascen_pred_33_40)
                                for mut_33_40 in list_33_40:
                                    mut_37_44_start = mut_33_40[4:]
                                    worst_mutation += mut_37_44_start
                                    if mut_37_44_start not in checked_low_mutations['37-44']:
                                        checked_low_mutations['37-44'].append(mut_37_44_start)
                                        list_37_44 = find_max(mut_37_44_start, ascen_pred_37_44)
                                        for mut_37_44 in list_37_44:
                                            mut_41_48_start = mut_37_44[4:]
                                            worst_mutation += mut_41_48_start
                                            if mut_41_48_start not in checked_low_mutations['41-48']:
                                                checked_low_mutations['41-48'].append(mut_41_48_start)
                                                list_41_48 = find_max(mut_41_48_start, ascen_pred_41_48)
                                                for mut_41_48 in list_41_48:
                                                    mut_45_52_start = mut_41_48[4:]
                                                    worst_mutation += mut_45_52_start
                                                    if mut_45_52_start not in checked_low_mutations['45-52']:
                                                        checked_low_mutations['45-52'].append(mut_45_52_start)
                                                        list_45_52 = find_max(mut_45_52_start, ascen_pred_45_52)
                                                        for mut_45_52 in list_45_52:
                                                            worst_mutations.append(worst_mutation)
                                                            worst_mutation = worst_mutation[:-4]
                                                        else:
                                                            worst_mutation = worst_mutation[:-4]
                                                    else:
                                                        short_worst_mutations.append(worst_mutation)
                                                        worst_mutation = worst_mutation[:-4]
                                                else:
                                                    worst_mutation = worst_mutation[:-4]
                                            else:
                                                short_worst_mutations.append(worst_mutation)
                                                worst_mutation = worst_mutation[:-4]
                                        else:
                                            worst_mutation = worst_mutation[:-4]
                                    else:
                                        short_worst_mutations.append(worst_mutation)
                                        worst_mutation = worst_mutation[:-4]
                                else:
                                    worst_mutation = worst_mutation[:-4]
                            else:
                                short_worst_mutations.append(worst_mutation)
                                worst_mutation = worst_mutation[:-4]
                        else:
                            worst_mutation = worst_mutation[:-4]
                    else:
                        short_worst_mutations.append(worst_mutation)
                        worst_mutation = worst_mutation[:-4]
                else:
                    worst_mutation = worst_mutation[:-4]
            else:
                short_worst_mutations.append(worst_mutation)
                worst_mutation = worst_mutation[:-4]
        else:
            worst_mutation = worst_mutation[:-4]
            
            
elongated_mutations = []
for short_mutation in short_worst_mutations:
    checked_motive = short_mutation[-4:]
    short_mot_len = len(short_mutation)
    for worst_mutation in worst_mutations:
        if worst_mutation[short_mot_len - 4: short_mot_len] == checked_motive:
            elongated_mutations.append(short_mutation + worst_mutation[short_mot_len:])
            
worst_mutations += elongated_mutations


worst_expr = 0        
for worst_mut in worst_mutations:
    mut_to_pred = worst_mut
    if len(mut_to_pred) == 40:
        
        mut_17_24 = mut_to_pred[:8]
        mut_21_28 = mut_to_pred[4:12]
        mut_25_32 = mut_to_pred[8:16]
        mut_29_36 = mut_to_pred[12:20]
        mut_33_40 = mut_to_pred[16:24]
        mut_37_44 = mut_to_pred[20:28]
        mut_41_48 = mut_to_pred[24:32]
        mut_45_52 = mut_to_pred[28:36]
        
        expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
        expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
        expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
        expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
        expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
        expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
        expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
        expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
        
        pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
        if pred_expr < worst_expr:
            worst_expr = pred_expr
            print(pred_expr, mut_to_pred)
            
            
#Prediction startes from 29-36 fragment. Set len(mut_list) >= 47!!!
checked_mutations = {'33-40': [], '37-44': [], '41-48': [], '45-52': []}
best_mutations = []
short_best_mutations = []

for mutation in pred_29_36.Mutation[0:245]:
    best_mutation = mutation
    mutation_end = mutation[4:]
    if mutation_end not in checked_mutations['33-40']:
        checked_mutations['33-40'].append(mutation_end) 
        list_33_40 = find_max(mutation_end, pred_33_40)
        for mut_33_40 in list_33_40:
            mut_37_44_start = mut_33_40[4:]
            best_mutation += mut_37_44_start
            if mut_37_44_start not in checked_mutations['37-44']:
                checked_mutations['37-44'].append(mut_37_44_start)
                list_37_44 = find_max(mut_37_44_start, pred_37_44)
                for mut_37_44 in list_37_44:
                    mut_41_48_start = mut_37_44[4:]
                    best_mutation += mut_41_48_start
                    if mut_41_48_start not in checked_mutations['41-48']:
                        checked_mutations['41-48'].append(mut_41_48_start)
                        list_41_48 = find_max(mut_41_48_start, pred_41_48)
                        for mut_41_48 in list_41_48:
                            mut_45_52_start = mut_41_48[4:]
                            best_mutation += mut_45_52_start
                            if mut_45_52_start not in checked_mutations['45-52']:
                                checked_mutations['45-52'].append(mut_45_52_start)
                                list_45_52 = find_max(mut_45_52_start, pred_45_52)
                                for mut_45_52 in list_45_52:
                                    best_mutation += mut_45_52[4:]
                                    best_mutations.append(best_mutation)
                                    best_mutation = best_mutation[:-4]
                                else:
                                    best_mutation = best_mutation[:-4]
                            else:
                                short_best_mutations.append(best_mutation)
                                best_mutation = best_mutation[:-4]
                        else:
                            best_mutation = best_mutation[:-4]
                    else:
                        short_best_mutations.append(best_mutation)
                        best_mutation = best_mutation[:-4]
                else:
                    best_mutation = best_mutation[:-4]
            else:
                short_best_mutations.append(best_mutation)
                best_mutation = best_mutation[:-4]
        else:
            best_mutation = best_mutation[:-4]


elongated_mutations = []
for short_mutation in short_best_mutations:
    checked_motive = short_mutation[-4:]
    short_mot_len = len(short_mutation)
    for best_mutation in best_mutations:
        if best_mutation[short_mot_len - 4: short_mot_len] == checked_motive:
            elongated_mutations.append(short_mutation + best_mutation[short_mot_len:])
            
best_mutations += elongated_mutations


checked_mutations = {'17-24': [], '21-28': [], '25-32': []}
full_best_mutations = []
short_best_mutations = []

for mutation in best_mutations:
    best_mutation = mutation
    mutation_start = mutation[:4]
    if mutation_start not in checked_mutations['25-32']:
        checked_mutations['25-32'].append(mutation_start) 
        list_25_32 = find_max(mutation_start, pred_25_32, 'reverse')
        for mut_25_32 in list_25_32:
            mut_21_28_end = mut_25_32[:4]
            best_mutation = mut_21_28_end + best_mutation
            if mut_21_28_end not in checked_mutations['21-28']:
                checked_mutations['21-28'].append(mut_21_28_end)
                list_21_28 = find_max(mut_21_28_end, pred_21_28, 'reverse')
                for mut_21_28 in list_21_28:
                    mut_17_24_end = mut_21_28[:4]
                    best_mutation = mut_17_24_end + best_mutation
                    if mut_17_24_end not in checked_mutations['17-24']:
                        checked_mutations['17-24'].append(mut_17_24_end)
                        list_17_24 = find_max(mut_17_24_end, pred_17_24, 'reverse')
                        for mut_17_24 in list_17_24:
                            best_mutation = mut_17_24[:4] + best_mutation
                            full_best_mutations.append(best_mutation)
                            best_mutation = best_mutation[4:]
                        else:
                            best_mutation = best_mutation[4:]
                    else:
                        short_best_mutations.append(best_mutation)
                        best_mutation = best_mutation[4:]
                else:
                    best_mutation = best_mutation[4:]
            else:
                short_best_mutations.append(best_mutation)
                best_mutation = best_mutation[4:]
        else:
            best_mutation = best_mutation[4:]
    else:
        short_best_mutations.append(best_mutation)


elongated_mutations = []
for short_mutation in short_best_mutations:
    checked_motive = short_mutation[:4]
    short_mot_len = len(short_mutation)
    for full_best_mutation in full_best_mutations:
        full_best_mut_len = len(full_best_mutation)
        if full_best_mutation[full_best_mut_len-short_mot_len:full_best_mut_len-short_mot_len+4] == checked_motive:
            elongated_mutations.append(full_best_mutation[:full_best_mut_len-short_mot_len] + short_mutation)
            
full_best_mutations += elongated_mutations

best_expr = 0        
for best_mut in full_best_mutations:
    mut_to_pred = best_mut
    if len(mut_to_pred) == 40:
        
        mut_17_24 = mut_to_pred[:8]
        mut_21_28 = mut_to_pred[4:12]
        mut_25_32 = mut_to_pred[8:16]
        mut_29_36 = mut_to_pred[12:20]
        mut_33_40 = mut_to_pred[16:24]
        mut_37_44 = mut_to_pred[20:28]
        mut_41_48 = mut_to_pred[24:32]
        mut_45_52 = mut_to_pred[28:36]
        
        expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
        expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
        expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
        expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
        expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
        expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
        expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
        expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
        
        pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
        if pred_expr >= best_expr:
            best_expr = pred_expr
            print(pred_expr, mut_to_pred)


checked_low_mutations = {'33-40': [], '37-44': [], '41-48': [], '45-52': []}
worst_mutations = []
short_worst_mutations = []

for mutation in ascen_pred_29_36.Mutation[0:100]:
    worst_mutation = mutation
    mutation_end = mutation[4:]
    if mutation_end not in checked_low_mutations['33-40']:
        checked_low_mutations['33-40'].append(mutation_end) 
        list_33_40 = find_max(mutation_end, ascen_pred_33_40)
        for mut_33_40 in list_33_40:
            mut_37_44_start = mut_33_40[4:]
            worst_mutation += mut_37_44_start
            if mut_37_44_start not in checked_low_mutations['37-44']:
                checked_low_mutations['37-44'].append(mut_37_44_start)
                list_37_44 = find_max(mut_37_44_start, ascen_pred_37_44)
                for mut_37_44 in list_37_44:
                    mut_41_48_start = mut_37_44[4:]
                    worst_mutation += mut_41_48_start
                    if mut_41_48_start not in checked_low_mutations['41-48']:
                        checked_low_mutations['41-48'].append(mut_41_48_start)
                        list_41_48 = find_max(mut_41_48_start,ascen_pred_41_48)
                        for mut_41_48 in list_41_48:
                            mut_45_52_start = mut_41_48[4:]
                            worst_mutation += mut_45_52_start
                            if mut_45_52_start not in checked_low_mutations['45-52']:
                                checked_low_mutations['45-52'].append(mut_45_52_start)
                                list_45_52 = find_max(mut_45_52_start, ascen_pred_45_52)
                                for mut_45_52 in list_45_52:
                                    worst_mutation += mut_45_52[4:]
                                    worst_mutations.append(worst_mutation)
                                    worst_mutation = worst_mutation[:-4]
                                else:
                                    worst_mutation = worst_mutation[:-4]
                            else:
                                short_worst_mutations.append(worst_mutation)
                                worst_mutation = worst_mutation[:-4]
                        else:
                            worst_mutation = worst_mutation[:-4]
                    else:
                        short_worst_mutations.append(worst_mutation)
                        worst_mutation = worst_mutation[:-4]
                else:
                    worst_mutation = worst_mutation[:-4]
            else:
                short_worst_mutations.append(worst_mutation)
                worst_mutation = worst_mutation[:-4]
        else:
            worst_mutation = worst_mutation[:-4]
            
            
elongated_mutations = []
for short_mutation in short_worst_mutations:
    checked_motive = short_mutation[-4:]
    short_mot_len = len(short_mutation)
    for worst_mutation in worst_mutations:
        if worst_mutation[short_mot_len - 4: short_mot_len] == checked_motive:
            elongated_mutations.append(short_mutation + worst_mutation[short_mot_len:])
            
worst_mutations += elongated_mutations


checked_worst_mutations = {'17-24': [], '21-28': [], '25-32': []}
full_worst_mutations = []
short_worst_mutations = []

for mutation in worst_mutations:
    worst_mutation = mutation
    mutation_start = mutation[:4]
    if mutation_start not in checked_worst_mutations['25-32']:
        checked_worst_mutations['25-32'].append(mutation_start) 
        list_25_32 = find_max(mutation_start, ascen_pred_25_32, 'reverse')
        for mut_25_32 in list_25_32:
            mut_21_28_end = mut_25_32[:4]
            worst_mutation = mut_21_28_end + worst_mutation
            if mut_21_28_end not in checked_worst_mutations['21-28']:
                checked_worst_mutations['21-28'].append(mut_21_28_end)
                list_21_28 = find_max(mut_21_28_end, ascen_pred_21_28, 'reverse')
                for mut_21_28 in list_21_28:
                    mut_17_24_end = mut_21_28[:4]
                    worst_mutation = mut_17_24_end + worst_mutation
                    if mut_17_24_end not in checked_worst_mutations['17-24']:
                        checked_worst_mutations['17-24'].append(mut_17_24_end)
                        list_17_24 = find_max(mut_17_24_end, ascen_pred_17_24, 'reverse')
                        for mut_17_24 in list_17_24:
                            worst_mutation = mut_17_24[:4] + worst_mutation
                            full_worst_mutations.append(worst_mutation)
                            worst_mutation = worst_mutation[4:]
                        else:
                            worst_mutation = worst_mutation[4:]
                    else:
                        short_worst_mutations.append(worst_mutation)
                        worst_mutation = worst_mutation[4:]
                else:
                    worst_mutation = worst_mutation[4:]
            else:
                short_worst_mutations.append(worst_mutation)
                worst_mutation = worst_mutation[4:]
        else:
            worst_mutation = worst_mutation[4:]


elongated_mutations = []
for short_mutation in short_worst_mutations:
    checked_motive = short_mutation[:4]
    short_mot_len = len(short_mutation)
    for full_worst_mutation in full_worst_mutations:
        full_worst_mut_len = len(full_worst_mutation)
        if full_worst_mutation[full_worst_mut_len-short_mot_len:full_worst_mut_len-short_mot_len+4] == checked_motive:
            elongated_mutations.append(full_worst_mutation[:full_worst_mut_len-short_mot_len] + short_mutation)
            
full_worst_mutations += elongated_mutations


worst_expr = 0        
for worst_mut in full_worst_mutations:
    mut_to_pred = worst_mut
    if len(mut_to_pred) == 36:
        
        mut_17_24 = mut_to_pred[:8]
        mut_21_28 = mut_to_pred[4:12]
        mut_25_32 = mut_to_pred[8:16]
        mut_29_36 = mut_to_pred[12:20]
        mut_33_40 = mut_to_pred[16:24]
        mut_37_44 = mut_to_pred[20:28]
        mut_41_48 = mut_to_pred[24:32]
        mut_45_52 = mut_to_pred[28:36]
        
        expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
        expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
        expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
        expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
        expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
        expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
        expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
        expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
        
        pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
        if pred_expr < worst_expr:
            worst_expr = pred_expr
            print(pred_expr, mut_to_pred)


#Prediction startes from 49-56 fragment.
checked_mutations = {'17-24': [], '21-28': [], '25-32': [], '29-36': [], '33-40': [], '37-44': [], '41-48': []}
best_mutations = []
short_best_mutations = []

for mutation in pred_45_52.Mutation[0:100]:
    best_mutation = mutation
    mutation_start = mutation[:4]
    if mutation_start not in checked_mutations['41-48']:
        checked_mutations['41-48'].append(mutation_start)
        list_41_48 = find_max(mutation_start, pred_41_48, 'reverse')
        for mut_41_48 in list_41_48:
            mut_37_44_end = mut_41_48[:4]
            best_mutation = mut_37_44_end + best_mutation
            if mut_37_44_end not in checked_mutations['37-44']:
                checked_mutations['37-44'].append(mut_37_44_end)
                list_37_44 = find_max(mut_37_44_end, pred_37_44, 'reverse')
                for mut_37_44 in list_37_44:
                    mut_33_40_end = mut_37_44[:4]
                    best_mutation = mut_33_40_end + best_mutation
                    if mut_33_40_end not in checked_mutations['33-40']:
                        checked_mutations['33-40'].append(mut_33_40_end)
                        list_33_40 = find_max(mut_33_40_end, pred_33_40, 'reverse')
                        for mut_33_40 in list_33_40:
                            mut_29_36_end = mut_33_40[:4]
                            best_mutation = mut_29_36_end + best_mutation
                            if mut_29_36_end not in checked_mutations['29-36']:
                                checked_mutations['29-36'].append(mut_29_36_end)
                                list_29_36 = find_max(mut_29_36_end, pred_29_36, 'reverse')
                                for mut_29_36 in list_29_36:
                                    mut_25_32_end = mut_29_36[:4]
                                    best_mutation = mut_25_32_end + best_mutation
                                    if mut_25_32_end not in checked_mutations['25-32']:
                                        checked_mutations['25-32'].append(mut_25_32_end)
                                        list_25_32 = find_max(mut_25_32_end, pred_25_32, 'reverse')
                                        for mut_25_32 in list_25_32:
                                            mut_21_28_end = mut_25_32[:4]
                                            best_mutation = mut_21_28_end + best_mutation
                                            if mut_21_28_end not in checked_mutations['21-28']:
                                                checked_mutations['21-28'].append(mut_21_28_end)
                                                list_21_28 = find_max(mut_21_28_end, pred_21_28, 'reverse')
                                                for mut_21_28 in list_21_28:
                                                    mut_17_24_end = mut_21_28[:4]
                                                    best_mutation = mut_17_24_end + best_mutation
                                                    if mut_17_24_end not in checked_mutations['17-24']:
                                                        checked_mutations['17-24'].append(mut_17_24_end)
                                                        list_17_24 = find_max(mut_17_24_end, pred_17_24, 'reverse')
                                                        for mut_17_24 in list_17_24:
                                                            best_mutation = mut_17_24[:4] + best_mutation
                                                            full_best_mutations.append(best_mutation)
                                                            best_mutation = best_mutation[4:]
                                                        else:
                                                            best_mutation = best_mutation[4:]
                                                    else:
                                                        short_best_mutations.append(best_mutation)
                                                        best_mutation = best_mutation[4:]
                                                else:
                                                    best_mutation = best_mutation[4:]
                                            else:
                                                short_best_mutations.append(best_mutation)
                                                best_mutation = best_mutation[4:]
                                        else:
                                            best_mutation = best_mutation[4:]
                                    else:
                                        short_best_mutations.append(best_mutation)
                                        best_mutation = best_mutation[4:]
                                else:
                                    best_mutation = best_mutation[4:]
                            else:
                                short_best_mutations.append(best_mutation)
                                best_mutation = best_mutation[4:]
                        else:
                            best_mutation = best_mutation[4:]
                    else:
                        short_best_mutations.append(best_mutation)
                        best_mutation = best_mutation[4:]
                else:
                    best_mutation = best_mutation[4:]
            else:
                short_best_mutations.append(best_mutation)
                best_mutation = best_mutation[4:]
        else:
            best_mutation = best_mutation[4:]

elongated_mutations = []
for short_mutation in short_best_mutations:
    checked_motive = short_mutation[:4]
    short_mot_len = len(short_mutation)
    for best_mutation in best_mutations:
        best_mut_len = len(best_mutation)
        if best_mutation[best_mut_len-short_mot_len:best_mut_len-short_mot_len+4] == checked_motive:
            elongated_mutations.append(best_mutation[:best_mut_len-short_mot_len] + short_mutation)
            
            
best_mutations += elongated_mutations


best_expr = 0        
for best_mut in best_mutations:
    mut_to_pred = best_mut
    if len(mut_to_pred) == 40:
        
        mut_17_24 = mut_to_pred[:8]
        mut_21_28 = mut_to_pred[4:12]
        mut_25_32 = mut_to_pred[8:16]
        mut_29_36 = mut_to_pred[12:20]
        mut_33_40 = mut_to_pred[16:24]
        mut_37_44 = mut_to_pred[20:28]
        mut_41_48 = mut_to_pred[24:32]
        mut_45_52 = mut_to_pred[28:36]
        
        expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
        expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
        expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
        expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
        expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
        expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
        expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
        expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
        
        pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
        if pred_expr >= best_expr:
            best_expr = pred_expr
            print(pred_expr, mut_to_pred)


checked_low_mutations = {'17-24': [], '21-28': [], '25-32': [], '29-36': [], '33-40': [], '37-44': [], '41-48': []}
worst_mutations = []
short_worst_mutations = []

for mutation in ascen_pred_45_52.Mutation[0:100]:
    worst_mutation = mutation
    mutation_end = mutation[:4]
    if mutation_end not in checked_low_mutations['41-48']:
        checked_low_mutations['41-48'].append(mutation_end)
        list_41_48 = find_max(mutation_end, ascen_pred_41_48, 'reverse')
        for mut_41_48 in list_41_48:
            mut_37_44_end = mut_41_48[:4]
            worst_mutation = mut_37_44_end + worst_mutation
            if mut_37_44_end not in checked_low_mutations['37-44']:
                checked_low_mutations['37-44'].append(mut_37_44_end)
                list_37_44 = find_max(mut_37_44_end, ascen_pred_37_44, 'reverse')
                for mut_37_44 in list_37_44:
                    mut_33_40_end = mut_37_44[:4]
                    worst_mutation = mut_33_40_end + worst_mutation
                    if mut_33_40_end not in checked_low_mutations['33-40']:
                        checked_low_mutations['33-40'].append(mut_33_40_end)
                        list_33_40 = find_max(mut_33_40_end, ascen_pred_33_40, 'reverse')
                        for mut_33_40 in list_33_40:
                            mut_29_36_end = mut_33_40[:4]
                            worst_mutation = mut_29_36_end + worst_mutation
                            if mut_29_36_end not in checked_low_mutations['29-36']:
                                checked_low_mutations['29-36'].append(mut_29_36_end)
                                list_29_36 = find_max(mut_29_36_end, ascen_pred_29_36, 'reverse')
                                for mut_29_36 in list_29_36:
                                    mut_25_32_end = mut_29_36[:4]
                                    worst_mutation = mut_25_32_end + worst_mutation
                                    if mut_25_32_end not in checked_low_mutations['25-32']:
                                        checked_low_mutations['25-32'].append(mut_25_32_end)
                                        list_25_32 = find_max(mut_25_32_end, ascen_pred_25_32, 'reverse')
                                        for mut_25_32 in list_25_32:
                                            mut_21_28_end = mut_25_32[:4]
                                            worst_mutation = mut_21_28_end + worst_mutation
                                            if mut_21_28_end not in checked_low_mutations['21-28']:
                                                checked_low_mutations['21-28'].append(mut_21_28_end)
                                                list_21_28 = find_max(mut_21_28_end, ascen_pred_21_28, 'reverse')
                                                for mut_21_28 in list_21_28:
                                                    mut_17_24_end = mut_21_28[:4]
                                                    worst_mutation = mut_17_24_end + worst_mutation
                                                    if mut_17_24_end not in checked_low_mutations['17-24']:
                                                        checked_low_mutations['17-24'].append(mut_17_24_end)
                                                        list_17_248 = find_max(mut_17_24_end, ascen_pred_17_24, 'reverse')
                                                        for mut_17_24 in list_17_24:
                                                            worst_mutation = mut_17_24[:4] + worst_mutation
                                                            worst_mutations.append(worst_mutation)
                                                            worst_mutation = worst_mutation[4:]
                                                        else:
                                                            worst_mutation = worst_mutation[4:]
                                                    else:
                                                        short_worst_mutations.append(worst_mutation)
                                                        worst_mutation = worst_mutation[4:]
                                                else:
                                                    worst_mutation = worst_mutation[4:]
                                            else:
                                                short_worst_mutations.append(worst_mutation)
                                                worst_mutation = worst_mutation[4:]
                                        else:
                                            worst_mutation = worst_mutation[4:]
                                    else:
                                        short_worst_mutations.append(worst_mutation)
                                        worst_mutation = worst_mutation[4:]
                                else:
                                    worst_mutation = worst_mutation[4:]
                            else:
                                short_worst_mutations.append(worst_mutation)
                                worst_mutation = worst_mutation[4:]
                        else:
                            worst_mutation = worst_mutation[4:]
                    else:
                        short_worst_mutations.append(worst_mutation)
                        worst_mutation = worst_mutation[4:]
                else:
                    worst_mutation = worst_mutation[4:]
            else:
                short_worst_mutations.append(worst_mutation)
                worst_mutation = worst_mutation[4:]
        else:
            worst_mutation = worst_mutation[4:]
            
            
elongated_mutations = []
for short_mutation in short_worst_mutations:
    checked_motive = short_mutation[:4]
    short_mot_len = len(short_mutation)
    for worst_mutation in worst_mutations:
        worst_mut_len = len(worst_mutation)
        if worst_mutation[worst_mut_len-short_mot_len:worst_mut_len-short_mot_len+4] == checked_motive:
            elongated_mutations.append(worst_mutation[:worst_mut_len-short_mot_len] + short_mutation)
            
worst_mutations += elongated_mutations


worst_expr = 0        
for worst_mut in worst_mutations:
    mut_to_pred = worst_mut
    if len(mut_to_pred) == 40:
        
        mut_17_24 = mut_to_pred[:8]
        mut_21_28 = mut_to_pred[4:12]
        mut_25_32 = mut_to_pred[8:16]
        mut_29_36 = mut_to_pred[12:20]
        mut_33_40 = mut_to_pred[16:24]
        mut_37_44 = mut_to_pred[20:28]
        mut_41_48 = mut_to_pred[24:32]
        mut_45_52 = mut_to_pred[28:36]
        
        expr_17_24 = pred_17_24.prediction_label[pred_17_24.Mutation == mut_17_24].tolist()
        expr_21_28 = pred_21_28.prediction_label[pred_21_28.Mutation == mut_21_28].tolist()
        expr_25_32 = pred_25_32.prediction_label[pred_25_32.Mutation == mut_25_32].tolist()
        expr_29_36 = pred_29_36.prediction_label[pred_29_36.Mutation == mut_29_36].tolist()
        expr_33_40 = pred_33_40.prediction_label[pred_33_40.Mutation == mut_33_40].tolist()
        expr_37_44 = pred_37_44.prediction_label[pred_37_44.Mutation == mut_37_44].tolist()
        expr_41_48 = pred_41_48.prediction_label[pred_41_48.Mutation == mut_41_48].tolist()
        expr_45_52 = pred_45_52.prediction_label[pred_45_52.Mutation == mut_45_52].tolist()
        
        pred_expr = (expr_17_24[0]*r2_17_24 + expr_21_28[0]*r2_21_28 + expr_25_32[0]*r2_25_32 + expr_29_36[0]*r2_29_36 + expr_33_40[0]*r2_33_40 + expr_37_44[0]*r2_37_44 + expr_41_48[0]*r2_41_48 + expr_45_52[0]*r2_45_52) / r2_sum
        if pred_expr < worst_expr:
            worst_expr = pred_expr
            print(pred_expr, mut_to_pred)
