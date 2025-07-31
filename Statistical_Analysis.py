# Title: This script include the functions to run the statistical analysis, especially the Tukey's HSD post-hoc. 
#
# Author: Farhad Abdollahi (farhad.abdollahi.ctr@dot.gov)
# Date: 05/19/2025
# ======================================================================================================================

# Importing the required libraries.
import os
import numpy as np
import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd, MultiComparison


def Tukey_Grouping(data):
    """
    Define a function to perform the Tukey's HSD (Honestly Significant Difference) post-hoc test for ANOVA comparisons.

    :param data: input data as a DataFrame with two columns: 1st column called "value" and include all data values, and 
        second column called "group" and include the correspnding group name. 
    :return NewDF, df: two variable, "NewDF" which is the pairwise group comparisons and the second one "df" which 
        include the grouping based on Tukey's HSD post-hpc test. 
    """
    # Perform the analysis. 
    mc = MultiComparison(data['value'], data['group'])
    tukey_result = mc.tukeyhsd()
    # print(tukey_result.summary())
    # ----------------------------------------------------
    # Covert the results into a dataframe. 
    NewDF = {'Group1': [], 'Group2': [], 'MeanDiff': [], 'P-adj': [], 'Reject': []}
    GroupsUnique = tukey_result.groupsunique
    # def natural_key(s):
    #     return [int(text) if text.isdigit() else text for text in re.split('(\d+)', s)]
    # try:
    #     GroupsUnique = sorted(GroupsUnique, key=natural_key)
    # except: pass
    MeanDiff     = tukey_result.meandiffs
    Pvalues      = tukey_result.pvalues
    Rejects      = tukey_result.reject
    counter = 0
    for i, g1 in enumerate(GroupsUnique):
        for g2 in GroupsUnique[i+1:]:
            NewDF['Group1'].append(g1)
            NewDF['Group2'].append(g2)
            NewDF['MeanDiff'].append(MeanDiff[counter])
            NewDF['P-adj'].append(Pvalues[counter])
            NewDF['Reject'].append(Rejects[counter])
            counter += 1
    NewDF = pd.DataFrame(NewDF)
    # NewDF
    # ----------------------------------------------------
    # Perform the grouping. 
    LetterCodes = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    counter = 0
    Classes = [[] for i in range(len(GroupsUnique))]
    for i in range(len(GroupsUnique) - 1):
        letter = LetterCodes[counter]
        Used = False
        for j in range(i+1, len(GroupsUnique)):
            g1 = GroupsUnique[i]
            g2 = GroupsUnique[j]
            df = NewDF[(NewDF['Group1'] == g1) & (NewDF['Group2'] == g2)]
            if not df.iloc[0]['Reject']:        # If it is False, we need to have a letter for it. 
                # Check if they already have similar letter. 
                IsAvailable = False
                for letter1 in Classes[i]:
                    if letter1 in Classes[j]:
                        IsAvailable = True
                        break
                # Otherwise, add new letter. 
                if not IsAvailable:
                    if letter not in Classes[i]:
                        Classes[i].append(letter)
                    Classes[j].append(letter)
                    Used = True
        if Used:
            counter += 1
        else:
            if len(Classes[i]) == 0:
                Classes[i].append(letter)
                counter += 1
    if len(Classes[-1]) == 0:
        Classes[-1].append(LetterCodes[counter])
    # ----------------------------------------------------
    # Add the classes to the NewDF. 
    df = {'Group': [], 'Grouping': []}
    for i in range(len(GroupsUnique)):
        df['Group'].append(GroupsUnique[i])
        df['Grouping'].append(','.join(Classes[i]))
    df = pd.DataFrame(df)
    # ----------------------------------------------------
    # Return the results. 
    return NewDF, df
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


if __name__ == '__main__':
    pass
