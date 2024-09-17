import pandas as pd
import scipy.stats as sps
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sbn
import itertools as  itr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import sys

if __name__ == "__main__":
    if len(sys.argv) == 4:
        #args
        prj = sys.argv[1]
        output = sys.argv[2]
        delta = sys.argv[3]
        #data preprocessing
        data = pd.read_table(f"{prj}")
        def find_indices(words):
            indices = []
            for index, word in enumerate(words):
                if 'IJC' in word or 'SJC' in word:
                    indices.append(index)
            return indices
        indices = find_indices(data.columns)
        lst = []
        lst_ = []
        for i in (indices):
            lst.append(np.array(data[data.columns[i]]))
    
        for i in (indices):
            lst_.append(data.columns[i])
    

        IJC_SJC = [[]] * len(data.index)

        for i in data.index:
            for j in range(len(lst)):
                if type(lst[j][i]) == float:
                    print(j, i, lst[j][i])
                lst[j][i] = lst[j][i].split(',')
                lst[j][i] = list(map(int, lst[j][i]))
                if j % 2:
                    IJC_SJC[i] = lst[j][i] + IJC_SJC[i]
                else:
                    IJC_SJC[i] = IJC_SJC[i] + lst[j][i]
        
        #data filthering
        for i in range(len(lst)):
            data[lst_[i]] = lst[i]

        data['IJC_SJC'] = IJC_SJC

        lst = []
        ln = len(data['IJC_SJC'][0])
        for i in data.index:
            lst.append(np.abs(np.log(sum(data['IJC_SJC'][i][:ln // 2]) /
                    sum(data['IJC_SJC'][i][ln // 2:]))))
    
        data['flag'] = lst
        data = data.loc[data['flag'] < 5]

        lst = []
        for i in data.index:
            lst.append(min(list(map(sum, [data[lst_[0]][i], data[lst_[1]][i], data[lst_[2]][i], data[lst_[3]][i]]))))

        data['flag2'] = lst
        data = data.loc[data['flag2'] >= 3]

        #spearmanr main funcion
        def corr_pirs_gene(data, GeneID = 'Ablim2', delta = 0.0001, f = False):
            s = 0
            corr = []
            df = (data.loc[data['GeneID'] == GeneID])
            if len(df) < 2:
                return []
            for i in range(1, len(df.index)):
                for j in range(i):
                    k = min(len(df["IJC_SJC"][df.index[i]]), len(df["IJC_SJC"][df.index[j]]))
                    VAL = np.vstack((df["IJC_SJC"][df.index[i]][:k], df["IJC_SJC"][df.index[j]][:k]))
                    l = 0
                    psi1, psi2 = sum(VAL[0][k // 2:]) / sum(VAL[0]), sum(VAL[1][k // 2:]) / sum(VAL[1])
                    VAL = np.vstack((VAL[0][k // 2:], VAL[1][k // 2:], VAL[0][:k // 2], VAL[1][:k // 2]))
                    #en = -((psi1 + psi2)/2) * np.log2((psi1 + psi2)/2) * (2-abs(np.corrcoef(VAL[0], VAL[1])[0][1]))
                    pvalue = sps.spearmanr(VAL[0], VAL[1]).pvalue
                    corr = (sps.spearmanr(VAL[0], VAL[1]).statistic * np.sum((VAL[0], VAL[1])) + 
                            sps.spearmanr(VAL[2], VAL[3]).statistic * np.sum((VAL[2], VAL[3]))) / np.sum(VAL)
                    if pvalue < delta:
                        s += 1
                        if f:
                            with open(f, 'a') as out:
                                print(GeneID, df['GeneSymbol'][df.index[i]],
                                      df['ID'][df.index[i]], df['ID'][df.index[j]], 
                                      pvalue,
                                      corr,
                                      #df['flag'][df.index[i]], df['flag'][df.index[j]],
                                      df['exonStart_0base'][df.index[i]], df['exonStart_0base'][df.index[j]], 
                                      df['exonEnd'][df.index[i]], df['exonEnd'][df.index[j]],
                                      #df['num'][df.index[i]], df['num'][df.index[j]], 
                                      psi1,
                                      psi2,
                                      #en,
                                      #sum(VAL[0][k // 2:]) * sum(VAL[1]) / (sum(VAL[1][k // 2:]) * sum(VAL[0]))
                                      sep = '\t', file = out)
                        else:
                            print(GeneID, df['geneSymbol'][df.index[i]],
                                  df['ID'][df.index[i]], df['ID'][df.index[j]], 
                                  pvalue,
                                  corr,
                                  #df['flag'][df.index[i]], df['flag'][df.index[j]],
                                  df['exonStart_0base'][df.index[i]], df['exonStart_0base'][df.index[j]], 
                                  df['exonEnd'][df.index[i]], df['exonEnd'][df.index[j]],
                                  #df['num'][df.index[i]], df['num'][df.index[j]], 
                                  psi1,
                                  psi2,
                                  sps.spearmanr(VAL[0], VAL[1]).statistic,
                                  #en,
                                  #sum(VAL[0][k // 2:]) * sum(VAL[1]) / (sum(VAL[1][k // 2:]) * sum(VAL[0])),  
                                  sep = '\t')
            return corr
        sum_, corr, f = 0, [], output
        with open(f, 'w') as out:
                print('GeneID', 'GeneSymbol', 'id1', 'id2', 'pvalue', 'spearmanr', 'exonstart1', 'exonstart2', "exonend1", "exonend2", "psi1", "psi2", sep = '	', file = out)
        print('GeneID', 'GeneSymbol', 'id1', 'id2', 'pvalue', 'spearmanr', 'log1', 'log2', 'exonstart1', 'exonstart2', "exonend1", "exonend2", "psi1", "psi2", sep = '	')
        for i in sorted(list(set(data['GeneID'].tolist()))):
            corr += corr_pirs_gene(data, delta = delta, GeneID = i, f = f)
        
        mxe = pd.read_table(f, sep = "	")
        #filthered correlated data
        p = np.zeros(len(mxe), dtype = int)+1
        for i in mxe.index:
            st_end = [mxe["exonstart1"][i], mxe["exonstart2"][i], mxe["exonend1"][i], mxe["exonend2"][i]]
            for j in range(4):
                if st_end[(j-1)%4] <= st_end[(j)%4] <= st_end[(j+1)%4]:
                    p[i] = 0
                    break
        mxe = mxe.loc[p == 1]
        mxe = mxe.drop_duplicates(subset=['exonstart1', 'exonstart2', 'exonend1', 'exonend2'], keep=False)
        mxe.to_excel(f"{output[:-4]}_spearmanr_filthered.xlsx", index = False)
    else:
        print(f"Expected 3 arguments, give {len(sys.argv)-1}")