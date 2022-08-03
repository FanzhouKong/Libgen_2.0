import pandas as pd
from tqdm import tqdm
import itertools
import numpy as np
def readin_MSDIAL(file):
    df = pd.read_csv(file, sep = '\t', header=[4])
    return(df)
def precursor_matching(df, sample_list, mzmmu, adducts, comments, CE = None,threshold = 0.3, ifppm = False, iffill = True):
    # print("i am updated!")
    # mz_mmu = mzmmu
    data_raw = pd.DataFrame(columns = ['NAME','key','PRECURSORMZ','InChIKey','Formula',
                       'ExactMass','Adduct','Spectrum_type','RETENTIONTIME','Average_mz',
                       'Comment', 'Alignment_ID','ms1','msms','Collision_energy','intensity','mix_label'])
    if iffill == True:
        ms2df = df.loc[(df['MS/MS assigned']) & (df['Fill %'] < threshold)]
    else:
        ms2df = df.loc[(df['MS/MS assigned'])]
    # parameters
    # this is already matching the spectrum in the same data file
    spectrum_type = 'MS2'
    # adducts = ['[M+H]+', '[M+NH4]+', '[M+Na]+']
    # comments = ' EAD spectrum'
    for mix in tqdm(ms2df.columns[32:-2]):
        for index, df_row in ms2df[ms2df['Spectrum reference file name'] == mix].iterrows():
            msms = df_row['MS/MS spectrum'].replace(' ', '\n').replace(':', '\t')
            ms1 = df_row['MS1 isotopic spectrum'].replace(' ', '\n').replace(':', '\t')
            nlines = msms.count('\n')+1
            if ifppm:
                if df_row['Average Mz']<400:
                    mz_lower= df_row['Average Mz'] - 0.004
                    mz_upper = df_row['Average Mz'] + 0.004
                else:
                    mz_lower = df_row['Average Mz'] - df_row['Average Mz']*mzmmu/(1E6)
                    mz_upper = df_row['Average Mz'] + df_row['Average Mz']*mzmmu/(1E6)
            else:
                mz_lower= df_row['Average Mz'] - mzmmu
                mz_upper = df_row['Average Mz'] + mzmmu

            for index, ls_row in sample_list[sample_list['Mix label'] == mix].iterrows():
                for adduct in adducts:
                    if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
                        if adduct in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
                            ion_mode = 'P'
                        else:
                            ion_mode = 'N'
                    if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
                        data_raw = data_raw.append({
                                                        'NAME':ls_row['Name'],
                                                       'key':str(ls_row['InChIKey'])+adduct,
                                                        'Ion_mode':ion_mode,
                                                       'PRECURSORMZ':(float(ls_row[adduct])),
                                                       'InChIKey':ls_row['InChIKey'],
                                                        'Formula':ls_row['Formula'],
                                                        'ExactMass':(float(ls_row['Exact Mass'])),
                                                        'Adduct':adduct,
                                                        'Collision_energy':CE,
                                                        'Spectrum_type':spectrum_type,
                                                        'RETENTIONTIME':(float(df_row['Average Rt(min)'])),
                                                        'Average_mz':(float(df_row['Average Mz'])),
                                                        'Comment':str(df_row['Alignment ID']) + comments + ' intensity ' + str(df_row[mix]),
                                                        'Alignment_ID':(df_row['Alignment ID']),
                                                        # 'Num_Peaks':int(nlines),
                                                        'ms1':ms1,
                                                        'msms':msms,
                                                        'intensity':(df_row[mix]),
                                                        'mix_label':mix,
                                                       # 'count':str('na')
                            }, ignore_index=True)
    # occ = data_raw.key.value_counts(normalize=False, sort = False)


    # data_raw_count = pd.DataFrame()
    # for i in range(len(occ)):
    #     temp_df = data_raw.loc[data_raw.key == occ.index[i],:]
    #     temp_df['count'] = temp_df['count'].replace(['na'], int(occ[i]))
    #     data_raw_count=pd.concat([data_raw_count, temp_df], ignore_index = True, axis = 0)
    # data_raw_count['row_num'] = np.arange(len(data_raw_count))
    data_raw.round({'PRECURSORMZ':6})
    # data_raw['PRECURSORMZ']=pd.to_numeric(data_raw['PRECURSORMZ'])
    # data_raw['ExactMass']=pd.to_numeric(data_raw['ExactMass'])
    # data_raw['Average_mz']=pd.to_numeric(data_raw['Average_mz'])
    # data_raw['intensity']=pd.to_numeric(data_raw['intensity'])
    # data_raw['Num_Peaks']=pd.to_numeric(data_raw['Num_Peaks'])
    # data_raw['RETENTIONTIME']=pd.to_numeric(data_raw['RETENTIONTIME'])
    data_raw.reset_index(inplace=True, drop=True)
    return(data_raw)

def precursor_matching_rt(df, sample_list, mzmmu, adducts, comments, CE = None,threshold = 0.3, ppm = False):
    # print("i am updated!")
    # mz_mmu = mzmmu
    data_raw = pd.DataFrame(columns = ['NAME','key','PRECURSORMZ','InChIKey','Formula',
                       'ExactMass','Adduct','Spectrum_type','RETENTIONTIME','Average_mz',
                       'Comment', 'Alignment_ID','ms1','msms','Collision_energy','intensity','mix_label'])
    ms2df = df.loc[(df['MS/MS assigned']) & (df['Fill %'] < threshold)]
    # parameters
    # this is already matching the spectrum in the same data file
    spectrum_type = 'MS2'
    # adducts = ['[M+H]+', '[M+NH4]+', '[M+Na]+']
    # comments = ' EAD spectrum'
    for mix in tqdm(ms2df.columns[32:-2]):
        for index, df_row in ms2df[ms2df['Spectrum reference file name'] == mix].iterrows():
            msms = df_row['MS/MS spectrum'].replace(' ', '\n').replace(':', '\t')
            ms1 = df_row['MS1 isotopic spectrum'].replace(' ', '\n').replace(':', '\t')
            nlines = msms.count('\n')+1
            if ppm:
                if df_row['Average Mz']<400:
                    mz_lower= df_row['Average Mz'] - 0.004
                    mz_upper = df_row['Average Mz'] + 0.004
                else:
                    mz_lower = df_row['Average Mz'] - df_row['Average Mz']*mzmmu/(1E6)
                    mz_upper = df_row['Average Mz'] + df_row['Average Mz']*mzmmu/(1E6)
            else:
                mz_lower= df_row['Average Mz'] - mzmmu
                mz_upper = df_row['Average Mz'] + mzmmu
            rt = df_row['Average Rt(min)']
            for index, ls_row in sample_list[sample_list['Mix label'] == mix].iterrows():
                for adduct in adducts:
                    if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
                        if adduct in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
                            ion_mode = 'P'
                        else:
                            ion_mode = 'N'
                    if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
                        if abs(ls_row['RT [min]']-rt)<0.5:
                            data_raw = data_raw.append({
                                                        'NAME':ls_row['Name'],
                                                       'key':ls_row['InChIKey']+adduct,
                                                        'Ion_mode':ion_mode,
                                                       'PRECURSORMZ':(float(ls_row[adduct])),
                                                       'InChIKey':ls_row['InChIKey'],
                                                        'Formula':ls_row['Formula'],
                                                        'ExactMass':(float(ls_row['Exact Mass'])),
                                                        'Adduct':adduct,
                                                        'Collision_energy':CE,
                                                        'Spectrum_type':spectrum_type,
                                                        'RETENTIONTIME':(float(df_row['Average Rt(min)'])),
                                                        'Average_mz':(float(df_row['Average Mz'])),
                                                        'Comment':str(df_row['Alignment ID']) + comments + ' intensity ' + str(df_row[mix]),
                                                        'Alignment_ID':(df_row['Alignment ID']),
                                                        # 'Num_Peaks':int(nlines),
                                                        'ms1':ms1,
                                                        'msms':msms,
                                                        'intensity':(df_row[mix]),
                                                        'mix_label':mix,
                                                       # 'count':str('na')
                            }, ignore_index=True)

    # occ = data_raw.key.value_counts(normalize=False, sort = False)


    # data_raw_count = pd.DataFrame()
    # for i in range(len(occ)):
    #     temp_df = data_raw.loc[data_raw.key == occ.index[i],:]
    #     temp_df['count'] = temp_df['count'].replace(['na'], int(occ[i]))
    #     data_raw_count=pd.concat([data_raw_count, temp_df], ignore_index = True, axis = 0)
    # data_raw_count['row_num'] = np.arange(len(data_raw_count))
    data_raw.round({'PRECURSORMZ':6})
    # data_raw['PRECURSORMZ']=pd.to_numeric(data_raw['PRECURSORMZ'])
    # data_raw['ExactMass']=pd.to_numeric(data_raw['ExactMass'])
    # data_raw['Average_mz']=pd.to_numeric(data_raw['Average_mz'])
    # data_raw['intensity']=pd.to_numeric(data_raw['intensity'])
    # data_raw['Num_Peaks']=pd.to_numeric(data_raw['Num_Peaks'])
    # data_raw['RETENTIONTIME']=pd.to_numeric(data_raw['RETENTIONTIME'])
    return(data_raw)
print("I am precursor matching!")