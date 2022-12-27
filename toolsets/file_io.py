import pandas as pd 
from toolsets.spectra_operations import clean_spectrum, break_spectra, pack_spectra, convert_nist_to_string, normalize_spectrum
from toolsets.search import num_search, string_search
import re
import os
from toolsets.msp_file import read_one_spectrum
from toolsets.filename import smart_io
import re
from tqdm import tqdm
def parse_file_name(filepath):
    import ntpath
    return(ntpath.basename(filepath))
def readin_MSDIAL(file):
    df = pd.read_csv(file,
                     sep = '\t',
                     header=[4]
                     )
    return(df)
def read_msp_files(msp_path, clean = True):
    msp_file = smart_io(msp_path, mode = 'r')
    specs = []
    for spec in read_one_spectrum(msp_file, include_raw=False):
        specs.append(spec)
    msp = pd.DataFrame.from_dict(specs)
    msp_with_spectrum = msp[msp['spectrum'].map(lambda d: len(d)) > 0]
    if clean == False:
        return(msp_with_spectrum)
    
    msp_with_spectrum['PRECURSORMZ']=pd.to_numeric(msp_with_spectrum['PRECURSORMZ'])
    msp_with_spectrum['RETENTIONTIME']=pd.to_numeric(msp_with_spectrum['RETENTIONTIME'])
    peaks = []
    print('cleaning incoming spectrum')
    for index, row in tqdm(msp_with_spectrum.iterrows(), total = len(msp_with_spectrum)):
        peaks.append(clean_spectrum(convert_nist_to_string(row['spectrum']), max_mz = row['PRECURSORMZ'],
            tolerance = 0.05, ifppm = False, noise_level = 0.00))
    msp_with_spectrum['peaks']=peaks
    msp_with_spectrum = msp_with_spectrum[msp_with_spectrum['peaks'].map(lambda d: len(d)) > 0]
    return (msp_with_spectrum)

def find_files(base, pattern):
    '''Return list of files matching pattern in base folder.'''
    return [filename for filename in os.listdir(base) if re.search(pattern, filename, re.IGNORECASE)]
def check_missing_files(file_name_list, tail, source_dir):
    duplicate_files = []
    missing_files = []
    good_files = []
    for mix in (file_name_list):
        found_file = find_files(source_dir, f"{str(mix)+tail}")
        if len(found_file)==0:
            missing_files.append(mix)
        elif len(found_file)==1:
            good_files.append(mix)
        elif len(found_file)>1:
            duplicate_files.append(mix)
    return(good_files, duplicate_files, missing_files)
def read_in_alphapept(path_to_features, peak_purity = 0.001):
    features = pd.read_csv(path_to_features)
    for index, row in features.iterrows():
        try:
            features.at[index, "peaks"]=parse_feature_peaks(row['peaks'])
        except ValueError:
            print("the scannumber is ", row['scan_number'])
            print("the problematic feature is ", path_to_features.split("/")[-1])
            # break
    features.drop_duplicates(subset = ['query_idx'], keep = "first", inplace = True, ignore_index = True)
    features = num_search(features, 'ms1_precursor_intensity', 0, ">", inclusion = False)
    features = num_search(features, 'ms1_intensity_ratio', peak_purity, ">", inclusion = False)
    vc = features['charge'].value_counts().rename_axis('charge').reset_index(name='counts')
    vc.sort_values(by=['counts'], inplace=True, ascending=False)
    if vc.iloc[0]['charge']>0:
        features = string_search(features, 'charge', 1)
    elif vc.iloc[0]['charge']<0:
        features = string_search(features, 'charge', -1)
    # features_plus = string_search(features, 'charge', 1)
    # features_minus = string_search(features, 'charge', -1)
    # features['peaks'] = features.apply(parse_feature_peaks, axis = 1)
    return features

def parse_feature_peaks(peaks):
    peaks_split = peaks.replace('[','').replace(']','').split('\n')
    mass = []
    intensity = []
    for peak in peaks_split:
        peak = peak.strip()
        mass_temp = float(peak.split(" ")[0])
        intensity_temp = float(peak.split(" ")[-1])
        mass.append(mass_temp)
        intensity.append(intensity_temp)

    return(pack_spectra(mass, intensity))





from tqdm import tqdm
from toolsets.spectra_operations import num_peaks
def export_library_msp(data,output_location, typeofmsms='peaks', ifcollision_energy = False):
    entry = ''
    for index, row in tqdm(data.iterrows(), total = len(data)):
        entry = entry + 'Name: ' + row['reference_name'] + '\n'
        # entry = entry + 'RETENTIONTIME: ' + row['retention_time_wa'] + '\n'
        entry = entry +'Spectrum_type: '+'MS2'+ '\n'
        entry = entry + 'PrecursorMZ: ' + str(row['reference_precursor_mz']) + '\n'
        entry = entry + 'InChIKey: ' + str(row['reference_inchikey']) + '\n'
        entry = entry + 'Formula: ' + row['reference_formula'] + '\n'
        entry = entry + 'ExactMass: ' + str(row['reference_mono_mass']) + '\n'
        entry = entry + 'Precursor_type: ' + row['reference_adduct'] + '\n'
        if ifcollision_energy:
            entry = entry + 'Collision_enerty: ' + str(row['Collision_energy']) + '\n'
        entry = entry + 'RETENTIONTIME: ' + str(row['retention_time_wa']) + '\n'
        if row['reference_adduct'][-1]=='+':
            # charge = '1+'
            ionmode = 'P'
        else:
            ionmode = 'N'
        entry = entry+'Ion_mode: '+ionmode+ '\n'
        entry = (entry + 'Comment: ' + 'method_'+str(row['comment'])+'ms1intensity'+'_'+str(row['ms1_precursor_intensity'])+"_"
                +'ei_'+str(row['ei']) + '\n')
        entry = entry + 'Spectrum_entropy: ' +str((row['spectrum_entropy'])) + '\n'
        entry = entry + 'Normalized_entropy: ' + str((row['normalized_entropy'])) + '\n'
        entry = entry + 'Num peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        entry = entry + row[typeofmsms]
        # entry = entry +str(row['count'])
        entry = entry + '\n'
        entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')
     
    #write string to file
    text_file.write(entry)
     
    #close file
    text_file.close()

def export_ms_sirius(row, output):
    # if row['Adduct'][-1]=='+':
    #         charge = '1+'
    #     else:
    #         charge = '1-'
    mass_1, intensity_1 = break_spectra(row['ms1'])
    pep_mass =de.find_parention(mass_1,intensity_1, row['PRECURSORMZ'])
    entry = ''
    entry = entry + '>compound '+str(row['NAME'])+'\n'
    entry = entry + '>parentmass '+str((pep_mass))+'\n'
    entry = entry + '>ionization '+str((row['Adduct']))+'\n'
    entry = entry +'\n'
    entry = entry + '>collision 35' + '\n'
    entry = entry+(row['msms'])+'\n'
    entry = entry +'\n'
    entry = entry +'\n'
    entry = entry + '>ms1peaks' + '\n'
    entry = entry+(row['ms1'])+'\n'
    entry = entry +'\n'
    text_file = open(output+'.ms', "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()



def export_mgf_sirius(inputfile, output):

    entry = ''
    for index, row in inputfile.iterrows():

        mass_1, intensity_1 = break_spectra(row['ms1'])
        if row['Adduct'][-1]=='+':
            charge = '1+'
        else:
            charge = '1-'
        pep_mass =de.find_parention(mass_1,intensity_1, row['PRECURSORMZ'])
        # output = os.path.join(output_dir, row['NAME']+'.mgf')
        # ms1
        entry = entry + 'BEGIN IONS'+'\n'
        entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
        entry = entry + 'MSLEVEL=1'+'\n'
        entry = entry+'CHARGE=' + charge +'\n'
#         entry = entry+'Adduct=' +str(row['adduct']) +'\n'
        entry = entry+(row['ms1'])+'\n'
        entry = entry + 'END IONS'+'\n'
        entry = entry +'\n'
    #     ms2
        entry = entry + 'BEGIN IONS'+'\n'
        entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
        entry = entry + 'MSLEVEL=2'+'\n'
        entry = entry+'CHARGE=' + charge +'\n'
        entry = entry+(row['msms'])+'\n'
        entry = entry + 'END IONS'+'\n'
    text_file = open(output+'.mgf', "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()
        # break

def export_mat(data,output_location, typeofms1='ms1',typeofmsms = 'msms', ifcollision_energy = True):
    entry = ''
    for index, row in data.iterrows():
        entry = entry + 'NAME: ' + row['NAME'] + '\n'
        entry = entry + 'PrecursorMZ: ' + str(row['PRECURSORMZ']) + '\n'
        entry = entry + 'PRECURSORTYPE: ' + row['Adduct'] + '\n'
        entry = entry + 'INSTRUMENTTYPE: ' +'\n' #empty line
        entry = entry + 'INSTRUMENT: ' +'\n'
        entry = entry + 'Authors: ' +'Arpana, Parker and Fanzhou'+'\n'
        entry = entry + 'License: '+'\n'
        entry = entry + 'FORMULA: ' + str(row['Formula']) + '\n'
        entry = entry + 'ONTOLOGY: ' +'\n'
        entry = entry + 'SMILES: ' +'\n'
        entry = entry + 'INCHIKEY: ' + row['InChIKey'] + '\n'
        entry = entry + 'INCHI: ' +'\n'
        entry = entry+'IONMODE: '+row['Ion_mode']+ '\n'
        if ifcollision_energy:
            entry = entry + 'Collision_enerty: ' + str(row['Collision_energy']) +'eV'+ '\n'
        entry = entry+'SPECTRUMTYPE: Centroid and composite'+ '\n'
        entry = entry + 'METABOLITENAME: ' + '\n'
        entry = entry + 'SCANNUMBER: Alignment ID '+str(row['Alignment_ID']) + '\n'
        entry = entry + 'RETENTIONTIME: ' + str(row['RETENTIONTIME']) + '\n'
        entry = entry + 'RETENTIONINDEX: N/A' +'\n'
        entry = entry + 'CCS: ' +'\n'
        entry = entry + 'INTENSITY: ' +str(row['intensity'])+'\n'
        entry = entry + '#Specific field for labeled experiment' +'\n'
        entry = entry + 'IsMarked: False' +'\n'
        entry = entry + 'Comment: ' + str(row['Comment']) + '\n'
        entry = entry + 'MSTYPE: MS1'+ '\n'
        entry = entry + 'Num Peaks: ' + str(num_peaks(row[typeofms1])) + '\n'
        mass1, intensity1= break_spectra(row[typeofms1])
        for i in range(0, len(mass1)):
            entry = entry + str(mass1[i]) +'\t'+ str(intensity1[i]) + '\t'+'"'+str(mass1[i])+'"'+'\n'
        entry = entry + 'MSTYPE: MS2'+ '\n'
        entry = entry + 'Num Peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        mass2, intensity2= break_spectra(row[typeofmsms])
        for i in range(0, len(mass2)):
            entry = entry + str(mass2[i]) +'\t'+ str(intensity2[i]) + '\t'+'"'+str(mass2[i])+'"'+'\n'
        entry = entry + '\n'
        #
        #
        # entry = entry +'Spectrum_type: '+row['Spectrum_type']+ '\n'
        #
        # entry = entry + 'InChIKey: ' + row['InChIKey'] + '\n'
        #
        # entry = entry + 'ExactMass: ' + str(row['ExactMass']) + '\n'
        #
        #
        #
        # entry = entry + 'Num peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        # entry = entry + row[typeofmsms]
        # # entry = entry +str(row['count'])

        # entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')

    #write string to file
    text_file.write(entry)

    #close file
    text_file.close()

#%%
