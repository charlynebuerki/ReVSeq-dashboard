import pandas as pd
import numpy as np

df=pd.read_csv('dashboard/static/example_data.csv')
strains=df['strain_name'].unique()

coverage=pd.read_csv('dashboard/static/pileup_all.csv')
coverage.set_index('idx')
coverage=coverage.drop(columns='idx')

for strain in strains:
    print(strain)
    df_strain=df[df['strain_name']==strain]

    match strain:
        case 'Rhino- / Enterovirus':
            strain_name = 'Enterovirus'
        case 'RSV - A/B':
            strain_name = 'RSV'
        case _:
            strain_name= strain.replace(' ' ,'_')

    match strain:
        case 'Influenza A':
            segments = {'PB2': ['KU509700.1', 'KJ942623.1'], 'PB1': ['KU509701.1', 'KJ942622.1'], 'PA': ['KU509702.1', 'KJ942621.1'], 'HA': ['KU509703.1', 'KJ942616.1'],'NP': ['KU509704.1', 'KJ942619.1'], 'NA': ['KU509705.1', 'KJ942618.1'],'MP': ['KU509706.1', 'KJ942617.1'], 'NS': ['KU509707.1', 'KJ942620.1']}
            for seg in segments.keys():
                cols=[s for s in coverage.columns if (s.split('_')[0] in df_strain['pseudonymized_id'].tolist()) &  (s.split('_')[1] in segments[seg])] 
                coverage_strain=coverage[cols].fillna(0)
                coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                pileup.to_csv('dashboard/static/pileup/Influenza_A_'+seg+'_all.csv', index=False)
            print('\tH1N1')
            for seg in segments.keys():
                cols=[s for s in coverage.columns if (s.split('_')[0] in df_strain['pseudonymized_id'].tolist()) &  (s.split('_')[1] == segments[seg][0])] 
                coverage_strain=coverage[cols].fillna(0)
                coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                pileup.to_csv('dashboard/static/pileup/Influenza_A_'+seg+'_H1N1.csv', index=False)
            print('\tH3N2')
            for seg in segments.keys():
                cols=[s for s in coverage.columns if (s.split('_')[0] in df_strain['pseudonymized_id'].tolist()) &  (s.split('_')[1] == segments[seg][1])] 
                coverage_strain=coverage[cols].fillna(0)
                coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                pileup.to_csv('dashboard/static/pileup/Influenza_A_'+seg+'_H3N2.csv', index=False)

        case 'Influenza B':
            segments = {'PB2': ['KC866604.1'], 'PB1': ['KC866603.1'], 'HA': ['KX058884.1'],'NP': ['KC866605.1'], 'NA': ['FJ766839.1'],'MP': ['KC866607.1'], 'NS': ['KC866606.1']} #'PA': ['MN325123.1'],
            for seg in segments.keys():
                cols=[s for s in coverage.columns if (s.split('_')[0] in df_strain['pseudonymized_id'].tolist()) and  (s.split('_')[1] in segments[seg])] 
                coverage_strain=coverage[cols].fillna(0)
                coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                pileup.to_csv('dashboard/static/pileup/Influenza_B_'+seg+'_all.csv', index=False)
        case _:
            cols=[s for s in coverage.columns if s.split('_')[0] in df_strain['pseudonymized_id'].tolist()] 
            if not len(cols):
                print('\t'+'Missing coverage info')
            else:
                coverage_strain=coverage[cols].fillna(0)
                coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                pileup.to_csv('dashboard/static/pileup/'+strain_name+'_all.csv', index=False)

                substrains = df_strain['substrain_name'].unique()
                if len(substrains) > 1:
                    for substrain in substrains:
                        print('\t'+substrain)
                        df_strain=df[df['substrain_name']==substrain]

                        cols=[s for s in coverage.columns if s.split('_')[0] in df_strain['pseudonymized_id'].tolist()] 
                        coverage_substrain=coverage[cols].fillna(0)
                        coverage_substrain=coverage_substrain[:np.where(coverage_substrain.sum(axis=1))[0][-1]]
                        pileup= pd.DataFrame({'idx': coverage_substrain.index+1, 'mean': coverage_substrain.mean(axis=1), 'ci_lower': coverage_substrain.mean(axis=1) - 1.96 *coverage_substrain.std(axis=1)/ np.sqrt(len(coverage_substrain.columns)), 'ci_upper': coverage_substrain.mean(axis=1) + 1.96 *coverage_substrain.std(axis=1)/ np.sqrt(len(coverage_substrain.columns))})
                        pileup.to_csv('dashboard/static/pileup/'+strain_name+'_'+substrain.replace(' ', '_')+ '.csv', index=False)
