import pandas as pd
import numpy as np

def error_of_vectors(df, original_df):
    error_dfs = []
    CFL_lens = df['CFL'].unique()
    for i in range(len(CFL_lens)-1):
    
        CFL_len_0 = CFL_lens[i]
        CFL_len_1 = CFL_lens[i+1]
        
        print("Calculating error: " + str(CFL_len_0) + " vs. " + str(CFL_len_1))
        
        df_CFL_0 = df[df['CFL'] == CFL_len_0]
        df_CFL_1 = df[df['CFL'] == CFL_len_1]
        
        error_height = (df_CFL_1['Height'] - df_CFL_0['Height'])
        error_x_velocity = (df_CFL_1['X_velocity'] - df_CFL_0['X_velocity'])
        error_y_velocity = (df_CFL_1['Y_velocity'] - df_CFL_0['Y_velocity'])
        
        CFL_len_comparison = f"{CFL_len_0} vs. {CFL_len_1}"
        
        error = pd.DataFrame({
            'CFL': CFL_len_comparison,
            'X': df_CFL_0['X'],
            'Y': df_CFL_0['Y'],
            'Height_error': error_height,
            'X_velocity_error': error_x_velocity,
            'Y_velocity_error': error_y_velocity
        })
        
        error_dfs.append(error)
        
    df_error_all = pd.concat(error_dfs)
    return df_error_all