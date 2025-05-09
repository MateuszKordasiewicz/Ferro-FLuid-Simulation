import pandas as pd
import numpy as np

def error_of_vectors(df, original_df):
    error_dfs = []
    char_lens = df['CharLen'].unique()
    for i in range(len(char_lens)-1):
        char_len_0 = char_lens[i]
        char_len_1 = char_lens[i+1]
        
        char_len_0_len = original_df['CharLen'][original_df['CharLen'] == char_len_0].size
        char_len_1_len = original_df['CharLen'][original_df['CharLen'] == char_len_1].size
        
        print("Calculating error: " + str(char_len_0) + " vs. " + str(char_len_1))
        
        df_char_0 = df[df['CharLen'] == char_len_0]
        df_char_1 = df[df['CharLen'] == char_len_1]
        
        error_height = (df_char_1['Height'] - df_char_0['Height'])/(char_len_1_len - char_len_0_len)
        error_x_velocity = (df_char_1['X_velocity'] - df_char_0['X_velocity'])/(char_len_1_len - char_len_0_len)
        error_y_velocity = (df_char_1['Y_velocity'] - df_char_0['Y_velocity'])/(char_len_1_len - char_len_0_len)
        
        char_len_comparison = f"{char_len_0} vs. {char_len_1}"
        
        error = pd.DataFrame({
            'CharLen': char_len_comparison,
            'X': df_char_0['X'],
            'Y': df_char_0['Y'],
            'Height_error': error_height,
            'X_velocity_error': error_x_velocity,
            'Y_velocity_error': error_y_velocity
        })
        
        error_dfs.append(error)
        
    df_error_all = pd.concat(error_dfs)
    return df_error_all