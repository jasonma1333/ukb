#!/usr/bin/env python3
"""
UK Biobank Protein Data Processing Scr        for file in sorted(olink1_files):
         for file in sorted(olink2_files):
            print(f"  Processing {file.name}...")
            df = pd.read_csv(file)
            
            if len(df.columns) > 1:
                id_col = df.columns[1]
                df = df.rename(columns={id_col: 'Participant_ID'})
                
                if 'olink_instance' in str(df.columns[0]).lower():
                    df = df.drop(columns=[df.columns[0]])
                
                df = df.set_index('Participant_ID')
                
                # Add instance number to avoid column conflicts
                instance_num = file.name.split('_')[-1].split('.')[0]
                df.columns = [f"olink2_inst{instance_num}_{col}" for col in df.columns]
                
                if combined_olink is None:
                    combined_olink = df.copy()
                else:
                    combined_olink = combined_olink.join(df, how='outer')"  Processing {file.name}...")
            df = pd.read_csv(file)
            
            if len(df.columns) > 1:
                id_col = df.columns[1]
                df = df.rename(columns={id_col: 'Participant_ID'})
                
                if 'olink_instance' in str(df.columns[0]).lower():
                    df = df.drop(columns=[df.columns[0]])
                
                df = df.set_index('Participant_ID')
                
                # Add instance number to avoid column conflicts
                instance_num = file.name.split('_')[-1].split('.')[0]
                df.columns = [f"olink1_inst{instance_num}_{col}" for col in df.columns]
                
                if combined_olink is None:
                    combined_olink = df.copy()
                else:
                    combined_olink = combined_olink.join(df, how='outer')cleans multiple datasets based on patient ID
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import warnings
import pickle
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore')

def load_sample_format():
    """Load the sample.csv to understand target format structure."""
    print("Loading sample format...")
    sample_df = pd.read_csv('sample.csv', index_col=0)
    print(f"Sample format: {sample_df.shape[0]} rows, {sample_df.shape[1]} columns")
    return sample_df

def combine_olink_data():
    """Combine all olink protein datasets."""
    print("Processing Olink protein data...")
    
    # Initialize combined dataframe
    combined_olink = None
    
    # Process olink-0 data (multiple instances)
    olink0_path = Path('olink-0')
    if olink0_path.exists():
        print("Processing olink-0 data...")
        olink0_files = list(olink0_path.glob('data_olink_instance_*.csv'))
        
        for file in sorted(olink0_files):
            print(f"  Processing {file.name}...")
            df = pd.read_csv(file)
            
            # Use the second column as participant ID (standard format)
            if len(df.columns) > 1:
                id_col = df.columns[1]  # Usually 'Participant ID(participant - eid)'
                df = df.rename(columns={id_col: 'Participant_ID'})
                
                # Drop the first column if it's a duplicate ID column
                if 'olink_instance' in str(df.columns[0]).lower():
                    df = df.drop(columns=[df.columns[0]])
                
                # Set participant ID as index
                df = df.set_index('Participant_ID')
                
                # Add suffix to column names to indicate olink-0 source and instance
                instance_num = file.name.split('_')[-1].split('.')[0]  # Extract instance number
                df.columns = [f"olink0_inst{instance_num}_{col}" for col in df.columns]
                
                # Combine with existing data
                if combined_olink is None:
                    combined_olink = df.copy()
                else:
                    combined_olink = combined_olink.join(df, how='outer')
    
    # Process olink-1 data
    olink1_path = Path('olink-1')
    if olink1_path.exists():
        print("Processing olink-1 data...")
        olink1_files = list(olink1_path.glob('data_olink_instance_*.csv'))
        
        for file in sorted(olink1_files):
            print(f"  Processing {file.name}...")
            df = pd.read_csv(file)
            
            if len(df.columns) > 1:
                id_col = df.columns[1]
                df = df.rename(columns={id_col: 'Participant_ID'})
                
                if 'olink_instance' in str(df.columns[0]).lower():
                    df = df.drop(columns=[df.columns[0]])
                
                df = df.set_index('Participant_ID')
                df.columns = [f"olink1_{col}" for col in df.columns]
                
                if combined_olink is None:
                    combined_olink = df.copy()
                else:
                    combined_olink = combined_olink.join(df, how='outer', rsuffix='_dup')
    
    # Process olink2 data
    olink2_path = Path('olink2')
    if olink2_path.exists():
        print("Processing olink2 data...")
        olink2_files = list(olink2_path.glob('data_olink_instance_*.csv'))
        
        for file in sorted(olink2_files):
            print(f"  Processing {file.name}...")
            df = pd.read_csv(file)
            
            if len(df.columns) > 1:
                id_col = df.columns[1]
                df = df.rename(columns={id_col: 'Participant_ID'})
                
                if 'olink_instance' in str(df.columns[0]).lower():
                    df = df.drop(columns=[df.columns[0]])
                
                df = df.set_index('Participant_ID')
                df.columns = [f"olink2_{col}" for col in df.columns]
                
                if combined_olink is None:
                    combined_olink = df.copy()
                else:
                    combined_olink = combined_olink.join(df, how='outer', rsuffix='_dup')
    
    if combined_olink is not None:
        print(f"Combined Olink data: {combined_olink.shape[0]} participants, {combined_olink.shape[1]} protein measurements")
    else:
        print("No Olink data found")
        combined_olink = pd.DataFrame()
    
    return combined_olink

def load_cancer_data():
    """Load and combine cancer diagnosis data."""
    print("Processing cancer diagnosis data...")
    
    cancer_path = Path('ICD10 -癌症，三个指标连着用')
    combined_cancer = None
    
    if cancer_path.exists():
        # Load age at cancer diagnosis
        age_file = cancer_path / 'Age at cancer diagnosis 40008.csv'
        if age_file.exists():
            print("  Loading age at cancer diagnosis...")
            age_df = pd.read_csv(age_file)
            age_df = age_df.set_index('Participant ID')
            age_df.columns = [f"cancer_age_{col}" if 'Instance' in col else col for col in age_df.columns]
            combined_cancer = age_df
        
        # Load date of cancer diagnosis
        date_file = cancer_path / 'Date of cancer diagnosis 40005.csv'
        if date_file.exists():
            print("  Loading date of cancer diagnosis...")
            date_df = pd.read_csv(date_file)
            date_df = date_df.set_index('Participant ID')
            date_df.columns = [f"cancer_date_{col}" if 'Instance' in col else col for col in date_df.columns]
            
            if combined_cancer is None:
                combined_cancer = date_df
            else:
                combined_cancer = combined_cancer.join(date_df, how='outer')
        
        # Load type of cancer (Excel file)
        type_file = cancer_path / 'type of cancer ICD10 40006.xlsx'
        if type_file.exists():
            print("  Loading type of cancer diagnosis...")
            try:
                type_df = pd.read_excel(type_file)
                if 'Participant ID' in type_df.columns:
                    type_df = type_df.set_index('Participant ID')
                    type_df.columns = [f"cancer_type_{col}" if 'Instance' in col else col for col in type_df.columns]
                    
                    if combined_cancer is None:
                        combined_cancer = type_df
                    else:
                        combined_cancer = combined_cancer.join(type_df, how='outer')
            except Exception as e:
                print(f"    Error loading Excel file: {e}")
    
    if combined_cancer is not None:
        print(f"Combined cancer data: {combined_cancer.shape[0]} participants, {combined_cancer.shape[1]} features")
        # Reset index to have Participant ID as a column for consistency
        combined_cancer = combined_cancer.reset_index()
        combined_cancer = combined_cancer.rename(columns={'Participant ID': 'Participant_ID'})
        combined_cancer = combined_cancer.set_index('Participant_ID')
    else:
        print("No cancer data found")
        combined_cancer = pd.DataFrame()
    
    return combined_cancer

def load_basic_data():
    """Load basic demographic and examination data."""
    print("Processing basic data...")
    
    basic_path = Path('基本')
    combined_basic = None
    
    if basic_path.exists():
        excel_files = list(basic_path.glob('*.xlsx'))
        
        for file in excel_files:
            print(f"  Processing {file.name}...")
            try:
                df = pd.read_excel(file)
                
                # Look for participant ID column (various possible names)
                id_columns = [col for col in df.columns if 'participant' in col.lower() or 'id' in col.lower() or 'eid' in col.lower()]
                
                if id_columns:
                    id_col = id_columns[0]
                    df = df.rename(columns={id_col: 'Participant_ID'})
                    df = df.set_index('Participant_ID')
                    
                    # Add prefix to avoid column name conflicts
                    file_prefix = file.stem.replace(' ', '_')
                    df.columns = [f"{file_prefix}_{col}" for col in df.columns]
                    
                    if combined_basic is None:
                        combined_basic = df.copy()
                    else:
                        combined_basic = combined_basic.join(df, how='outer')
                else:
                    print(f"    Warning: No participant ID column found in {file.name}")
                    
            except Exception as e:
                print(f"    Error loading {file.name}: {e}")
    
    if combined_basic is not None:
        print(f"Combined basic data: {combined_basic.shape[0]} participants, {combined_basic.shape[1]} features")
    else:
        print("No basic data could be loaded")
        combined_basic = pd.DataFrame()
    
    return combined_basic

def clean_and_standardize_data(df):
    """Clean and standardize the combined dataset."""
    print("Cleaning and standardizing data...")
    
    # Convert all numeric columns to float
    numeric_columns = []
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            numeric_columns.append(col)
        except:
            pass
    
    print(f"  Converted {len(numeric_columns)} columns to numeric")
    
    # Remove columns that are entirely NaN
    initial_cols = df.shape[1]
    df = df.dropna(axis=1, how='all')
    print(f"  Removed {initial_cols - df.shape[1]} empty columns")
    
    # Remove rows that are entirely NaN
    initial_rows = df.shape[0]
    df = df.dropna(axis=0, how='all')
    print(f"  Removed {initial_rows - df.shape[0]} empty rows")
    
    # Calculate missing data statistics
    missing_stats = df.isnull().sum().sort_values(ascending=False)
    missing_percent = (missing_stats / len(df) * 100).round(2)
    
    print(f"  Final dataset shape: {df.shape}")
    print(f"  Columns with >50% missing data: {sum(missing_percent > 50)}")
    
    return df, missing_stats

def handle_missing_data(df, strategy='conservative', missing_threshold=0.5):
    """
    Handle missing data using various strategies.
    
    Parameters:
    - df: DataFrame with missing data
    - strategy: 'conservative', 'moderate', 'aggressive', or 'imputation'
    - missing_threshold: Threshold for dropping columns (e.g., 0.5 = 50% missing)
    
    Returns:
    - DataFrame with missing data handled
    - Dictionary with processing statistics
    """
    print(f"Handling missing data using '{strategy}' strategy...")
    print("=" * 50)
    
    df_processed = df.copy()
    initial_shape = df_processed.shape
    stats = {}
    
    # Calculate initial missing data statistics
    missing_stats = df_processed.isnull().sum()
    missing_percent = (missing_stats / len(df_processed) * 100)
    
    print(f"Initial shape: {initial_shape}")
    print(f"Initial missing data: {missing_stats.sum():,} cells ({(missing_stats.sum() / (initial_shape[0] * initial_shape[1]) * 100):.2f}%)")
    
    # Identify olink columns for special handling
    olink_columns = [col for col in df_processed.columns if 'olink' in col.lower()]
    other_columns = [col for col in df_processed.columns if 'olink' not in col.lower()]
    
    print(f"Olink protein columns: {len(olink_columns)}")
    print(f"Other columns: {len(other_columns)}")
    
    if strategy == 'conservative':
        print("\n🔒 CONSERVATIVE STRATEGY")
        print("- Remove columns with >50% missing data")
        print("- Keep only participants with complete protein data")
        
        # Remove high-missing columns
        high_missing_cols = missing_percent[missing_percent > 50].index.tolist()
        df_processed = df_processed.drop(columns=high_missing_cols)
        print(f"  Removed {len(high_missing_cols)} columns with >50% missing data")
        
        # For olink columns, keep only participants with complete data
        olink_cols_remaining = [col for col in df_processed.columns if 'olink' in col.lower()]
        if olink_cols_remaining:
            complete_olink_mask = df_processed[olink_cols_remaining].notna().all(axis=1)
            df_processed = df_processed[complete_olink_mask]
            print(f"  Kept {complete_olink_mask.sum():,} participants with complete protein data")
    
    elif strategy == 'moderate':
        print("\n⚖️ MODERATE STRATEGY")
        print("- Remove columns with >70% missing data")
        print("- Keep participants with >80% complete protein data")
        
        # Remove very high-missing columns
        high_missing_cols = missing_percent[missing_percent > 70].index.tolist()
        df_processed = df_processed.drop(columns=high_missing_cols)
        print(f"  Removed {len(high_missing_cols)} columns with >70% missing data")
        
        # Keep participants with mostly complete olink data
        olink_cols_remaining = [col for col in df_processed.columns if 'olink' in col.lower()]
        if olink_cols_remaining:
            olink_completeness = df_processed[olink_cols_remaining].notna().mean(axis=1)
            good_participants_mask = olink_completeness >= 0.8
            df_processed = df_processed[good_participants_mask]
            print(f"  Kept {good_participants_mask.sum():,} participants with >80% complete protein data")
    
    elif strategy == 'aggressive':
        print("\n🚀 AGGRESSIVE STRATEGY")
        print("- Remove columns with >80% missing data")
        print("- Keep participants with >60% complete protein data")
        
        # Remove extremely high-missing columns
        high_missing_cols = missing_percent[missing_percent > 80].index.tolist()
        df_processed = df_processed.drop(columns=high_missing_cols)
        print(f"  Removed {len(high_missing_cols)} columns with >80% missing data")
        
        # More lenient participant filtering
        olink_cols_remaining = [col for col in df_processed.columns if 'olink' in col.lower()]
        if olink_cols_remaining:
            olink_completeness = df_processed[olink_cols_remaining].notna().mean(axis=1)
            good_participants_mask = olink_completeness >= 0.6
            df_processed = df_processed[good_participants_mask]
            print(f"  Kept {good_participants_mask.sum():,} participants with >60% complete protein data")
    
    elif strategy == 'imputation':
        print("\n🔄 IMPUTATION STRATEGY")
        print("- Remove columns with >90% missing data")
        print("- Impute missing values using median/mode")
        
        # Remove extremely sparse columns
        high_missing_cols = missing_percent[missing_percent > 90].index.tolist()
        df_processed = df_processed.drop(columns=high_missing_cols)
        print(f"  Removed {len(high_missing_cols)} columns with >90% missing data")
        
        # Impute missing values
        olink_cols_remaining = [col for col in df_processed.columns if 'olink' in col.lower()]
        other_cols_remaining = [col for col in df_processed.columns if 'olink' not in col.lower()]
        
        # Impute olink columns with median
        if olink_cols_remaining:
            for col in olink_cols_remaining:
                if df_processed[col].isnull().any():
                    median_val = df_processed[col].median()
                    df_processed[col].fillna(median_val, inplace=True)
            print(f"  Imputed {len(olink_cols_remaining)} protein columns with median values")
        
        # Impute other columns with median for numeric, mode for categorical
        imputed_other = 0
        for col in other_cols_remaining:
            if df_processed[col].isnull().any():
                if df_processed[col].dtype in ['float64', 'float32', 'int64', 'int32']:
                    median_val = df_processed[col].median()
                    if pd.notna(median_val):
                        df_processed[col].fillna(median_val, inplace=True)
                        imputed_other += 1
                else:
                    mode_val = df_processed[col].mode()
                    if len(mode_val) > 0:
                        df_processed[col].fillna(mode_val[0], inplace=True)
                        imputed_other += 1
        print(f"  Imputed {imputed_other} other columns")
    
    # Calculate final statistics
    final_shape = df_processed.shape
    final_missing = df_processed.isnull().sum().sum()
    final_missing_percent = (final_missing / (final_shape[0] * final_shape[1]) * 100)
    
    stats = {
        'initial_shape': initial_shape,
        'final_shape': final_shape,
        'participants_removed': initial_shape[0] - final_shape[0],
        'columns_removed': initial_shape[1] - final_shape[1],
        'initial_missing_percent': (missing_stats.sum() / (initial_shape[0] * initial_shape[1]) * 100),
        'final_missing_percent': final_missing_percent,
        'strategy_used': strategy
    }
    
    print(f"\n📊 MISSING DATA HANDLING RESULTS:")
    print(f"  Initial shape: {initial_shape}")
    print(f"  Final shape: {final_shape}")
    print(f"  Participants removed: {stats['participants_removed']:,}")
    print(f"  Columns removed: {stats['columns_removed']}")
    print(f"  Missing data reduction: {stats['initial_missing_percent']:.2f}% → {stats['final_missing_percent']:.2f}%")
    
    return df_processed, stats

def apply_log_transformation(df, olink_columns=None):
    """
    Apply log(1+x) transformation to protein biomarker data.
    
    Parameters:
    - df: DataFrame with protein data
    - olink_columns: List of olink columns to transform. If None, auto-detect.
    
    Returns:
    - DataFrame with log-transformed data
    """
    print("Applying log(1+x) transformation to protein biomarkers...")
    
    df_transformed = df.copy()
    
    # Auto-detect olink columns if not provided
    if olink_columns is None:
        olink_columns = [col for col in df.columns if 'olink' in col.lower()]
    
    print(f"  Found {len(olink_columns)} olink protein columns")
    
    # Apply log(1+x) transformation to olink columns
    log_transformed_count = 0
    for col in olink_columns:
        if col in df_transformed.columns:
            # Only transform if the column contains numeric data
            if df_transformed[col].dtype in ['float64', 'float32', 'int64', 'int32']:
                # Ensure non-negative values for log transformation
                min_val = df_transformed[col].min()
                if pd.notna(min_val):
                    if min_val < 0:
                        print(f"    Warning: Column {col} has negative values. Shifting by {abs(min_val)}")
                        df_transformed[col] = df_transformed[col] - min_val
                    
                    # Apply log(1+x) transformation
                    df_transformed[col] = np.log1p(df_transformed[col])
                    log_transformed_count += 1
    
    print(f"  Successfully applied log transformation to {log_transformed_count} columns")
    return df_transformed

def apply_standardization(df, scaler_dict_path='scaler_dict.p', olink_columns=None):
    """
    Apply standardization using pre-trained StandardScaler objects.
    
    Parameters:
    - df: DataFrame with data to standardize
    - scaler_dict_path: Path to the saved scaler dictionary
    - olink_columns: List of olink columns to standardize. If None, auto-detect.
    
    Returns:
    - DataFrame with standardized data
    """
    print("Applying standardization using pre-trained scalers...")
    
    df_standardized = df.copy()
    
    # Auto-detect olink columns if not provided
    if olink_columns is None:
        olink_columns = [col for col in df.columns if 'olink' in col.lower()]
    
    # Try to load the scaler dictionary
    scaler_dict_file = Path(scaler_dict_path)
    scaler_dict = None
    
    if scaler_dict_file.exists():
        try:
            with open(scaler_dict_file, 'rb') as f:
                scaler_dict = pickle.load(f)
            print(f"  Loaded scaler dictionary from {scaler_dict_path}")
            print(f"  Available scalers: {len(scaler_dict) if scaler_dict else 0}")
        except Exception as e:
            print(f"  Warning: Could not load scaler dictionary: {e}")
            scaler_dict = None
    else:
        print(f"  Warning: Scaler dictionary file {scaler_dict_path} not found")
    
    # Apply standardization
    if scaler_dict is not None:
        standardized_count = 0
        for col in olink_columns:
            if col in df_standardized.columns and col in scaler_dict:
                try:
                    # Get the scaler for this column
                    scaler = scaler_dict[col]
                    
                    # Apply standardization to non-null values
                    mask = df_standardized[col].notna()
                    if mask.sum() > 0:
                        values = df_standardized.loc[mask, col].values.reshape(-1, 1)
                        scaled_values = scaler.transform(values)
                        df_standardized.loc[mask, col] = scaled_values.flatten()
                        standardized_count += 1
                except Exception as e:
                    print(f"    Warning: Could not standardize column {col}: {e}")
        
        print(f"  Successfully standardized {standardized_count} columns using pre-trained scalers")
    else:
        # Fallback: Create new StandardScalers if no pre-trained scalers available
        print("  No pre-trained scalers available. Creating new StandardScalers...")
        new_scaler_dict = {}
        standardized_count = 0
        
        for col in olink_columns:
            if col in df_standardized.columns:
                if df_standardized[col].dtype in ['float64', 'float32', 'int64', 'int32']:
                    try:
                        # Create and fit new scaler
                        scaler = StandardScaler()
                        mask = df_standardized[col].notna()
                        
                        if mask.sum() > 0:
                            values = df_standardized.loc[mask, col].values.reshape(-1, 1)
                            scaled_values = scaler.fit_transform(values)
                            df_standardized.loc[mask, col] = scaled_values.flatten()
                            new_scaler_dict[col] = scaler
                            standardized_count += 1
                    except Exception as e:
                        print(f"    Warning: Could not standardize column {col}: {e}")
        
        print(f"  Created and applied {standardized_count} new StandardScalers")
        
        # Save the new scaler dictionary
        try:
            with open('new_scaler_dict.p', 'wb') as f:
                pickle.dump(new_scaler_dict, f)
            print(f"  Saved new scaler dictionary to new_scaler_dict.p")
        except Exception as e:
            print(f"  Warning: Could not save new scaler dictionary: {e}")
    
    return df_standardized

def main():
    """Main data processing function."""
    print("Starting UK Biobank Protein Data Processing...")
    print("=" * 60)
    
    # Change to the data directory
    os.chdir(r'c:\Users\BMS1A-315\Desktop\ukb_protein')
    
    # Load sample format for reference
    sample_df = load_sample_format()
    
    # Combine all datasets
    print("\n" + "=" * 60)
    
    # Load Olink protein data
    olink_data = combine_olink_data()
    
    # Load cancer data
    cancer_data = load_cancer_data()
    
    # Load basic demographic data
    basic_data = load_basic_data()
    
    print("\n" + "=" * 60)
    print("Combining all datasets...")
    
    # Start with olink data as the base (most comprehensive)
    if not olink_data.empty:
        combined_df = olink_data.copy()
        print(f"Base dataset (Olink): {combined_df.shape}")
    else:
        combined_df = pd.DataFrame()
    
    # Add cancer data
    if not cancer_data.empty and not combined_df.empty:
        combined_df = combined_df.join(cancer_data, how='outer')
        print(f"After adding cancer data: {combined_df.shape}")
    elif not cancer_data.empty:
        combined_df = cancer_data.copy()
        print(f"Base dataset (Cancer): {combined_df.shape}")
    
    # Add basic data
    if not basic_data.empty and not combined_df.empty:
        combined_df = combined_df.join(basic_data, how='outer')
        print(f"After adding basic data: {combined_df.shape}")
    elif not basic_data.empty:
        combined_df = basic_data.copy()
        print(f"Base dataset (Basic): {combined_df.shape}")
    
    if combined_df.empty:
        print("Error: No data could be loaded!")
        return
    
    # Clean and standardize the data
    print("\n" + "=" * 60)
    final_df, missing_stats = clean_and_standardize_data(combined_df)
    
    # Handle missing data
    print("\n" + "=" * 60)
    final_df, missing_data_stats = handle_missing_data(final_df, strategy='moderate')
    
    # Apply log transformation to protein biomarkers
    print("\n" + "=" * 60)
    final_df = apply_log_transformation(final_df)
    
    # Apply standardization using pre-trained scalers
    print("\n" + "=" * 60)
    final_df = apply_standardization(final_df, scaler_dict_path='scaler_dict.p')
    
    # Save the processed data
    print("\n" + "=" * 60)
    print("Saving processed data...")
    
    # Save main combined dataset (with log transformation and standardization)
    output_file = 'combined_ukb_protein_data_processed.csv'
    final_df.to_csv(output_file)
    print(f"Saved processed dataset (log-transformed + standardized): {output_file}")
    print(f"Final dataset shape: {final_df.shape}")
    
    # Save protein biomarkers only (processed)
    olink_columns = [col for col in final_df.columns if 'olink' in col.lower()]
    if olink_columns:
        biomarkers_only = final_df[olink_columns]
        biomarkers_file = 'protein_biomarkers_processed.csv'
        biomarkers_only.to_csv(biomarkers_file)
        print(f"Saved processed protein biomarkers only: {biomarkers_file}")
        print(f"Protein biomarkers shape: {biomarkers_only.shape}")
    
    # Save missing data handling report
    missing_data_report = pd.DataFrame({
        'Strategy': [missing_data_stats['strategy_used']],
        'Initial_Participants': [missing_data_stats['initial_shape'][0]],
        'Final_Participants': [missing_data_stats['final_shape'][0]],
        'Initial_Features': [missing_data_stats['initial_shape'][1]],
        'Final_Features': [missing_data_stats['final_shape'][1]],
        'Participants_Removed': [missing_data_stats['participants_removed']],
        'Features_Removed': [missing_data_stats['columns_removed']],
        'Initial_Missing_Percent': [missing_data_stats['initial_missing_percent']],
        'Final_Missing_Percent': [missing_data_stats['final_missing_percent']]
    })
    missing_data_report.to_csv('missing_data_handling_report.csv', index=False)
    print(f"Saved missing data handling report: missing_data_handling_report.csv")
    
    # Save detailed missing data statistics
    final_missing_stats = final_df.isnull().sum().sort_values(ascending=False)
    missing_report = pd.DataFrame({
        'Column': final_missing_stats.index,
        'Missing_Count': final_missing_stats.values,
        'Missing_Percent': (final_missing_stats / len(final_df) * 100).round(2)
    })
    missing_report.to_csv('detailed_missing_data_report.csv', index=False)
    print(f"Saved detailed missing data report: detailed_missing_data_report.csv")
    
    # Create alternative datasets with different missing data strategies
    print(f"\nCreating alternative datasets with different missing data strategies...")
    
    # Conservative approach
    try:
        conservative_df, conservative_stats = handle_missing_data(combined_df, strategy='conservative')
        conservative_df = apply_log_transformation(conservative_df)
        conservative_df = apply_standardization(conservative_df, scaler_dict_path='scaler_dict.p')
        conservative_df.to_csv('conservative_processed_data.csv')
        print(f"Saved conservative dataset: conservative_processed_data.csv ({conservative_df.shape})")
    except Exception as e:
        print(f"Warning: Could not create conservative dataset: {e}")
    
    # Aggressive approach
    try:
        aggressive_df, aggressive_stats = handle_missing_data(combined_df, strategy='aggressive')
        aggressive_df = apply_log_transformation(aggressive_df)
        aggressive_df = apply_standardization(aggressive_df, scaler_dict_path='scaler_dict.p')
        aggressive_df.to_csv('aggressive_processed_data.csv')
        print(f"Saved aggressive dataset: aggressive_processed_data.csv ({aggressive_df.shape})")
    except Exception as e:
        print(f"Warning: Could not create aggressive dataset: {e}")
    
    # Imputation approach
    try:
        imputed_df, imputed_stats = handle_missing_data(combined_df, strategy='imputation')
        imputed_df = apply_log_transformation(imputed_df)
        imputed_df = apply_standardization(imputed_df, scaler_dict_path='scaler_dict.p')
        imputed_df.to_csv('imputed_processed_data.csv')
        print(f"Saved imputed dataset: imputed_processed_data.csv ({imputed_df.shape})")
    except Exception as e:
        print(f"Warning: Could not create imputed dataset: {e}")
    
    # Save a sample of the data for inspection
    sample_size = min(100, len(final_df))
    sample_data = final_df.head(sample_size)
    sample_data.to_csv('sample_combined_data_processed.csv')
    print(f"Saved sample data ({sample_size} rows): sample_combined_data_processed.csv")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("PROCESSING SUMMARY")
    print("=" * 60)
    print(f"Total participants: {len(final_df):,}")
    print(f"Total features: {len(final_df.columns)}")
    print(f"Data completeness: {((1 - final_df.isnull().sum().sum() / (final_df.shape[0] * final_df.shape[1])) * 100):.2f}%")
    
    # Count protein biomarkers
    olink_columns = [col for col in final_df.columns if 'olink' in col.lower()]
    print(f"Protein biomarkers: {len(olink_columns)}")
    
    # Show missing data handling results
    print(f"\nMissing data handling results:")
    print(f"  Strategy used: {missing_data_stats['strategy_used']}")
    print(f"  Participants: {missing_data_stats['initial_shape'][0]:,} → {missing_data_stats['final_shape'][0]:,}")
    print(f"  Features: {missing_data_stats['initial_shape'][1]} → {missing_data_stats['final_shape'][1]}")
    print(f"  Missing data: {missing_data_stats['initial_missing_percent']:.2f}% → {missing_data_stats['final_missing_percent']:.2f}%")
    
    # Show processing steps applied
    print(f"\nProcessing steps applied:")
    print(f"  ✓ Data cleaning and standardization")
    print(f"  ✓ Missing data handling ({missing_data_stats['strategy_used']} strategy)")
    print(f"  ✓ Log(1+x) transformation on protein biomarkers")
    print(f"  ✓ Standardization using pre-trained scalers")
    
    # Show protein biomarker statistics
    if olink_columns:
        print(f"\nProtein biomarker statistics (after processing):")
        biomarker_stats = final_df[olink_columns].describe()
        print(f"  Mean range: {biomarker_stats.loc['mean'].min():.3f} to {biomarker_stats.loc['mean'].max():.3f}")
        print(f"  Std range: {biomarker_stats.loc['std'].min():.3f} to {biomarker_stats.loc['std'].max():.3f}")
        print(f"  Missing values: {final_df[olink_columns].isnull().sum().sum():,}")
    
    # Show files created
    print(f"\nFiles created:")
    print(f"  📊 Main dataset: {output_file}")
    print(f"  🔬 Biomarkers only: {'protein_biomarkers_processed.csv' if olink_columns else 'N/A'}")
    print(f"  📋 Missing data reports: missing_data_handling_report.csv, detailed_missing_data_report.csv")
    print(f"  🎯 Alternative datasets: conservative_processed_data.csv, aggressive_processed_data.csv, imputed_processed_data.csv")
    
    # Show top columns by completeness
    print(f"\nTop 10 most complete features:")
    completeness = (1 - final_df.isnull().sum() / len(final_df)) * 100
    top_complete = completeness.sort_values(ascending=False).head(10)
    for feature, pct in top_complete.items():
        print(f"  {feature}: {pct:.1f}% complete")
    
    print(f"\nProcessing completed successfully!")
    print(f"Main output file: {output_file}")
    print(f"Protein biomarkers file: {'protein_biomarkers_processed.csv' if olink_columns else 'N/A'}")

if __name__ == "__main__":
    main()
