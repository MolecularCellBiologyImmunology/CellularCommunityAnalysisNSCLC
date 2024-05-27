import sys
import pandas as pd
import numpy as np
from scipy.stats import weibull_min
from scipy.optimize import minimize
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

def fit_weibull(data):
    # Function to fit Weibull distribution and return shape and scale parameters
    def negative_log_likelihood(params, data):
        shape, scale = params
        return -np.sum(weibull_min.logpdf(data, shape, scale=scale))

    initial_params = [1, 1]
    result = minimize(negative_log_likelihood, initial_params, args=(data,), method='Nelder-Mead')
    if result.success:
        return result.x
    else:
        return [np.nan, np.nan]

def main(input_file, output_file):
    # Load 1-NN X/Y histogram coordinates dataframe (output from script 1__get1NNdistances.R)
    all_distances_data = pd.read_csv(input_file, sep='\t').dropna()
    all_distances_data['distance_window'] = all_distances_data['WinMean']
    all_distances_data['phenotype_combo'] = all_distances_data['phenotype_from'] + '_to_' + all_distances_data['phenotype_to']
    all_distances_data = all_distances_data[['Patient_ID', 'phenotype_combo', 'count_scaled', 'distance_window']]
    all_distances_data['new'] = all_distances_data['count_scaled'] * 1000
    print(all_distances_data[all_distances_data.isna().any(axis=1)])
    # print("Max and nans in 'new'",all_distances_data['new'].max(), all_distances_data['new'].isna().sum())
    all_distances_data['new'] = all_distances_data['new'].round().astype(int)

    # Recreate 1-NN histogram from coordinates (this is required for function "fitdistrplus")
    all_combos_dists = pd.DataFrame()
    for tn in all_distances_data['Patient_ID'].unique():
        for combo in all_distances_data['phenotype_combo'].unique():
            df_filtered = all_distances_data[(all_distances_data['Patient_ID'] == tn) & (all_distances_data['phenotype_combo'] == combo)]
            dists = [1] + df_filtered['distance_window'].tolist()
            times = [1] + df_filtered['new'].tolist()
            
            if len(times) == 2:
                dists = list(range(1, 299))
                times = [0] * 298
            
            all_dists_tnum = []
            for i in range(len(times)):
                all_dists_tnum.extend([dists[i]] * times[i])
            
            aldistsdf = pd.DataFrame({
                'dists': all_dists_tnum,
                'Patient_ID': [tn] * len(all_dists_tnum),
                'combo': [combo] * len(all_dists_tnum)
            })
            all_combos_dists = pd.concat([all_combos_dists, aldistsdf])

    # Fit a Weibull distribution by Maximum likelihood MLE [this is an initial estimation that will be optimized in step 3)]
    initial_params = pd.DataFrame(columns=['term', 'estimate', 'std.error', 'Patient_ID', 'combo'])

    for combo in all_distances_data['phenotype_combo'].unique():
        print(combo)
        for tn in all_combos_dists['Patient_ID'].unique():
            if 'PanCK+' in combo or 'negative' in combo or 'Cancer' in combo or 'Negative' in combo:
                data = all_combos_dists[(all_combos_dists['combo'] == combo) & (all_combos_dists['Patient_ID'] == tn)]['dists']
            else:
                data = all_combos_dists[(all_combos_dists['combo'] == combo) & (all_combos_dists['Patient_ID'] == tn) & (all_combos_dists['dists'] < 100)]['dists']
            
            if len(data) > 0:
                shape, scale = fit_weibull(data)
                params = pd.DataFrame({
                    'term': ['shape', 'scale'],
                    'estimate': [shape, scale],
                    'std.error': [np.nan, np.nan],  # Standard error is not calculated here
                    'Patient_ID': [tn] * 2,
                    'combo': [combo] * 2
                })
                initial_params = pd.concat([initial_params, params])

    # Save the estimated parameters to a file
    initial_params.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("You must parse an input file")
        sys.exit(1)
    elif len(sys.argv) == 2:
        print("Output file not specified. Using default: ./results/step2_weibull_initial_params.tsv")
        sys.argv.append("./results/step2_weibull_initial_params.tsv")

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
