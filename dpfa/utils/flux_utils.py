import os
import logging
import pandas as pd
import cobra
import numpy as np


FLUX_ON = 1000
FLUX_OFF = 2000
FLUX_REVERSED = 3000


def process_models(slow_model_path: str, fast_model_path: str, tissue: str, output_dir: str):
    """
    Load models, solve optimization, save fluxes.
    Returns paths to saved CSV files with fluxes.
    """
    os.makedirs(output_dir, exist_ok=True)
    slow_model = cobra.io.load_json_model(slow_model_path)
    fast_model = cobra.io.load_json_model(fast_model_path)

    if tissue.lower() == 'liver':
        for m in [slow_model, fast_model]:
            try:
                m.reactions.get_by_id("EX_lac_L[e]").upper_bound = 0
            except Exception:
                pass

    slow_solution = slow_model.optimize()
    fast_solution = fast_model.optimize()

    slow_flux = slow_solution.fluxes.reset_index()
    fast_flux = fast_solution.fluxes.reset_index()
    slow_flux.columns = ['reaction_id', 'flux']
    fast_flux.columns = ['reaction_id', 'flux']
    slow_flux['flux'] = slow_flux['flux'].where(abs(slow_flux['flux']) >= 1e-5, 0)
    fast_flux['flux'] = fast_flux['flux'].where(abs(fast_flux['flux']) >= 1e-5, 0)

    slow_path = os.path.join(output_dir, f'slow_model_{tissue}.csv')
    fast_path = os.path.join(output_dir, f'fast_model_{tissue}.csv')
    slow_flux.to_csv(slow_path, index=False)
    fast_flux.to_csv(fast_path, index=False)
    logging.info(f"Fluxes saved to {slow_path} and {fast_path}")
    return slow_path, fast_path


def ratio_fluxes(flux_slow_path: str, flux_fast_path: str,
                 output_path: str) -> pd.DataFrame:
    """
    Save both flux columns and calculate ratio
    """
    flux_slow = pd.read_csv(flux_slow_path)
    flux_fast = pd.read_csv(flux_fast_path)

    merged = pd.merge(flux_slow, flux_fast, on='reaction_id',
                     how='outer', suffixes=('_slow', '_fast'))

    merged['flux_slow'] = merged['flux_slow'].fillna(0)
    merged['flux_fast'] = merged['flux_fast'].fillna(0)

    ratios = []
    for _, row in merged.iterrows():
        x, y = row['flux_slow'], row['flux_fast']

        if x == 0 and y != 0:
            ratios.append(FLUX_ON)
        elif x != 0 and y == 0:
            ratios.append(FLUX_OFF)
        elif x != 0 and y != 0 and ((x < 0 < y) or (y < 0 < x)):
            ratios.append(FLUX_REVERSED)
        elif x != 0 and y != 0:
            val = y / x
            ratios.append(val if abs(val) >= 1e-5 else 0)
        else:
            ratios.append(0)

    merged['flux_ratio'] = ratios
    merged = merged[(merged['flux_slow'] != 0) | (merged['flux_fast'] != 0)]

    merged.to_csv(output_path, index=False)
    logging.info(f"Flux data saved to {output_path}")

    return merged


def create_DRF(df: pd.DataFrame, threshold: float = 0.30) -> pd.DataFrame:
    """
    Categorize flux changes using absolute difference

    """
    slow = df['flux_slow'].abs()
    fast = df['flux_fast'].abs()

    on_mask = (slow == 0) & (fast > 0)
    off_mask = (slow > 0) & (fast == 0)
    reversed_mask = ((df['flux_slow'] < 0) & (df['flux_fast'] > 0)) | \
                    ((df['flux_slow'] > 0) & (df['flux_fast'] < 0))

    diff = fast - slow
    change_threshold = threshold * slow

    conditions = [
        off_mask,
        on_mask,
        reversed_mask,
        diff < -change_threshold,
        diff > change_threshold,
    ]

    choices = ['Off', 'On', 'Reversed', 'Decreased', 'Increased']
    df['DRF_category'] = np.select(conditions, choices, default='Unchanged')

    return df
