This post is a rough guide to implementing real-time adjustments to CHLA values in Argo NetCDF files. The real-time adjustments are relatively new and reflect a balance between simplicity and completeness, with the understanding that delayed-mode QC will be more robust and incorporate more float-specific considerations. In case it's not clear, nothing in this post should be construed as a definitive QC implementation and should be subject to considerable testing prior to deployment as actual live QC.

To implement the tests on existing ArgoNetCDF files in Python, we'll use the [argopandas](https://github.com/ArgoCanada/argopandas) package (to list float profiles and do some basic interaction with Argo NetCDF files) and [gsw](https://teos-10.github.io/GSW-Python/) (to perform seawater calculations required for some of the checks). We'll also use [pandas](https://pandas.pydata.org/) data frames and [numpy](https://numpy.org/) arrays to marshal the generic data handling and [matplotlib](https://matplotlib.org/) to plot as we go along.


```python
# pip install argopandas gsw matplotlib
# conda install -c conda-forge argopandas gsw matplotlib
import argopandas as argo
import gsw
import numpy as np
import matplotlib.pyplot as plt
```

### Example data

We'll need some data to practice with, too. As an example, I'll use all the profiles from float 6904117. It's probably best to use a fresh cache for every Python session (the default), but for the purposes of rendering this post I'll use a local cache to avoid downloading the files everytime I render it.


```python
profiles = argo.prof \
    .subset_float(6904117) \
    .subset_direction('ascending') 
bgc_profiles = argo.bio_prof \
    .subset_float(6904117) \
    .subset_direction('ascending')
```

    Downloading 'https://data-argo.ifremer.fr/ar_index_global_prof.txt.gz'
    Downloading 'https://data-argo.ifremer.fr/argo_bio-profile_index.txt.gz'
    


```python
profiles
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>file</th>
      <th>date</th>
      <th>latitude</th>
      <th>longitude</th>
      <th>ocean</th>
      <th>profiler_type</th>
      <th>institution</th>
      <th>date_update</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1826644</th>
      <td>coriolis/6904117/profiles/R6904117_001.nc</td>
      <td>2021-02-03 12:42:43+00:00</td>
      <td>57.295</td>
      <td>19.963</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-09-14 15:09:29+00:00</td>
    </tr>
    <tr>
      <th>1826646</th>
      <td>coriolis/6904117/profiles/R6904117_002.nc</td>
      <td>2021-02-05 00:55:14+00:00</td>
      <td>57.251</td>
      <td>19.854</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-09-14 15:09:42+00:00</td>
    </tr>
    <tr>
      <th>1826648</th>
      <td>coriolis/6904117/profiles/R6904117_003.nc</td>
      <td>2021-02-06 12:42:44+00:00</td>
      <td>57.165</td>
      <td>19.849</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-09-14 15:09:52+00:00</td>
    </tr>
    <tr>
      <th>1826650</th>
      <td>coriolis/6904117/profiles/R6904117_004.nc</td>
      <td>2021-02-08 00:59:14+00:00</td>
      <td>57.136</td>
      <td>19.973</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-09-14 15:10:04+00:00</td>
    </tr>
    <tr>
      <th>1826652</th>
      <td>coriolis/6904117/profiles/R6904117_005.nc</td>
      <td>2021-02-09 12:44:44+00:00</td>
      <td>57.163</td>
      <td>20.046</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-09-14 15:10:15+00:00</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1826978</th>
      <td>coriolis/6904117/profiles/R6904117_168.nc</td>
      <td>2021-10-12 12:26:21+00:00</td>
      <td>57.789</td>
      <td>19.844</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-10-12 16:43:24+00:00</td>
    </tr>
    <tr>
      <th>1826980</th>
      <td>coriolis/6904117/profiles/R6904117_169.nc</td>
      <td>2021-10-14 00:28:17+00:00</td>
      <td>57.828</td>
      <td>19.870</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-10-14 04:40:19+00:00</td>
    </tr>
    <tr>
      <th>1826982</th>
      <td>coriolis/6904117/profiles/R6904117_170.nc</td>
      <td>2021-10-15 12:30:17+00:00</td>
      <td>57.816</td>
      <td>19.850</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-10-15 17:48:45+00:00</td>
    </tr>
    <tr>
      <th>1826984</th>
      <td>coriolis/6904117/profiles/R6904117_171.nc</td>
      <td>2021-10-17 00:29:12+00:00</td>
      <td>57.831</td>
      <td>19.861</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-10-17 04:40:26+00:00</td>
    </tr>
    <tr>
      <th>1826986</th>
      <td>coriolis/6904117/profiles/R6904117_172.nc</td>
      <td>2021-10-18 13:12:23+00:00</td>
      <td>57.825</td>
      <td>19.852</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>2021-10-18 16:39:57+00:00</td>
    </tr>
  </tbody>
</table>
<p>172 rows × 8 columns</p>
</div>




```python
bgc_profiles
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>file</th>
      <th>date</th>
      <th>latitude</th>
      <th>longitude</th>
      <th>ocean</th>
      <th>profiler_type</th>
      <th>institution</th>
      <th>parameters</th>
      <th>parameter_data_mode</th>
      <th>date_update</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>172309</th>
      <td>coriolis/6904117/profiles/BR6904117_001.nc</td>
      <td>2021-02-03 12:42:43+00:00</td>
      <td>57.295</td>
      <td>19.963</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-09-14 15:09:29+00:00</td>
    </tr>
    <tr>
      <th>172311</th>
      <td>coriolis/6904117/profiles/BR6904117_002.nc</td>
      <td>2021-02-05 00:55:14+00:00</td>
      <td>57.251</td>
      <td>19.854</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-09-14 15:09:42+00:00</td>
    </tr>
    <tr>
      <th>172313</th>
      <td>coriolis/6904117/profiles/BR6904117_003.nc</td>
      <td>2021-02-06 12:42:44+00:00</td>
      <td>57.165</td>
      <td>19.849</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-09-14 15:09:52+00:00</td>
    </tr>
    <tr>
      <th>172315</th>
      <td>coriolis/6904117/profiles/BR6904117_004.nc</td>
      <td>2021-02-08 00:59:14+00:00</td>
      <td>57.136</td>
      <td>19.973</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-09-14 15:10:04+00:00</td>
    </tr>
    <tr>
      <th>172317</th>
      <td>coriolis/6904117/profiles/BR6904117_005.nc</td>
      <td>2021-02-09 12:44:44+00:00</td>
      <td>57.163</td>
      <td>20.046</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-09-14 15:10:15+00:00</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>172643</th>
      <td>coriolis/6904117/profiles/BR6904117_168.nc</td>
      <td>2021-10-12 12:26:21+00:00</td>
      <td>57.789</td>
      <td>19.844</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-10-12 16:43:24+00:00</td>
    </tr>
    <tr>
      <th>172645</th>
      <td>coriolis/6904117/profiles/BR6904117_169.nc</td>
      <td>2021-10-14 00:28:17+00:00</td>
      <td>57.828</td>
      <td>19.870</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-10-14 04:40:19+00:00</td>
    </tr>
    <tr>
      <th>172647</th>
      <td>coriolis/6904117/profiles/BR6904117_170.nc</td>
      <td>2021-10-15 12:30:17+00:00</td>
      <td>57.816</td>
      <td>19.850</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-10-15 17:48:45+00:00</td>
    </tr>
    <tr>
      <th>172649</th>
      <td>coriolis/6904117/profiles/BR6904117_171.nc</td>
      <td>2021-10-17 00:29:12+00:00</td>
      <td>57.831</td>
      <td>19.861</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-10-17 04:40:26+00:00</td>
    </tr>
    <tr>
      <th>172651</th>
      <td>coriolis/6904117/profiles/BR6904117_172.nc</td>
      <td>2021-10-18 13:12:23+00:00</td>
      <td>57.825</td>
      <td>19.852</td>
      <td>A</td>
      <td>834</td>
      <td>IF</td>
      <td>MTIME PRES C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY...</td>
      <td>RRRRRARRRRRRRRRRRARRRRRRRRAR</td>
      <td>2021-10-18 16:39:57+00:00</td>
    </tr>
  </tbody>
</table>
<p>172 rows × 10 columns</p>
</div>



### The recipe

As I understand it, this is what has to happen to the real-time measurements during the CHLA processing:

- `CHLA` values whose flags are unset are flagged as `Flag.PROBABLY_BAD`
- `CHLA_ADJUSTED` values whose flags are unset are flagged as `Flag.GOOD`
- Apply global range test (flag for `CHLA` and `CHLA_ADJUSTED` set to `Flag.BAD` for CHLA values outside the range -0.1 to 100)
- Apply the "dark correction" to `CHLA` values:
    - Read the `SCIENTIFIC_CALIB_COEFFICIENT` variable for the CHLA profile and parse the `FLOAT_DARK`, and `FLOAT_DARK_QC` values. It is likely that `FLOAT_DARK_QC` doesn't exist in any files yet, so default to `Flag.NO_QC`. For profiles that haven't gone through 5 cycles deeper than the "mixed layer depth", there should be a `PRELIM_DARK` variable in `SCIENTIFIC_CALIB_COEFFICIENT`
    - If `FLOAT_DARK` exists, apply the equation `CHLA_ADJUSTED = ((FLUORESCENCE_CHLA-FLOAT_DARK)*SCALE_CHLA)/2`.
    - If it doesn't, calculate try to calculate `PRELIM_DARK` and look for 4 previous instances of `PRELIM_DARK` with which `FLOAT_DARK/_QC` should be calculated. If there are 5 `PRELIM_DARK` values, calculate `FLOAT_DARK/_QC` and use that. If there aren't calculate the adjusted value using `PRELIM_DARK`.
- Apply the non-photochemical quenching correction (NPC) to the value calculated by the dark correction. This is the correction for the upper portion of the ocean whereby phytoplankton have a reduced response to UV light because they have already been exposed to it from the sun.

A few of these steps have details that need to be expanded upon, so let's do those first

### QC Flags

In addition to computing adjusted realtime values, the QC checks for CHLA also flag some values using the Argo flag scheme. One of the challeneges in doing this is that you don't want to make a data point look "better" by assigning a QC flag (e.g., flag a value as "probably bad" when it's already been flagged as "bad"). This logic was implemented by the experimental [argortqcpy](https://github.com/euroargodev/argortqcpy) package ([here](https://github.com/euroargodev/argortqcpy/blame/main/argortqcpy/checks.py)) and I've modified it slightly to get all of the flag logic in one happy Python class:


```python
class Flag:
    """
    Flags for check output. These values are valid values of the
    ``qc`` and ``adjusted_qc`` attributes of a
    :class:`~medsrtqc.core.Trace` object. Utility functions are
    provided as static methods to get the name or value of a flag
    or to update flag values ensuring that values that are already
    marked at a "worse" QC level are not inadvertently changed.
    """

    @staticmethod
    def label(flag):
        """Return the label of a QC flag"""
        return Flag._names[flag]

    @staticmethod
    def value(label):
        """Return the value of a QC flag"""
        for value, lab in Flag._names.items():
            if label == lab:
                return value
        raise KeyError(f"'{label}' is not the name of a QC flag")

    @staticmethod
    def update_safely(qc, to, where=None):
        """
        Safely update ``qc`` to the value ``to``. Values that are
        already marked at a "worse" QC level are not modified.
        """
        where = slice(None) if where is None else where
        flags = qc[where]
        for overridable_flag in Flag._precedence[to]:
            flags[flags == overridable_flag] = to
        qc[where] = flags

    NO_QC = b'0'
    GOOD = b'1'
    PROBABLY_GOOD = b'2'
    PROBABLY_BAD = b'3'
    BAD = b'4'
    CHANGED = b'5'
    # '6' not used
    # '7' not used
    ESTIMATED = b'8'
    MISSING = b'9'
    FILL_VALUE = b''

    _names = {
        NO_QC: 'NO_QC',
        GOOD: 'GOOD',
        PROBABLY_GOOD: 'PROBABLY_GOOD',
        PROBABLY_BAD: 'PROBABLY_BAD',
        BAD: 'BAD',
        CHANGED: 'CHANGED',
        ESTIMATED: 'ESTIMATED',
        MISSING: 'MISSING',
        FILL_VALUE: 'FILL_VALUE'
    }

    _precedence = {
        NO_QC: set(),
        GOOD: {
            NO_QC,
        },
        PROBABLY_GOOD: {
            NO_QC,
            GOOD,
            CHANGED,
        },
        PROBABLY_BAD: {
            NO_QC,
            GOOD,
            PROBABLY_GOOD,
            CHANGED,
        },
        BAD: {
            NO_QC,
            GOOD,
            PROBABLY_GOOD,
            CHANGED,
            PROBABLY_BAD,
        },
        CHANGED: {
            NO_QC,
        },
        ESTIMATED: {
            NO_QC,
            GOOD,
            PROBABLY_GOOD,
        },
        MISSING: {
            NO_QC,
        },
    }
```

### Parsing SCIENTIFIC_CALIB_COEFFICIENT

This variable in the NetCDF files looks like this:

```
DARK_CHLA = 47, SCALE_CHLA = 0.0073
```

and is in string form. It can be parsed by splitting a few times to return a `dict()` and reconstituted by pasting together a few times.


```python
def parse_scientific_calib_coefficient(val):
    parts = val.split(',')
    vals =  [part.split('=') for part in parts]
    return {k.strip(): float(v.strip()) for k, v in vals}

def unparse_calib_coefficient(coefs):
    return ', '.join(f'{k} = {v}' for k, v in coefs.items())

coefs = parse_scientific_calib_coefficient('DARK_CHLA = 47, SCALE_CHLA = 0.0073')
print(coefs)
coefs['SOMETHING_ELSE'] = 1.23
unparse_calib_coefficient(coefs)
```

    {'DARK_CHLA': 47.0, 'SCALE_CHLA': 0.0073}
    




    'DARK_CHLA = 47.0, SCALE_CHLA = 0.0073, SOMETHING_ELSE = 1.23'



### The mixed layer depth

As I understand it, there are as many ways to calculate the mixed layer depth as there are oceanographers. The following is an implementation based on the guidance provided in the draft CHLA QC document. I'm using `plt` as an argument here so that you can debug this stuff interactively (but skip installing matplotlib in production). The gist is that it looks for density changes > 0.03 below 10 dbar.


```python
import gsw

def calc_mixed_layer_depth(pres, temp, psal, longitude=0, latitude=0, plt=None):
    abs_salinity = gsw.SA_from_SP(psal, pres, longitude, latitude)
    conservative_temp = gsw.CT_from_t(abs_salinity, temp, pres)
    density = gsw.sigma0(abs_salinity, conservative_temp)

    if plt:
        plt.plot(density, pres)
        plt.gca().invert_yaxis()
        plt.gca().set_xlabel('sigma0')

    mixed_layer_start = (np.diff(density) > 0.03) & (pres[:-1] > 10)
    if not np.any(mixed_layer_start):
        # Can't determine mixed layer depth (no density changes > 0.03 below 10 dbar)
        return None

    mixed_layer_start_index = np.where(mixed_layer_start)[0][0]
    mixed_layer_depth = pres[mixed_layer_start_index]

    if plt:
        plt.gca().axhline(y = mixed_layer_depth, linestyle='--')

    return mixed_layer_depth
```

It's a little easier to parameterize this in terms of a `NetCDFWrapper` object since that's what we get handed.


```python
def calc_mixed_layer_depth_nc(nc_core, plt=None):
    if 'PRES' not in nc_core.variables or \
        'TEMP' not in nc_core.variables or \
        'PSAL' not in nc_core.variables:
        return None
    
    pres = nc_core['PRES'][:][0]
    temp = nc_core['TEMP'][:][0]
    psal = nc_core['PSAL'][:][0]

    if 'LONGITUDE' in nc_core.variables and 'LATITUDE' in nc_core.variables:
        longitude = nc_core['LONGITUDE'][:][0]
        latitude = nc_core['LATITUDE'][:][0]
    else:
        longitude = 0
        latitude = 0

    return calc_mixed_layer_depth(
        pres, temp, psal, 
        longitude=longitude,
        latitude=latitude, 
        plt=plt)

calc_mixed_layer_depth_nc(
    argo.nc('dac/' + profiles.reset_index().file[0]), 
    plt=plt)
```

    Downloading 'https://data-argo.ifremer.fr/dac/coriolis/6904117/profiles/R6904117_001.nc'
    




    60.5




    
![png](post_files/post_14_2.png)
    


### Calculating FLOAT_DARK

The guidance on calculating `PRELIM_DARK` is that it's the minimum `CHLA` value below the mixed layer depth boundary. In Python:


```python
def calc_prelim_dark(pres, chla, mixed_layer_depth):
    if mixed_layer_depth is None:
        return None
    
    chla_filter = chla[(pres > mixed_layer_depth) & ~np.isnan(chla)]
    if len(chla_filter) > 0:
        return chla_filter.min()
    else:
        return None
```

Because no `PRELIM_DARK` values have been calculated since the spec is new, we need something to do this too. Again, it's a little easier below if we parameterize this in terms of a `NetCDFWrapper`.


```python
def calc_prelim_dark_nc(nc_bgc, nc_core):
    if 'CHLA' not in nc_bgc.variables or 'PRES' not in nc_bgc.variables:
        raise KeyError(f"'CHLA' or 'PRES' not found")
    
    # we need a mixed layer depth for this calculation
    mixed_layer_depth = calc_mixed_layer_depth_nc(nc_core)
    if mixed_layer_depth is None:
        return None

    chla_prof_i = nc_bgc.param[nc_bgc.param.STATION_PARAMETERS.str.strip() == 'CHLA'] \
        .reset_index() \
        .iloc[0] \
        .N_PROF

    pres = nc_bgc['PRES'][:][chla_prof_i]
    chla = nc_bgc['CHLA'][:][chla_prof_i]
    
    return calc_prelim_dark(pres, chla, mixed_layer_depth)

calc_prelim_dark_nc(
    argo.nc('dac/' + bgc_profiles.reset_index().file[0]),
    argo.nc('dac/' + profiles.reset_index().file[0]))
```

    Downloading 'https://data-argo.ifremer.fr/dac/coriolis/6904117/profiles/BR6904117_001.nc'
    




    0.2263



Calculating `FLOAT_DARK` is a little more involved and it's where we'll make use of our float/cycle index that we created earlier. We need to start from the first NetCDF file for the float and collect any `PRELIM_DARK` values from `SCIENTIFIC_CALIB_COEFFICIENT` that might have been calculated. If there were 3 or fewer, we do nothing. If there are exactly 4, we add our `PRELIM_DARK` that we just calculated, calculate the mean of the values, and use that. The `FLOAT_DARK_QC` gets added based on the standard deviation of the `PRELIM_DARK` values.


```python
def accumulate_prelim_dark(bio_prof_files):
    prelim_dark = []
    for nc in argo.nc(bio_prof_files):
        coefs_df = nc.calib[nc.calib.PARAMETER.str.strip() == 'CHLA'] \
            .filter(['SCIENTIFIC_CALIB_COEFFICIENT'])
        if len(coefs_df) == 0 or 'SCIENTIFIC_CALIB_COEFFICIENT' not in coefs_df:
            continue

        coefs = parse_scientific_calib_coefficient(coefs_df.iloc[0][0])
        if 'PRELIM_DARK' not in coefs:
            continue
        
        prelim_dark.append(coefs['PRELIM_DARK'])
    
    return np.array(prelim_dark)

accumulate_prelim_dark(['dac/' + f for f in bgc_profiles.head(10).file])
```

    Downloading 9 files from 'https://data-argo.ifremer.fr/dac/coriolis/6904117/profiles'
    




    array([], dtype=float64)



An empty result here is expected since the scheme was just invented.

### The whole game

The rest is not very oceanographic but does involve navigating the structure of BGC variables' encoding in Argo profile NetCDFs. First, we get a `NetCDFWrapper` handle and a local filename to the .nc file. In argopandas, a `NetCDFWrapper` is a thin wrapper around a `netCDF4.Dataset` that implements a few common accessors in addition to the most common accessors for the dataset (notably, `obj.variables` and `obj['variable name']`).


```python
nc_filename = argo.filename('dac/coriolis/6904117/profiles/BR6904117_171.nc')
nc_core_filename = argo.filename('dac/coriolis/6904117/profiles/R6904117_171.nc')
nc = argo.nc(nc_filename)
nc_core = argo.nc(nc_core_filename)

# make sure we've got a CHLA variable
if 'CHLA' not in nc.variables or 'PRES' not in nc.variables:
    raise KeyError(f"'CHLA' or 'PRES' not found in '{nc_filename}'")

# find the profile index associated with CHLA
chla_prof_i = nc.param[nc.param.STATION_PARAMETERS.str.strip() == 'CHLA'] \
    .reset_index() \
    .iloc[0] \
    .N_PROF

# get the PRES/CHLA/CHLA_QC series we'll work with
# (note that there will be trailing NaN values here
# but we want to keep those because we'll have to keep the
# size the same to reassign values to a copy later)
pres = nc['PRES'][:][chla_prof_i]
chla = nc['CHLA'][:][chla_prof_i]
chla_qc_original = nc['CHLA_QC'][:][chla_prof_i]
chla_qc = chla_qc_original.copy()  # keep original so we can compare!

# create the chla_adjusted and chla_qc variables from originals
chla_adjusted = chla.copy()
chla_adjusted_qc = chla_qc.copy()

# reset chla_qc to Flag.NO_QC
chla_qc[:] = Flag.NO_QC

# plot to verify!
plt.plot(chla, pres)
plt.gca().invert_yaxis()
```

    Downloading 'https://data-argo.ifremer.fr/dac/coriolis/6904117/profiles/BR6904117_171.nc'
    Downloading 'https://data-argo.ifremer.fr/dac/coriolis/6904117/profiles/R6904117_171.nc'
    


    
![png](post_files/post_23_1.png)
    


The first step is to set the initial QC values to `Flag.PROBABLY_BAD`.


```python
chla_qc[:] = Flag.PROBABLY_BAD
chla_adjusted_qc[:] = Flag.GOOD
```

Then we apply the global range test:


```python
Flag.update_safely(chla_qc, Flag.BAD, where=(chla < 0.1) | (chla > 100))
Flag.update_safely(chla_adjusted_qc, Flag.BAD, where=(chla < 0.1) | (chla > 100))
```

Next, we go through the steps for the dark correction. This profile was collected before any of this QC was implemented so it doesn't contain any `PRELIM_DARK` or `FLOAT_DARK`. I'll write out the logic anyway.


```python
coefs_df = nc.calib[nc.calib.PARAMETER.str.strip() == 'CHLA'] \
    .filter(['SCIENTIFIC_CALIB_COEFFICIENT'])
if len(coefs_df) == 0 or 'SCIENTIFIC_CALIB_COEFFICIENT' not in coefs_df:
    raise ValueError(f"Can't find 'SCIENTIFIC_CALIB_COEFFICIENT' for 'CHLA' in file '{nc_filename}'")

# keep original and modified coefs so we know if they need updating
coefs = parse_scientific_calib_coefficient(coefs_df.iloc[0][0])
coefs_mod = coefs.copy()

if 'FLOAT_DARK' in coefs:
    float_dark_qc = coefs['FLOAT_DARK_QC'] if 'FLOAT_DARK_QC' in coefs else Flag.NO_QC
    chla_adjusted[:] = chla - coefs['FLOAT_DARK']
    Flag.update_safely(chla_adjusted_qc, float_dark_qc)
else:
    prelim_dark_acc = accumulate_prelim_dark(['dac/' + f for f in bgc_profiles.file])
    prelim_dark_this = calc_prelim_dark_nc(nc, nc_core)
    if prelim_dark_this is None:
        prelim_dark_this = np.array([])
    else:
        coefs_mod['PRELIM_DARK'] = prelim_dark_this
        prelim_dark_this = np.array([prelim_dark_this])
    
    prelim_dark_acc = np.concatenate([prelim_dark_acc, prelim_dark_this])

    if len(prelim_dark_acc) > 0:
        float_dark = prelim_dark_acc.mean()
        float_dark_qc = Flag.PROBABLY_GOOD  # not sure how this is actually calculated

        chla_adjusted[:] = chla - float_dark
        Flag.update_safely(chla_adjusted_qc, float_dark_qc)

        # potentially store in coefs_mod if there are enough PRELIM_DARK values
        if len(prelim_dark_acc) == 5:
            coefs_mod['FLOAT_DARK'] = float_dark
            coefs_mod['FLOAT_DARK_QC'] = float_dark_qc
    else:
        # we can't calculate the adjusted value
        chla_adjusted[:] = np.nan

# should show the adjustment
print(chla[::100])
print(chla_adjusted[::100])
```

    Downloading 161 files from 'https://data-argo.ifremer.fr/dac/coriolis/6904117/profiles'
    

    [1.1241999864578247 1.1461000442504883 1.292099952697754
     1.0292999744415283 0.23360000550746918 0.20440000295639038
     0.21170000731945038 0.21170000731945038 0.23360000550746918
     0.2773999869823456 0.2847000062465668 0.2847000062465668
     0.2919999957084656]
    [0.9343999624252319 0.9563000202178955 1.1022999286651611
     0.8394999504089355 0.04380001127719879 0.014600008726119995
     0.021900013089179993 0.021900013089179993 0.04380001127719879
     0.0875999927520752 0.09490001201629639 0.09490001201629639
     0.10220000147819519]
    

There are no previous `PRELIM_DARK` values and no `FLOAT_DARK` value, so this applies the calculation based on the current profile.

From my reading of the document, that's the operation! What's left is to figure out what changed and modify the NetCDF. There is probably an object diff library in Python that could do a good job of reporting what changed and I'd encourage readers to implement it! The objects that might have changed are `coefs_mod` (check differences with `coefs`) and `chla_qc` (check differences via `chla_qc_original`). `chla_adjusted` and `chla_adjusted_qc` can be written without checking diffs since in theory a brand new file won't have those variables anyway.

I imagine modification would be done by making a copy of the NetCDF file and opening via `netCDF4.Dataset()` in write mode.


```python
if coefs_mod != coefs:
    coefs_serialized = unparse_calib_coefficient(coefs_mod)
    # write coefs_searialized to the appropriate location in 
    # 'SCIENTIFIC_CALIB_COEFFICIENT'

chla_qc_updated = chla_qc != chla_qc_original
```

### Tests

In the course of writing all of this I also wrote some unit tests and it's worth putting them here so they don't get lost!


```python
import unittest

class TestFlag(unittest.TestCase):

    def test_value(self):
        self.assertEqual(Flag.NO_QC, Flag('NO_QC'))
        self.assertEqual(Flag.label(Flag.NO_QC), 'NO_QC')
        with self.assertRaises(KeyError):
            Flag('not a QC key')
        with self.assertRaises(KeyError):
            Flag.label(b'def not a flag')

    def test_update(self):
        qc = np.array([Flag.GOOD, Flag.PROBABLY_BAD, Flag.MISSING])
        Flag.update_safely(qc, to=Flag.BAD)
        self.assertTrue(np.all(qc == np.array([Flag.BAD, Flag.BAD, Flag.MISSING])))

        qc = np.array([Flag.GOOD, Flag.PROBABLY_BAD, Flag.MISSING])
        Flag.update_safely(qc, to=Flag.BAD, where=np.array([False, True, False]))
        self.assertTrue(np.all(qc == np.array([Flag.GOOD, Flag.BAD, Flag.MISSING])))

```
