
## TESTING ENVIRONMENT FOR PYTHON CODE TO BE EMBEDDED IN RMARKDOWN

import pandas as pd

df = pd.read_csv('deployment/canada_deployments.csv')
df['DEPLOYMENT DATE'] = df['DEPLOYMENT DATE'].apply(pd.Timestamp, args=(None, 'utc'))
df = df.sort_values('DEPLOYMENT DATE')
df['DEPLOYMENT DATE'] = [d.strftime('%d %b, %Y') for d in df['DEPLOYMENT DATE']]
df['WMO'] = [str(int(w)) if pd.notna(w) else '' for w in df.WMO]
df['IMEI'] = [str(int(i)) if pd.notna(i) else '' for i in df.IMEI]
df['SERIAL NUMBER'] = [s if pd.notna(s) else '' for s in df['SERIAL NUMBER']]

mapper = {
    	'PROGRAM':'Program',
        'STATUS':'Status',
        'MODEL':'Model',
        'DEPLOYMENT DATE':'Date',
        'DEPLOYMENT LAT':'Latitude',
        'DEPLOYMENT LON':'Longitude',
        'DEPLOYMENT SHIP':'Ship',
        'SERIAL NUMBER':'Serial No.',
}

df = df.rename(columns=mapper)
df = df.reset_index()