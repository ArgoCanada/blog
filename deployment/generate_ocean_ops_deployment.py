
import pandas as pd

df = pd.read_csv('canada_deployments.csv')
serial_map = pd.read_csv('serial_info.csv').set_index('serial_prefix')
program_map = pd.read_csv('program_info.csv').set_index('institute')
df['DATE'] = df['DEPLOYMENT DATE'].apply(pd.Timestamp)

recent_deployment = pd.Timestamp('2026-07-01')

df = df.loc[(df.STATUS == 'OPERATIONAL') & (df.DATE > recent_deployment)]
df['SERIAL NUMBER'] = [sn if len(sn.split('-')) > 1 else f'TWR-{sn}' for sn in df['SERIAL NUMBER']]
df = df.drop(['STATUS', 'MODEL', 'MODEL_DETAIL', 'BASIN'], axis=1)
df['MODEL'] = [serial_map.loc[s.split('-')[0], 'model'] for s in df['SERIAL NUMBER']]
df['PROGRAM'] = [program_map.loc[inst, 'program'] for inst in df['INSTITUTE']]
df['SERIAL NUMBER'] = [sn if sn.split('-')[0] != 'TWR' else sn.split('-')[1] for sn in df['SERIAL NUMBER']]

df = df.drop('DATE', axis=1)
now = pd.Timestamp('now', tz='utc')
df.to_csv(f'OceanOps/canada_recent_deployment_{now.year}-{now.month:02d}-{now.day:02d}T{now.hour:02d}{now.minute:02d}{now.second:02d}.csv', index=False, float_format='%.16g')