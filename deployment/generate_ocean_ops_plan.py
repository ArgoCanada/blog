
import pandas as pd

df = pd.read_csv('canada_deployments.csv')
serial_map = pd.read_csv('nke_serial_info.csv').set_index('serial_prefix')

# NKE floats only, write APEX interpreter later
df = df.loc[(df.STATUS == 'CONFIRMED') & (df.MODEL != 'APEX')]
df = df.drop(['STATUS', 'MODEL', 'MODEL_DETAIL', 'BASIN'], axis=1)
df['MODEL'] = [serial_map.loc[s.split('-')[0], 'model'] for s in df['SERIAL NUMBER']]
df['SENSOR_MODEL'] = [serial_map.loc[s.split('-')[0], 'sensor_model'] for s in df['SERIAL NUMBER']]

now = pd.Timestamp('now', tz='utc')
df.to_csv(f'OceanOps/canada_deployment_plan_{now.year}-{now.month:02d}-{now.day:02d}T{now.hour:02d}{now.minute:02d}{now.second:02d}.csv', index=False, float_format='%.16g')