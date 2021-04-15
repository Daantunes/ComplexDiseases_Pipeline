import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

data=pd.read_csv('t2d.relatedness',sep='\t')
print(data)

ind=[]
for i in data['INDV1']:
	ind.append(i)
	
indv=list(dict.fromkeys(ind))
#indv=list(set(ind))

d = np.zeros(shape=(len(indv),len(indv)))
d = pd.DataFrame(d, columns=indv, index=indv)
print(d)

for i in range(len(data['INDV1'])):
	a=data.loc[i,'INDV1']
	b=data.loc[i,'INDV2']
	d.loc[a,b]=data.loc[i,'RELATEDNESS_AJK']
	d.loc[b,a]=data.loc[i,'RELATEDNESS_AJK']
	
print(d)
d=pd.DataFrame(MinMaxScaler().fit_transform(d), columns=d.columns, index=d.index)
print(d)

import plotly.graph_objects as go
x=list(d.columns)
y=list(d.index)
z=d.values


print(x)
print(z)
print(d.loc['Ex41','Ex51'])

fig = go.Figure(data=go.Heatmap(
        z=z,
        x=x,
        y=y,
        colorscale='Greys'))
fig.write_image("t2d_relatedness.png", width = 1500, height = 1500)
