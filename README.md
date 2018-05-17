
# Pymaceuticals Observed Trends
* Trend 1
* Trend 2
* Trend 3


```python
# Import dependencies
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
```


```python
# Read in the clinical trials data into data frames
csv_path = 'data/clinicaltrial_data.csv'

ct_df = pd.read_csv(csv_path)

ct_df.head()
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
      <th>Mouse ID</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>b128</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>f932</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>g107</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>a457</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>c819</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Read in the clinical trials data into data frames
csv_path = 'data/mouse_drug_data.csv'

mouse_df = pd.read_csv(csv_path)

mouse_df.head()
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
      <th>Mouse ID</th>
      <th>Drug</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>f234</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>1</th>
      <td>x402</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>2</th>
      <td>a492</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>3</th>
      <td>w540</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>4</th>
      <td>v764</td>
      <td>Stelasyn</td>
    </tr>
  </tbody>
</table>
</div>




```python
# In order to analyze the data by treatment, we need to merge the data
merged_data = pd.merge(ct_df, mouse_df, on="Mouse ID", how="inner")
```


```python
# Then we select the columns we care about for the next scatter plot
vol_vs_time = merged_data.groupby(["Drug","Timepoint"])[["Tumor Volume (mm3)"]]

vol_vs_time = vol_vs_time.mean()

scattered_df = vol_vs_time.reset_index().pivot(index="Timepoint", columns="Drug", values = "Tumor Volume (mm3)")[["Capomulin","Ketapril","Placebo","Infubinol"]]

scattered_df.head()
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
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Ketapril</th>
      <th>Placebo</th>
      <th>Infubinol</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>45.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>44.266086</td>
      <td>47.389175</td>
      <td>47.125589</td>
      <td>47.062001</td>
    </tr>
    <tr>
      <th>10</th>
      <td>43.084291</td>
      <td>49.582269</td>
      <td>49.423329</td>
      <td>49.403909</td>
    </tr>
    <tr>
      <th>15</th>
      <td>42.064317</td>
      <td>52.399974</td>
      <td>51.359742</td>
      <td>51.296397</td>
    </tr>
    <tr>
      <th>20</th>
      <td>40.716325</td>
      <td>54.920935</td>
      <td>54.364417</td>
      <td>53.197691</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Compute error bar by looping through each timepoint and computing the sem at each value, for each treatment
drugs = ["Capomulin", "Ketapril", "Placebo", "Infubinol"]

# We will use these dictionaries to reference the relative means and sems
sems = {}
means = {}

for drug in drugs: 
    subset_df = merged_data.loc[merged_data["Drug"]==drug, ["Timepoint", "Tumor Volume (mm3)"]]
    sems[drug] = [stats.sem(subset_df.loc[subset_df["Timepoint"] == x, ["Tumor Volume (mm3)"]]) for x in np.arange(0,50,5)]
    means[drug] = scattered_df[drug]
```

## Tumor Response to Treatment


```python
# This scatter plot shows how the tumor volume changes over time for each treatment.
x_axis = np.arange(0,50,5)

plt.errorbar(np.arange(0, 50, 5), means["Capomulin"], yerr = sems["Capomulin"], fmt = '*--', capsize = 2.5)
plt.errorbar(np.arange(0, 50, 5), means["Ketapril"], yerr = sems["Ketapril"], fmt = '^--', capsize = 2.5)
plt.errorbar(np.arange(0, 50, 5), means["Placebo"], yerr = sems["Placebo"], fmt = 's--', capsize = 2.5)
plt.errorbar(np.arange(0, 50, 5), means["Infubinol"], yerr = sems["Infubinol"], fmt = 'o--', capsize = 2.5)

# Add legend
plt.legend(loc="best")

# Add gridlines
plt.grid(alpha = 0.5)

# Add labels
plt.title('Tumor Response to Treatment')
plt.xlabel('Time (days)')
plt.ylabel('Tumor Volume (mm3)')

# Plot the graph
plt.show()
```


![png](output_8_0.png)



```python
scattered_df.keys
```




    <bound method NDFrame.keys of Drug       Capomulin   Ketapril    Placebo  Infubinol
    Timepoint                                            
    0          45.000000  45.000000  45.000000  45.000000
    5          44.266086  47.389175  47.125589  47.062001
    10         43.084291  49.582269  49.423329  49.403909
    15         42.064317  52.399974  51.359742  51.296397
    20         40.716325  54.920935  54.364417  53.197691
    25         39.939528  57.678982  57.482574  55.715252
    30         38.769339  60.994507  59.809063  58.299397
    35         37.816839  63.371686  62.420615  60.742461
    40         36.958001  66.068580  65.052675  63.162824
    45         36.236114  70.662958  68.084082  65.755562>


