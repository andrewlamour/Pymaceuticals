#!/usr/bin/env python
# coding: utf-8

# Analysis:
# 1. Capomulin does a good job reducing the overall size of the tumor in each mouse. 
# 2. The last graph shows there seems to be a correlation between the average tumor volume and overall weight of the mouse
# 3. Ketapril seems to be the least effective drug as the average tumor size amongst the mice on that drug is the highest according this sample size. 

# In[2]:


import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st


# In[3]:


# Read CSV
mouse_metadata_path = "Data/mousedata.csv"
study_results_path= "Data/studyresults.csv"


# In[4]:


# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

study_data = pd.merge(study_results, mouse_metadata, how="left", on="Mouse ID")
study_data.head(5)


# In[5]:


means = study_data.groupby('Drug Regimen').mean()['Tumor Volume (mm3)']
medians = study_data.groupby('Drug Regimen').median()['Tumor Volume (mm3)']
variances = study_data.groupby('Drug Regimen').var()['Tumor Volume (mm3)']
sds = study_data.groupby('Drug Regimen').std()['Tumor Volume (mm3)']
sems = study_data.groupby('Drug Regimen').sem()['Tumor Volume (mm3)']
summary_table = pd.DataFrame({"Mean Tumor Volume":means,
                              "Median Tumor Volume":medians,
                              "Tumor Volume Variance":variances,
                              "Tumor Volume Std. Dev.":sds,
                              "Tumor Volume Std. Err.":sems})
summary_table


# In[6]:


#bar charts using pandas
counts = study_data['Drug Regimen'].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Data Points")
plt.show()


# In[7]:


#bar chart using Pyplot
counts = study_data['Drug Regimen'].value_counts()
plt.bar(counts.index.values,counts.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Data Points")
plt.show()


# In[8]:


#pie chart using pandas
counts = study_data.Sex.value_counts()
counts.plot(kind="pie",autopct='%1.1f%%')
plt.show()


# In[9]:


#pie chart using pyplot
counts = study_data.Sex.value_counts()
plt.pie(counts.values,labels=counts.index.values,autopct='%1.1f%%')
plt.ylabel("Sex")
plt.show()


# In[12]:


#Tumor volume for each mouse with each treatment
max_tumor = study_data.groupby(["Mouse ID"])['Timepoint'].max()
max_tumor = max_tumor.reset_index()


merged_data = max_tumor.merge(study_data,on=['Mouse ID','Timepoint'],how="left")

treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]


tumor_vol_list = []


for drug in treatment_list:
    
 
    final_tumor_vol = merged_data.loc[merged_data["Drug Regimen"] == drug, 'Tumor Volume (mm3)']
    
 
    tumor_vol_list.append(final_tumor_vol)
    
    
    quartiles = final_tumor_vol.quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
    outliers = final_tumor_vol.loc[(final_tumor_vol < lower_bound) | (final_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers: {outliers}")


# In[13]:


#Box blot for tummer volume for each treatment
orange_out = dict(markerfacecolor='red',markersize=12)
plt.boxplot(tumor_vol_list, labels = treatment_list,flierprops=orange_out)
plt.ylabel('Final Tumor Volume (mm3)')
plt.show()


# In[19]:


#line plot for tumor volume over time with Capomulin
capomulin_table = study_data.loc[study_data['Drug Regimen'] == "Capomulin"]
mousedata = capomulin_table.loc[capomulin_table['Mouse ID']== 'b128']
plt.plot(mousedata['Timepoint'],mousedata['Tumor Volume (mm3)'])
plt.xlabel('Timepoint (days)')
plt.ylabel('Tumor Volume (mm3)')
plt.title('Capomulin treatment of mouse b128')
plt.show()


# In[21]:


#scatter plot of mouse weight vs avg tumor volume
capomulin_table = study_data.loc[study_data['Drug Regimen'] == "Capomulin"]
capomulin_average = capomulin_table.groupby(['Mouse ID']).mean()
plt.scatter(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()


# In[22]:


#Correlation Coefficient and linear regression modeal 
corr=round(st.pearsonr(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])[0],2)
print(f"The correlation between mouse weight and the average tumor volume is {corr}")
model = st.linregress(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])

y_values = capomulin_average['Weight (g)']*model[0]+model[1]
plt.scatter(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])
plt.plot(capomulin_average['Weight (g)'],y_values,color="red")
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()


# In[ ]:




