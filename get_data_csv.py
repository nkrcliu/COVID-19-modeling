import csv
import os
import pandas as pd
import numpy as np


# 获取当前目录下的CSV文件名
filepaths = []
for info in os.listdir('E:\\references\\COVID-19-master\\COVID-19-master\\csse_covid_19_data\\csse_covid_19_daily_reports_us'):
    domain = os.path.abspath(r'E:\\references\\COVID-19-master\\COVID-19-master\\csse_covid_19_data\\csse_covid_19_daily_reports_us')
    #获取文件夹的路径
    filepath = os.path.join(domain,info) #将路径与文件名结合起来就是每个文件的完整路径
    filepaths.append(filepath)
print(filepaths)

needed = []
for path in filepaths:
    df = pd.read_csv(path,header=0,index_col=0)
    print(df['People_Hospitalized']['New York'])
    needed.append(df['People_Hospitalized']['New York'])
print(needed)
needed = pd.DataFrame(needed)
needed.to_csv('E:\\references\\COVID-19-master\\COVID-19-master\\csse_covid_19_data\\csse_covid_19_daily_reports_us'+'//'+'hospital.csv',index=False)


