import matplotlib.pyplot as plt
from numpy.random import rand
import pandas as pd
import numpy as np

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
#read in csv file
infile = '/projects/trans_scratch/validations/workspace/szong/AML_capture/all_patients/itd/AML_SPC_FLT3_ITD_results.csv'
df = pd.read_csv(infile, delimiter=r"\t")

print df
# rows
print df[2:5]
# columns
print df[['supported_reads1(strand+)','supported_reads2(strand+)', 'length(ITD)']]
x = df[['supported_reads1(strand+)']]
y = df[['supported_reads2(strand+)']]
z = df[['length(ITD)']]
#color = ['red', 'green']
#legend = 'read2/read1'


ax1.scatter(x, y, alpha=1, edgecolors='none',label='support reads')
#ax1.scatter(x,y, label=legend, alpha=1, edgecolors='none')
#ax1.set_aspect(1./ax1.get_data_ratio()) # make axes square`
ax1.set_ylim([0,50])
ax1.set_xlim([0,50])
ax1.legend(loc=0, frameon=False)

ax2.scatter(x, z, color='red', marker= 'v', label='ITD length')
ax2.legend(loc=2, frameon=False)

ax1.set_xlabel('support read1')
ax1.set_ylabel('support read2')
ax2.set_ylabel('ITD length')
ax2.set_ylim([0,300])
ax1.grid(True)
fig.suptitle('ITD support reads and length', fontsize=20)
plt.tight_layout()
plt.savefig('aml_support_reads.pdf')
plt.show()
