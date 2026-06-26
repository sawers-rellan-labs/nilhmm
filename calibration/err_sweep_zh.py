import sys, numpy as np, pandas as pd
from scipy.stats import ks_2samp
sys.path.insert(0,"agent/nilhmm")
from nilhmm.io import read_vcf_counts
from nilhmm.core import introgression_hmm_counts
import logging; logging.basicConfig(level=logging.ERROR)

sim=pd.read_csv("agent/bc2s2_segments.csv")["mb"].values
ref,alt,md,samp,mi=read_vcf_counts("data/nilhmm/Zh_counts.vcf.gz")
rb,ab,mdb,sb,mib=read_vcf_counts("data/nilhmm/B73_counts.vcf.gz")
chrom=mi["CHROM"].values; pos=mi["POS"].values.astype(float)

def blocks(c):
    out=[]
    for i in range(c.shape[0]):
        for ch in np.unique(chrom):
            m=chrom==ch; st=c[i,m]>0; p=pos[m]
            if not st.any(): continue
            mid=(p[:-1]+p[1:])/2; lo=np.r_[p[0],mid]; hi=np.r_[mid,p[-1]]
            d=np.diff(np.r_[0,st.astype(int),0]); s=np.where(d==1)[0]; e=np.where(d==-1)[0]-1
            out+=list((hi[e]-lo[s])/1e6)
    return np.array(out)

R=3e-5
print(f"{'err':>5} {'D':>7} {'HET%':>6} {'ALT%':>7} {'dosage':>7} {'HET:ALT':>8} {'B73max':>7}", flush=True)
for err in [0.01,0.02,0.05,0.1,0.2,0.3]:
    c=introgression_hmm_counts(ref,alt,md,err=err,conc=20,r=R,f_1=0.0625,f_2=0.0938)
    het=(c==1).mean()*100; al=(c==2).mean()*100; dos=(het+2*al)/200
    D=ks_2samp(blocks(c),sim).statistic
    cb=introgression_hmm_counts(rb,ab,mdb,err=err,conc=20,r=R,f_1=0.0625,f_2=0.0938)
    bdos=((cb==1).sum(1)+2*(cb==2).sum(1))/(2*cb.shape[1])
    ratio=het/al if al>0 else float('inf')
    print(f"{err:>5} {D:>7.4f} {het:>6.2f} {al:>7.3f} {dos:>7.4f} {ratio:>8.1f} {bdos.max():>7.4f}", flush=True)
