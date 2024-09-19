from numpy import mean, median
from statistics import mode
import os

files=os.listdir('.')
files=[i for i in files if i.endswith('cov.tsv')]


for f in files:
	x=open(f).readlines()
	nuc_co={}
	nuc_ep={}
	nuc_li={}
	nuc_or={}
	mt={}
	pt={}
	for l in x:
		if l.startswith('mt'):
			gene=l.split()[0]
			try:mt[gene].append(l.split()[2])
			except KeyError:mt[gene]=[l.split()[2]]
		elif l.startswith('pt'):
			gene=l.split()[0]
			try:pt[gene].append(l.split()[2])
			except KeyError:pt[gene]=[l.split()[2]]
		else:
			gene=l.split()[0]
			if gene.endswith('americana'):
				gene=gene.split('_')[0]
				try:nuc_co[gene].append(l.split()[2])
				except KeyError:nuc_co[gene]=[l.split()[2]]
			elif gene.endswith('virginiana'):
				gene=gene.split('_')[0]
				try:nuc_ep[gene].append(l.split()[2])
				except KeyError:nuc_ep[gene]=[l.split()[2]]
			elif gene.endswith('philippensis'):
				gene=gene.split('_')[0]
				try:nuc_li[gene].append(l.split()[2])
				except KeyError:nuc_li[gene]=[l.split()[2]]
			elif gene.endswith('fasciculata'):
				gene=gene.split('_')[0]
				try:nuc_or[gene].append(l.split()[2])
				except KeyError:nuc_or[gene]=[l.split()[2]]			
	mt_site_cov=[]
	for k in mt.keys():mt_site_cov=mt_site_cov+mt[k][150:-150]
	mt_site_cov=[int(i) for i in mt_site_cov]
	#mt_site_cov_mean=mean(mt_site_cov)
	mt_site_cov_median=median(mt_site_cov)
	pt_site_cov=[]
	for k in pt.keys():pt_site_cov=pt_site_cov+pt[k][150:-150]
	pt_site_cov=[int(i) for i in pt_site_cov]
	#pt_site_cov_mean=mean(pt_site_cov)
	pt_site_cov_median=median(pt_site_cov)
	nuc_site_cov=[]
	for k in nuc_co.keys():
		nuc_site_cov=nuc_site_cov+nuc_co[k][150:-150]
		try:nuc_site_cov=nuc_site_cov+nuc_ep[k][150:-150]
		except KeyError:pass
		try:nuc_site_cov=nuc_site_cov+nuc_li[k][150:-150]
		except KeyError:pass
		try:nuc_site_cov=nuc_site_cov+nuc_or[k][150:-150]
		except KeyError:pass
	nuc_site_cov=[int(i) for i in nuc_site_cov if int(i)<50]
	nuc_site_cov_mode=mode(nuc_site_cov)
	print(f,str(mt_site_cov_median),str(pt_site_cov_median),str(nuc_site_cov_mode))


#output and plot in r
ref=open('mt_pt.average_cov.tsv')
mt_ref={}
pt_ref={}
for l in ref:
	mt_ref[l.split()[0]]=int(l.split()[1])
	pt_ref[l.split()[0]]=int(l.split()[2])


for f in files:
	x=open(f).readlines()
	f=f.split('.')[0]
	mt={}
	pt={}
	for l in x:
		if l.startswith('mt') or l.startswith('NODE'):
			gene=l.split()[0]
			try:mt[gene].append(l.split()[2])
			except KeyError:mt[gene]=[l.split()[2]]
		elif l.startswith('pt'):
			gene=l.split()[0]
			try:pt[gene].append(l.split()[2])
			except KeyError:pt[gene]=[l.split()[2]]	
	mt_site_cov=[]
	for k in mt.keys():mt_site_cov=mt_site_cov+mt[k][150:-150]
	mt_site_cov=[int(i) for i in mt_site_cov if int(i)<3*mt_ref[f]]
	out=open(f+'.mt_cov.tsv','a')
	out.write('\n'.join([str(i) for i in mt_site_cov]))
	out.close()
	pt_site_cov=[]
	for k in pt.keys():pt_site_cov=pt_site_cov+pt[k][150:-150]
	pt_site_cov=[int(i) for i in pt_site_cov if int(i)<3*pt_ref[f]]
	#pt_site_cov_mean=mean(pt_site_cov)
	out=open(f+'.pt_cov.tsv','a')
	out.write('\n'.join([str(i) for i in pt_site_cov]))
	out.close()