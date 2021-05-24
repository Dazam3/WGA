#!/usr/env python3

'''
Example Command-line Initiation:
python3 Multi_Alignment.py Species_List.txt
'''

import sys
import subprocess
import os
import time
import datetime

Species = sys.argv[1]

GENOME = '/home/gareeves/WGA/1_genomes/'
TWOB = '/home/gareeves/WGA/2_2bit/'
LAV = '/home/gareeves/WGA/3_lav/'
PSL = '/home/gareeves/WGA/4_psl/'
CHAINS = '/home/gareeves/WGA/5_chain/'
MERCHAINS = '/home/gareeves/WGA/6_mchain/'
PNETS = '/home/gareeves/WGA/7_pnet/'
NETS = '/home/gareeves/WGA/8_net/'
AXTS = '/home/gareeves/WGA/9_axt/'
MAFS = '/home/gareeves/WGA/10_maf/'
ROAST = '/home/gareeves/WGA/11_Final/'
REFER = 'Nfur'

#####################################
### Define Functions
#####################################

def Submit_Jobs(group1,group2):

	'''
	'''
	
	for scaf1 in os.listdir(TWOB + group1[:4] + '/'):
		for scaf2 in os.listdir(TWOB + group2[:4] + '/'):
			if (scaf1[:-5] + '_' + scaf2[:-5] + '.lav') not in os.listdir(LAV + group1[:4] + '_' + group2[:4] + '/'):

				var1 = (TWOB + group1[:4] + '/' + scaf1)
				var2 = (TWOB + group2[:4] + '/' + scaf2)
				var3 = (LAV + group1 + '_' + group2 + '/' + scaf1[:-5] + '_' + scaf2[:-5] + '.lav')
				var4 = (scaf1[:-5] + '_' + scaf2[:-5])
				subprocess.call(['sbatch', '--job-name=' + var4, '-o', LAV + group1 + '_' + group2 + '/out/out.%j', '-e', LAV + group1 + '_' + group2 + '/err/err.%j', '--export=ONE=' + var1 + ',TWO=' + var2 + ',THREE=' + var3, '/home/gareeves/WGA/Scripts/Align_Submission.sh'])

	jobber = 0
	while jobber == 0:			    
		time.sleep(600)
		with open('/home/gareeves/WGA/Scripts/Tacker.txt', 'wt') as jobby:
			subprocess.call(['squeue', '-u', 'gareeves'], stdout=jobby)
		with open('/home/gareeves/WGA/Scripts/Tacker.txt', 'r') as jobby:
			follower = 0
			for line in jobby:
				follower = follower + 1
		subprocess.call(['rm', '/home/gareeves/WGA/Scripts/Tacker.txt'])
		if follower > 1:
			pass
		else:
			jobber = 1

	NEW = len(os.listdir(LAV + group1[:4] + '_' + group2[:4])) - 2
	OLD1 = len(os.listdir(TWOB + group1[:4] + '/'))
	OLD2 = len(os.listdir(TWOB + group2[:4] + '/'))
	
	if NEW == (OLD1 * OLD2):
		return(1)
	elif NEW < (OLD1 * OLD2):
		return(0)
	else:
		return(2)
		
def Submit_Jobs2(groupx):

	'''
	'''
	
	for align in os.listdir(PSL + groupx + '/'):
		if align[:-4] + '.chain' not in os.listdir(CHAINS + groupx + '/'):
			tmp_name = align[:-4].split('_')
			var1 = TWOB + groupx[:4] + '/' + tmp_name[0] + '.2bit'
			var2 = TWOB + groupx[5:] + '/' + tmp_name[1] + '.2bit'
			var3 = PSL + groupx + '/' + align
			var4 = CHAINS + groupx + '/' + align[:-4] + '.chain'
			var5 = align[:-4]
			subprocess.call(['sbatch', '--job-name=' + var5, '-o', CHAINS + groupx + '/out/out.%j', '-e', CHAINS + groupx + '/err/err.%j', '--export=THREE=' + var3 + ',ONE=' + var1 + ',TWO=' + var2 + ',FOUR=' + var4, '/home/gareeves/WGA/Scripts/Chain_Submission.sh'])
	        
	        
	jobber2 = 0
	while jobber2 == 0:
                time.sleep(600)
                with open('/home/gareeves/WGA/Scripts/Tacker2.txt', 'wt') as jobby2:
                        subprocess.call(['squeue', '-u', 'gareeves'], stdout=jobby2)
                with open('/home/gareeves/WGA/Scripts/Tacker2.txt', 'r') as jobby2:
                        follower2 = 0
                        for line in jobby2:
                                follower2 = follower2 + 1
                subprocess.call(['rm', '/home/gareeves/WGA/Scripts/Tacker2.txt'])
                if follower2 > 1:
                        pass
                else:
                        jobber2 = 1

	NEW = len(os.listdir(CHAINS + groupx + '/')) - 2
	OLD = len(os.listdir(PSL + groupx + '/'))
	
	if NEW == OLD:
		return(1)
	elif NEW < OLD:
		return(0)
	else:
		Return(2)

#####################################
### Step 1 Split Fastas and Reference
#####################################

tracker = 1
with open(Species, 'r') as inner:
	for line in inner:
		if tracker == 1:
			users = line.rstrip('\n').split(',')
		if tracker == 2:
			line = line.rstrip('\n').split('_')
			TYPE = line[0]
			TREE = line[1]
		tracker = tracker + 1


for sample in os.listdir(GENOME + 'whole/'):
	if sample[:4] in users:
	
		with open(GENOME + 'sizes/' + sample[:4] + '.sizes', 'wt') as outer1:
			subprocess.call(['/home/gareeves/WGA/packages/faSize', GENOME + 'whole/' + sample, '-detailed',], stdout=outer1)
			
		subprocess.call(['mkdir', GENOME + 'split/' + sample[:4]])
		subprocess.call(['faToTwoBit', GENOME + 'whole/' + sample, GENOME + 'w2b/' + sample[:-3] + '.2bit'])
		subprocess.call(['faSplit', 'byName', GENOME + 'whole/' + sample, GENOME + 'split/' + sample[:4] + '/'])
		
		Sp_Group = {'Nfur':'NW', 'Aaus':'ZZ', 'Alim':'ZZ', 'Olat':'len', 'Gacu':'scaffold', 'Drer':'len', 'Blan':'xfSc', 'Test':'NW', 'Exam':'NW', 'Samp':'NW'}
		subprocess.call(['mkdir', GENOME + 'split/' + sample[:4] + '/temp'])
		for scaf in os.listdir(GENOME + 'split/' + sample[:4] + '/'):
		    if Sp_Group[sample[:4]] == 'len':
		        if len(scaf[:-3]) >= 3:
		            subprocess.call(['mv', GENOME + 'split/' + sample[:4] + '/' + scaf, GENOME + 'split/' + sample[:4] + '/temp/'])
		    else:
		       if Sp_Group[sample[:4]] in scaf:
		            subprocess.call(['mv', GENOME + 'split/' + sample[:4] + '/' + scaf, GENOME + 'split/' + sample[:4] + '/temp/'])
		
		if len(os.listdir(GENOME + 'split/' + sample[:4] + '/temp/')) != 0:
		    with open(GENOME + 'split/' + sample[:4] + '/NonChr.fa', 'wt') as joiner:
		        for cont in os.listdir(GENOME + 'split/' + sample[:4] + '/temp/'):
		            with open(GENOME + 'split/' + sample[:4] + '/temp/' + cont, 'r') as temper:
		                for line in temper:
		                    joiner.write(line)
		    subprocess.call(['rm', '-rf', GENOME + 'split/' + sample[:4] + '/NonChr.fa'])
		
		subprocess.call(['rm', '-rf', GENOME + 'split/' + sample[:4] + '/temp'])
		print(sample[:4] + ' split by chromosome')

print('\n\nStep 1 Complete\n\n')


#####################################
### Step 2 Generate 2bit files
#####################################		

		
for folder in os.listdir(GENOME + 'split/'):
	if folder[:4] in users:
		subprocess.call(['mkdir', TWOB + folder[:4]])
		for scaffold in os.listdir(GENOME + 'split/' + folder[:4] + '/'):
			#if (scaffold[0:2] == 'NC') or (len(scaffold) <= 5 and scaffold[0:2] != 'MT') or (len(scaffold) <= 12 and scaffold[0:2] == 'AA'):
			#for test run, comment line below for real run
			if 'yes' == 'yes':
				subprocess.call(['faToTwoBit', GENOME + 'split/' + folder[:4] + '/' + scaffold, TWOB + folder[:4] + '/' + scaffold[:-3] + '.2bit'])
	
		print('Converted ' + folder[:4] + 'to 2bit format' )
		
print('\n\nStep 2 Complete\n\n')


#####################################
### Step 3 Pairwise Alignment
#####################################


for folder1 in os.listdir(TWOB):
	for folder2 in os.listdir(TWOB):
		if (folder1[:4] == REFER) and (folder2[:4] in users) and folder1[:4] != folder2[:4]:
			if (folder1[:4] + '_' + folder2[:4] + '/') not in os.listdir(LAV):
				print('Alinging ' + folder1[:4] + ' and ' + folder2[:4])
				subprocess.call(['mkdir', LAV + folder1[:4] + '_' + folder2[:4]])
				subprocess.call(['mkdir', LAV + folder1[:4] + '_' + folder2[:4] + '/out'])
				subprocess.call(['mkdir', LAV + folder1[:4] + '_' + folder2[:4] + '/err'])

				counter = 0
				while counter == 0:
					counter = Submit_Jobs(folder1,folder2)
				if counter == 2:
					print('Overuse Error in Processing' + folder1[:4] + ' and ' + folder2[:4])
			
print('\n\nStep 3 Complete\n\n')


#####################################
### Step 4 PSL Conversion
#####################################	
	
for folder in os.listdir(LAV):
	if len(folder) > 3 and folder[:4] in users and folder[5:] in users:
		subprocess.call(['mkdir', PSL + folder])
		
		for align in os.listdir(LAV + folder + '/'):
			if '.lav' in align:
				subprocess.call(['lavToPsl', LAV + folder + '/' + align, PSL + folder + '/' + align[:-4] + '.psl'])
		print('Converted ' + folder + ' to psl format')
			
print('\n\nStep 4 Complete\n\n')


#####################################
### Step 5 Chain Generation
#####################################
	
for folder in os.listdir(PSL):
	if folder[:4] in users and folder[5:] in users:
		print('Chaining ' + folder[:4] + ' and ' + folder[5:])
		subprocess.call(['mkdir', CHAINS + folder])
		subprocess.call(['mkdir', CHAINS + folder + '/out'])
		subprocess.call(['mkdir', CHAINS + folder + '/err'])
		counter = 0
		while counter == 0:
			counter = Submit_Jobs2(folder)
		if counter == 2:
			print('Overuse Error in Processing' + folder)

print('\n\nStep 5 Complete\n\n')


#####################################
### Step 6 Chain Merging and Sorting
#####################################

for folder in os.listdir(CHAINS):
	if len(folder) > 3 and folder[:4] in users and folder[5:] in users:
		with open(MERCHAINS + folder + '.chain', 'wt') as outer2:
			crystal = []
			for shard in os.listdir(CHAINS + folder + '/'):
				crystal.append(CHAINS + folder + '/' + shard)
			with open('./chain_list.txt', 'wt') as temper1:
				temper1.write(str('\n'.join(crystal)))
			subprocess.call(['chainMergeSort', '-inputList=./chain_list.txt'], stdout=outer2)
			subprocess.call(['rm', './chain_list.txt'])
			print('Merged and Sorted Chains of ' + folder[:4] + ' and ' + folder[5:])
			
print('\n\nStep 6 Complete\n\n')


#####################################
### Step 7 Pre-Netting Chains
#####################################

for chain in os.listdir(MERCHAINS):
	if chain[:4] in users and chain[5:-6] in users:
		back1 = GENOME + 'sizes/' + chain[:4] + '.sizes'
		back2 = GENOME + 'sizes/' + chain[5:-6] + '.sizes'
		subprocess.call(['chainPreNet', MERCHAINS + chain, back1, back2, PNETS + 'pn_' + chain])
		print('Pre-net Generated for ' + chain[:4] + ' and ' + chain[5:-6])
		

print('\n\nStep 7 Complete\n\n')


#####################################
### Step 8 Net Generation
#####################################

for prenet in os.listdir(PNETS):
	if prenet[3:7] in users and prenet[8:-6] in users:
		back1 = GENOME + 'sizes/' + prenet[3:7] + '.sizes'
		back2 = GENOME + 'sizes/' + prenet[8:-6] + '.sizes'
		subprocess.call(['chainNet', PNETS + prenet, '-minSpace=1', back1, back2, NETS + 'temp.net', './dump.net'])
		subprocess.call(['netSyntenic', NETS + 'temp.net', NETS + prenet[3:7] + '_' + prenet[8:-6] + '.net'])
		subprocess.call(['rm', NETS + 'temp.net'])
		subprocess.call(['rm', './dump.net'])
		print('Net Generated for ' + prenet[3:7] + ' and ' + prenet[8:-6])
		
print('\n\nStep 8 Complete\n\n')


#####################################
### Step 9 AXT Conversion and Sorting
#####################################		

for net in os.listdir(NETS):
	if net[:4] in users and net[5:-4] in users:
		ref1 = GENOME + 'w2b/' + net[:4] + '.2bit'
		ref2 = GENOME + 'w2b/' + net[5:-4] + '.2bit'
		subprocess.call(['netToAxt', NETS + net, PNETS + 'pn_' + net[:-4] + '.chain', ref1, ref2, AXTS + 'temp.axt'])
		subprocess.call(['axtSort', AXTS + 'temp.axt', AXTS + net[:-4] + '.axt'])
		subprocess.call(['rm', AXTS + 'temp.axt'])
		print('Converted and Sorted ' + net[:4] + ' and ' + net[5:-4] + ' to axt format')
		
		
print('\n\nStep 9 Complete\n\n')


#####################################
### Step 10 Maf Generation
#####################################		
		
for axt in os.listdir(AXTS):
	if axt[:4] in users and axt[5:-4] in users:
		back1 = GENOME + 'sizes/' + axt[:4] + '.sizes'
		back2 = GENOME + 'sizes/' + axt[5:-4] + '.sizes'
		subprocess.call(['axtToMaf', AXTS + axt, back1, back2, MAFS + axt[:4] + '.' + axt[5:-4] + '.sing.maf', '-tPrefix=' + axt[:4] + '.', '-qPrefix=' + axt[5:-4] + '.'])
		print('Maf Generated for' + axt[:4] + ' and ' + axt[5:-4])	
		
print('\n\nStep 10 Complete\n\n')


#####################################
### Step 11 Multi-Alignment
#####################################

now = datetime.datetime.now()
current = now.strftime("%Y-%m-%d")
Namer = '_'.join(users)

if TYPE == 'REFERENCE-BASED':
	subprocess.call(['mkdir', MAFS + 'RB-' + current + '/'])
	for file in os.listdir(MAFS):
		if file[:4] == REFER and file[5:-9] in users:
			subprocess.call(['cp', MAFS + file, MAFS + 'RB-' + current + '/'])
	
	os.chdir(MAFS + 'RB-' + current + '/')	
	subprocess.call(['roast', '+', 'E=' + REFER, TREE, '*.*.sing.maf', '../' + ROAST + 'RB_' + Namer + '.maf'])
	os.chdir('../../Scripts')
	subprocess.call(['rm', '-r', MAFS + 'RB-' + current])		

elif TYPE == 'ALL-WAY':
	subprocess.call(['mkdir', MAFS + 'AW-' + current + '/'])
	for file in os.listdir(MAFS):
		if file[:4] in users and file[5:-9] in users:
			subprocess.call(['cp', MAFS + file, MAFS + 'AW-' + current + '/'])
	
	os.chdir(MAFS + 'AW-' + current + '/')
	subprocess.call(['tba', '+', TREE, MAFS + 'AW-' + current + '/' + '*.*.sing.maf', '../' + ROAST + 'AW_' + Namer + '.maf'])
	os.chdir('../../Scripts')
	subprocess.call(['rm', '-r', MAFS + 'AW-' + current])		

elif TYPE == 'BOTH':
	subprocess.call(['mkdir', MAFS + 'AWR-' + current + '/'])
	for file in os.listdir(MAFS):
		if file[:4] in users and file[5:-9] in users:
			subprocess.call(['cp', MAFS + file, MAFS + 'AWR-' + current + '/'])
				
	os.chdir(MAFS + 'AWR-' + current + '/')
	subprocess.call(['tba', '+', 'E=' + REFER, TREE, MAFS + 'AWR-' + current + '/' + '*.*.sing.maf', '../' + ROAST + 'AWR_' + Namer + '.maf'])
	os.chdir('../../Scripts')
	subprocess.call(['rm', '-r', MAFS + 'AWR-' + current])
	
else:
	print('Multi-Alingment Type Not Recognized')		
		
print('\n\nProcess Complete\n\n')		
	

#####################################
### Optional File Cleanup
#####################################
'''
for file in os.listdir('../'):
	if file != '1_genomes' or file != 'Scripts' or file != 'packages' or file != '11_Final':
		subprocess.call(['rm', '-r', '../' + file + '/*'])
subprocess.call(['mkdir', '../3_lav/err'])
subprocess.call(['mkdir', '../3_lav/out'])
subprocess.call(['mkdir', '../5_chain/err'])
subprocess.call(['mkdir', '../5_chain/out'])
print('\n\nOptional File Clean-up Complete\n\n')
'''
