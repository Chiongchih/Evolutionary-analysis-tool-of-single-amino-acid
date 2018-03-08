datas=input('input the dataset name:')
ds=open(datas,'r')
ds2=open(datas,'r')
specie=input('input the specie name:(human/rat/mouse)')
aas=input('input the amino acid you would like to analyze(K/R.../all):')
se=open('1-sites.txt','w')
print('All files are ready now, start to run STEP1 1-sites.txt')
dl=ds2.readlines()
dll=len(dl)-1
i=0
ds.readline()
ds.readline()
ds.readline()
ds.readline()
while i<=dll-3:
    line=ds.readline().split('	')
    if len(line)>1:
        #print(line)
        # se.write(line[2]+' ')
        u = line[4]
        ii = 2
        spi = line[6]
        if str(aas)!='all':
            if str(u[0])==aas:
                if str(spi) == specie:
                    while ii < len(u):
                        # print(u[ii])
                        if u[ii] == '-':
                            iii = 1
                            sp = ''
                            while iii < ii:
                                sp = sp + u[iii]
                                # print(sp)
                                iii = iii + 1
                            break
                        ii = ii + 1
                    # print(sp)
                    if len(line[2]) == 6:
                        se.write(line[2] + ' ' + sp)
                        se.write('\n')
        if str(aas)=='all':
            if str(spi) == specie:
                while ii < len(u):
                    # print(u[ii])
                    if u[ii] == '-':
                        iii = 1
                        sp = ''
                        while iii < ii:
                            sp = sp + u[iii]
                            # print(sp)
                            iii = iii + 1
                        break
                    ii = ii + 1
                # print(sp)
                if len(line[2]) == 6:
                    se.write(line[2] + ' ' + sp)
                    se.write('\n')


    i=i+1
se.close()

print('1-sites.txt complete now, now running for 1-sites.txt')
ds=open('1-sites.txt','r')
ds2=open('1-sites.txt','r')
se2=open('1-set.txt','w')
dl=ds2.readlines()
dll=len(dl)
lib=dict()
#print(dll)
i=0
while i<dll:
    line=ds.readline().split(' ')
    name=line[0]
    sp=line[1]
    if lib.get(name) is None:
        lib[name]=sp.strip()
    else:
        lib[name]=lib[name]+';'+sp.strip()

    i=i+1
#print(lib)
for key in lib:
    se2.write(key+' '+lib[key])
    se2.write('\n')

se2.close()

import requests
datas='1-set.txt'
ds=open(datas,'r')
ds2=open(datas,'r')
output='2-seq-'
print('1-set.txt complete now, now running for STEP2 file name 2-seq-x.txt ')
dl=ds2.readlines()
dll=len(dl)
i=0
p=1
se=open(output+str(p)+'.txt','w')
nu=1+(dll//4500)
while i<dll:
    if i < p * 4500:
        line = ds.readline().split(' ')
        urll = "http://www.uniprot.org/uniprot/" + line[0] + ".fasta"
        rr = requests.get(urll)
        fasta = (rr.text).split('\n')
        yu = 1
        ser = ''
        while yu < len(fasta) - 1:
            ser = ser + fasta[yu]
            yu = yu + 1
        if len(ser)!=0:
            #print(line[0])
            se.write('>' + line[0])
            se.write('\n')
            #print(ser)
            se.write(ser)
            se.write('\n')
    else:
        line = ds.readline().split(' ')
        urll = "http://www.uniprot.org/uniprot/" + line[0] + ".fasta"
        rr = requests.get(urll)
        fasta = (rr.text).split('\n')
        yu = 1
        ser = ''
        while yu < len(fasta) - 1:
            ser = ser + fasta[yu]
            yu = yu + 1
        if len(ser)!=0:
            se.write('>' + line[0])
            se.write('\n')
            se.write(ser)
            se.write('\n')
        se.close()
        p = p + 1
        se = open(output + str(p) + '.txt', 'w')
    if i==dll//10:
        print('1/10 part of job done, keep running...')
    if i==(dll//10)*2:
        print('2/10 part of job done, keep running...')
    if i==(dll//10)*3:
        print('3/10 part of job done, keep running...')
    if i==(dll//10)*4:
        print('4/10 part of job done, keep running...')
    if i==(dll//10)*5:
        print('5/10 part of job done, keep running...')
    if i==(dll//10)*6:
        print('6/10 part of job done, keep running...')
    if i==(dll//10)*7:
        print('7/10 part of job done, keep running...')
    if i==(dll//10)*8:
        print('8/10 part of job done, keep running...')
    if i==(dll//10)*9:
        print('9/10 part of job done, keep running...')
    i=i+1

print('2-seqs.txt complete now, please upload them on eggnog for further analysis')