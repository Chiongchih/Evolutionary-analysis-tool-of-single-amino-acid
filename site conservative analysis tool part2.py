po = input('input the total number of datasets:')
fastaf='2-seq'
#print('the eggNOG output file name is ' + fastaf)
jj=0
seqinf = open("3-OG information.txt", 'w')
print('All files are ready now, running for STEP3 file name 3-OG information.txt')
while jj<int(po):
    file = open(fastaf + '-' + str(int(jj)+1) + '.txt.emapper.annotations', 'r')
    hs = len(file.readlines()) - 1
    #print(hs)
    files = open(fastaf + '-' + str(int(jj)+1) + '.txt.emapper.annotations', 'r')
    ii = 0
    while ii <= hs:
        line = files.readline().split('\t')
        mm = line[10].split('|')
        seqinf.write(line[0])
        seqinf.write(' ')
        seqinf.write(mm[0])
        seqinf.write(' ')
        seqinf.write(line[1])
        seqinf.write(' ')
        seqinf.write('\n')
        ii = ii + 1
    jj=jj+1

seqinf.close()
print('3-OG information.txt complete now, now running for STEP4 file name 4-OG (human only).txt')

fff = '3-OG information.txt'
file = open(fff,'r')
file2 = open(fff,'r')
iii='1-set.txt'
file3 = open(iii,'r')
fl2=len(file3.readlines())
fl=len(file2.readlines())
speice=input('input the specie number(for human, input: 9606):')
se=open('4-OG (human only).txt','w')
i = 0
import re
while i<fl:
    #print(i)
    line = file.readline().strip()
    ls=line.split(' ')
    #print(line)
    #print(ls)
    id = ls[0]
    og = ls[1]
    #print(str(ls[2].split('.')[0]))
    if ls[2].split('.')[0]==speice:
        se.write(id)
        se.write(' ')
        if re.match('KOG', ls[1]):
            se.write(ls[1])
            se.write(' ')
        else:
            se.write('ENOG41' + str(ls[1]))
            se.write(' ')
        se.write(ls[2])
        se.write(' ')
        file3 = open(iii, 'r')
        j=0
        while j<fl2:
            set = file3.readline().split(' ')
            if set[0]==id:
                #print(set[1])
                se.write(set[1])
                break
            j=j+1
        #se.write('\n')

    i=i+1
se.close()

print('4-OG (human only).txt complete now, now running for STEP5 file name 5-checked OG (human only).txt')

import requests
import re
import os
#PART 1Ace-k  meth K-R sumo K  Phos
path = os.path.split(os.path.realpath(__file__))[0]
seqf='4-OG (human only).txt'
sequences=open(seqf,'r')
seqs=sequences.readlines()
seqlen = len(seqs)-1
seqr = open(seqf,'r')
#aaa='K'
seqinf = open('5-checked OG (human only).txt','w')
#ty=open('not K2.txt','w')
i=0
print('Running...')
while i <= seqlen:
    seq = seqr.readline().strip().split(' ')#读取本地文件信息
    url = "http://eggnogapi.embl.de/nog_data/text/fasta/" + seq[1]  # URL搜索出相应orthology group的raw_alignment文件
    #print(url)
    r = requests.get(url)  # 读取raw_alignment文件
    ralgn =(r.text).strip().split()
    #print(ralgn)
    mm=len(ralgn)-1
    y=0
    yy=0
    ss=''
    while y<=mm:
        if re.match('>',ralgn[y]):
            if re.search(seq[2],ralgn[y]):
                yy=y+1
                while re.search('>',ralgn[yy])is None:
                    ss=ss+ralgn[yy]
                    if (yy+1)< mm:
                        yy=yy+1
                    else:
                        break

                #print(ss)
        y = y + 1
    ss = ''.join(ss.split())
    seqq=''
    #print('1/3 part of job has done.')
    urll = "http://www.uniprot.org/uniprot/" + seq[0] + ".fasta"
    rr = requests.get(urll)
    fasta = (rr.text).split('\n')
    yu = 1
    ser = ''
    while yu < len(fasta) - 1:
        ser = ser + fasta[yu]
        yu = yu + 1
    #print(ser)
    #print(ss)
    if ser==ss:
        seqinf.write(seq[0])
        seqinf.write(' ')
        seqinf.write(seq[1])
        seqinf.write(' ')
        seqinf.write(seq[2])
        seqinf.write(' ')
        seqinf.write(seq[3])
        seqinf.write('\n')
        #print('fuck')
    #else:
        #print('yes')
    if i==seqlen//3:
        print('1/3 sequences checked, keep running...')
    if i==(seqlen//3)*2:
        print('2/3 sequences checked, keep running...')
    i=i+1
seqinf.close()
print('5-checked OG (human only).txt complete now, now running for STEP6 file name 6-checked OG (human only).txt')

seqf='5-checked OG (human only).txt'
sequences=open(seqf,'r')
seqs=sequences.readlines()
seqlen = len(seqs)-1
seqr = open(seqf,'r')
aaa=input('input the amino acid you want to analyze(e.g. Lysine, input:K):')
seqinf = open("6-aligned sites("+aaa+").txt",'w')
seqinf2 = open("6-original sites("+aaa+").txt",'w')
#ty=open('not  '+aaa+'phoss.txt','w')
i=0
print('Running...')
while i < seqlen:
    seq = seqr.readline().strip().split(' ')
    seqinf.write(seq[0])
    seqinf.write(' ')
    seqinf.write(seq[1])
    seqinf.write(' ')
    seqinf.write(seq[2])
    seqinf.write(' ')
    seqinf2.write(seq[0])
    seqinf2.write(' ')
    seqinf2.write(seq[1])
    seqinf2.write(' ')
    seqinf2.write(seq[2])
    seqinf2.write(' ')
    url = "http://eggnogapi.embl.de/nog_data/text/raw_alg/" + seq[1]  # URL搜索出相应orthology group的raw_alignment文件
    #print(url)
    r = requests.get(url)  # 读取raw_alignment文件
    ralgn =(r.text).strip().split()
    #print(ralgn)
    mm=len(ralgn)-1
    y=0
    yy=0
    ss=''
    while y<=mm:
        if re.match('>',ralgn[y]):
            if re.search(seq[2],ralgn[y]):
                yy=y+1
                while re.search('>',ralgn[yy])is None:
                    ss=ss+ralgn[yy]
                    if (yy+1)< mm:
                        yy=yy+1
                    else:
                        break

                #print(ss)
        y = y + 1
    ss = ''.join(ss.split())
    seqq=''
    if re.search(';',seq[3]):
        seqq=seq[3].split(';')
        shn=len(seqq)-1
        pp=0
        while pp <= shn:
            #print(shn)
            site2=int(seqq[pp])
            ssn = len(ss) -1#raw_alignment序列长度
            tt=0
            ttr=0
            #print('mimi')
            while tt<=ssn:#进入在线序列搜索
                if ss[tt] != '-':
                    ttr = ttr+1#取到site正确的位置
                if ttr == site2:
                    if ss[tt] == aaa:
                        if pp==0:
                            seqinf.write(str(tt))
                            seqinf2.write(str(site2))
                        else:
                            seqinf.write(';')
                            seqinf2.write(';')
                            seqinf.write(str(tt))
                            seqinf2.write(str(site2))
                    ttr=100000
                tt=tt+1

            pp=pp+1
            #print(pp)

    else:
        site = int(seq[3])
        ssn = len(ss) - 1
        tt = 0
        ttr = 0
        while tt <= ssn:
            if ss[tt] != '-':
                ttr = ttr + 1
            if ttr == site:
                if ss[tt] == aaa:
                    seqinf.write(str(tt))
                    seqinf2.write(str(site))
                ttr=ttr+1
            tt = tt + 1
    seqinf.write('\n')
    seqinf2.write('\n')
    if i==seqlen//3:
        print('1/3 part of job done, keep running...')
    if i==(seqlen//3)*2:
        print('2/3 part of job done, keep running...')
    i=i+1
seqinf.close()
seqinf2.close()
#ty.close()
print('6-aligned sites('+str(aaa)+').txt and 6-original sites('+str(aaa)+').txt complete now, now running for STEP7 file name 7-modify aligned sites('+str(aaa)+').txt')

seqf='6-aligned sites('+str(aaa)+').txt'
file = open(seqf,'r')
lines= open(seqf,'r')
whole= open("7-modify aligned sites("+str(aaa)+').txt','w')
line = len(lines.readlines())-1
i=0
print('Running...')
while i< line:
    read = file.readline().strip().split(' ')  # sites文件中的一个蛋白讯息读取
    if len(read)==4:
        sites=read[3]
        whole.write(read[0])
        whole.write(' ')
        whole.write(read[1])
        whole.write(' ')
        whole.write(read[2])
        whole.write(' ')
        if sites[0]==';':
            n=len(sites)
            whole.write(sites[1:n-1])
            whole.write('\n')
        else:
            whole.write(sites)
            whole.write('\n')
    i=i+1

whole.close()

print('7-modify aligned sites('+str(aaa)+').txt complete now, now running for STEP8 file name 8-control sites of sites('+str(aaa)+').txt')

seqf2 = "7-modify aligned sites("+str(aaa)+').txt'
sites = open(seqf2,'r')
lines= open(seqf2,'r')
whole=open("8-control sites of sites("+str(aaa)+').txt','w')
line = len(lines.readlines())-1
import re
import requests
i=0
uu = 0
tuu = 0
print('Running...')
while i<= line:
    read = sites.readline().strip().split(' ')  # sites文件中的一个蛋白讯息读取
    id = read[0]
    whole.write(id)
    whole.write(' ')
    og = read[1]
    whole.write(og)
    whole.write(' ')
    name = read[2]
    whole.write(name)
    whole.write(' ')
    url = "http://eggnogapi.embl.de/nog_data/text/raw_alg/" + read[1]  # URL搜索出相应orthology group的raw_alignment文件
    r = requests.get(url)  # 读取raw alignment文件
    fasta = (r.text).split()
    pp = 0
    ppp = 0
    nn = 0
    while pp < (len(fasta) - 1):  # OG里面搜有无对应number的蛋白序列
        ss = ''
        if re.search(name, fasta[pp]):
            nn = nn + 1
            ppp = pp
            while re.search('>', fasta[ppp + 1]) is None:
                ss = ss + fasta[ppp + 1]
                if (ppp + 1) < len(fasta)-1:
                    ppp = ppp + 1  # get到raw alignment中的每一个序列
                else:
                    break

            ss = ''.join(ss.split())
            ssl=len(ss)-1
            tt=0
            while tt<=ssl:
                if ss[tt]==aaa:
                    whole.write(str(tt))
                    whole.write(' ')

                tt=tt+1

        pp = pp + 1
    whole.write('\n')
    if i==line//3:
        print('1/3 part of job done, keep running...')
    if i==(line//3)*2:
        print('2/3 part of job done, keep running...')
    i=i+1

whole.close()
print('8-control sites('+str(aaa)+').txt complete now, now running for STEP9 file name 9-checked modify aligned sites('+str(aaa)+').txt')

seqf='7-modify aligned sites('+str(aaa)+').txt'
filename =seqf
file = open(filename,'r')
file2 =  open(filename,'r')
fl=len(file2.readlines())
se = open('9-checked modify aligned sites('+str(aaa)+').txt','w')
i=0
import re
while i<fl:
    line = file.readline().split(' ')
    if re.search(';',line[3]):
        se.write(line[0]+' '+line[1]+' '+line[2]+' ')
        pos=line[3].strip().split(';')
        pl=len(pos)
        j=0
        while j<pl-1:
            se.write(pos[j]+' ')
            j=j+1
        se.write(pos[pl-1])
        se.write('\n')
    else:
        se.write(line[0]+' '+line[1]+' '+line[2]+' '+line[3])

    i=i+1
se.close()
print('9-checked modify aligned sites('+str(aaa)+').txt complete now, now running for STEP10 file name 10-rank information('+str(aaa)+').txt')

import re
import requests
import math
import os
import sys
from decimal import Decimal, getcontext

def pcal(aa,bb,cc,dd):
    ee = aa + bb + cc + dd
    a = math.factorial(aa + bb)
    b = math.factorial(cc + dd)
    c = math.factorial(aa + cc)
    d = math.factorial(bb + dd)
    e = math.factorial(aa)
    f = math.factorial(bb)
    g = math.factorial(cc)
    h = math.factorial(dd)
    i = math.factorial(ee)
    p1 = a * b * c * d
    p2 = e * f * g * h * i

    p = Decimal(p1)/Decimal(p2)
    #print(p)
    #pp = 1/p
    return p

def flen(fasta):
    pp=0
    nn=0
    while pp < (len(fasta) - 1):  # OG里面搜有无对应number的蛋白序列
        if re.match('>', fasta[pp]):
            nn = nn + 1
        pp=pp+1
    return nn

def tcal(tread,tll,fasta,aaa):
    #print(tread)
    #print(tread)
    pp = 0
    ppp = 0
    uu = 0
    nn = 0
    ii=3

    while pp < (len(fasta) - 1):  # OG里面搜有无对应number的蛋白序列
        ss = ''
        if re.match('>', fasta[pp]):
            nn = nn + 1
            ppp = pp
            while re.search('>', fasta[ppp + 1]) is None:
                ss = ss + fasta[ppp + 1]
                if (ppp + 1) < (len(fasta) - 1):
                    ppp = ppp + 1  # get到raw alignment中序列
                else:
                    break

            ss = ''.join(ss.split())

            # print(ss[number])
            while ii <= tll:  # 一个蛋白的每个位点搜索
                number = int(tread[ii])

                if ss[number] == aaa:
                    uu = uu + 1

                ii = ii + 1
            ii = 3

        pp = pp + 1
    return uu

def pepg(name, n):
    import requests
    urll = "http://www.uniprot.org/uniprot/"+name+".fasta"
    rr = requests.get(urll)
    seq = rr.text.split('\n')
    sel = len(seq)
    u = 1
    feq = ''
    while u < sel:
        # print(seq[u])
        feq = feq + seq[u]
        u = u + 1
    ly = len(feq)
    if n - 15 > 0:
        if n + 15 > ly:
            yu = feq[(n - 16):ly]
        else:
            yu = feq[(n - 16):(n + 15)]

    if n - 15 <= 0:
        if n + 15 > ly:
            yu = feq[0:ly]
        else:
            yu = feq[0:n + 15]
    print(yu)
    return yu

def acal(aread, all, fasta, a,b,c,og,name,aaa):

    d = 0
    nn = 0
    ii = 3

    while ii <= all:  # 一个蛋白的每个位点搜索
        number = int(aread[ii])
        pp = 0
        ppp = 0
        while pp < (len(fasta) - 1):  # OG里面搜有无对应number的蛋白序列
            ss = ''
            if re.match('>', fasta[pp]):
                nn = nn + 1
                ppp = pp
                while re.search('>', fasta[ppp + 1]) is None:
                    ss = ss + fasta[ppp + 1]
                    if (ppp + 1) < (len(fasta) - 1):
                        ppp = ppp + 1  # get到raw alignment中序列
                    else:
                        break

                ss = ''.join(ss.split())
                if ss[number] == aaa:
                    d = d + 1
            pp=pp+1
        #print(a - d, b - c, a, b, c, d,d/a, c/b)
        #pv = pcal(a - d, b - c, c, d)
       # d=0
        #if d/a>c/b:
        rank.write(aread[0])
        rank.write(' ')
        rank.write(og)
        rank.write(' ')
        rank.write(name)
        rank.write(' ')
        rank.write(aread[ii])
        rank.write(' ')
        rank.write(str(d / a))
        rank.write(' ')
        rank.write(str(d))
        rank.write(' ')
        rank.write(str(a))
        rank.write(' ')
        rank.write(str(c / b))
        rank.write(' ')
        rank.write(str(c))
        rank.write(' ')
        rank.write(str(b))
            #rank.write(' ')
            #rank.write(str(pcal(a,b,c,d)))
        rank.write('\n')

        d=0

        #print(pv)
        ii = ii + 1

afile = '9-checked modify aligned sites('+str(aaa)+').txt'
tfile='8-control sites of sites('+str(aaa)+').txt'
acetyl_files = open(afile,'r')
acetyl_file = open(afile,'r')
total_files = open(tfile,'r')
total_file= open(tfile,'r')
rank=open('10-rank information('+aaa+').txt','w')
alen = len(acetyl_files.readlines())
tlen = len(total_files.readlines())
i = 0
print('Files are ready, job start to run.')
print('Running...')
while i < alen:
    aread = acetyl_file.readline().strip().split(' ')# acetyl sites文件中的一个蛋白讯息读取
    tread = total_file.readline().strip().split(' ')# total sites文件中的一个蛋白讯息读取
    if aread[0]==tread[0]:
        id = aread[0]
        og = aread[1]
        #print(tread)
        name = aread[2]
        all = len(aread) - 1
        tll = len(tread)-1
        url = "http://eggnogapi.embl.de/nog_data/text/raw_alg/" + aread[1]  # URL搜索出相应orthology group的raw_alignment文件
        r = requests.get(url)  # 读取fasta文件
        fasta = (r.text).split()
        #print(flen(fasta))
        fl=flen(fasta)
        a=fl
        b=(tll-2)*fl
        #print(b)
        pp=0
        c=tcal(tread,tll,fasta,aaa)
        #print(c)
        acal(aread, all, fasta,a,b,c,og,name,aaa)
    if i==alen//3:
        print('1/3 part of job done, keep running...')
    if i==(alen//3)*2:
        print('2/3 part of job done, keep running...')
    i = i + 1

rank.close()
print('10-rank information('+str(aaa)+').txt complete now, now running for STEP11 file name 11-position return.txt')

rankn='10-rank information('+str(aaa)+').txt'
ranking = open(rankn,'r')
rl=ranking.readlines()
rank = open(rankn,'r')
opo='6-original sites('+aaa+').txt'
apo='6-aligned sites('+aaa+').txt'
opof=open(opo,'r')
apof=open(apo,'r')
lopo=len(opof.readlines())
apof=len(apof.readlines())
result=open('11-position return.txt','w')
rlen=len(rl)-1
i=0
import re
import requests
print('files are ready, job start to run')
print('Running...')
while i<=rlen:
    read = rank.readline().strip().split(' ')
    #print(read)
    name=read[0]
    og=read[1]
    seed=read[2]
    site=read[3]
    pv1=read[4]
    dd=read[5]
    aa=read[6]
    pv2=read[7]
    cc=read[8]
    bb=read[9]
    #pvv=read[10]
    opof = open(opo, 'r')
    apof = open(apo, 'r')
    ii=0
    oss=0
    while ii<lopo:
        lol = opof.readline().strip().split(' ')
        lal = apof.readline().strip().split(' ')
        if str(lol[0])==str(name):
            #print(lol[0])
            #print(lol[0])
            slal=lal[3].split(';')
            slol=lol[3].split(';')
            #print(lal)
            #print(slal)
            lslal=len(slal)
            iii=0
            aos = 0
            #print(site)
            while iii<lslal:

                if str(slal[iii])==str(int(site)):
                   # print('yes')
                    aos=iii
                    oss=int(slol[aos])
                iii=iii+1
            #print(slal)
        ii=ii+1
    result.write(name)
    result.write(' ')
    result.write(og)
    result.write(' ')
    result.write(str(oss))
    result.write(' ')
    result.write(pv1)
    result.write(' ')
    result.write(dd)
    result.write(' ')
    result.write(aa)
    result.write(' ')
    result.write(pv2)
    result.write(' ')
    result.write(cc)
    result.write(' ')
    result.write(bb)
    result.write('\n')
    i=i+1

result.close()

print('11-position return.txt complete now, now running for STEP12 file name 12-result.txt')

import math
import os
import sys
from decimal import Decimal, getcontext

def pepg(name, n):
    import requests
    urll = "http://www.uniprot.org/uniprot/"+name+".fasta"
    rr = requests.get(urll)
    seq = rr.text.split('\n')
    sel = len(seq)
    u = 1
    feq = ''
    while u < sel:
        # print(seq[u])
        feq = feq + seq[u]
        u = u + 1
    ly = len(feq)
    if n - 15 > 0:
        if n + 15 > ly:
            yu = feq[(n - 16):ly]
        else:
            yu = feq[(n - 16):(n + 15)]

    if n - 15 <= 0:
        if n + 15 > ly:
            yu = feq[0:ly]
        else:
            yu = feq[0:n + 15]
    #print(yu)
    return yu

def pcal(aa,bb,cc,dd):
    #print(aa,bb,cc,dd)
    ee = aa + bb + cc + dd
    a = math.factorial(aa + dd)
    b = math.factorial(cc + bb)
    c = math.factorial(aa + bb)
    d = math.factorial(cc + dd)
    e = math.factorial(aa)
    f = math.factorial(bb)
    g = math.factorial(cc)
    h = math.factorial(dd)
    i = math.factorial(ee)
    p1 = a * b * c * d
    p2 = e * f * g * h * i

    p = Decimal(p1)/Decimal(p2)
    #print(p)
    #pp = 1/p
    return p

jn = '11-position return.txt'
frank2 = open(jn, 'r')
rankl = frank2.readlines()
rl = len(rankl)
result = open('12-result.txt', 'w')
i = 0
frank = open(jn, 'r')

while i < rl:
    line = frank.readline().split(' ')
    if i<rl:
        #print(str(i) + '//' + str(rl))
        id = line[0]
        og = line[1]
        site = line[2]
        pvv1 = line[3]
        dd = line[4]
        aa = line[5]
        pvv2 = line[6]
        cc = line[7]
        bb = line[8]
        if pvv1==pvv1:
            result.write(id)
            result.write('-')
            result.write(site)
            result.write('	')
            result.write(pepg(id, int(site)))
            result.write('	')
            result.write(og)
            result.write('	')
            result.write(pvv1)
            result.write('	')
            result.write(dd)
            result.write('	')
            result.write(aa)
            result.write('	')
            result.write(pvv2)
            result.write('	')
            result.write(cc)
            result.write('	')
            result.write(bb)
            result.write('	')
            result.write(str(pcal(int(aa), int(bb), int(cc), int(dd))))
            result.write('\n')
    if i==rl//3:
        print('1/3 part of job done, keep running...')
    if i==(rl//3)*2:
        print('2/3 part of job done, keep running...')

    i = i + 1

result.close()
print('12-result.txt complete now')
print('Congratulation! This job is finished successfully!')