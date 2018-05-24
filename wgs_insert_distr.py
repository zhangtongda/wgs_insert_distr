#!/usr/bin/env python
#  zhangtongda
#  zhangtongda@genomics.cn

import sys
import numpy 
#from operator import itemgetter
from optparse import OptionParser
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

# some constants for sam/bam field ids
SAM_FLAG = 1
SAM_REFNAME = 2
SAM_MATE_REFNAME = 6
SAM_ISIZE = 8

parser = OptionParser()
parser.add_option("-N",
    dest="N",
    type="int",
    help="Number of aligment items to read. default:100000",
    default=1000000)

parser.add_option("-o",
    dest="output_file",
    help="Output file")

parser.add_option("-b",
     dest="bam_file",
     help="bam file")

parser.add_option("-p",
     dest="png_file",
     help="save png file")

(options, args) = parser.parse_args()


#if not options.N:
#    parser.error('N not given')

if not options.output_file:
    parser.error('Output file not given')

if not options.bam_file:
    parser.error('Bam file not given')
if not options.png_file:
    parser.error('Output png file not given')

required = 97
restricted = 3484
flag_mask = required | restricted
max_isize_analysis = 1000
L = [0 for n in range(max_isize_analysis + 1)]
c = 0
L2 = []

bf = pysam.AlignmentFile(options.bam_file, 'rb')
for l in bf:
    if c >= options.N:
        break
    flag = l.flag
    refname = l.reference_name
    isize = l.isize

    want = (isize >= 0 and isize <= max_isize_analysis)
    if want:
        c += 1
        L[isize] += 1
        L2.append(isize)
        #print(str(isize) + '\n')
bf.close()

#L = numpy.array(L)

def mean_std(L):
    s = sum(L)
    mean = s / float(len(L))
    sq_sum = 0.0
    for v in L:
        sq_sum += (v - mean)**2.0
    var = sq_sum / float(len(L))
    return mean, var**0.5


mean, stdev = mean_std(L2)

#mean = numpy.mean(L)
#stdev = numpy.std(L)

f = open(options.output_file, 'w')
for i in range(1,1000):
    o = str(i) + "\t" + str(float(L[i])/float(options.N)) + "\n"
    f.write(o)
f.close()
print('mean:' + str(mean) + '\tstdev:' + str(stdev))

#matplotlib.pyplot.figure(figsize=(1000, 500)) 
x, y = numpy.loadtxt(options.output_file, delimiter='\t', unpack=True)
matplotlib.pyplot.plot(x,y)
matplotlib.pyplot.xlabel('Length of insert fragment')  
matplotlib.pyplot.ylabel('density distribution')  
matplotlib.pyplot.title('Insert fragment ')  
matplotlib.pyplot.savefig(options.png_file, format='png', dpi=300) 

