import numpy as np
import sys

if len(sys.argv)<4:
	print("python ychi2bsin.py par_in yorp_chi2 outpar")
	print("par_in: the input_parameter file of the inversion (yorp)")
	print("yorp_chi2: the results file of inversion (yorp)")
	print("outpar: the output file of this program, which is required by yorp_e pragram")
	exit(1)

par=open(sys.argv[1]).readlines()
fflambda=par[0].split()[2]
ffbeta=par[1].split()[2]
ffp=par[2].split()[2]
ffnu=par[3].split()[2]
conreg=par[6].split()[0]
lml=par[7].split()[0]
lmm=par[7].split()[1]
rows=par[8].split()[0]
ffa=par[9].split()[2]
ffd=par[10].split()[2]
ffk=par[11].split()[2]
ffb=par[12].split()[2]
ffc=par[13].split()[2]
stopcondition=par[14]

f=np.loadtxt(sys.argv[2],comments='#')

f1=open(sys.argv[3],'w')

f1.write("%.15g %s inital lambda [deg] (0/1 - fixed/free) \n" % (f[0][0], fflambda))
f1.write("%.15g %s initial beta [deg] (0/1 - fixed/free)\n" % (f[0][1], ffbeta))
f1.write("%.15g %s inital period [hours] (0/1 - fixed/free)\n" % (f[0][2], ffp))
f1.write("%.15g %s inital YORP strength [rad/day^2] (0/1 - fixed/free)\n" % (f[0][3], ffnu))
f1.write("%.15g zero time [JD]\n" % f[0][5])
f1.write("%.15g initial rotation angle [deg]\n" % f[0][6])
f1.write("%s convexity regularization\n" % conreg)
f1.write("%s %s degree and order of spherical harmonics expansion\n" % (lml, lmm))
f1.write("%s number of rows\n" % rows)
f1.write("%.15g %s phase funct. param. 'a' (0/1 - fixed/free)\n" % (f[0][7], ffa))
f1.write("%.15g %s phase funct. param. 'd' (0/1 - fixed/free)\n" % (f[0][8], ffd))
f1.write("%.15g %s phase funct. param. 'k' (0/1 - fixed/free)\n" % (f[0][9], ffk))
f1.write("%.15g %s phase funct. param. 'b' (0/1 - fixed/free)\n" % (f[0][10], ffb))
f1.write("%.15g %s phase funct. param. 'c' (0/1 - fixed/free)\n" % (f[0][11], ffc))
f1.write("%s" % stopcondition)

f1.close()
