import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# plt.style.use('default')

if len(sys.argv)<3:
	print("python fig1st.py path yorp_chi2")
	print("path: the path where you want to save the results")
	print("yorp_chi2: the path of the results file of yorp_e")
	exit(1)


os.system("rm -rf %s" % sys.argv[1])
os.system("mkdir %s" % sys.argv[1])

outf=open('%s/meanstd.txt' % sys.argv[1],'w')
outf.write('# parameter mean 3*std\n')

f=np.loadtxt(sys.argv[2],comments='#')

plt.figure()
plt.hist(f[:,0],density=1,bins=20,color='skyblue')
plt.xlabel(r'$\lambda$ (deg)')
plt.ylabel(r'density')
plt.savefig('%s/lambda.eps' % sys.argv[1])
plt.savefig('%s/lambda.pdf' % sys.argv[1])
outf.write('lambda %.15g %.15g\n' % (np.mean(f[:,0]), 3*np.std(f[:,0])))

plt.figure()
plt.hist(f[:,1],density=1,bins=20,color='skyblue')
plt.xlabel(r'$\beta$ (deg)')
plt.ylabel(r'density')
plt.savefig('%s/beta.eps' % sys.argv[1])
plt.savefig('%s/beta.pdf' % sys.argv[1])
outf.write('beta %.15g %.15g\n' % (np.mean(f[:,1]), 3*np.std(f[:,1])))

plt.figure()
plt.hist(f[:,2],density=1,bins=20,color='skyblue')
plt.xlabel(r'$P$ (h)')
plt.ylabel(r'density')
plt.savefig('%s/P.eps' % sys.argv[1])
plt.savefig('%s/P.pdf' % sys.argv[1])
outf.write('P %.15g %.15g\n' % (np.mean(f[:,2]), 3*np.std(f[:,2])))

plt.figure()
plt.hist(f[:,3],density=1,bins=20,color='skyblue')
plt.xlabel(r'$\upsilon$ (rad day$^{-2}$)')
plt.ylabel(r'density')
plt.savefig('%s/yorp.eps' % sys.argv[1])
plt.savefig('%s/yorp.pdf' % sys.argv[1])
outf.write('yorp %.15g %.15g\n' % (np.mean(f[:,3]), 3*np.std(f[:,3])))

plt.figure()
plt.hist(f[:,7],density=1,bins=20,color='skyblue')
plt.xlabel(r'$a$')
plt.ylabel(r'density')
plt.savefig('%s/fa.eps' % sys.argv[1])
plt.savefig('%s/fa.pdf' % sys.argv[1])
outf.write('fa %.15g %.15g\n' % (np.mean(f[:,7]), 3*np.std(f[:,7])))

plt.figure()
plt.hist(f[:,8],density=1,bins=20,color='skyblue')
plt.xlabel(r'$d$')
plt.ylabel(r'density')
plt.savefig('%s/fd.eps' % sys.argv[1])
plt.savefig('%s/fd.pdf' % sys.argv[1])
outf.write('fd %.15g %.15g\n' % (np.mean(f[:,8]), 3*np.std(f[:,8])))

plt.figure()
plt.hist(f[:,9],density=1,bins=20,color='skyblue')
plt.xlabel(r'$k$')
plt.ylabel(r'density')
plt.savefig('%s/fk.eps' % sys.argv[1])
plt.savefig('%s/fk.pdf' % sys.argv[1])
outf.write('fk %.15g %.15g\n' % (np.mean(f[:,9]), 3*np.std(f[:,9])))

plt.figure()
plt.hist(f[:,10],density=1,bins=20,color='skyblue')
plt.xlabel(r'$b$')
plt.ylabel(r'density')
plt.savefig('%s/fb.eps' % sys.argv[1])
plt.savefig('%s/fb.pdf' % sys.argv[1])
outf.write('fb %.15g %.15g\n' % (np.mean(f[:,10]), 3*np.std(f[:,10])))

plt.figure()
plt.hist(f[:,11][f[:,11]<10.0],density=1,bins=20,color='skyblue')
plt.xlabel(r'$c$')
plt.ylabel(r'density')
plt.savefig('%s/fc.eps' % sys.argv[1])
plt.savefig('%s/fc.pdf' % sys.argv[1])
outf.write('fc %.15g %.15g\n' % (np.mean(f[:,11][f[:,11]<10.0]), 3*np.std(f[:,11][f[:,11]<10.0])))
