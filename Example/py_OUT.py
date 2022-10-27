#--------------------------------------------------Import function
from datetime import datetime
#------------------------------------------------------------------
start=datetime.now()
import matplotlib.pyplot as plt
import numpy as np

print('\n************************************ Pyhon ')
#Structure of data 
data = np.fromfile('data_structure.dat', dtype=np.int32, count=6)
#print(data)
EDO_len=data[0]
TxD_len=data[1]
N_Step=data[2]
N_POP=data[3]
Number_EDOTest=data[4]
Number_TxDTest=data[5]
print('EDO_len:',EDO_len)
print('TxD_len:',TxD_len)
print('N_Step:',N_Step)
print('N_POP:',N_POP)
print('Number_EDOTest:',Number_EDOTest)
print('Number_TxDTest:',Number_TxDTest)


data_EDO=np.zeros((EDO_len,3))
data = np.fromfile('X_EDO.dat', dtype=np.float64, count=EDO_len*3)
data_EDO[:,0] = data.data[0:EDO_len]
data_EDO[:,1] = data.data[EDO_len:2*EDO_len]
data_EDO[:,2] = data.data[2*EDO_len:3*EDO_len]

data_TxD=np.zeros((TxD_len,4))
data = np.fromfile('X_TxD.dat', dtype=np.float64, count=TxD_len*4)
data_TxD[:,0] = data.data[0:TxD_len]
data_TxD[:,1] = data.data[TxD_len:2*TxD_len]
data_TxD[:,2] = data.data[2*TxD_len:3*TxD_len]
data_TxD[:,3] = data.data[3*TxD_len:4*TxD_len]
#print('data_TxD\n',data_TxD)

data_HX_EDO=np.zeros((N_Step,3,N_POP,Number_EDOTest))
data= np.fromfile('HX_EDO.dat', dtype=np.float64, count=N_Step*3*N_POP*Number_EDOTest)
q=0
w=1
e=2
for k in range(0,Number_EDOTest):
	for i in range(0,N_POP):
		data_HX_EDO[:,0,i,k] =data[q*N_Step:(q+1)*N_Step]
		data_HX_EDO[:,1,i,k] =data[w*N_Step:(w+1)*N_Step]
		data_HX_EDO[:,2,i,k] =data[e*N_Step:(e+1)*N_Step]
		q=q+3
		w=w+3
		e=e+3

data_RX_TxD=np.zeros((N_Step,3,N_POP,Number_TxDTest))
data= np.fromfile('RX_TxD.dat', dtype=np.float64, count=N_Step*3*N_POP*Number_TxDTest)
q=0
w=1
e=2
for k in range(0,Number_TxDTest):
	for i in range(0,N_POP):
		data_RX_TxD[:,0,i,k] =data[q*N_Step:(q+1)*N_Step]
		data_RX_TxD[:,1,i,k] =data[w*N_Step:(w+1)*N_Step]
		data_RX_TxD[:,2,i,k] =data[e*N_Step:(e+1)*N_Step]
		q=q+3
		w=w+3
		e=e+3

print('Read data 1 *********************************1/2')

# =============================================================================
# #----------------------------------------------------------------- PLOT    
## =============================================================================

N_c_plot=1 # 1 - N_POP numbers of output curve plot

# parameters for plotting
# aspect ration with the golden number looks good
phi = 1.55
# parameters for the size and font 
szFig = 6
ftSze = 12
lnWdth = 2
plt.rcParams.update({'font.size': ftSze-2,
                     'legend.fontsize': ftSze})
plt.rcParams["text.usetex"] =True
plt.rcParams["mathtext.fontset"] = "cm"

#EDO----------------------------------------
plt.figure(num=1,figsize=(szFig,szFig/phi))
for k in range(0,Number_EDOTest):
	for i in range(0,N_c_plot): # range(0,N_POP)
	   plt.plot(-data_HX_EDO[:,0,i,k],data_HX_EDO[:,2,i,k],'--k')
for k in range(1,int(max(data_EDO[:,2])+1)):
	mask=data_EDO[:,2]==k
	mask_data_EDO=data_EDO[mask]
	plt.scatter(mask_data_EDO[:,1],mask_data_EDO[:,0],s=25,label=k) 
	# legend, labels etc
#plt.xlim(0, 1.25)
#plt.ylim(0.62, 0.75)
plt.legend()
plt.grid(linestyle='dotted')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.ylabel('$ e $  (-)',fontsize=ftSze)
plt.xlabel('$\sigma_v$  (MPa)',fontsize=ftSze)
# save the figure
name = 'OE.pdf'
plt.savefig(name)


# TxD---------------------------------------
# plotting
plt.figure(num=2,figsize=(szFig,szFig/phi))
for k in range(0,Number_TxDTest):
	for i in range(0,N_c_plot): # range(0,N_POP)
	   plt.plot(data_RX_TxD[:,2,i,k]*100,data_RX_TxD[:,0,i,k],'--k')
for k in range(1,int(max(data_TxD[:,3])+1)):
	mask=data_TxD[:,3]==k
	mask_data_TxD=data_TxD[mask]
	#print(mask_data_EDO)
	plt.scatter(mask_data_TxD[:,0],mask_data_TxD[:,2],s=25,label=k) 
	# legend, labels etc
#plt.xlim(0, 1.25)
#plt.ylim(0.62, 0.75)
plt.legend()
plt.grid(linestyle='dotted')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.xlabel(r'$\varepsilon_a$  (\%)',fontsize=ftSze)
plt.ylabel('$q$  (MPa)',fontsize=ftSze)
# save the figure
name = 'TD_q.pdf'
plt.savefig(name)

# plotting
plt.figure(num=3,figsize=(szFig,szFig/phi))
for k in range(0,Number_TxDTest):
	for i in range(0,N_c_plot): # range(0,N_POP)
	   plt.plot(data_RX_TxD[:,2,i,k]*100,data_RX_TxD[:,1,i,k]*100,'--k')
for k in range(1,int(max(data_TxD[:,3])+1)):
	mask=data_TxD[:,3]==k
	mask_data_TxD=data_TxD[mask]
	#print(mask_data_EDO)
	plt.scatter(mask_data_TxD[:,0],mask_data_TxD[:,1],s=25,label=k) 
	# legend, labels etc
#plt.xlim(0, 1.25)
#plt.ylim(0.62, 0.75)
plt.legend()
plt.grid(linestyle='dotted')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel(r'$\varepsilon_a$  (\%)',fontsize=ftSze)
plt.ylabel(r'$\varepsilon_v$  (\%)',fontsize=ftSze)
# save the figure
name = 'TD_epv.pdf'
plt.savefig(name)

print('Plot data  *********************************2/2')
print('Execution time:',datetime.now()-start)

