import sys,os,numpy as np
sys.path.insert(1,'/home/dsaha3/Phx-hoomd/hoomd-blue/hoomd')
import hoomd, hoomd.md as md
from hoomd import azplugins
import gsd, gsd.hoomd

def extend(s,chain_length,spm_length):
	boxdim = s.configuration.box[:3]
	zmin,zmax,dz = -boxdim[2]/2., boxdim[2]/2., boxdim[2]
	pos1 =  s.particles.position
	pos = pos1.copy()
	natoms=len(s.particles.typeid)
	spm_ind=s.particles.types.index('SPM')
	i=0
	while i<natoms:
		if s.particles.typeid[i]==spm_ind: # this is spm
			c_length=spm_length
		else: # this is protein
			c_length=chain_length
		mol_coord=pos[i:i+c_length,2]
		for j in range(1,c_length):
			dist2=(mol_coord[j]-mol_coord[j-1])**2
			if dist2>8:
				excess = np.sign(mol_coord[j] - mol_coord[j-1])*dz
				mol_coord[j] = mol_coord[j] - excess
				com = np.mean(mol_coord)
				if com < zmin:
					mol_coord += dz
				elif com > zmax:
					mol_coord -= dz
		pos[i:i+c_length,2]=mol_coord
		i+=c_length
	
	return pos

import sys,os,numpy as np
sys.path.insert(1,'/home/dsaha3/Phx-hoomd/hoomd-blue/hoomd')
import hoomd, hoomd.md as md
from hoomd import azplugins
import gsd, gsd.hoomd

# input parameters
# protein name
fileroot=sys.argv[1]
# ion concentration in M
ion=float(sys.argv[2])
# epsilon in kcal/mol
epsilon_lj=4.184*float(sys.argv[3])
# temperature in K
temp = int(sys.argv[4])

n_spm  = int(sys.argv[5])                                                       
                                                                                
# spermine parameter                                                            
spm_b_side=0.48                                                                 
spm_b_center=0.57 
slab_z_length=280

# HPS-salt: hpssalt_a = 1.5; hpssalt_b = 0.06
hpssalt_a=1.5; hpssalt_b=0.06

# T-HPS: T-HPS_on=1 or 0
#thps_on=1
thps_on=1

kappa= 1/((0.1/ion)**0.5)
r_coul=35*0.1/kappa

seq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP','B':'ASPX','Z':'GLUX','J':'SPM'}

## get epsilon, salting out parameters from predetermined values
file_base = os.getcwd()

# 1.1 Read sequence

filein=fileroot
fileout='%s.seq3'%(fileroot)
nline=1
count=0
fout=open(fileout,'w')
with open(filein,'r') as fid:
    for i in fid:
        if i[0]!='#':
            for j in i:
                if j in seq:
                    fout.write(' %s'%seq[j])
                    count+=1
                    if count==nline:
                        fout.write('\n')
                        count=0
fout.close()

# #### 1.2 Read force field parameters

# Input parameters for all the amino acids (force field)
ff_para='''#AA     Mass    Charge  Sigma   Lambda  ks
Ala     71.08   0.00    5.040   0.730   -0.01000
Arg     156.20  1.00    6.560   0.000   -0.25528
Asn     114.10  0.00    5.680   0.432   -0.11500
Asp     115.10  -1.00   5.580   0.378   -0.09000
Cys     103.10  0.00    5.480   0.595   -0.14000
Gln     128.10  0.00    6.020   0.514   -0.09500
Glu     129.10  -1.00   5.920   0.459   -0.07000
Gly     57.05   0.00    4.500   0.649   -0.04000
His     137.10  0.50    6.080   0.514   -0.08050
Ile     113.20  0.00    6.180   0.973   0.08000
Leu     113.20  0.00    6.180   0.973   0.08000
Lys     128.20  1.00    6.360   0.514   -0.08050
Met     131.20  0.00    6.180   0.838   0.02968
Phe     147.20  0.00    6.360   1.000   0.07000
Pro     97.12   0.00    5.560   1.000   0.08477
Ser     87.08   0.00    5.180   0.595   -0.05295
Thr     101.10  0.00    5.620   0.676   -0.02541
Trp     186.20  0.00    6.780   0.946   0.07000
Tyr     163.20  0.00    6.460   0.865   0.07000
Val     99.07   0.00    5.860   0.892   0.06000
Aspx     115.10  -0.69   5.580   0.378   -0.09000
Glux     129.10  -0.36   5.920   0.459   -0.07000
Spm     51.50   1.00   5.040   0.514   0.00000'''

aalist={}
for i in ff_para.split('\n'):
    if i[0]!='#':
        tmp=i.rsplit()
        aalist[tmp[0].upper()]=np.loadtxt(tmp[1:],dtype=float)
aakeys=list(aalist.keys())
## This translates each amino acid type into a number, which will be used in HOOMD
## For example, GLY is with an ID of 10
#print('Gly has a number code of',aakeys.index('Gly'))

# dihedral force field K-B
ff_dih={}
with open('karanicolas_dihe_parm.dat', 'r') as fid:
	for i in fid:
		ls=i.split()
		aa1, aa2, k, m, d = ls[0],ls[1],float(ls[2]),int(ls[3]),float(ls[4])
		if aa1+aa2 not in ff_dih.keys():
			ff_dih[aa1+aa2]=[]
		ff_dih[aa1+aa2].append((k,m,d))
#dih_types=list(ff_dih.keys())
dih_types=list(ff_dih.keys()) + ['SPM_dihed'] 
	
# define temperature adjustments
def ExT(a, b, c, T):
    return a*(T*T) + (b*T) + c
def temp_adjust(aa, T):
    type = {'THR': 'P', 'GLU': 'C', 'ASP': 'C', 'LYS': 'C', 'ILE': 'H',
            'SER': 'P', 'ARG': 'C', 'ALA': 'H', 'VAL': 'H', 'GLY': 'O',
            'PRO': 'O', 'GLN': 'P', 'PHE': 'A', 'TYR': 'A', 'LEU': 'H',
            'HIS': 'A', 'ASN': 'P', 'CYS': 'O', 'MET': 'H', 'TRP': 'A', 'ASPX':'C','GLUX':'C',
            'SPM':'C'}
    parabola = {'H': (-0.00025597, 0.15379, -22.657), 'A': (-0.00026696, 0.15876, -23.364),
                'O': (0.000026, -0.015064, 2.1607), 'P': (0.0001201, -0.071482, 10.475),
                'C': (0.000093317, -0.057676, 8.5997)}
    a, b, c = parabola[type[aa]]
    alpha = 0.7836
    Tref = 296.7
    Tshift = 61.97
    ex1 = ExT(a, b, c, (T-Tshift))
    ex2 = ExT(a, b, c, (Tref-Tshift))
    return alpha*(ex1-ex2)

# adjust parameters
aamass=[]
aacharge=[]
aaradius=[]
aahps=[]
aasalt=[]
for i in aakeys:
    aamass.append(aalist[i][0])
    aacharge.append(aalist[i][1])
    aaradius.append(aalist[i][2])
    hps=aalist[i][3]
    aahps.append(hps)
    # grab amino acid salting out constant
    salt=aalist[i][4]
    aasalt.append(salt)
    # adjust hydrophobicity for salting out effect
    hps_sadjust= hpssalt_a * (hpssalt_b + salt) * (ion - 0.1)
    # adjust hydrophobicity for temperature
    hps_tadjust = temp_adjust(i, temp)*thps_on
    aalist[i][3] += hps_sadjust + hps_tadjust
print('Name of amino acids:',aakeys)
print('Mass of amino acids:',aamass)
print('Charge of amino acids:',aacharge)
print('Hydrophobicity of amino acids:',aahps)
print('Salting out constant of amino acids:', aasalt)

chain_id=[]
chain_mass=[]
chain_charge=[]

with open(fileout,'r') as fid:
    for i in fid:
        iname=i.rsplit()[0]
        chain_id.append(aakeys.index(iname))
        chain_mass.append(aalist[iname][0])
        chain_charge.append(aalist[iname][1])
print('Sequence coded in numbers:',chain_id)
print('Mass of sequence:',chain_mass)
print('Charge of sequence:',chain_charge)

bond_length=0.38  # nm
chain_length=len(chain_id)
spm_length=4
print('Chain length=',chain_length)
print('SPM length=',spm_length)
#atmnum = chain_length + n_spm * 4
f = gsd.hoomd.open(name='resize3.gsd', mode='rb')
s=f[0]
# WZ: I took out these two below, so now the box is not extended, but the Ashbaugh energy still explodes, which means it's not equilibrated after resizing
s.particles.position=extend(s,chain_length,spm_length)
s.configuration.box=[s.configuration.box[0],s.configuration.box[1],slab_z_length,0,0,0]
s.configuration.step = 0
outfile = gsd.hoomd.open('for-ext.gsd','wb')
outfile.append(s)
outfile.close()

# ## 2. Run HOOMD

hoomd.context.initialize("--mode=gpu")
system = hoomd.init.read_gsd('for-ext.gsd',time_step=0)

#Bonds
harmonic=hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond',k=8368,r0=bond_length)
# spm bonds                                                                     
harmonic.bond_coeff.set('SPM_bond_s',k=8368,r0=spm_b_side)                      
harmonic.bond_coeff.set('SPM_bond_c',k=8368,r0=spm_b_center)  

#Angle potentials
#ang_coeffs=hoomd.md.angle.coeff()
angle_pot=hoomd.md.angle.table(width=1000)
coeff_list_angle = (0.1, 106.4, 1.60, 4.3, 26.3, 2.27)
def angle_table_func(theta, coeff_list):
    cal2j=4.184
    (gamma, k_alpha, theta_a, epsilon_a, k_beta, theta_b) = coeff_list
    k_alpha *=cal2j
    k_beta *=cal2j
    gamma  /= cal2j
    epsilon_a *=cal2j
    alpha_value = k_alpha * (theta - theta_a) ** 2 + epsilon_a
    beta_value = k_beta * (theta - theta_b) ** 2
    alpha_exp = np.exp(-gamma * alpha_value)
    beta_exp = np.exp(-gamma * beta_value)
    sumval = alpha_exp + beta_exp
    E_theta = -np.log(sumval) / gamma
    # getting analytical derivative
    factor1=-2*gamma*k_alpha*(theta-theta_a)
    factor2=-2*gamma*k_beta*(theta-theta_b)
    num = factor1*alpha_exp + factor2*beta_exp
    deriv = (1/gamma) * (num/sumval)
    return E_theta, deriv
angle_pot.angle_coeff.set('AA_angle', func=angle_table_func, coeff=dict(coeff_list=coeff_list_angle))

# SPM angle potential                                                           
# new angle potential                                                           
coeff_list_angle2 = (200,1.6,1.2,2.8,3)                                         
def angle_table_func2(theta,coeff_list):                                        
    cal2j=4.184                                                                 
    (k1,k2,theta1,theta2,theta3)=coeff_list                                     
    k1*=cal2j                                                                   
    k2*=cal2j                                                                   
    if theta <0.5:                                                              
        E_theta=k1 * (theta+np.pi-theta3)**2 + k2*(theta3-theta2)**2            
        deriv= -2*k1*(theta+np.pi-theta3)                                       
    elif theta>=0.5 and theta<theta1:                                           
        E_theta=k1 * (theta-theta1)**2 + k2*(theta1-theta2)**2                  
        deriv= -2*k1*(theta-theta1)                                             
    elif theta>=theta1 and theta<theta3:                                        
        E_theta=k2 * (theta-theta2)**2                                          
        deriv= -2*k2*(theta-theta2)                                             
    elif theta>=theta3:                                                         
        E_theta=k1 * (theta-theta3)**2 + k2*(theta3-theta2)**2                  
        deriv= -2*k1*(theta-theta3)                                             
    return E_theta,deriv                                                        
angle_pot.angle_coeff.set('SPM_angle', func=angle_table_func2, coeff=dict(coeff_list=coeff_list_angle2))

#Dihedral potentials
dihedral1=hoomd.md.dihedral.harmonic()
dihedral2=hoomd.md.dihedral.harmonic()
dihedral3=hoomd.md.dihedral.harmonic()
dihedral4=hoomd.md.dihedral.harmonic()
dihedral5=hoomd.md.dihedral.harmonic() 
for i in dih_types[:-1]:
	para1,para2,para3,para4=ff_dih[i[:6]]
	dihedral1.dihedral_coeff.set(i,k=2*4.184*para1[0],d=1,n=para1[1],phi_0=para1[2]/180*np.pi)
	dihedral2.dihedral_coeff.set(i,k=2*4.184*para2[0],d=1,n=para2[1],phi_0=para2[2]/180*np.pi)
	dihedral3.dihedral_coeff.set(i,k=2*4.184*para3[0],d=1,n=para3[1],phi_0=para3[2]/180*np.pi)
	dihedral4.dihedral_coeff.set(i,k=2*4.184*para4[0],d=1,n=para4[1],phi_0=para4[2]/180*np.pi)
	dihedral5.dihedral_coeff.set(i,k=2*4.184*(-1.16),d=1,n=1,phi_0=297.35/180*np.pi)
dihedral1.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=1,phi_0=0)                   
dihedral2.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=2,phi_0=0)                   
dihedral3.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=3,phi_0=0)                   
dihedral4.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=4,phi_0=0)                   
dihedral5.dihedral_coeff.set('SPM_dihed',k=2*4.184*(-0.123),d=1,n=1,phi_0=0) 

# Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', '1-3', '1-4', 'body'])
#nl.reset_exclusions(exclusions=['1-2', 'body'])

# Pairwise interactions
# This is LJ Lambda potential
nb = azplugins.pair.ashbaugh(r_cut=1.5, nlist=nl)
for i in aakeys:
    for j in aakeys:
        nb.pair_coeff.set(i,j,lam=(aalist[i][3]+aalist[j][3])/2.,
                          epsilon=epsilon_lj, sigma=(aalist[i][2]+aalist[j][2])/10./2.,r_cut=(aalist[i][2]+aalist[j][2])/10./2.*3)    
        
# Electrostatics
# 8.98755e9*1e9*1.6e-19**2*6.02e23/1000./80.=1.73136
prefactor_yukawa=1.73136
yukawa = hoomd.md.pair.yukawa(r_cut=r_coul, nlist=nl)
for i in aakeys:
	for j in aakeys:
		charge=aalist[i][1]*aalist[j][1]
		if charge!=0.:
			print(i,j,charge)
			yukawa.pair_coeff.set(i,j, epsilon=prefactor_yukawa*aalist[i][1]*aalist[j][1], kappa=kappa, r_cut=r_coul)
		else:
			yukawa.pair_coeff.set(i,j, epsilon=0, kappa=kappa, r_cut=0)

## Group Particles
all = hoomd.group.all()

## Set up integrator
hoomd.md.integrate.mode_standard(dt=0.01)

#Temperature=298  # K
kTinput=temp * 8.3144598/1000.
#integrator = hoomd.md.integrate.nph(group=all, tauP=1.0,P=2,gamma=0.1,x=True,y=True,z=True,couple='none')
integrator = hoomd.md.integrate.langevin(group=all, kT=kTinput, seed=63535)
gamma=0.01
for i in aakeys:
    integrator.set_gamma(i, gamma=aalist[i][0]*gamma)

## Outputs
hoomd.analyze.log('extend.log', quantities=['bond_harmonic_energy','pair_ashbaugh_energy','pair_yukawa_energy','potential_energy','kinetic_energy','temperature','lx','ly','lz'], period=1000, overwrite=True, header_prefix='#')
hoomd.dump.gsd('extend.gsd', period=10000, group=all, truncate=True)
#old.dump.xml(all,fileroot+'.xml', period=1000000)
hoomd.dump.dcd('extend.dcd', period=10000, group=all, overwrite=True)

## Run simulation
hoomd.run(tsteps=2e7)
