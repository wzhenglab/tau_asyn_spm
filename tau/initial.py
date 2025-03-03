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
# number of spermine molecules
n_spm  = int(sys.argv[5])

# spermine parameter
spm_b_side=0.48
spm_b_center=0.57

# region with different epsilon
#region=[244,372] # start counting at 1
region=[244,372] # start counting at 1
epsilon_lj_region=4.184*float(sys.argv[6])
aaregion=['ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']

# HPS-salt: hpssalt_a = 1.5; hpssalt_b = 0.06
hpssalt_a=1.5; hpssalt_b=0.06

# T-HPS: T-HPS_on=1 or 0
thps_on=1

kappa= 1/((0.1/ion)**0.5)
r_coul=35*0.1/kappa
if r_coul>=5:
	r_coul=5

# added 'SPM' here
seq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP','J':'SPM'}

## get epsilon, salting out parameters from predetermined values
file_base = os.getcwd()

# 1.1 Read sequence

filein=fileroot
fileout='%s.seq3'%(fileroot)
nline=1
count=0
fout=open(fileout,'w')
seq3=[]
with open(filein,'r') as fid:
    for i in fid:
        if i[0]!='#':
            for j in i:
                if j in seq:
                    fout.write(' %s'%seq[j])
                    seq3.append(seq[j])
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
Spm     51.50   1.00   5.040   0.514   0.00000'''
# need to think about SPM parameters

aalist={}
for i in ff_para.split('\n'):
    if i[0]!='#':
        tmp=i.rsplit()
        aalist[tmp[0].upper()]=np.loadtxt(tmp[1:],dtype=float)

#aakeys=list(aalist.keys())
aakeys=list(set(seq3))
aakeys.append('SPM')

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
	
# define temperature adjustments
def ExT(a, b, c, T):
    return a*(T*T) + (b*T) + c
def temp_adjust(aa, T):
    type = {'THR': 'P', 'GLU': 'C', 'ASP': 'C', 'LYS': 'C', 'ILE': 'H',
            'SER': 'P', 'ARG': 'C', 'ALA': 'H', 'VAL': 'H', 'GLY': 'O',
            'PRO': 'O', 'GLN': 'P', 'PHE': 'A', 'TYR': 'A', 'LEU': 'H',
            'HIS': 'A', 'ASN': 'P', 'CYS': 'O', 'MET': 'H', 'TRP': 'A',
			'SPM': 'C'}
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

# #### 1.3 Build input gsd file for initial conformation and coordinate

## Read sequence
## Now we can translate the entire sequence into a number code according to the order in 'aakeys'
chain_id=[]
chain_mass=[]
chain_charge=[]

#with open(fileout,'r') as fid:
#    for i in fid:
for iname in seq3:
        #iname=i.rsplit()[0]
        chain_id.append(aakeys.index(iname))
        chain_mass.append(aalist[iname][0])
        chain_charge.append(aalist[iname][1])

chain_length=len(chain_id)
print('Protein chain length=',chain_length)

# add spermine into chain_id, chain_mass and chain_charge
iname='SPM'
for i in range(n_spm):
	for j in range(4):
		chain_id.append(aakeys.index(iname))
		chain_mass.append(aalist[iname][0])
		chain_charge.append(aalist[iname][1])	

print('Sequence coded in numbers:',chain_id)
print('Mass of sequence:',chain_mass)
print('Charge of sequence:',chain_charge)

bond_length=0.38  # nm
box_length=bond_length*chain_length+50

# Now we can build HOOMD structure
# data structure for one single frame
s=gsd.hoomd.Snapshot()

# add number of spermines
s.particles.N = chain_length + 4 * n_spm
print('Number of particles:',s.particles.N )

# Claim the name of amino acids
s.particles.types = aakeys
s.particles.typeid = chain_id
s.particles.mass = chain_mass
s.particles.charge = chain_charge

# Build initial position as a linear chain
pos=[]
cy=-(box_length/2)+0.01
cz=-(box_length/2)+0.01
for i in range(chain_length):
# Change the z-coordinate to have a #linear chain# zig zag pattern
    pos.append((-(box_length/2)+0.01,cy,cz))
    cy+=bond_length*3**0.5/2
    cz+=bond_length/2*(-1)**i
#pos=[]
#cy=-box_length/2+0.1
#cz=0
#for i in range(chain_length):
#    # Change the z-coordinate to have a #linear chain# zig zag pattern
#    pos.append((0,cy,cz))
#    cy+=bond_length*3**0.5/2
#    cz+=bond_length/2*(-1)**i

# build spermine initial coordinates
spm_template_y=[0,spm_b_side,spm_b_side+spm_b_center,spm_b_side*2+spm_b_center]
spm_gap=2
cx=5 # far away from protein
n_spm_y=10
for i in range(n_spm):
	cy=(i%n_spm_y)*spm_gap	
	cz=(i//n_spm_y)*spm_gap
	pos.append((cx,cy+spm_template_y[0],cz))
	pos.append((cx,cy+spm_template_y[1],cz))
	pos.append((cx,cy+spm_template_y[2],cz))
	pos.append((cx,cy+spm_template_y[3],cz))
pos=np.array(pos)
s.particles.position= pos

# Initialize bond
# add spermine bonds, each with three bonds
nbonds_protein=max(chain_length-1,0)
nbonds_spm=n_spm *3
nbonds=nbonds_protein+nbonds_spm
s.bonds.N = nbonds
# add new bond types, s means side bond and c means center bond
s.bonds.types = ['AA_bond','SPM_bond_s','SPM_bond_c']
bond_typeid=[0]*nbonds_protein
for i in range(n_spm):
    bond_typeid.append(1)
    bond_typeid.append(2)
    bond_typeid.append(1)
s.bonds.typeid = bond_typeid
bond_pairs=np.zeros((nbonds,2),dtype=int)
# protein bond pairs
for i in range(nbonds_protein):
	bond_pairs[i,:]=np.array([i,i+1])
for i in range(n_spm):
    start_i=chain_length + i*4
    bond_pairs[nbonds_protein+i*3,:]= np.array([start_i,start_i+1])
    bond_pairs[nbonds_protein+i*3+1,:]= np.array([start_i+1,start_i+2])
    bond_pairs[nbonds_protein+i*3+2,:]= np.array([start_i+2,start_i+3])
s.bonds.group = bond_pairs

# Initialize angles
nangle_protein=max(chain_length-2,0)
nangle_spm=n_spm*2
nangle = nangle_protein+nangle_spm
s.angles.N = nangle
s.angles.types = ['AA_angle','SPM_angle']
s.angles.typeid = [0]*(nangle_protein) + [1]*(nangle_spm)
angle_trios = np.zeros((nangle,3),dtype=int)
for i in range(nangle):
    angle_trios[i,:]=np.array([i, i+1, i+2])
for i in range(n_spm):
	start_i=chain_length + i*4
	angle_trios[nangle_protein+i*2,:]= np.array([start_i,start_i+1,start_i+2])
	angle_trios[nangle_protein+i*2+1,:]= np.array([start_i+1,start_i+2,start_i+3])
s.angles.group = angle_trios

# Initialize dihedral quads
ndihed_protein = max(chain_length-3,0)
ndihed_spm= n_spm
ndihed = ndihed_protein + ndihed_spm
s.dihedrals.N = ndihed
dih_types=list(ff_dih.keys()) + ['SPM_dihed']
dih_typeid=np.zeros(ndihed)
dihed_quads = np.zeros((ndihed,4),dtype=int)
for i in range(ndihed_protein):
        dihed_quads[i,:] = np.array([i, i+1, i+2, i+3])
        ckey=aakeys[s.particles.typeid[i+1]][:3]+aakeys[s.particles.typeid[i+2]][:3]
        dih_typeid[i]=dih_types.index(ckey)
for i in range(ndihed_spm):
        start_i=chain_length + i*4
        dihed_quads[ndihed_protein+i,:]=np.array([start_i,start_i+1,start_i+2,start_i+3])
        dih_typeid[ndihed_protein+i]=dih_types.index('SPM_dihed')
s.dihedrals.types=dih_types
s.dihedrals.typeid=dih_typeid
s.dihedrals.group = dihed_quads

# dual-epsilon

for i in range(chain_length):
	if i>=region[0]-1 and i<=region[1]-1:
		oldkey=aakeys[chain_id[i]]
		if oldkey in aaregion:
			newkey=oldkey+'R'
			if newkey not in aakeys:
				aakeys.append(newkey)
			chain_id[i]=aakeys.index(newkey)
s.particles.typeid = chain_id

print('Number of types of particles:',len(aakeys))

# Box size
s.configuration.dimensions=3
s.configuration.box=[box_length,box_length,box_length,0,0,0]
s.configuration.step=0
# Write gsd file
f = gsd.hoomd.open(name='start.gsd', mode='wb')
f.append(s)
f.close()

# write pdb
fout=open('%s.pdb'%(fileroot),'w')
for i in range(s.particles.N):
    fout.write('ATOM  %5i  CA  %s %5i    %8.3f%8.3f%8.3f  1.00  0.00\n'%(i+1,s.particles.types[s.particles.typeid[i]],i+1,s.particles.position[i,0]*10,s.particles.position[i,
1]*10,s.particles.position[i,2]*10))
fout.write('END')

# ## 2. Run HOOMD

hoomd.context.initialize("--mode=gpu")
system = hoomd.init.read_gsd('start.gsd',time_step=0)
system.replicate(nx=2,ny=2,nz=2)

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
# we turn off dihedral 1 to 4 and only use dihedral 5 for SPM_dihed
dihedral1.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=1,phi_0=0)
dihedral2.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=2,phi_0=0)
dihedral3.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=3,phi_0=0)
dihedral4.dihedral_coeff.set('SPM_dihed',k=0,d=1,n=4,phi_0=0)
dihedral5.dihedral_coeff.set('SPM_dihed',k=2*4.184*(-0.123),d=1,n=1,phi_0=0)

# Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', '1-3', '1-4','body'])

# Pairwise interactions
# This is LJ Lambda potential
nb = azplugins.pair.ashbaugh(r_cut=1.5, nlist=nl)
for i in aakeys:
	for j in aakeys:
		if len(i)==4 and len(j)==4 and i[3]=='R' and j[3]=='R':
			epsilon_current=epsilon_lj_region
		else:
			epsilon_current=epsilon_lj
		nb.pair_coeff.set(i,j,lam=(aalist[i[:3]][3]+aalist[j[:3]][3])/2.,epsilon=epsilon_current, sigma=(aalist[i[:3]][2]+aalist[j[:3]][2])/10./2.,r_cut=(aalist[i[:3]][2]+aalist[j[:3]][2])/10./2.*3)    

# Electrostatics
# 8.98755e9*1e9*1.6e-19**2*6.02e23/1000./80.=1.73136
prefactor_yukawa=1.73136
yukawa = hoomd.md.pair.yukawa(r_cut=r_coul, nlist=nl)
yukawa.pair_coeff.set(aakeys,aakeys, epsilon=0.0,kappa=kappa,r_cut=0)
for i in aakeys:
	for j in aakeys:
		charge=aalist[i[:3]][1]*aalist[j[:3]][1]
		if charge!=0.:
			yukawa.pair_coeff.set(i,j, epsilon=prefactor_yukawa*aalist[i[:3]][1]*aalist[j[:3]][1], kappa=kappa, r_cut=r_coul)
		#else:
		#	yukawa.pair_coeff.set(i,j, epsilon=0, kappa=kappa, r_cut=0)

## Group Particles
all = hoomd.group.all()

## Set up integrator
hoomd.md.integrate.mode_standard(dt=0.005) # Time units in ps
#Temperature=298  # K
kTinput=temp * 8.3144598/1000.
integrator = hoomd.md.integrate.langevin(group=all, kT=kTinput, seed=61535)
gamma=0.01
for i in aakeys:
	integrator.set_gamma(i, gamma=aalist[i[:3]][0]*gamma)

## Outputs
hoomd.analyze.log('init.log', quantities=['bond_harmonic_energy','pair_ashbaugh_energy','pair_yukawa_energy','potential_energy','kinetic_energy','temperature','lx','ly','lz'], period=10000, overwrite=True, header_prefix='#')
hoomd.dump.gsd('restart_init1.gsd', period=1000000, group=all, truncate=True)    
hoomd.dump.gsd('restart_init2.gsd', period=1000000, group=all, truncate=True, phase=500000)
resize=hoomd.update.box_resize(Lx=hoomd.variant.linear_interp([(0, system.box.Lx), (2e7, 40)]), Ly=hoomd.variant.linear_interp([(0, system.box.Ly), (2e7, 40)]),Lz=hoomd.variant.linear_interp([(0, system.box.Lz), (2e7, 40)]) )

hoomd.dump.dcd('init.dcd', period=10000, group=all, overwrite=True)
hoomd.dump.gsd('init.gsd', period=10000, group=all, truncate=True, overwrite=True)

hoomd.run(tsteps=2e7)
