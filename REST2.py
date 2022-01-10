#!/usr/bin/env python3
import os
import numpy as np
import pathlib
import re
import glob




def run_md(pdbname, scale, loop=0):
	"""
	Let molecule in pdb files go into the equilibrium state through em, nvt and npt simulations. The boxes information, solvent, ions are added too.
	All initial structures and walkers have the exact same solvent number, ion number, ion type and box size. Concentration of saline is set as 0.15M.
	For this purpose, we record the information of the first structure as the tamplate.
	"""
	global num_sol, box_size, num_Na, num_Cl
	print('echo -e "1\n1\n" | gmx pdb2gmx -f %s -o processed.gro -ignh	> grompp.log 2>&1' % pdbname)
	os.system('echo -e "1\n1\n" | gmx pdb2gmx -f %s -o processed.gro -ignh	> grompp.log 2>&1' % pdbname)

	if loop == 0:
		print('gmx editconf -f processed.gro -o newbox.gro -d 1 -c -bt dodecahedron')
		os.system(
			'gmx editconf -f processed.gro -o newbox.gro -d 1 -c -bt dodecahedron')
		print('gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top > sol.log 2>&1')
		os.system(
			'gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top > sol.log 2>&1')
		with open('solv.gro', 'r') as sol_gro:
			line =sol_gro.readlines()
			info=line[-1].split()
			box_size = [float(k)+0.10000 for k in info[0:3]]
		print(box_size)
		print('=================')
		with open('topol.top', 'r') as top:
			for line in top.readlines():
				line_sp = line.split()
				if line_sp == []:
					continue
				if line.split()[0] == 'SOL' and line_sp[1].isdigit():
					num_sol = line_sp[1]
		print('Max number of solvents is:', num_sol)
		os.system(
			'gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 > grompp_ion.log 2>&1')
		os.system(
			'echo -e "13\n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15')
		with open('topol.top', 'r') as top:
			for line in top.readlines():
				line_sp = line.split()
				if line_sp == []:
					continue
				if line.split()[0] == 'NA':
					num_Na = line_sp[1]
				if line.split()[0] == 'CL':
					num_Cl = line_sp[1]
		with open('../box_information.txt', 'w') as box_info:
			box_info.write('num_sol={}\nbox_size={},{},{}\nnum_Na={}\nnum_Cl={}'.format(
				num_sol, box_size[0], box_size[1], box_size[2], num_Na, num_Cl))
		print(box_size)
	else:
		print('gmx editconf -f processed.gro -o newbox.gro -box {} {} {} -c -bt dodecahedron'.format(
			box_size[0], box_size[1], box_size[2]))
		os.system('gmx editconf -f processed.gro -o newbox.gro -box {} {} {} -c -bt dodecahedron'.format(
			box_size[0], box_size[1], box_size[2]))
		print('gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top > sol.log 2>&1')
		os.system(
			'gmx solvate -cp newbox.gro -cs spc216.gro -maxsol {} -o solv.gro -p topol.top > sol.log 2>&1'.format(int(num_sol)))

		with open('topol.top', 'r') as top:
			for line in top.readlines():
				line_sp = line.split()
				if line_sp == []:
					continue
				if line.split()[0] == 'SOL' and line_sp[1].isdigit():
					print('Max number of solvents is:', line_sp[1])

		os.system(
			'gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 > grompp_ion.log 2>&1')
		os.system('echo -e "13\n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -np {} -nn {}'.format(num_Na, num_Cl))

	os.system('gmx grompp -f mini.mdp -c solv_ions.gro -p topol.top -pp	 > grompp_em.log 2>&1')

	count=0
	with open('processed.top', 'r') as top:
		for line in top.readlines():
			count+=1
			line_sp = line.split()
			if "residue 355" in line:
				line_seek=count
				break

	with open('processed_hotatom.top', 'w') as fw:
		count_line=0
		with open('processed.top', 'r') as top:
			for line in top.readlines():
				count_line+=1
				line_sp = line.split()
				if (line_seek < count_line <(line_seek+1567)):
					line_hotatom=line.replace(line_sp[1],line_sp[1]+'_',1)
				elif ('Include Position restraint file' in line):
					line_hotatom=line.strip('\n')+'''
#ifdef POSRES
#include "posre.itp"
#endif'''
				else:line_hotatom=line
				print(line_hotatom.strip('\n'),file=fw)
				
	os.system('plumed partial_tempering {} < processed_hotatom.top > topol_scaled.top'.format(scale))
	
	os.system('gmx grompp -f mini.mdp -c solv_ions.gro -p topol_scaled.top -o em.tpr -maxwarn 1')
	os.system('gmx mdrun -v -deffnm em')
	
	os.system('gmx grompp -f nvt.mdp -c em.gro -p topol_scaled.top -o nvt.tpr -r em.gro -maxwarn 1')
	os.system('gmx mdrun -deffnm nvt -ntmpi 1 -v -nt 4')

	os.system('gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol_scaled.top -o npt.tpr -r nvt.gro -maxwarn 1')
	os.system('gmx mdrun -deffnm npt -ntmpi 1 -v -nt 4')
	
	os.system('gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol_scaled.top -o md.tpr -r npt.gro -maxwarn 1')


if __name__ == '__main__':
	
	temp=np.loadtxt('temperature.dat',delimiter=',')
	num_sol = None
	box_size = []
	num_Na, num_Cl = None, None
	ccc=''
	for i,t in enumerate(temp):
		pp_dir="replica%d" %i
		pdbname='replica%d.pdb' %i
		pathlib.Path(pp_dir).mkdir(parents=True, exist_ok=True)
		os.chdir(pp_dir)
		os.system('cp ../%s ./' %pdbname )
		os.system('cp ../*.mdp ./' )
		os.system('cp -r ../charmm36.ff ./')
		os.system('cp ../plumed.dat ./' )
		scaling=300/float(t)
		print('temperature:{} ;scale:{}'.format(t,scaling))
		
		run_md(pdbname,scaling,loop=i)
		os.chdir('..')
		ccc+='c'+str(i)+' '
	print(ccc)
