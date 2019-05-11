'''
Created on 2018 12 8 

@author: Andrew
'''
import numpy as np
import re
#import quickSort as qSort
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

class TBmodel(object):
    Ham=[]
    Ham_k=[]
    lattice=[]
    norbital=0
    orbital_coor=[]
    onsite_energy=[]
    
    nhoppings=0
    hopping=[]
    
    kpath=[]
    kmesh=[]
    
    wcc_path=[]
    wcc_polar_direction=-1
    nelectrons=-1
    
    def __init__(self):
        self.readInputs()
    
    def readInputs(self):
        try:
            f=open('POSCAR','r')
        except IOError as e:
            print e
            exit()
        data=f.read()
        f.close()
        data=[line for line in data.split('\n') if line]
        
        # read the position of each tag
        positionOfOrbital=-1
        positionOfHopping=-1
        positionOfKPath=-1
        positionOfKMesh=-1
        positionOfWCC=-1
        for index, line in enumerate(data):
            if(re.findall('[a-zA-Z]+',line)):
                if(re.findall('[a-zA-Z]+',line)[0]=='Direct'): positionOfOrbital=index+1
                if(re.findall('[a-zA-Z]+',line)[0]=='Hopping'):positionOfHopping=index+1
                if(re.findall('[a-zA-Z]+',line)[0]=='KPath'):positionOfKPath=index+1
                if(re.findall('[a-zA-Z]+',line)[0]=='KMesh'):positionOfKMesh=index+1
                if(re.findall('[a-zA-Z]+',line)[0]=='WCC'):positionOfWCC=index+1
        
        # read lattice
        length_factor=float(data[1])
        for i in range(2,5):
            ax,ay,az=re.findall('[\-0-9.eE]+',data[i])
            self.lattice.append([float(ax),float(ay),float(az)])
        self.lattice=length_factor*np.array(self.lattice)
        
        # read number of orbitals
        for norb in re.findall('[\-0-9.eE]+',data[6]):
            self.norbital+=int(norb)
        
        # read orbital coordinates and onsite energies
        for i in range(positionOfOrbital,positionOfOrbital+self.norbital):
            ox,oy,oz,oe=re.findall('[\-0-9.eE]+',data[i])
            self.orbital_coor.append([float(ox),float(oy),float(oz)])
            self.onsite_energy.append(float(oe))
        self.orbital_coor=np.array(self.orbital_coor)
        self.onsite_energy=np.array(self.onsite_energy)
            
        # read hoppings
        # firstly, construct the architecture
        for origin_orbital_index in range(self.norbital):
            hopping_tmp=[]
            for destiny_orbital_index in range(self.norbital):
                hopping_tmp.append([])
            self.Ham.append(hopping_tmp)
        
        # then read all hoppings
        self.nhoppings=int(re.findall('[0-9]+',data[positionOfHopping])[0])
        #print self.nhoppings
        
        for ihopping in range(positionOfHopping+2,positionOfHopping+2+self.nhoppings):
            iorigin,idestiny,Rx,Ry,Rz,amplify_str=re.findall('[\-0-9.eEi+]+',data[ihopping])
            R=np.array([int(Rx),int(Ry),int(Rz)])
            # convert amplify
            amplify_list=[str for str in re.findall('[\-0-9.eEi]+',amplify_str)]
            amplify_real=0.
            amplify_imag=0.
            if(len(amplify_list)==2):
                amplify_real=float(amplify_list[0])
                amplify_imag=float(re.findall('[\-0-9.eE]+',amplify_list[1])[0])
            else:
                amplify_imag_check=[str for str in re.findall('[i]',amplify_str)]
                if(amplify_imag_check):
                    amplify_imag=float(re.findall('[\-0-9.eE]+',amplify_list[0])[0])
                else:
                    amplify_real=float(amplify_list[0])
                    
            amplify=amplify_real+1j*amplify_imag
            self.Ham[int(iorigin)][int(idestiny)].append([R,amplify])
            
        # read KPath
        if(positionOfKPath>0):
            nHighSymKpts=int(re.findall('[0-9]+',data[positionOfKPath])[0])
            highSymKpts_list=[]
            # read high symmetry k points
            for ihsk in range(positionOfKPath+1,positionOfKPath+1+nHighSymKpts):
                kx,ky,kz,ninter=re.findall('[\-0-9.eE]+',data[ihsk])
                highSymKpts_list.append([[float(kx),float(ky),float(kz)],int(ninter)])
                
            # generate k path
            for ihsk in range(1,nHighSymKpts):
                kx_list=np.linspace(highSymKpts_list[ihsk-1][0][0],
                                    highSymKpts_list[ihsk][0][0],
                                    highSymKpts_list[ihsk-1][1])
                ky_list=np.linspace(highSymKpts_list[ihsk-1][0][1],
                                    highSymKpts_list[ihsk][0][1],
                                    highSymKpts_list[ihsk-1][1])
                kz_list=np.linspace(highSymKpts_list[ihsk-1][0][2],
                                    highSymKpts_list[ihsk][0][2],
                                    highSymKpts_list[ihsk-1][1])
                k_list=zip(kx_list,ky_list,kz_list)
                # delete the end point
                k_list.pop(highSymKpts_list[ihsk-1][1]-1)
                self.kpath.extend(k_list)
                
        # read Kmesh
        if(positionOfKMesh>0):
            nkx, nky, nkz=[int(i) for i in re.findall('[0-9]+',data[positionOfKMesh])]
                
        # read WCC
        if(positionOfWCC>0):
            wccPolarDirection, integrationStep, occ=re.findall('[0-9]+',data[positionOfWCC])
            nHighSymKpts=int(re.findall('[0-9]+',data[positionOfWCC+2])[0])
            self.wcc_polar_direction=int(wccPolarDirection)
            self.nelectrons=int(occ)
            integrationStep=int(integrationStep)
            stepLength=1./integrationStep
            
            highSymKpts_list=[]
            for ihsk in range(positionOfWCC+3,positionOfWCC+3+nHighSymKpts):
                #print data[ihsk]
                kx,ky,kz,ninter=re.findall('[\-0-9.eE]+',data[ihsk])
                highSymKpts_list.append([[float(kx),float(ky),float(kz)],int(ninter)])
                
            kpath=[]
            for ihsk in range(1,nHighSymKpts):
                kx_list=np.linspace(highSymKpts_list[ihsk-1][0][0],
                                    highSymKpts_list[ihsk][0][0],
                                    highSymKpts_list[ihsk-1][1])
                ky_list=np.linspace(highSymKpts_list[ihsk-1][0][1],
                                    highSymKpts_list[ihsk][0][1],
                                    highSymKpts_list[ihsk-1][1])
                kz_list=np.linspace(highSymKpts_list[ihsk-1][0][2],
                                    highSymKpts_list[ihsk][0][2],
                                    highSymKpts_list[ihsk-1][1])
                k_list=zip(kx_list,ky_list,kz_list)
                k_list.pop(highSymKpts_list[ihsk-1][1]-1)
                kpath.extend(k_list)
            kpath.append(kpath[0])
            
            for ik in kpath:
                #print ik
                polar_k_line=[]
                for ik_online in range(integrationStep):
                    kpt_tmp=[ik[0],ik[1],ik[2]]
                    kpt_tmp[self.wcc_polar_direction]=ik_online*stepLength
                    polar_k_line.append(np.array(kpt_tmp))
                    #print polar_k_line
                self.wcc_path.append(polar_k_line)
            
    def constructHk(self,kpt=[0.,0.,0.]):
        kpt=np.array(kpt)
        Ham_k=[]
        for ibloch in range(self.norbital):
            Hk_tmp=[]
            for jbloch in range(self.norbital):
                matrix_element=0.
                for hopping in self.Ham[ibloch][jbloch]:
                    matrix_element+=hopping[1]*np.exp(-2j*np.pi*np.dot(kpt,hopping[0]),dtype=complex)
                    #print ibloch, jbloch, hopping
                    #print kpt, hopping[0], np.exp(2j*np.pi*np.dot(kpt,hopping[0]),dtype=complex)
                Hk_tmp.append(matrix_element)
            Ham_k.append(Hk_tmp)
        Ham_k=np.array(Ham_k)+self.onsite_energy*np.eye(self.norbital,dtype=complex)
        
        #print Ham_k
        return Ham_k
    
    def solveHk(self,kpt=[0.,0.,0.],return_orb=False):
        Ham_k=self.constructHk(kpt=kpt)
        #print Ham_k
        eig, vec = np.linalg.eigh(Ham_k,'U')
        
        if return_orb:
            return eig, vec
        else:
            return eig
    
    def plotbands(self):
        fig = plt.figure()
        for ikpt, kpt in enumerate(self.kpath):
            eig_list = self.solveHk(kpt=kpt)
            #print kpt, eig_list
            for eig in eig_list:
                plt.scatter(ikpt,eig,color='black')
        plt.show()
    
    def calcOverlap(self,vec_k0,vec_k1,b,dim=-1):
        '''
        Core algorithm for overlap calculation
        Ref: Z2Pack: Numerical implementation of hybridWannier centers for identifying topological materials
             PHYSICAL REVIEW B 95, 075146 (2017)
        page 18, equation (D3)
        '''
        dim = self.nelectrons if dim < 0 else dim
        br=np.array([np.exp(-2j*np.pi*np.dot(b,r),dtype=complex) for r in self.orbital_coor])
        
        #print 'v0', vec_k0
        #print 'v1', vec_k1
        OverlapMatrix=np.eye(dim,dtype=complex)
        for m in range(dim):
            for n in range(dim):
                #print vec_k0[m]
                #print m, n
                OverlapMatrix[m][n]=sum(vec_k0[m].conjugate()*vec_k1[n]*br)
        
        #print OverlapMatrix
        return OverlapMatrix
        
    
    def calcWilsonLoop(self,kpts=[],bandset=-1):
        '''
        Wilson loop algorithm
        '''
        D_overlapMatrix=self.nelectrons if bandset < 0 else 1
        #print D_overlapMatrix
        vec_list=[]
        for ikpt in kpts:
            eig, vec = self.solveHk(ikpt,True)
            idx=np.real(eig).argsort()
            idx=idx[:self.nelectrons] if bandset < 0 else [idx[bandset]]
            vec=vec[:, idx].transpose()
            vec_list.append(np.array(vec))
            
        vec_list.append(vec_list[0])
        
        Wilson_loop=np.eye(D_overlapMatrix,dtype=complex)
        for ikpt in range(len(kpts)):
            b=kpts[ikpt+1]-kpts[ikpt] if ikpt<(len(kpts)-1) else kpts[0]-kpts[ikpt]+1
            overlapMatrix=self.calcOverlap(vec_list[ikpt],vec_list[ikpt+1],b,dim=D_overlapMatrix)
            Wilson_loop=np.dot(Wilson_loop,overlapMatrix)
        
        eig, _ = np.linalg.eig(Wilson_loop)
        return [(1j*np.log(z,dtype=complex)).real/(2*np.pi) %1 for z in eig]
    
    def calcParallelTransport(self,kpts=[],bandset=-1):
        '''
        Parallel transport algorithm
        Ref: Maximally localized generalized Wannier functions for composite energy bands
             PHYSICAL REVIEW B 56, 20 (1997)
        page 12853, the highlight sentences
        '''
        D_overlapMatrix=self.nelectrons if bandset < 0 else 1
        vec_list=[]
        for ikpt in kpts:
            eig, vec =self.solveHk(ikpt,True)
            sortIndex=eig.argsort()
            sortIndex=sortIndex[:self.nelectrons] if bandset < 0 else [sortIndex[bandset]]
            vec=np.array(vec)
            vec=vec[:,sortIndex].transpose()
            vec_list.append(vec)
        vec_list.append(vec_list[0])
        
        VW_stack=np.eye(D_overlapMatrix,dtype=complex)
        for ikpt in range(len(kpts)):
            b=kpts[ikpt+1]-kpts[ikpt] if ikpt<(len(kpts)-1) else 1+kpts[0]-kpts[ikpt]
            overlap_matrix=self.calcOverlap(vec_list[ikpt],vec_list[ikpt+1],b,dim=D_overlapMatrix)
            #print overlap_matrix
            V,E,W=np.linalg.svd(overlap_matrix)
            VW_stack=np.dot(VW_stack,np.dot(V,W))
        
        eig, _=np.linalg.eig(VW_stack)
        wcc=np.array([(1j*np.log(z,dtype=complex)).real/(2*np.pi) %1 for z in eig])
        return wcc
    
    def calcGradualDecompose(self,kpts=[],bandset=-1):
        '''
        GradualDecompose algorithm
        the most expensive way, just for fun
        '''
        D_overlapMatrix=self.nelectrons if bandset < 0 else 1
        #print D_overlapMatrix
        vec_list=[]
        for ikpt in kpts:
            eig, vec =self.solveHk(ikpt,True)
            sortIndex=eig.argsort()
            sortIndex=sortIndex[:self.nelectrons] if bandset < 0 else [sortIndex[bandset]]
            vec=np.array(vec)
            vec=vec[:,sortIndex].transpose()
            vec_list.append(vec)
        vec_list.append(vec_list[0])
        
        wcc=np.linspace(0,0,D_overlapMatrix)
        for ikpt in range(len(kpts)):
            b=kpts[ikpt+1]-kpts[ikpt] if ikpt<(len(kpts)-1) else 1+kpts[0]-kpts[ikpt]
            #print vec_list[ikpt]
            #print vec_list[ikpt+1]
            overlap_matrix=self.calcOverlap(vec_list[ikpt],vec_list[ikpt+1],b,dim=D_overlapMatrix)
            #print overlap_matrix
            eig, _=np.linalg.eig(overlap_matrix.transpose().conjugate())
            wcc+=np.array([(1j*np.log(z,dtype=complex)).real for z in eig])
            #print wcc
        
        #eig, _=np.linalg.eig(VW_stack)
        #wcc=np.array([(1j*np.log(z,dtype=complex)).real/(2*np.pi) %1 for z in eig])
        wcc_normalized=np.array([iwcc/(2*np.pi) %1 for iwcc in wcc])
        return wcc_normalized
        
    
    def plotWCC(self, algorithm='WilsonLoop'):
        '''
        There I provide three algorithms to calculate the one-dimensional wccs
        algorithms=WilsonLoop | ParallelTransport | GradualDecompose
        '''
        fig = plt.figure()
        for ikline, kline in enumerate(self.wcc_path):
            wcc_sum=0.
            if (algorithm=='ParallelTransport'):
                for wcc in self.calcParallelTransport(kline):
                    wcc_sum+=wcc
                    plt.scatter(ikline,wcc,color='black')
            elif algorithm=='WilsonLoop':
                for wcc in self.calcWilsonLoop(kline):
                    wcc_sum+=wcc
                    plt.scatter(ikline,wcc,color='black')
            elif algorithm=='GradualDecompose':
                for wcc in self.calcGradualDecompose(kline):
                    wcc_sum+=wcc
                    plt.scatter(ikline,wcc,color='black')
            else:
                print 'there are three optional algorithms: WilsonLoop | ParallelTransport | GradualDecompose'
                exit()
                
            plt.scatter(ikline,wcc_sum%1,color='red')
            
        plt.show()
        
    def calcBerryCurvAtSingleKpt(self,k0,b1,b2,bandset=-1, algorithm='WilsonLoop'):
        kline=[k0,k0+b1,k0+b1+b2]
        wcc_sum=0.
        if (algorithm=='ParallelTransport'):
            for wcc in self.calcParallelTransport(kline,bandset=bandset):
                wcc_sum+=wcc
        elif algorithm=='WilsonLoop':
            for wcc in self.calcWilsonLoop(kline,bandset=bandset):
                wcc_sum+=wcc
        elif algorithm=='GradualDecompose':
            for wcc in self.calcGradualDecompose(kline,bandset=bandset):
                wcc_sum+=wcc
                
        return wcc_sum%1

class WannierKit(object):
    '''
    classdocs
    '''
    kpts_list=[]
    lattice=[]
    lattice_inv=[]
    reciprocal_lattice=[]
    overlap_matrix=[]  # [kpt][nnkpt][ibnd][jbnd]
    overlap_corr_nnkpts=[] # [startkp][inkpt]
    overlap_nkpts_index=[] # [startkp][inkpt]
    overlap_b=[] # [startkp][inkpt]
    
    b_list=[]
    wb_list=[]
    
    nkpts=-1
    nbands=-1
    nnkpts=-1
    
    def __init__(self):
        '''
        Constructor
        '''
    def calcWannierCenter(self,iband_range=[0,1],startkpt=0,endkpt=-1):
        '''
        <a(r)|r|a(r)>=Sum_k <u_k|id/dk|u_k>
                     =Sum_k_b wb*b*Im Ln<u_k|u_k+b>
        '''
        endkpt=self.nkpts if endkpt<0 else endkpt
        nband=iband_range[1]-iband_range[0]+1
        lambda_=np.eye(nband)
        
        kpath=self.kpts_list[endkpt]-self.kpts_list[startkpt]
        #print '',self.getM(ikpt=0,inkpt=0,iband_range=iband_range)
        Wilson_loop=np.eye(nband)
        VW_stack=np.eye(nband)
        for ikpt in range(startkpt,endkpt+1):
            for inkpt in range(self.nnkpts):
                # only include one direction
                if(np.dot(self.overlap_b[ikpt][inkpt],kpath)>0):
                    M=self.getM(ikpt=ikpt,inkpt=inkpt,iband_range=iband_range)
                    [V,E,W]=np.linalg.svd(M)
                    VW_stack=np.dot(VW_stack,np.dot(V,W))
                    
                    Wilson_loop=np.dot(Wilson_loop,M)
                    
 
        eig2,_=np.linalg.eig(Wilson_loop)
        print eig2
        print [(np.log(z)/(2*np.pi)).imag % 1 for z in eig2]
        
        [eig,_]=np.linalg.eig(VW_stack.conjugate().transpose())
        print eig
        wcc=[(np.log(z)/(2*np.pi)).imag % 1 for z in eig]
        return wcc
    
    def getM(self,ikpt=0,inkpt=0,iband_range=[0,1]):
        M=[]
        for iband in range(iband_range[0],iband_range[1]+1):
            M_tmp=[]
            for jband in range(iband_range[0],iband_range[1]+1):
                M_tmp.append(self.overlap_matrix[ikpt][inkpt][iband][jband])
            M.append(M_tmp)
        return np.array(M)
    
    def cartToFrac(self,cart_coor):
        #print 'cart:',cart_coor
        #print self.lattice_inv
        frac_coor=np.dot(cart_coor,self.lattice_inv)
        #print 'frac:',frac_coor
        frac_coor_to_home=[]
        for coor in frac_coor:
            #print 'former',coor
            frac_coor_to_home.append(coor%1)
            #print 'latter',coor
        return np.array(frac_coor_to_home)
    
    def readData(self):
        self.readKpts()
        self.readLattice()
        self.readOverlap()
        
    def readKpts(self):
        f=open('wannier90.win','r')
        data=f.read()
        f.close()
        data=[line for line in data.split('\n') if line]

        kpoints_block_met=False
        for line in data:
            if(line=='end kpoints'):break
            if(kpoints_block_met):
                kx, ky, kz =re.findall('[\-0-9.eE]+',line)
                self.kpts_list.append(np.array([float(kx),float(ky),float(kz)]))
            if(line=='begin kpoints'):kpoints_block_met=True
        
    def readLattice(self):
        f=open('wannier90.wout','r')
        data=f.read()
        f.close()
        data=[line for line in data.split('\n') if line]
        lattice_position=-1
        klattice_position=-1
        nnkpt_position=-1
        b_list_position=-1
        for line_index, line in enumerate(data):
            if(line=='                              Lattice Vectors (Ang)'):lattice_position=line_index+1
            if(line=='                        Reciprocal-Space Vectors (Ang^-1)'):klattice_position=line_index+1
            if(line==' |                        Shell   # Nearest-Neighbours                        |'):nnkpt_position=line_index+2
            if(line==' |                  b_k Vectors (Ang^-1) and Weights (Ang^2)                  |'):
                b_list_position=line_index+4
                break
            
        lattice=[]
        for line_index in range(lattice_position,lattice_position+3):
            #print data[line_index]
            _, ax,ay,az=re.findall('[\-0-9.eE]+',data[line_index])
            lattice.append([float(ax),float(ay),float(az)])
        self.lattice=np.array(lattice)
        self.lattice_inv=np.linalg.inv(self.lattice)
        
        reciprocal_lattice=[]
        for line_index in range(klattice_position,klattice_position+3):
            #print data[line_index]
            _, kx,ky,kz=re.findall('[\-0-9.eE]+',data[line_index])
            reciprocal_lattice.append([float(kx),float(ky),float(kz)])
        self.reciprocal_lattice=np.array(reciprocal_lattice)

        snn=0
        line_index=nnkpt_position
        while(True):
            try:
                _, nn=re.findall('[\-0-9.eE]+',data[line_index])
                line_index+=1
                snn+=int(nn)
            except ValueError as e:
                break   

        for line_index in range(b_list_position,b_list_position+snn):
            _, bkx,bky,bkz, wb=re.findall('[\-0-9.eE]+',data[line_index])
            self.b_list.append(np.array([float(bkx),float(bky),float(bkz)]))
            self.wb_list.append(float(wb))
            
    def readOverlap(self):
        f=open('wannier90.mmn','r')
        data=f.read()
        f.close()
        data=[line for line in data.split('\n') if line]
        nbands, nkpts, nnkpts=re.findall('[\-0-9.eE]+',data[1])
        nbands=int(nbands)
        nkpts=int(nkpts)
        nnkpts=int(nnkpts)
        self.nbands=nbands
        self.nkpts=nkpts
        self.nnkpts=nnkpts
        #Mmn=[] #with four index: k, nnk, ibnd, jbnd
        #Overlap_kpts=[]
        line_index=2
        for ikpt in range(nkpts):
            overlap_nnkpts_tmp=[]
            overlap_nnkpts_index_tmp=[]
            overlap_b_tmp=[]
            overlap_per_k=[]
            for innkpt in range(nnkpts):
                start_kp, end_kp, enhance_kpx, enhance_kpy, enhance_kpz = re.findall('[\-0-9.eE]+',data[line_index])
                end_kp_coor=self.kpts_list[int(end_kp)-1]+np.array([int(enhance_kpx),int(enhance_kpy),int(enhance_kpz)])
                overlap_nnkpts_index_tmp.append(int(end_kp)-1)
                overlap_nnkpts_tmp.append(end_kp_coor)
                overlap_b_tmp.append(end_kp_coor-self.kpts_list[int(start_kp)-1])
                line_index+=1
        
                overlap_per_k_per_nnk=[]
                for iband in range(nbands):
                    overlap_per_k_per_nnk_per_ibnd=[]
                    for jband in range(nbands):
                        relpart, imgpart = re.findall('[\-0-9.eE]+',data[line_index])
                        line_index+=1
                        overlap_per_k_per_nnk_per_ibnd.append(float(relpart)+float(imgpart)*1j)
                    overlap_per_k_per_nnk.append(overlap_per_k_per_nnk_per_ibnd)
                overlap_per_k.append(overlap_per_k_per_nnk)
            self.overlap_corr_nnkpts.append(overlap_nnkpts_tmp)
            self.overlap_b.append(overlap_b_tmp)
            self.overlap_nkpts_index.append(overlap_nnkpts_index_tmp)
            self.overlap_matrix.append(overlap_per_k)
    