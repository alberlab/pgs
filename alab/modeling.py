#!/usr/bin/env python

# Copyright (C) 2015 University of Southern California and
#                          Nan Hua
# 
# Authors: Nan Hua
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Prerequests:
# IMP 2.5 is required for this module

__author__  = "Nan Hua"
__credits__ = ["Nan Hua","Ke Gong","Harianto Tjong"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"

import matrix as alabmatrix #alabmatrix
import utils  as alabutils #alabutils
import time
import numpy as np
import IMP
import IMP.core
import IMP.container
import IMP.algebra
import IMP.atom
import random
import h5py
import logging
import io #system io
import re

class tadmodel(object):
    """
    A wrapper to do IMP modeling in TAD level beads
    
    Parameters
    ----------
    
    probfile : alab.matrix.contactmatrix instant
        probability matrix at TAD level, alab.matrix.contactmatrix hdf5 format is required
    nucleusRadius : float
        radius of nucleus, default 5000(nm)
    contactRange : int
        folds for surface-surface contact coefficient
    level : loglevel, 
        default:None will record everything during caculation
        debug, info, warning, error, critical are supported
    """
    def __init__(self,probfile,nucleusRadius=5000.0,contactRange=1,level=None):
        self.probmat = alabmatrix.contactmatrix(probfile)
        self.nbead   = len(self.probmat)
        #setup log
        LEVELS={'debug':logging.DEBUG,'info':logging.INFO,'warning':logging.WARNING,'error':logging.ERROR,'critical':logging.CRITICAL}
        loglevel = LEVELS.get(level,logging.NOTSET)
        self.logger = logging.getLogger()
        self.logger.setLevel(loglevel)
        self._log_capture_string = io.StringIO()
        chhandler = logging.StreamHandler(self._log_capture_string)
        chhandler.setLevel(loglevel)
        self.logger.addHandler(chhandler)
        self.logger.setLevel(loglevel)
        #CONST
        rscale               = 1.38                  # 20% occupancy
        self.nucleusRadius   = nucleusRadius         # nm
        cdensity             = 107.45                # bp/nm assuming 197 bp/nucleosomes and 6 nucleosome/11 nm
        kscale               = (0.75*15**2)**(1.0/3.0) # 3/4*r**2 where r=15nm
        self.contactRange    = contactRange          # surface to surface distance scale of (r1+r2)
                                                     # for which 2 beads are considered as contact 
        #get radius of each bead
        self.beadRadius = [rscale * kscale * ((index['end'] - index['start'])/cdensity) ** (1.0/3.0) for index in self.probmat.idx]
        #calculate the total volumn of DNA (diploid) and nucleus
        dnavol   = sum(4. * 3.1415/3. * np.array(self.beadRadius)**3) * 2 
        nucvol   = (4*3.1415/3)*self.nucleusRadius**3
        #And chromosome occupancy
        dnaocc   = dnavol / nucvol
        self.logger.debug(u'occupancy: %.2f with Rnuc %d'%(dnaocc,self.nucleusRadius))
        #diploid Rb; 2xtotal haploid beads 
        self.beadRadius = self.beadRadius + self.beadRadius
        # Chromosome territory apply
        self.genome = alabutils.genome(self.probmat.genome)
        cscale=1.0
        chrvol = nucvol * self.genome.info['length']/sum(self.genome.info['length'])/2
        self.chromRadius=cscale*((chrvol/4*3/3.1415)**(1./3.))
        
        #record starting time
        self.model      = IMP.Model()
        self.chain      = IMP.container.ListSingletonContainer(self.model)
        self.restraints = IMP.RestraintSet(self.model)
        #IMP.set_check_level(IMP.USAGE)
        IMP.set_check_level(IMP.NONE)
        IMP.set_log_level(IMP.SILENT)
        #setup nucleus envelope
        self.center = IMP.algebra.Vector3D(0,0,0)
    #=========================
    def updateScoringFunction(self,restraintset=None):
        if restraintset is None:
            self.scoringFunction = IMP.core.RestraintsScoringFunction(self.restraints)
        else:
            self.scoringFunction = IMP.core.RestraintsScoringFunction(restraintset)
        return self.scoringFunction
    
    def cache_coordinates(self):
        i = -1
        for p in self.chain.get_particles():
            i += 1
            pattr = IMP.core.XYZR(p)
            self.xyz[i,0] = pattr.get_x()
            self.xyz[i,1] = pattr.get_y()
            self.xyz[i,2] = pattr.get_z()
            self.r[i]     = pattr.get_radius()
        #---
    def set_coordinates(self,coordinates=None):
        if coordinates is None:
            randomCoordinates = True
            random.seed()
            randomCap = IMP.algebra.Sphere3D(self.center, 1.5*self.nucleusRadius)
        else:
            randomCoordinates = False
        self.xyz = np.zeros((2*self.nbead,3))
        self.r   = np.zeros((2*self.nbead,1))
        for i in range(2*self.nbead):
            if (randomCoordinates):
                coor = IMP.algebra.get_random_vector_in(randomCap)
            else:
                coor = IMP.algebra.Vector3D(coordinates[i,0],coordinates[i,1],coordinates[i,2])
               
            sph = IMP.algebra.Sphere3D(coor,self.beadRadius[i])
            p0  = IMP.Particle(self.model)
            self.chain.add(p0)
            ma  = IMP.atom.Mass.setup_particle(p0,1)   #set mass = 1 g/mol
            sp  = IMP.core.XYZR.setup_particle(p0,sph)
            sp.set_coordinates_are_optimized(True)
        #-
        self.cache_coordinates()
    #--
    def set_excludedVolume(self,ksping=1,slack=10):
        # Set up excluded volume
        self.excludedVolumeRestraint = IMP.core.ExcludedVolumeRestraint(self.chain,ksping,slack)
        self.restraints.add_restraint(self.excludedVolumeRestraint) #1
        return self.excludedVolumeRestraint 
    
    def set_nucleusEnvelope(self,kspring):
        # Set up Nucleus cap  #
        ubnuc = IMP.core.HarmonicUpperBound(self.nucleusRadius,kspring)
        ssnuc = IMP.core.DistanceToSingletonScore(ubnuc,self.center) #center-to-center distance
        self.nucleusEnvelopeRestraint = IMP.container.SingletonsRestraint(ssnuc,self.chain)
        self.restraints.add_restraint(self.nucleusEnvelopeRestraint) #2
        return self.nucleusEnvelopeRestraint
    
    def set_consecutiveBeads(self,lowprob=0.1):
        # Setup consecutive bead upper bound constraints
        self.consecutiveBeadRestraints = self._get_consecutiveBeadRestraints(lowprob=lowprob,kspring=10)
        for rs in self.consecutiveBeadRestraints: #3
            self.restraints.add_restraint(rs) 
        self.logger.debug(u"Total consecutive bead restraints %d"%(len(self.consecutiveBeadRestraints)))
        print "Total consecutive bead restraints %d"%(len(self.consecutiveBeadRestraints))
        return self.consecutiveBeadRestraints
    
    def set_fmaxRestraints(self,kspring=1):
        self.fmaxRestraints = self._get_fmaxRestraints(kspring=kspring)
        for rs in self.fmaxRestraints: #4
            self.restraints.add_restraint(rs)
        self.logger.debug(u"Total fmax restraints %d"% (len(self.fmaxRestraints)))
        print "Total fmax restraints %d"% (len(self.fmaxRestraints))
        return self.fmaxRestraints
    
    def set_contactRestraints(self,actdist,kspring=1):
        self.intraContactRestraints, self.interContactRestraints = self._get_contactRestraints(actdist,kspring)
        for rs in self.intraContactRestraints:#4 add intra chromosome contacts
            self.restraints.add_restraint(rs)
        for rs in self.interContactRestraints:#4 add inter chromosome contacts
            self.restraints.add_restraint(rs)
        return self.intraContactRestraints, self.interContactRestraints
    
    #=========================chromosome territory functions
    def CondenseChromosome(self,rrange=0.5):
        """
        Collapse chains around centromere beads
        
        Parameters
        ----------
        rrange : float
            scale parameter in [0,1] for the radius limit
        """
        i = -1
        for chrom in self.genome.info['chrom']:
            i += 1
            #find centromere
            cenbead = np.flatnonzero((self.probmat.idx['chrom'] == chrom) & (self.probmat.idx['flag'] == 'CEN'))[0]
            p0A=IMP.core.XYZ(self.chain.get_particles()[cenbead]) #fetch indexes
            p0B=IMP.core.XYZ(self.chain.get_particles()[cenbead+self.nbead])
            coorA = p0A.get_coordinates()
            coorB = p0B.get_coordinates()
            rlimit = self.chromRadius[i]*rrange
            for j in np.flatnonzero(self.probmat.idx['chrom'] == chrom):
                p1A=IMP.core.XYZ(self.chain.get_particles()[j])
                p1B=IMP.core.XYZ(self.chain.get_particles()[j+self.nbead])
                dx=(2*random.random()-1)*rlimit
                dy=(2*random.random()-1)*rlimit
                dz=(2*random.random()-1)*rlimit
                randA = coorA
                randA[0] += dx
                randA[1] += dy
                randA[2] += dz
                dx=(2*random.random()-1)*rlimit
                dy=(2*random.random()-1)*rlimit
                dz=(2*random.random()-1)*rlimit
                randB = coorB
                randB[0] += dx
                randB[1] += dy
                randB[2] += dz
                p1A.set_coordinates(randA) #placed nearby cen
                p1B.set_coordinates(randB) #placed nearby cen
            #--
        #--
    #=============================end chromosome terriory
    #=============================restraint set functions
    def _get_beadDistanceRestraint(self,bead1,bead2,dist,kspring=1):
        """
        get distance upper bound restraint to bead1 and bead2
        
        
        Parameters
        -----------
        
        bead1,bead2 : int 
            bead id
        dist : int
            distance upper bound
        kspring : int
            harmonic constant k
        
        Return
        ------
        
        restraint 
        """
        restraintName = "Bead (%d,%d) : %f k = %.1f" % (bead1,bead2,dist,kspring)
        ds = IMP.core.SphereDistancePairScore(IMP.core.HarmonicUpperBound(dist,kspring))
        pr = IMP.core.PairRestraint(self.model,ds,IMP.ParticlePair(self.chain.get_indexes()[bead1],self.chain.get_indexes()[bead2]),restraintName)
        return pr
    
    def _get_consecutiveBeadRestraints(self,lowprob=0.1,kspring=10):
        """
        calculate distance constraints to consecutive beads
        
        Parameters
        -----------
        lowprob : Min probility for consecutive beads
        """
        consecRestraints = []
        for i in range(self.nbead-1):
            if self.probmat.idx[i]['chrom'] != self.probmat.idx[i+1]['chrom']:
                continue
            p = max(self.probmat.matrix[i,i+1],lowprob)
            b1 = i;b2 = i+1
            b3 = b1 + self.nbead
            b4 = b2 + self.nbead
            #calculate upper bound for consecutive domains
            consecDist = consecutiveDistanceByProbability(self.beadRadius[b1],self.beadRadius[b2],p,self.contactRange+1)
            
            # set k = 10 to be strong interaction
            rs1 = self._get_beadDistanceRestraint(b1,b2,consecDist,kspring) 
            rs2 = self._get_beadDistanceRestraint(b3,b4,consecDist,kspring) 
            
            #push restraint into list
            consecRestraints.append(rs1)
            consecRestraints.append(rs2)
            
            if i>0 and self.probmat.idx[i]['chrom'] == self.probmat.idx[i-1]['chrom'] and self.probmat.idx[i]['flag']!="domain" and self.probmat.idx[i-1]!="gap":
                p = max(self.probmat.matrix[i-1,i+1],lowprob)
                b1 = i-1;b2 = i+1
                b3 = b1 + self.nbead
                b4 = b2 + self.nbead
                consecDist = consecutiveDistanceByProbability(self.beadRadius[b1],self.beadRadius[b2],p,self.contactRange+1)
                rs1 = self._get_beadDistanceRestraint(b1,b2,consecDist,kspring) 
                rs2 = self._get_beadDistanceRestraint(b3,b4,consecDist,kspring) 
                consecRestraints.append(rs1)
                consecRestraints.append(rs2)
            #---
        #-------
        return consecRestraints
    #--------------probmat restraints
    def _get_minPairRestraints(self,bpair,dist,minnum,kspring = 1):
        """
        Return restraint decorater of min pair restraints
        for minnum out of bpairs are satisfied 
        
        Parameters
        -----------
        
        bpair : tuple list
            contact pair candidates
        dist : int
            distance upperbound for contact
        minnum : int
            minimun number of pairs required to satisify
        contactRange : int 
            scale of (r1+r2) where a contact is defined   
        """
        ambi = IMP.container.ListPairContainer(self.model)
        restraintName = "Bead [ "
        for p in bpair:
            restraintName += "(%d,%d) " % (p[0],p[1])
            p0 = self.chain.get_particles()[p[0]]
            p1 = self.chain.get_particles()[p[1]]
            pair = IMP.ParticlePair(p0,p1)
            ambi.add(pair)
        restraintName += "] : %f k = %.1f" % (dist,kspring)
        ds = IMP.core.SphereDistancePairScore(IMP.core.HarmonicUpperBound(dist,kspring))
        minpr = IMP.container.MinimumPairRestraint(ds,ambi,minnum,restraintName)
        return minpr
    
    def _get_fmaxRestraints(self,kspring=1):
        """
        return restraints list for prob=1.0
        """
        fmaxrs = []
        for i in range(self.nbead):
            for j in range(i+1,self.nbead):
                if self.probmat.matrix[i,j] <= 0.999:
                    continue
                if self.probmat.idx[i]['chrom'] == self.probmat.idx[j]['chrom']: #intra
                    if j-i > 1:
                        rs1 = self._get_beadDistanceRestraint(i,j,self.contactRange*(self.beadRadius[i]+self.beadRadius[j]),kspring)
                        rs2 = self._get_beadDistanceRestraint(i+self.nbead,j+self.nbead,self.contactRange*(self.beadRadius[i]+self.beadRadius[j]),kspring)
                        fmaxrs.append(rs1)
                        fmaxrs.append(rs2)
                else: #inter
                    bpair = [(i,j),(i,j+self.nbead),(i+self.nbead,j),(i+self.nbead,j+self.nbead)] #bead pair
                    minprrs = self._get_minPairRestraints(bpair,self.contactRange*(self.beadRadius[i]+self.beadRadius[j]),minnum=2,kspring=kspring)
                    fmaxrs.append(minprrs)
                #--
            #--
        #--
        return fmaxrs
    def _get_contactRestraints(self,actdist,kspring=1):
        """
        return restraints list given actdist list
        
        Parameters
        ----------
        actdist : array like
            Activation distanct array [i,j,dist]
        """
        intracontactrs = []
        intercontactrs = []
        for (b1,b2,dcutoff) in actdist:
            lastdist1 = surfaceDistance(self.chain,(b1,b2))
            lastdist2 = surfaceDistance(self.chain,(b1+self.nbead,b2+self.nbead))
            if self.probmat.idx[b1]['chrom'] == self.probmat.idx[b2]['chrom']: #intra chromosome
                if max(lastdist1,lastdist2) <= dcutoff:
                    rs1 = self._get_beadDistanceRestraint(b1,b2,self.contactRange*(self.beadRadius[b1]+self.beadRadius[b2]),kspring)
                    rs2 = self._get_beadDistanceRestraint(b1+self.nbead,b2+self.nbead,self.contactRange*(self.beadRadius[b1]+self.beadRadius[b2]),kspring)
                    intracontactrs.append(rs1)
                    intracontactrs.append(rs2)
                elif min(lastdist1,lastdist2) <= dcutoff:
                    bpair = [(b1,b2),(b1+self.nbead,b2+self.nbead)] #bead pair
                    minprrs = self._get_minPairRestraints(bpair,self.contactRange*(self.beadRadius[b1]+self.beadRadius[b2]),minnum=1)
                    intracontactrs.append(minprrs)
                #else none contact enfored
            else:
                lastdist3 = surfaceDistance(self.chain,(b1,b2+self.nbead))
                lastdist4 = surfaceDistance(self.chain,(b1+self.nbead,b2))
                sdists    = sorted([lastdist1,lastdist2,lastdist3,lastdist4])
                if sdists[1] <= dcutoff:
                    nchoose = 2
                elif sdists[0] <= dcutoff:
                    nchoose = 1
                else: #non enforced
                    continue
                
                bpair = [(b1,b2),(b1,b2+self.nbead),(b1+self.nbead,b2),(b1+self.nbead,b2+self.nbead)] #bead pair
                minprrs = self._get_minPairRestraints(bpair,self.contactRange*(self.beadRadius[b1]+self.beadRadius[b2]),minnum=nchoose)
                intercontactrs.append(minprrs)
        return intracontactrs, intercontactrs
    #==========================end probmat restraint
    #-------------------------modeling
    def cgstep(self,step,silent=False):
        """
            perform conjugate gradient on model using scoring function sf
        """
        t0 = time.time()
        self.cache_coordinates()
        o = IMP.core.ConjugateGradients(self.model)
        o.set_check_level(IMP.USAGE)
        o.set_scoring_function(self.scoringFunction)
        try:
            s = o.optimize(step)
        except:
            #tt = str(np.random.random_integers(1000))+'.txt'
            #np.savetxt(tt,self.xyz,fmt='%.4f')
            #print tt
            #self.cache_coordinates()
            #np.savetxt(tt+'.af',self.xyz,fmt='%.4f')
            s = o.optimize(step)
        if not silent:
            self.logger.info(u'CG %d steps done @ %.1fs score = %f'%(step,alabutils.timespend(t0),s))
            print 'CG %d steps done @ %.1fs score = %f'%(step,alabutils.timespend(t0),s)
        return s

    def mdstep(self,t,step,gamma=0.1,silent=False):
        t0 = time.time()
        #self.cache_coordinates()
        xyzr = self.chain.get_particles()
        o    = IMP.atom.MolecularDynamics(self.model)
        o.set_scoring_function(self.scoringFunction)
        md   = IMP.atom.VelocityScalingOptimizerState(self.model,xyzr,t)
        #md   = IMP.atom.LangevinThermostatOptimizerState(self.model,xyzr,t,gamma)
        o.add_optimizer_state(md)
        s    = o.optimize(step)
        o.remove_optimizer_state(md)
        if not silent:
            self.logger.info(u'MD %d steps done @ %.1fs score = %f'%(step,alabutils.timespend(t0),s))
            print 'MD %d steps done @ %.1fs score = %f'%(step,alabutils.timespend(t0),s)
        return s
    def mdstep_withChromosomeTerritory(self,t,step):
        """
        perform an mdstep with chromosome terriory restraint
        
        Parameters
        -----------
        t : int
            temperature
        step : int
            optimization steps
        """
        t0 = time.time()
        chrContainers=[]
        for chrom in self.genome.info['chrom']:
            chromContainer1=IMP.container.ListSingletonContainer(self.model,'Container %s s1'%chrom)
            chromContainer2=IMP.container.ListSingletonContainer(self.model,'Container %s s2'%chrom)
            for j in np.flatnonzero(self.probmat.idx['chrom'] == chrom):
                p=self.chain.get_particle(j)
                chromContainer1.add(p)
                p=self.chain.get_particle(j+self.nbead)
                chromContainer2.add(p)
            chrContainers.append(chromContainer1)
            chrContainers.append(chromContainer2)
        # set each chromosome to different containers
        for st in range(step/10):
            ctRestraintSet = IMP.RestraintSet(self.model)
            for i in range(len(self.genome.info['chrom'])):
                comxyz = centerOfMass(chrContainers[2*i])
                comc   = IMP.algebra.Vector3D(comxyz[0],comxyz[1],comxyz[2])
                ub     = IMP.core.HarmonicUpperBound(self.chromRadius[i],0.2)
                ss     = IMP.core.DistanceToSingletonScore(ub,comc)
                ct_rs  = IMP.container.SingletonsRestraint(ss,chrContainers[2*i])
                
                ctRestraintSet.add_restraint(ct_rs)
                #another one
                comxyz = centerOfMass(chrContainers[2*i+1])
                comc   = IMP.algebra.Vector3D(comxyz[0],comxyz[1],comxyz[2])
                ub     = IMP.core.HarmonicUpperBound(self.chromRadius[i],0.2)
                ss     = IMP.core.DistanceToSingletonScore(ub,comc)
                ct_rs  = IMP.container.SingletonsRestraint(ss,chrContainers[2*i+1])
                
                ctRestraintSet.add_restraint(ct_rs)
            self.updateScoringFunction([self.restraints,ctRestraintSet])
            s = self.mdstep(t,10,silent=True)
        #---
        #s = mdstep(model,chain,sf,t,1000)
        self.updateScoringFunction()
        self.logger.info(u'CT-MD %d steps done @ %.1fs score = %f'%(step,alabutils.timespend(t0),s))
        print 'CT-MD %d steps done @ %.1fs score = %f'%(step,alabutils.timespend(t0),s)
        return s
    def SimulatedAnnealing(self,hot,cold,nc=10,nstep=500):
        """
            perform a cycle of simulated annealing from hot to cold
        """
        t0 = time.time()
        dt = (hot-cold)/nc
        for i in range(nc):
            t = hot-dt*i
            s = self.mdstep(t,nstep,silent=True)
            self.logger.info(u"      Temp=%d Step=%d Time=%.1fs Score = %.8f"%(t,nstep,alabutils.timespend(t0),s))
            print "      Temp=%d Step=%d Time=%.1fs Score = %.8f"%(t,nstep,alabutils.timespend(t0),s)
        s= self.mdstep(cold,nstep,silent=True)
        self.logger.info(u"      Temp=%d Step=%d Time=%.1fs Score = %.8f"%(cold,nstep,alabutils.timespend(t0),s))
        print "      Temp=%d Step=%d Time=%.1fs Score = %.8f"%(cold,nstep,alabutils.timespend(t0),s)
        self.cgstep(100,silent=True)

    def SimulatedAnnealing_Scored(self,hot,cold,nc=10,nstep=500,lowscore=10):
        """perform a cycle of simulated annealing but stop if reach low score"""
        t0 = time.time()
        dt = (hot-cold)/nc
        for i in range(nc):
            t = hot-dt*i
            self.mdstep(t,nstep,silent=True)
            score = self.cgstep(100,silent=True)
            self.logger.info(u"      Temp=%s Step=%s Time=%.1fs Score=%.8f"%(t,nstep,alabutils.timespend(t0),score))
            print "      Temp=%s Step=%s Time=%.1fs Score=%.8f"%(t,nstep,alabutils.timespend(t0),score)
            if score < lowscore:
                return t,score
                break
        self.mdstep(cold,300,silent=True)
        score = self.cgstep(100,silent=True)
        self.logger.info(u"      Temp=%s Step=%s Time=%.1fs Score=%.8f"%(cold,nstep,alabutils.timespend(t0),score))
        print "      Temp=%s Step=%s Time=%.1fs Score=%.8f"%(cold,nstep,alabutils.timespend(t0),score)
        return t,score
    
    def shrinkingOptimization(self,drange,shrinkScore,minscore,interScale):
        radStartShrink = int((1+drange)*self.nucleusRadius)
        radEndShrink = int((1-drange)*self.nucleusRadius)
        nucExpand=drange+1.3
        radExpand = int(nucExpand*self.nucleusRadius)
        incr = int(interScale*self.nucleusRadius)
        nucrads = range(radEndShrink,radStartShrink,incr)
        nucrads.append(radStartShrink)
        nucrads = sorted(nucrads, reverse=True)
        self.logger.debug(u"Optimization with decreasing NE %s"%(nucrads))
        print "Optimization with decreasing NE %s"%(nucrads)
        expanded=False
        self.logger.debug(u"\t--- Start shrinking ---")
        print "\t--- Start shrinking ---"
        self.restraints.remove_restraint(self.nucleusEnvelopeRestraint) #will be replaced by temporary
        for r_nuc in nucrads:
            ubnuc = IMP.core.HarmonicUpperBound(r_nuc,1.0)
            ssnuc = IMP.core.DistanceToSingletonScore(ubnuc,self.center)
            rnuc = IMP.container.SingletonsRestraint(ssnuc,self.chain)
            
            self.updateScoringFunction([self.restraints,rnuc])
            temp, score = self.SimulatedAnnealing_Scored(2000,300,nc=2, lowscore=shrinkScore)
            score = self.cgstep(500,silent=True)
            self.logger.debug(u"Score optimizing temporary NE %dmm done: %.8f"%(r_nuc,score))
            print "Score optimizing temporary NE %dmm done: %.8f"%(r_nuc,score)
            if score > minscore:
                if expanded == False: #haven't done expansion steps before
                    self.logger.debug(u"\t    --- Start expanding ---")
                    print "\t    --- Start expanding ---"
                    ubexpend = IMP.core.HarmonicUpperBound(radExpand,1.0)
                    ssexpend = IMP.core.DistanceToSingletonScore(ubexpend,self.center)
                    rexpend = IMP.container.SingletonsRestraint(ssexpend,self.chain)
                    expanded = True
                    self.updateScoringFunction([self.restraints,rexpend])
                    temp, score = self.SimulatedAnnealing_Scored(75000,500,nstep=1000, nc=2, lowscore=shrinkScore)
                    self.logger.debug(u"\t     ...back to shrinking...")
                    print "\t     ...back to shrinking..."
                else:
                    self.logger.debug(u'\t     Score is still high after expansion: %.8f' %(score))
                    print '\t     Score is still high after expansion: %.8f' %(score)
                    break
            #-
        #--
        # Putting back original radius cap  #
        self.restraints.add_restraint(self.nucleusEnvelopeRestraint) #2
        self.updateScoringFunction()
        temp, score = self.SimulatedAnnealing_Scored(2000,300, nc=2, nstep=1000, lowscore=1)
        self.logger.debug(u"Recover nucleus %.1f nm at T=%.1f K, score: %.8f"%(self.nucleusRadius,temp,score))
        print "Recover nucleus %.1f nm at T=%.1f K, score: %.8f"%(self.nucleusRadius,temp,score)
        return 0
    #==================utils
    def evaluateRestraints(self,restraintset,tolerance=0.05):
        total = 0
        for rs in restraintset:
            score = rs.get_score()
            rsstr = rs.get_name()
            dist,k = re.findall(': (\d+.\d+) k = (\d+.\d+)',rsstr)[0]
            dist = float(dist)
            k    = float(k)
            if (2*score/k)**0.5 > tolerance*dist:
                total += 1
                self.logger.warning(u"%s %f" % (rsstr,score))
                print "%s %f" % (rsstr,score)
        return total
    def savepym(self,filename):
        pymfile = IMP.display.PymolWriter(filename)
        g = IMP.core.XYZRsGeometry(self.chain)
        g.set_name('beads')
        g.set_color(IMP.display.Color(1,1,1))
        pymfile.add_geometry(g)

    def savepym_withChromosome(self,filename,s=1,v=1):
        import colorsys
        pymfile = IMP.display.PymolWriter(filename)
        nchrom  = len(self.genome.info['chrom'])
        h       = np.array(range(nchrom+1))*2.0/(nchrom+1)/3
        i       = -1
        for chrom in self.genome.info['chrom']:
            i += 1
            chromContainer1=IMP.container.ListSingletonContainer(self.model,'Container %s s1'%chrom)
            chromContainer2=IMP.container.ListSingletonContainer(self.model,'Container %s s2'%chrom)
            for j in np.flatnonzero(self.probmat.idx['chrom'] == chrom):
                p=self.chain.get_particle(j)
                chromContainer1.add(p)
                p=self.chain.get_particle(j+self.nbead)
                chromContainer2.add(p)
            chrColor = colorsys.hsv_to_rgb(h[i],s,v)
            color = IMP.display.Color(chrColor[0],chrColor[1],chrColor[2])
            g1 = IMP.core.XYZRsGeometry(chromContainer1)
            g1.set_name(chrom+' s1')
            g1.set_color(color)
            pymfile.add_geometry(g1)
            g2 = IMP.core.XYZRsGeometry(chromContainer2)
            g2.set_name(chrom+' s2')
            g2.set_color(color)
            pymfile.add_geometry(g2)

    def saveCoordinates(self,filename,prefix):
        import cPickle
        import cStringIO
        if (filename[-4:] != '.hms'):
            filename += '.hms'
        log_contents = self._log_capture_string.getvalue()
        #self._log_capture_string.close()
        #print log_contents
        pymhandler = cStringIO.StringIO()
        self.savepym_withChromosome(pymhandler)
        #xyz = np.zeros((len(self.chain.get_particles()),3))
        #r   = np.zeros((len(self.chain.get_particles()),1))
        #i = -1
        #for p in self.chain.get_particles():
            #i += 1
            #pattr = IMP.core.XYZR(p)
            #xyz[i,0] = pattr.get_x()
            #xyz[i,1] = pattr.get_y()
            #xyz[i,2] = pattr.get_z()
            #r[i] = pattr.get_radius()
        #---
        self.cache_coordinates()
        h5f = h5py.File(filename,'a')
        if not 'genome' in h5f.keys():
            h5f.create_dataset('genome',data = cPickle.dumps(self.probmat.genome))
        if not 'idx' in h5f.keys():
            h5f.create_dataset('idx',data = self.probmat.idx,compression='gzip')
        if prefix in h5f.keys():
            grp = h5f[prefix]
            if 'xyz' in grp.keys():
                grp['xyz'][...] = self.xyz
            else:
                grp.create_dataset('xyz',data=self.xyz,compression='gzip')
            if 'r' in grp.keys():
                grp['r'][...] = self.r
            else:
                grp.create_dataset('r',data=self.r,compression='gzip')
            if 'log' in grp.keys():
                grp['log'][...] = cPickle.dumps(log_contents)
            else:
                grp.create_dataset('log',data=cPickle.dumps(log_contents))
            if 'pym' in grp.keys():
                grp['pym'][...] = cPickle.dumps(pymhandler.getvalue())
            else:
                grp.create_dataset('pym',data=cPickle.dumps(pymhandler.getvalue()))
        else:
            grp = h5f.create_group(prefix)
            grp.create_dataset('xyz',data=self.xyz,compression='gzip')
            grp.create_dataset('r',data=self.r,compression='gzip')
            grp.create_dataset('log',data=cPickle.dumps(log_contents))
            grp.create_dataset('pym',data=cPickle.dumps(pymhandler.getvalue()))
        h5f.close()
#====================================end tadmodel

def consecutiveDistanceByProbability(r1,r2,p,xcontact=2):
    """
    Upper bound distance constraints for consecutive domains
    return surface to surface distance.
    
    Parameters
    -----------
    r1,r2 : int
        Radius for beads
    p : float Probability for contact
    xcontact : int
        scaling of (r1+r2) where a contact is defined. By default, 
        center to center distance D = 2*(r1+r2) is defined as contact.
    """
    if p > 0:
        d = (r1+r2)*(1. + (xcontact**3-1)/p)**(1./3.)
    else:
        d = 100*(r1+r2) # just a big number
    return d-r1-r2 # surface to surface distance

def centerOfMass(chain):
    xyzm = np.zeros((len(chain.get_particles()),4))
    i = -1
    for p in chain.get_particles():
        i += 1
        pattr = IMP.core.XYZR(p)
        xyzm[i,3] = pattr.get_radius()**3
        xyzm[i,0] = pattr.get_x()*xyzm[i,3]
        xyzm[i,1] = pattr.get_y()*xyzm[i,3]
        xyzm[i,2] = pattr.get_z()*xyzm[i,3]
    #---
    mass = sum(xyzm[:,3])
    return (sum(xyzm[:,0])/mass,sum(xyzm[:,1])/mass,sum(xyzm[:,2])/mass)

def surfaceDistance(chain,pair):
    '''
    calculate surface distance for a particle index pair in chain container 
    '''
    p1 = chain.get_particles()[pair[0]]
    p2 = chain.get_particles()[pair[1]]
    p1attr = IMP.core.XYZR(p1)
    p2attr = IMP.core.XYZR(p2)
    x1,y1,z1,r1 = [p1attr.get_x(),p1attr.get_y(),p1attr.get_z(),p1attr.get_radius()]
    x2,y2,z2,r2 = [p2attr.get_x(),p2attr.get_y(),p2attr.get_z(),p2attr.get_radius()]
    return np.linalg.norm([x1-x2,y1-y2,z1-z2])-r1-r2


#============================i/o
def readCoordinates(filename,prefix):
    h5f = h5py.File(filename,'r')
    try:
        grp = h5f[prefix]
    except:
        raise RuntimeError, "Group name %s doesn't exist!"%(prefix)
    xyz = h5f[prefix]['xyz'][:]
    r   = h5f[prefix]['r'][:]
    h5f.close()
    return xyz,r
    

