import meep as mp
import numpy as np
import argparse
import math,random
import cmath
import sys
import datetime

sys.path.insert(0,'/Users/wdai11/function')
import my_output as my
import meep.materials as mat

#########################################
# common parameters
#########################################
class common:
    pols=["Ex",]
    geoms=["Empty","Wg"]
    fcens=np.arange(0.26,0.27,0.02)
    dfs=np.ones(fcens.size)*0.02

    pol="Ex"
    geom="Wg"
    fcen=0.260
    df=0.020

    nDBR=5
    a=1.0  # lattice constant
    t=1.0  # waveguide width
    dx=0    # location of the dipole source 
    dy=0    # (0,0) is beteen grooves
    dz=0.5
    
    
    sig='{}_F{:0.3f}'.format(geom,fcen)
    if geom != 'Empty':
        sig='{}_{}_F{:0.3f}'.format(geom,pol,fcen)
        pstring='_L{:d}_A{:0.2f}T{:0.2f}'.format(nDBR,a,t)
        dstring='_DX{:0.2f}Y{:0.2f}Z{:0.2f}'.format(dx,dy,dz)
        sig+=pstring+dstring


######################################### 
# simulation function
#########################################
def simulation_fun():
    resolution=50
    geom=common.geom
    sig=common.sig
    
    # source parameters
    fcen = common.fcen
    df=common.df
    df2=df*1.4           # source frequency width
    nfreq=31             # number of frequency bins

    k=False              # ensure PEC boundary condition

    mat_env=mp.Medium(epsilon=1)
    mat_wg=mp.Medium(epsilon=14.44)
    mat_metal=mat.Al # mp.Medium(epsilon=20)
 
    nDBR=common.nDBR
    mat_high=mp.Medium(epsilon=14.44)
    mat_low=mp.Medium(epsilon=10.24)
    d_high=0.25
    d_low=0.30

    pol=common.pol
    src_cmpt=eval('mp.{}'.format(pol))

    ############################################
    # begin empty and actual simulation
    if geom=='Empty':
        dpml=1.0           
        sx=1               
        sy=1     
        sz=1
        sxx=sx+2*dpml
        syy=sy+2*dpml
        szz=sz+2*dpml
        cell_size = mp.Vector3(sxx,syy,szz)
        boundary_layers = [mp.PML(thickness=dpml)]

        # geometry
        geometry=[mp.Block(size=mp.Vector3(sxx,syy,szz),
                           material=mat_wg)]

        # define point
        pt_src=mp.Vector3(0,0,0)
        pt_field=mp.Vector3(sx/4,sy/4,0)
        pt_flux=mp.Vector3(0,0,0)
        l_flux=0.1
        frs_src=my.FluxRegions_Cube_Center(l_flux,pt_flux)

        # sources
        sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df2), 
                             component=src_cmpt, 
                             center=pt_src,
                             size=mp.Vector3(0,0,0))]

        # symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]
        if src_cmpt==mp.Ex:
            symmetries=[mp.Mirror(mp.X,-1),mp.Mirror(mp.Y,+1),
                        mp.Mirror(mp.Z,+1)]
        elif src_cmpt==mp.Ey:
            symmetries=[mp.Mirror(mp.X,+1),mp.Mirror(mp.Y,-1),
                        mp.Mirror(mp.Z,+1)]
        elif src_cmpt==mp.Ez:
            symmetries=[mp.Mirror(mp.X,+1),mp.Mirror(mp.Y,+1),
                        mp.Mirror(mp.Z,-1)]
        else:
            return None

    else:
        ### actually simulations
        # from bottom up
        t_front=0.1
        t_metal=0.2
        t_wg=common.t
        t_DBR=nDBR*(d_high+d_low)
        t_back=0.2
        
        
        a=common.a
        dx=common.dx
        dy=common.dy
        dz=common.dz

        dpmlxy=a           # PML thickness
        dpmlz=1
        sx=4
        sy=12  
        sz=t_front+t_metal+t_wg+t_DBR+t_back
        sxx=sx+2*dpmlxy
        syy=sy+2*dpmlxy
        szz=sz+dpmlz
        cell_size = mp.Vector3(sxx,syy,szz)

        # pml
        boundary_layers = [mp.Absorber(thickness=dpmlxy, 
                                       direction=mp.X),
                           mp.Absorber(thickness=dpmlxy, 
                                       direction=mp.Y),
                           mp.PML(thickness=dpmlz, 
                                  direction=mp.Z,side=mp.High)]

    
        # geometry 
        tts=[t_front,t_metal,t_wg]
        ccs=my.center_from_thickness(tts,-szz/2)

        geometry=[mp.Block(size=mp.Vector3(2*sxx,2*syy,2*szz),
                           center=mp.Vector3(0,0,0),
                           material=mat_wg),
                  mp.Block(size=mp.Vector3(sxx,syy,tts[1]),
                           center=mp.Vector3(0,0,ccs[1]),
                           material=mat_metal), 
                  mp.Block(size=mp.Vector3(sxx,syy,tts[2]),
                           center=mp.Vector3(0,0,ccs[2]),
                           material=mat_wg), ]

        z0=-szz/2+t_front+t_metal+t_wg # waveguide/DBR interface
        for i in range(nDBR):
            # add low
            geometry.append(mp.Block(size=mp.Vector3(sxx,syy,d_low),
                                     center=mp.Vector3(0,0,z0+d_low/2),
                                     material=mat_low))
            z0+=d_low
            # add high
            geometry.append(mp.Block(size=mp.Vector3(sxx,syy,d_high),
                                     center=mp.Vector3(0,0,z0+d_high/2),
                                     material=mat_high))
            z0+=d_high

        # sources and flux box
        z0=-szz/2+t_front+t_metal # metal surface
        z1=z0+t_wg+t_DBR  # DBR surface

        pt_src=mp.Vector3(dx,dy,z0+dz)
        pt_field=mp.Vector3(0.08*sx/2,0.07*sy/2,z0+0.77*t_wg)
        l_flux,pt_flux=my.FluxBox3D(size=0.1,src=pt_src,
                                    zmin=z0,zmax=z0+t_wg)
        frs_src=my.FluxRegions_Cube_Center(l_flux,pt_flux)
        assert dz<t_wg,"The source is inside the waveguide."

        # big box, the box include the whole GaSb area
        frs_big=my.FluxRegions_Box_Corner(-sx/2,sx/2,-sy/2,sy/2,
                                          z0, z1)
        # near2far region
        near_sizes=np.arange(sy,sy-6.05,-1.0)
        near_z=z1+0.1 # above DBR
        near_fluxes=[mp.FluxRegion(center=mp.Vector3(0,0,near_z), 
                                   size=mp.Vector3(sx,nsize,0))
                     for nsize in near_sizes]
        near_regions=[mp.Near2FarRegion(center=mp.Vector3(0,0,near_z), 
                               size=mp.Vector3(sx,nsize,0))
                      for nsize in near_sizes]

        sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df2), 
                             component=src_cmpt, 
                             center=pt_src,
                             size=mp.Vector3(0,0,0))]
        
        # symmetries 
        if src_cmpt==mp.Ex:
            symmetries=[mp.Mirror(mp.X,-1),mp.Mirror(mp.Y,+1)]
        elif src_cmpt==mp.Ey:
            symmetries=[mp.Mirror(mp.X,+1),mp.Mirror(mp.Y,-1)]
        elif src_cmpt==mp.Ez:
            symmetries=[mp.Mirror(mp.X,+1),mp.Mirror(mp.Y,+1)]
   
    # end of big if between emtpy and geom
    ####################################
    print('Real simulation area '+str(mp.Vector3(sx,sy,sz)))
    print('Total simulation area '+str(mp.Vector3(sxx,syy,szz)))
    print('Source at '+str(pt_src))
    print('Field is monitored at '+str(pt_field))
    print('Flux center at '+str(pt_flux))
    print('Flux region size is '+str(l_flux))

    # set up simulation
    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        geometry=geometry,
                        boundary_layers=boundary_layers,
                        dimensions=3,
                        sources=sources,
                        k_point=k,
                        force_complex_fields=True,
                        symmetries=symmetries,
                        default_material=mat_wg,
                        filename_prefix=sig )
                        
    # calculate flux
    if geom=='Empty':
        trans_src=[sim.add_flux(fcen,df,nfreq,fr) for fr in frs_src]
    else:
        trans_src=[sim.add_flux(fcen,df,nfreq,fr) for fr in frs_src]
        trans_big=[sim.add_flux(fcen,df,nfreq,fr) for fr in frs_big]

        near_trans=[sim.add_flux(fcen,df,nfreq,fr) for fr in near_fluxes]
        near_fields=[sim.add_near2far(fcen, df, nfreq, rgn) 
                    for rgn in near_regions]
        # nearfield=sim.add_near2far(fcen, df, nfreq, *frs)

 
    #### define two functions to organize output
    def get_flux_data_from_many(trans):
        # simulation data are saved in trans
        freqs=mp.get_flux_freqs(trans[0])
        fluxes=[mp.get_fluxes(tran) for tran in trans]
        data=np.column_stack((freqs, *fluxes,)) 
        return data

    def get_flux_data_from_one(ntran,nsize):
        # simulation data are saved in trans
        freqs=mp.get_flux_freqs(ntran)
        fluxes=mp.get_fluxes(ntran)
        data=np.column_stack((freqs, fluxes,))
        data=np.insert(data,1,nsize,axis=1)
        return data

    def get_far_data_from_one(nearfield,nsize):
        # simulation data are saved in nearfield
        freqs=mp.get_near2far_freqs(nearfield)
        r = 1000*3.8  # 1000 wavelengths out from the source

        thetas=np.concatenate([np.arange(0.5,20,1),
                               np.arange(22.5,89,3)])  # polar angle
        phis=[45,]  # azimuthal angle
        data=np.empty(shape=[0,16])
        for phi in phis:
            for theta in thetas:
                phi0=np.radians(phi)
                theta0=np.radians(theta)
                pt=r*mp.Vector3(math.sin(theta0)*math.cos(phi0),
                                math.sin(theta0)*math.sin(phi0),
                                math.cos(theta0))
                fields=sim.get_farfield(nearfield,pt)
                mat=np.vstack((np.real(fields),np.imag(fields)))
                fields=mat.T.reshape((1,-1)).reshape((-1,12))
                
                infos=np.vstack((freqs,nsize*np.ones(len(freqs)),
                                 theta*np.ones(len(freqs)),
                                 phi*np.ones(len(freqs))  ))

                data=np.vstack(( data, np.hstack((infos.T,fields)) ))
        return data

        
                            
    ####
    # define step function to display fluxes
    my_display_count=0
    step_time_flush=10
    step_time_print=100
    step_time_terminate=20
    def my_display_fluxes(sim):
        nonlocal my_display_count
        nonlocal step_time_print
        print('='*40)
        print('='*40)
        # power flow around the source
        data=get_flux_data_from_many(trans_src)
        my.matrix_output(None,data,"{:10.3e}","flux_src")

        if geom != "Empty":
            print('<'*40)
            data=get_flux_data_from_many(trans_big)
            my.matrix_output(None,data,"{:10.3e}","flux_big")
            print('<'*40)

            # power flow through near region
            for idx,ntran in enumerate(near_trans):
                nsize=near_sizes[idx]
                data=get_flux_data_from_one(ntran,nsize)
                my.matrix_output(None,data,"{:10.3e}",
                                 "flux{:0.2f}".format(nsize))
                print('<'*40)
                print('='*40)
            # far fields
            for idx,nfield in enumerate(near_fields):
                nsize=near_sizes[idx]
                data=get_far_data_from_one(nfield,nsize)
                my.matrix_output(None,data,"{:10.3e}",
                                 "farfield{:4.2f}".format(nsize))
                print('>'*40)
      
        my_display_count+=1
        print('='*40)
        print('No. {} display at t={}'.format(
            my_display_count,step_time_print))
        my.my_flush_step(sim)
    
    # run simulations
    sim.run(# mp.at_beginning(mp.output_epsilon),
        mp.at_every(step_time_flush, my.my_flush_step),
        mp.at_every(step_time_print, my_display_fluxes),
        until_after_sources=mp.stop_when_fields_decayed(
            step_time_terminate, src_cmpt, pt_field, 1e-9))
    sys.stdout.flush()


    # power flow around the source
    data=get_flux_data_from_many(trans_src)
    fname=sig+"_T0.dat"
    my.matrix_output(fname,data,"{:10.6e}","flux_src")
    if geom != "Empty":
        data=get_flux_data_from_many(trans_big)
        fname=sig+"_T1.dat"
        my.matrix_output(fname,data,"{:10.6e}","flux_big")

        # power flow through near region
        for idx,ntran in enumerate(near_trans):
            nsize=near_sizes[idx]
            data=get_flux_data_from_one(ntran,nsize)
            fname=sig+"_S{:0.2f}_T.dat".format(nsize)
            my.matrix_output(fname,data,"{:10.6e}",
                             "flux{:0.2f}".format(nsize))
            # far field
            for idx,nfield in enumerate(near_fields):
                nsize=near_sizes[idx]
                data=get_far_data_from_one(nfield,nsize)
                fname=sig+"_S{:0.2f}_FF.dat".format(nsize)
                my.matrix_output(fname,data,"{:10.3e}",
                                 "farfield{:4.2f}".format(nsize))


######################### 
# main function
#########################
# 3D simulaitons
time0=datetime.datetime.now()
print('-'*50)
print(common.sig)
simulation_fun()

time1=datetime.datetime.now()
dtime=time1-time0
print("Simulation used {:0.2f} seconds".format(
    dtime.total_seconds()))
sys.stdout.flush()

