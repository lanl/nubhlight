################################################################################
#                                                                              #
#  UTILITIES FOR PLOTTING                                                      #
#                                                                              #
# Authors: Jonah Miller (jonahm@lanl.gov) and Ben Ryan (brryan@lanl.gov)       #
#                                                                              #
# PURPOSE:                                                                     #
# To provide plotting routines in both Cartesian and Black Hole coordinates    #
# for nubhlight simulations. Handles re-ordering data when necessary.          #
#                                                                              #
################################################################################

from __future__ import print_function, division
import numpy as np
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt

# copy-pasted from stack overflow:
# https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

# GET XZ SLICE OF GRID DATA
def flatten_xz(array, hdr, flip=False):
    sign = 1.
    flat = np.zeros([2*hdr['N1'],hdr['N2']])
    for j in range(hdr['N2']):
        for i in range(hdr['N1']):
            flat[i,j] = sign*array[hdr['N1'] - 1 - i,j,hdr['N3']//2]
        for i in range(hdr['N1']):
            flat[i+hdr['N1'],j] = array[i,j,0]
    if flip:
        flat[:,0] = 0
        flat[:,-1] = 0
    return flat

# GET XY SLICE OF GRID DATA
def flatten_xy(array, hdr):
  if  hdr['stopx'][3] >= np.pi:
    return np.vstack((array.transpose(),array.transpose()[0])).transpose()
  else:
    return array.copy()

def plot_X1X2(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
  label=None, ticks=None, shading='gouraud', show_axes=False):
  X1 = geom['X1'][:,:,0]
  X2 = 1.-geom['X2'][:,:,0]
  mesh = ax.pcolormesh(X1, X2, var[:,:,0], cmap=cmap, vmin=vmin, vmax=vmax)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  if not show_axes:
      ax.set_xticklabels([]); ax.set_yticklabels([])

def plot_X1X3(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
  label=None, ticks=None, shading='gouraud',show_axes=False):
  j = dump['hdr']['N2']//2
  X1 = geom['X1'][:,j,:]
  X3 = geom['X3'][:,j,:]
  mesh = ax.pcolormesh(X1, X3, var[:,j,:], cmap=cmap, vmin=vmin, vmax=vmax)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  if not show_axes:
      ax.set_xticklabels([]); ax.set_yticklabels([])

def plot_xz(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True, 
            label=None, ticks=None, shading='gouraud',
            l_scale = None,
            reverse_x = False,
            reverse_z = False,
            fix_poles = False):
  x = geom['x']
  y = geom['y']
  z = geom['z']
  if dump['hdr']['N3'] > 1. and dump['hdr']['stopx'][3] >= np.pi:
    x = flatten_xz(x, dump['hdr'], flip=True)
    y = flatten_xz(y, dump['hdr'], flip=True)
    z = flatten_xz(z, dump['hdr'])
    var = flatten_xz(var, dump['hdr'])
    rcyl = np.sqrt(x**2 + y**2)
    rcyl[np.where(x<0)] *= -1
  else:
    x = x[:,:,0]
    y = y[:,:,0]
    z = z[:,:,0]
    var = var[:,:,0]
    rcyl = np.sqrt(x**2 + y**2)
  if reverse_x:
    rcyl *= -1.
  if reverse_z:
    z *= -1.
  if fix_poles:
    rcyl[:,0] = 0.
    rcyl[:,-1] = 0.
  ehr = dump['hdr']['Reh']
  if l_scale is not None:
    rcyl = l_scale*rcyl
    z = l_scale*z
    ehr = l_scale*ehr
  mesh = ax.pcolormesh(rcyl, z, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading=shading)
  circle1=plt.Circle((0,0),ehr,color='k'); 
  ax.add_artist(circle1)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  ax.set_xlabel('x/M'); ax.set_ylabel('z/M')
  return mesh
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)

def contour_xz(ax, geom, var, dump,
               levels = None,
               l_scale = None,
               reverse_x = False,
               reverse_z = False,
               fix_poles = False,
               **kwargs):
  x = geom['x']
  y = geom['y']
  z = geom['z']
  if dump['hdr']['N3'] > 1. and dump['hdr']['stopx'][3] >= np.pi:
    x = flatten_xz(x, dump['hdr'], flip=True)
    y = flatten_xz(y, dump['hdr'], flip=True)
    z = flatten_xz(z, dump['hdr'])
    var = flatten_xz(var, dump['hdr'])
    rcyl = np.sqrt(x**2 + y**2)
    rcyl[np.where(x<0)] *= -1
  else:
    x = x[:,:,0]
    y = y[:,:,0]
    z = z[:,:,0]
    var = var[:,:,0]
    rcyl = np.sqrt(x**2 + y**2)
  if reverse_x:
    rcyl *= -1.
  if reverse_z:
    z *= -1.
  if fix_poles:
    rcyl[:,0] = 0.
    rcyl[:,-1] = 0.
  if l_scale is not None:
    rcyl = l_scale*rcyl
    z = l_scale*z
    ehr = l_scale*erh
  if levels is not None:
    contour = ax.contour(rcyl,z,var,levels,**kwargs)
  else:
    contour = ax.contour(rcyl,z,var,**kwargs)
  ax.clabel(contour,inline=True)

  return contour

def plot_vxvz(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True, 
              label=None, ticks=None, shading='gouraud',
              l_scale = None,
              reverse_x = False,
              reverse_z = False,
              fix_poles = False):
  ucon_cart = dump['ucon_cart']
  x = geom['x']
  y = geom['y']
  z = geom['z']
  rcyl = np.sqrt(x**2 + y**2)
  vx = ucon_cart[:,:,:,1]/ucon_cart[:,:,:,0]
  vy = ucon_cart[:,:,:,2]/ucon_cart[:,:,:,0]
  vz = ucon_cart[:,:,:,3]/ucon_cart[:,:,:,0]
  vcyl = vx*x/rcyl + vy*y/rcyl
  if dump['hdr']['N3'] > 1.:
    vcyl = flatten_xz(vcyl, dump['hdr'], flip=True)
    vz = flatten_xz(vz, dump['hdr'])
    var = flatten_xz(var, dump['hdr'])
  else:
    vz = vz[:,:,0]
    var = var[:,:,0]
    vcyl = vcyl[:,:,0]
  if reverse_x:
    vcyl *= -1.
  if reverse_z:
    vz *= -1.
  if fix_poles:
    vcyl[:,0] = 0.
    vcyl[:,-1] = 0.
  ehr = dump['hdr']['Reh']
  if l_scale is not None:
    vcyl = l_scale*vcyl
    vz = l_scale*vz
    ehr = l_scale*ehr
  mesh = ax.pcolormesh(vcyl, vz, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading=shading)
  # circle1=plt.Circle((0,0),ehr,color='k'); 
  # ax.add_artist(circle1)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  #ax.set_aspect('equal')
  ax.set_xlabel(r'$v_r/c$'); ax.set_ylabel(r'$v_z/c$')
  return mesh
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)

# TODO: The quiver plots can be pretty difficult to read.
#       it would be good to have a generic mechanism for plotting
#       streamlines. Matplotlib can only natively do it on a Cartesian
#       grid, so we'd either need to generate our own, or interpolate.
#       either option is possible, but stinks. Think about this later.
#       ~JMM
def quiver_xz(ax, geom, dump, varx, varz, C=None,
              qk=None, qkpos = (0.8, 0.8), qklen = 1, qkloc='E',
              **kwargs):
    x = geom['x']
    y = geom['y']
    z = geom['z']
    if dump['hdr']['N3'] > 1.:
        x = flatten_xz(x, dump['hdr'], flip=True)
        y = flatten_xz(y, dump['hdr'], flip=True)
        z = flatten_xz(z, dump['hdr'])
        varx = flatten_xz(varx, dump['hdr'])
        varz = flatten_xz(varz, dump['hdr'])
        if C is not None:
            C = flatten_xz(C, dump['hdr'])
    else:
        x = x[:,:,0]
        y = y[:,:,0]
        z = z[:,:,0]
        varx = varx[:,:,0]
        varz = vary[:,:,0]
        if C is not None:
            C = C[:,:,0]
    rcyl = np.sqrt(x**2 + y**2)
    if C is not None:
        quiv = ax.quiver(rcyl,z,varx,varz,C,**kwargs)
    else:
        quiv = ax.quiver(rcyl,z,varx,varz,**kwargs)
    if qk is not None:
        qk = ax.quiverkey(quiv, qkpos[0], qkpos[1], qklen,
                          qk, labelpos=qkloc,
                          coordinates='figure')

def quiver_xy(ax, geom, dump, varx, vary, C=None,
              qk=None, qkpos = (0.8, 0.8), qklen = 1, qkloc='E',
              **kwargs):
    x = geom['x']
    y = geom['y']
    x = flatten_xy(x[:,dump['hdr']['N2']//2,:], dump['hdr'])
    y = flatten_xy(y[:,dump['hdr']['N2']//2,:], dump['hdr'])
    varx = flatten_xy(varx[:,dump['hdr']['N2']//2,:], dump['hdr'])
    vary = flatten_xy(vary[:,dump['hdr']['N2']//2,:], dump['hdr'])
    if C is not None:
        C = flatten_xy(color[:,dump['hdr']['N2']//2,:], dump['hdr'])
        quiv = ax.quiver(x,y,varx,vary,C,**kwargs)
    else:
        quiv = ax.quiver(x,y,varx,vary,**kwargs)
    if qk is not None:
        qk = ax.quiverkey(quiv, qkpos[0], qkpos[1], qklen,
                          qk, labelpos=qkloc,
                          coordinates='figure')

def overlay_field(ax, geom, dump, NLEV=20, linestyle='-', linewidth=1,
  linecolor='k'):
  from scipy.integrate import trapz
  hdr = dump['hdr']
  N1 = hdr['N1']; N2 = hdr['N2']
  x = flatten_xz(geom['x'], hdr).transpose()
  z = flatten_xz(geom['z'], hdr).transpose()
  A_phi = np.zeros([N2, 2*N1])
  gdet = geom['gdet'].transpose()
  B1 = dump['B1'].mean(axis=-1).transpose()
  B2 = dump['B2'].mean(axis=-1).transpose()
  print(gdet.shape)
  for j in range(N2):
    for i in range(N1):
      A_phi[j,N1-1-i] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx'][1]) - 
                         trapz(gdet[:j, i]*B1[:j, i], dx=hdr['dx'][2]))
      A_phi[j,i+N1] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx'][1]) - 
                         trapz(gdet[:j, i]*B1[:j, i], dx=hdr['dx'][2]))
  A_phi -= (A_phi[N2//2-1,-1] + A_phi[N2//2,-1])/2.
  Apm = np.fabs(A_phi).max()
  if np.fabs(A_phi.min()) > A_phi.max():
    A_phi *= -1.
  #NLEV = 20
  levels = np.concatenate((np.linspace(-Apm,0,NLEV)[:-1], 
                           np.linspace(0,Apm,NLEV)[1:]))
  ax.contour(x, z, A_phi, levels=levels, colors=linecolor, linestyles=linestyle,
    linewidths=linewidth)

def plot_xy(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
  label=None, ticks=None, shading='gouraud', fix_bounds=True):
  hdr = dump['hdr']
  x = geom['x']
  y = geom['y']
  x = flatten_xy(x[:,dump['hdr']['N2']//2,:], dump['hdr'])
  y = flatten_xy(y[:,dump['hdr']['N2']//2,:], dump['hdr'])
  if dump['hdr']['stopx'][3] < np.pi and fix_bounds:
      x[:,-1] = 0
      y[:,0] = 0
  var = flatten_xy(var[:,dump['hdr']['N2']//2,:], dump['hdr'])
  mesh = ax.pcolormesh(x, y, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading=shading)
  circle1=plt.Circle((0,0),dump['hdr']['Reh'],color='k'); 
  ax.add_artist(circle1)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  ax.set_xlabel('x/M'); ax.set_ylabel('y/M')
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)
  return mesh

def get_streamlines(geom,dump,nlines):
    """Integrates streamlines inward from outer boundary
    Assumes statistically average flow field.
    """
    from scipy import interpolate,integrate

    hdr = dump['hdr']
    X1,X2 = geom['X1'][:,0,0],geom['X2'][0,:,0]
    rcyl = geom['rcyl']
    z = geom['z']
    
    vcon = dump['ucon'][...,1:]/dump['ucon'][...,0,np.newaxis]
    vcon_interp = [interpolate.interp2d(X1,X2,vcon[:,:,0,i].T) \
                   for i in range(3)]
    
    rcyl_interp = interpolate.interp2d(X1,X2, rcyl[...,0].T)
    z_interp = interpolate.interp2d(X1, X2, z[...,0].T)
    
    def rhs_vec(t,x_vec):
        #Assumes x_vec is shape (nlines,2)
        vcon_vec = np.empty_like(x_vec)
        for l in range(nlines):
            for i in range(2):
                vcon_vec[l,i] = vcon_interp[i](x_vec[l,0],
                                               x_vec[l,1])
        vcon_vec *= -1.0
        return vcon_vec
    
    def rhs(t,x):
        x_vec = x.reshape((nlines,2))
        v_vec = rhs_vec(t,x_vec)
        v = v_vec.ravel()
        return v
    
    initial_data = np.array([[X1[-1],X2[(i+2)*hdr['N2']//(nlines+2)]]\
                             for i in range(0,nlines)])
    initial_data = initial_data.ravel()
    
    integrator = integrate.ode(rhs)
    integrator.set_initial_value(initial_data,0)
    
    tgrid = np.linspace(1,1e6,1e3)
    xygrid = np.empty((len(tgrid),nlines,2))
    for i,t in enumerate(tgrid):
        integrator.integrate(t)
        xygrid[i] = integrator.y.reshape((nlines,2))
    
    rcyl_interp = interpolate.interp2d(X1,X2,rcyl[...,0].T)
    z_interp = interpolate.interp2d(X1,X2,z[...,0].T)
    
    rcyl_path,z_path = np.empty_like(xygrid),np.empty_like(xygrid)
    for i in range(len(tgrid)):
        for l in range(nlines):
            rcyl_path[i,l] = rcyl_interp(xygrid[i,l,0],xygrid[i,l,1])
            z_path[i,l] = z_interp(xygrid[i,l,0],xygrid[i,l,1])
    
    return rcyl_path[...,0],z_path[...,0]
