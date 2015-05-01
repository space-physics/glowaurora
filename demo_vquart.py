from numpy import tile,roots,real,spacing,unique,sqrt
from numpy.testing import assert_allclose
#
from vquart import vquartmod
jmax = 120
#%% demo the quartic root finder
def demovquart(A):

    """
    A: N x 5 array of coefficients for polynomial a4*x**4 + a3*x**3 + a2*x**2 + a1*x**1 + a0
    """
    return vquartmod.vquart(tile(A,(jmax,1)),jmax)
    
def pyvquart(A):
    return roots(A)

if __name__ == '__main__':
    """
    find quartic roots 
    """
    A = [-1,0,0,0,1]
    #A =[-221552219390.86914, -237693.87553329580, 145299.46383845527, 0.93531342393315953,     2.5711612024742349E-010]
    #A =[-481611022949.21875, -1626383.7460428476, 8.2219608898981278, 1.4729792838469140E-004, 3.7933283123395760E-011]
    #A=[1,6*sqrt(10),2,0,1]    
    #A=[360,-42,-41,2,1]
    
    rglow = demovquart(A)[0]
    rpy = pyvquart(A[::-1])
    print('python roots: {}'.format(rpy))
    rpyreal = unique(real(rpy[real(rpy)>spacing(1)]))
    assert_allclose(rpyreal,rglow)
