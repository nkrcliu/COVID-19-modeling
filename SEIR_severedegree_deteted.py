import numpy as numpy
import scipy as scipy
import scipy.integrate

#定义动力学微分方程组
def system_defes(t, variables, beta1, beta2, beta3, beta1d, beta2d, beta3d, alpha, alphad, thetae, psie, gamma1, p1,
                 thetai1, psii1, p1d, gamma1d, gamma2, p2, thetai2, psii2, p2d, gamma2d, gamma3, mu, thetai3, psii3,
                 mud, gamma3d):
    #设置变量
    S, E, I1, I2, I3, D_E, D_I1, D_I2, D_I3, R, D = variables

    N =S+E+I1+I2+I3+D_E+D_I1+D_I2+D_I3+R+D

    #动力学方程组
    dS = -beta1*I1*S-beta2*I2*S-beta3*I3*S-beta1d*D_I1*S-beta2d*D_I2*S-beta3d*D_I3
    dE = beta1*I1*S+beta2*I2*S+beta3*I3*S+beta1d*D_I1*S+beta2d*D_I2*S+beta3d*D_I3-alpha*E-thetae*psie*E
    dD_E = thetae*psie*E-alphad*D_E
    dI1 = alpha*E-gamma1*I1-p1*I1-thetai1*psii1*I1
    dD_I1 = thetai1*psii1*I1+alphad*D_E-p1d*D_I1-gamma1d*D_I1
    dI2 = p1*I1-gamma2*I2-p2*I2-thetai2*psii2*I2
    dD_I2 = thetai2*psii2*I2+p1d*D_I1-p2d*D_I2-gamma2d*D_I2
    dI3 = p2*I2-gamma3*I3-mu*I3-thetai3*psii3*I3
    dD_I3 = thetai3*psii3*I3+p2d*D_I2-mud*D_I3-gamma3d*D_I3
    dR = gamma1*I1+gamma2*I2+gamma3*I3+gamma1d*D_I1+gamma2d*D_I2+gamma3d*D_I3
    dD = mu*I3+mud*D_I3

    #返回该变量，这种形式是适应solve_ivp的形式，这个方法比odeint要好一些
    return [dS, dE, dD_E, dI1, dD_I1, dI2, dD_I2, dI3, dD_I3, dR, dD]

#动力学方程数值解，totalstep是总步长，步长默认设置成1，可认为totalstep是天数，给step重新赋值可以改
def run(totalstep, odeset, initN, initS, initE, initD_E, initI1, initI2, initI3, initD_I1, initD_I2, initD_I3, initR, initD, step=1):

    #roundnum = totalstep/step
    #创建空组用于储存初值，未来考虑不同时期参数不同需要迭代的时候，可把这里定义全局变量（或者结构体）用于迭代
    Sset = []
    Eset = []
    D_Eset = []
    I1set = []
    D_I1set = []
    I2set = []
    D_I2set = []
    I3set = []
    D_I3set = []
    Rset = []
    Dset = []

    #空组中添加初值
    Sset.append(initS)
    Eset.append(initE)
    D_Eset.append(initD_E)
    I1set.append(initI1)
    D_I1set.append(initD_I1)
    I2set.append(initI2)
    D_I2set.append(initD_I2)
    I3set.append(initI3)
    D_I3set.append(initD_I3)
    Rset.append(initR)
    Dset.append(initD)

    t_eval = numpy.arange(start=0, stop=totalstep, step=step)

    #取数组（python里叫列表，一个东西）的最后一个值作为初始值，多期的时候迭代就是说这里，总取最后一个
    init_cond = [Sset[-1], Eset[-1], D_Eset[-1], I1set[-1], D_I1set[-1], I2set[-1], D_I2set[-1], I3set[-1],
                 D_I3set[-1], Rset[-1], Dset[-1]]
    #解初值问题
    solution = scipy.integrate.solve_ivp(odeset, t_span=[0, totalstep], y0=init_cond, t_eval=t_eval)
    print(solution)

    #拷贝数值解，转化为数组（列表）
    solutionset = numpy.array(solution['y'])
    print(solutionset)
    return solutionset

#计算每个时间点总的染病人数，算其他的可以加进去，画图也可以加
def caculatetotalnumber(solutionset):
    Sset = []
    Eset = []
    D_Eset = []
    I1set = []
    D_I1set = []
    I2set = []
    D_I2set = []
    I3set = []
    D_I3set = []
    Rset = []
    Dset = []

    Sset.append(solutionset[0])
    Eset.append(solutionset[1])
    D_Eset.append(solutionset[2])
    I1set.append(solutionset[3])
    D_I1set.append(solutionset[4])
    I2set.append(solutionset[5])
    D_I2set.append(solutionset[6])
    I3set.append(solutionset[7])
    D_I3set.append(solutionset[8])
    Rset.append(solutionset[9])
    Dset.append(solutionset[10])

    #计算总的时点上的感染人数总数
    currenttotalI = [I1set[i]+I2set[i]+I3set[i]+D_I1set[i]+D_I2set[i]+D_I3set[i] for i in range(0, len(Sset))]
    print(currenttotalI)
    return currenttotalI


#参数是瞎取的
odeset = lambda t ,X:system_defes(t, X, 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                                  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
solutionset = run(10, odeset, 10000, 10000, 0, 0, 1000, 0, 0, 0, 0, 0, 0, 0)
caculatetotalnumber(solutionset)
