import matplotlib.pyplot as plt
import numpy as np


'''
-----------------------------------------------------------------------------------------------------------------------------------
'''
###
fname='out_p.txt'
model='cox'
###
xT=[]
yT=[]
ypT=[]
yppT=[]
cT=[]
ybT=[]
dbT=[]
devT=[]
xdevT=[]
comp_dev=[]
log_goal=[]

with open(fname, 'r') as f:
    for line in f:
        if "START_NEW" in line:# and len(xT)==0:
            xT.append([])
            yT.append([])
            ypT.append([])
            yppT.append([])
            cT.append([])
            ybT.append([])
            dbT.append([])
            devT.append([])
            xdevT.append([])
        if 'df101' in line:
#            print(line + " 222")
            strline = line.replace("(","").replace(")","").replace(" ","").split(",")
            if True:
                xT[-1].append(1e-6*float(strline[1]))
                yT[-1].append(float(strline[3]))
                ypT[-1].append(float(strline[4]))
                yppT[-1].append(float(strline[5]))
                cT[-1].append(float(strline[2]))
                ybT[-1].append(float(strline[6]))
                dbT[-1].append(float(strline[7]))
                try:
                    devT[-1].append(float(strline[8]))
                    xdevT[-1].append(1e-6*float(strline[1]))
                except:
                    pass
#                else:
#                    print(float(strline[2]))
        if "[1]" in line and len(xT)>0:
            temp=line.split()
            if len(temp)==2:
                try:
                    comp_dev.append(float(temp[1]))
                except:
                    pass
            elif len(temp)==3:
                try:
                    log_goal.append(float(temp[2]))
                except:
                    pass
        pass


all_cmap='rainbow'
for j in range(len(xT)):
    fig = plt.figure(figsize=(17,10))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    x=xT[j]
    y=yT[j]
    yp=ypT[j]
    ypp=yppT[j]
    c=cT[j]
    yb=ybT[j]
    db=dbT[j]
    c_vmin=min(c)
    c_vmax=max(c)
    x,y,yp,ypp,c,yb,db = zip(*sorted(zip(x,y,yp,ypp,c,yb,db)))
    for ax in [ax1,ax2,ax3,ax4]:
        ax.grid(which='both',axis='both')
#        ax.set_xlim([0,100])
#        ax.hlines([0],min(x),max(x),color='black',alpha=.15)
    #
    if j==0:
        pass
#        ax1.set_ylim([-1e3,-300])
#        ax2.set_ylim([-2.5e2,2.5e2])
#        ax3.set_ylim([-1e3,1e3])
#        ax4.set_ylim([-1,1])
    elif j==1:
        pass
#        ax1.set_ylim([-84000,-82000])
#        ax2.set_ylim([-2.5e3,2.5e3])
#        ax3.set_ylim([-5e5,5e5])
#        ax4.set_ylim([-.01,.01])
    elif j==2:
        pass
#        ax1.set_ylim([-91000,-90000])
#        ax2.set_ylim([-2.5e2,2.5e2])
#        ax3.set_ylim([-5e2,5e2])
#        ax4.set_ylim([-.01,.01])
    im1=ax1.scatter(x[4:],y[4:],c=c[4:],s=10,cmap=all_cmap,vmin=c_vmin, vmax=c_vmax)
#    ax1.plot(x,y,color='black',alpha=.5)
#    ax1.plot(x,yb,color='red',alpha=.75)
    ax2.scatter(x,yp,c=c,s=10,cmap=all_cmap,vmin=c_vmin, vmax=c_vmax)
    ax2.plot(x,yp,color='black',alpha=.5)
    ax3.scatter(x,ypp,c=c,s=10,cmap=all_cmap,vmin=c_vmin, vmax=c_vmax)
    ax3.plot(x,ypp,color='black',alpha=.5)
    ax4.scatter(x,db,c=c,s=10,cmap=all_cmap,vmin=c_vmin, vmax=c_vmax)
    ax4.plot(x,db,color='black',alpha=.5)
    ##)
    ##

    ax1.set_xlabel("Time (sec)")
    ax2.set_xlabel("Time (sec)")
    ax3.set_xlabel("Time (sec)")
    ax4.set_xlabel("Time (sec)")

    ax1.set_ylabel("Log-Likelihood")
    ax2.set_ylabel("Log-Likelihood First Derivative")
    ax3.set_ylabel("Log-Likelihood Second Derivative")
    ax4.set_ylabel("Parameter Change")



    plt.tight_layout()
    plt.colorbar(im1,ax=[ax2,ax4],label="Step Number in Iteration")
    plt.savefig("PLOT/FULL_"+model+"_big_out_{}.png".format(j),bbox_inches='tight')
    plt.show()
    ##
    fig=plt.figure(figsize=(14,7))
    ax=fig.add_subplot(1,1,1)
    #
    dur=[]
    labels=['Initial\nPrep Matrices','Initial\nPrep Risks','Initial\nCalculation','Initial\nSteps','Initial\nRevert']+['Iteration\nPrep Matrices','Iteration\nUpdate Risk','Iteration\nCalculation','Iteration\nSteps','Iteration\nRevert']
    #
    z_counts=0
    #
    Ip=0.0
    Iu=0.0
    Ic=0.0
    Ic=0.0
    Ir=0.0
    In=0.0
    #
    s=0.0
    Iin=0
    #
    for i in range(len(x)):
        if c[i]==0:
            if z_counts==0:
                dur.append(x[0])
                z_counts+=1
            elif z_counts<3:
                dur.append(x[i]-x[i-1])
                z_counts+=1
            elif z_counts==3:
                if Iin>0:
                    dur.append(s/Iin)
                else:
                    dur.append(-.1)
                s=0.0
                Iin=0
                dur.append(x[i]-x[i-1])
                z_counts+=1
            else:
                if (z_counts-4)%4==0:
                    In+=1
                    Ip+=x[i]-x[i-1]
                elif (z_counts-4)%4==1:
                    Iu+=x[i]-x[i-1]
                elif (z_counts-4)%4==2:
                    Ic+=x[i]-x[i-1]
                else:
                    Ir+=x[i]-x[i-1]
                z_counts+=1
        else:
            if z_counts==3:
                s+=x[i]-x[i-1]
                Iin+=1
            else:
                s+=x[i]-x[i-1]
                Iin+=1
#    plt.plot(dur)
    #
    dur.append(Ip/In)
    dur.append(Iu/In)
    dur.append(Ic/In)
    if Iin>0:
        dur.append(s/Iin)
    else:
        dur.append(-.1)
    dur.append(Ir/In)
    #
#    print(Ip)
    #
    bars=ax.bar(range(len(dur)),dur,tick_label=labels,color='black')
    ax.bar_label(bars)
    ax.grid(which='both',axis='y',alpha=.5)
    ax.set_ylabel("Duration $(s)$")
    fig.suptitle("Average Section Duration\n{:.4f} seconds, {:.0f} Iterations".format((x[-1]-x[0]),In))
    plt.tight_layout()
    plt.savefig("PLOT/FULL_"+model+"_bar_{}.png".format(j),bbox_inches='tight')
    plt.show()
    try:
        dev=devT[j]
        xdev=xdevT[j]
        xdev,dev = zip(*sorted(zip(xdev,dev)))
        fig=plt.figure(figsize=(14,7))
        ax=fig.add_subplot(1,1,1)
        ax.plot(xdev,dev,'r:')
        try:
            cdev=comp_dev[j]
            ax.plot([xdev[0],xdev[-1]],[cdev,cdev],'k-')
        except:
            pass
        ax.semilogy()
        ##
        ax.set_xlabel("Time (seconds)")
        ax.set_ylabel("Standard Deviation of Error")
        plt.savefig("PLOT/FULL_"+model+"_conv_{}.png".format(j),bbox_inches='tight')
        ##
        plt.show()
        ##
    except:
#        clog=log_goal[j]
        fig=plt.figure(figsize=(14,7))
        ax=fig.add_subplot(1,1,1)
        ax.plot(x[3:],yb[3:],'r:')
        try:
            clog=log_goal[j]
            ax.plot([x[3],x[-1]],[clog,clog],'k-')
        except:
            pass
        ##
        ax.set_xlabel("Time (seconds)")
        ax.set_ylabel("Log-Likelihood")
        plt.savefig("PLOT/FULL_"+model+"_conv_{}.png".format(j),bbox_inches='tight')
        ##
        plt.show()
pass
