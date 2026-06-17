import numpy as np


ith_para=4
n_ref = 8
skipnum = 5000000

mean_n_ref = np.zeros(n_ref)
var_n_ref = np.zeros(n_ref)

for i in range(n_ref):
    i = i
    #
    # compare 1ref
    fname = "%s%s%s" % ("./results/chain", i, ".dat")
    # compare multi-ref
    #
    #print(fname)
    fname_ch = "%s" % (fname)
    print(fname_ch)
    #quit()
    # read chain 
    i_ch_tup = np.genfromtxt(fname_ch, skip_header=skipnum)
    #print(i_ch_tup.shape)
    #### change to array
    i_ch = np.asarray(i_ch_tup)
    #print(i_ch)
    #print(i_ch[0,:])
    #print(i_ch[0,0])
    #print(i_ch.shape)
    #omgsum = i_ch[:,2] + i_ch[:,3]
    mean_n_ref[i] = i_ch[:,ith_para].mean()
    var_n_ref[i] = i_ch[:,ith_para].var()
    print("mean: ", mean_n_ref[i])
    print("var: ", var_n_ref[i])
    print("sig: ", np.sqrt(var_n_ref[i]))

all_mean = mean_n_ref.mean()
print(all_mean)

L = i_ch_tup.shape[0]
#print(L)

B1 = L/(n_ref-1)
#print(B1)
B12 = (mean_n_ref-all_mean)
#print(B12)
B2 = (mean_n_ref-all_mean)**2.0
#print(B2)
B3 = B2.sum()
#print(B3)
B = B1*B3
print(B)


W = 1/n_ref*var_n_ref.sum()
print("W", W)

var_theta=(1-1/L)*W + 1/L*B
print("var_theta", var_theta)

Rc=np.sqrt(var_theta/W)
print(Rc)
