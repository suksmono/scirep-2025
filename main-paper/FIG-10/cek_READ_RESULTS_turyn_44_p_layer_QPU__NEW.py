########################################################
# Read results from IBMQ sites/queue
########################################################

from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_provider import IBMProvider
##
import matplotlib.pyplot as plt
from collections import Counter
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
from prb_qaoa import *
#--- insert path of used module
import sys
sys.path.append('../PRB')
#
#-------------------------------------
# convert binary {0,1} -> spin {-1, +1}
# -- variables
def b2s(x):
    if (x<1):
        s= 1.0
    if (x>0):
        s= -1.0
    return s
#
############################################################
# Evaluasi langsung terhadap fungsi objectiv
# masukan: bit string mis. '01110...'
############################################################
def obj_si_44(tv):
    # calculate the value of objective function
    # i.e. Hamiltonian's energy of a solution/bitstring
    # order 44--> 5 qubits
    # /////////// harusnya konversi [0,1] -> [-1, +1] ...??????
    ## >> revised
    ##-H44->> NQ=5 (non-hybrid)
    s0=b2s(int(tv[0]))
    s1=b2s(int(tv[1]))
    s2=b2s(int(tv[2]))
    s3=b2s(int(tv[3]))
    s4=b2s(int(tv[4]))
    ##################################### 
    # ------------------------------------------------------------------------
    obj= (s0*s1*s2 + s0*s3*s4 + s0*s3 - s0*s4 + s1*s2*s3*s4 + s1*s2*s3 
          + 2*s1*s2 + s1*s3*s4 + s1*s3 + s1*s4 + s1 + s2*s3*s4 + s2*s3 
          + s2*s4 + s2 + s4 + 5) 
    # ------------------------------------------------------------------------
    obj=abs(obj) ## since we seek for a minimum
    return obj
###
#mytoken='2ad695180e8492a69d5c74a7b1a402074af1513dd49c760fc87c5ae0d025894583dc16807fa7754d8a471b88f6a2221125a2ff2d33f28842ecf4af72d243bccd'
mytoken='73ad71c9f04043225a8b48c9896fe16a4df38bf5325d2dd39eeeeb235fc30184b444e3857c81b6fbc66e72713e8f4fb09b8df341bbdbe09a2369a8a2ab32aef7'
##
if __name__ == '__main__':
    #from qiskit.visualization import plot_histogram
    import matplotlib.pyplot as plt
    IBMProvider.save_account(token=mytoken, overwrite=True)
    # Initialize the account first.
    service = QiskitRuntimeService()
    ##supply job_id for old jobs already finished #####
    #job_id = 'cmcpg48pd0wg0081cbz0'
    job_id='ch2dp481n9ktavnkf4h0'
    job = service.job(job_id)
    qpu_counts=job.result().get_counts()
    #plot_histogram(qpu_counts)
    ##
    #
    
    qcntr=Counter(qpu_counts)
    #NSOL= int(len(qpu_counts)/2.0) #20*25 # 10 #int(2**NQ/10)
    NSOL= int(len(qpu_counts)/1.0) #20*25 # 10 #int(2**NQ/10)
    vbs=qcntr.most_common(NSOL) # get solution candidates
    # /// scan all solution candidates
    cnt_sol=0  # solution counter
    cnt_freq=0 # frequency of solution
    TOT_SOL=0
    min_obj=1000*1000
    TOT_ERR=0
    for m in range(0,NSOL): 
        tbs=vbs[m][0] # bitstring of solutions
        #---
        tfr=vbs[m][1] # frequency
        TOT_SOL=TOT_SOL+tfr
        vobj=obj_si_44(tbs) #     obj_si_bs
        TOT_ERR = TOT_ERR + vobj
        #---
        if vobj<min_obj:
          min_obj=vobj
          min_sol=tbs
        #
        if vobj<1:
          cnt_sol=+1
          cnt_freq=cnt_freq+tfr
          sol=tbs
    #--
    if cnt_sol<1:
      print('NO EXACT solution among the most frequent')
      sol=min_sol
    else:
      print('Number of correct solution found->', cnt_freq, 'among', TOT_SOL)
    #--
    bs=sol
    print('tbs->', bs)
    print('mean error->', TOT_ERR/NSOL)
    #// problem untuk order 68 
    tt=[]
    Nbs=len(bs)
    for nn in range(0, Nbs):
      sx=b2s(int(bs[nn])) # no reverse
      tt.append(sx)
    #--
    
    [X,Y,Z,W]=construct_XYZW_normal(tt)
    H=construct_hadamard(X,Y,Z,W)
    DG=np.matmul(H, H.transpose())
    #/// show image of Indicator-matrix
    plt.imshow(DG.astype(float), cmap="gray")
    plt.show()
    #/// show image of H-matrix
    plt.imshow(H.astype(float), cmap="gray")
    plt.show()
    
    #plot_histogram(qpu_counts)
    #######################################
    ## plot histogram
    # Total shots
    total = sum(qpu_counts.values())
    
    # Sort the labels
    sorted_labels = sorted(qpu_counts.keys())
    
    # Get probabilities in sorted order
    probs = [100*qpu_counts[label] / total for label in sorted_labels]
    
    # Plot
    #'''
    plt.figure(figsize=(10,6))
    plt.bar(sorted_labels, probs, width=0.8)
    
    # Label the axes
    plt.xlabel('Measurement Outcome (Bitstring)')
    plt.ylabel('Probability (%)')
    #plt.title('Measurement Results')
    plt.ylim(0, 1.25*max(probs))
    #plt.grid(axis='y', linestyle='--', alpha=0.7)
    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right')  # rotate 45 derajat dan rata kanan
    
    # Add text labels on top of bars
    for i in range(len(sorted_labels)):
        plt.text(sorted_labels[i], probs[i] + 0*0.02, f'{probs[i]:.2f}',
                 ha='center', va='bottom', fontsize=9, rotation= 60)
    plt.show()



###---## /////////////////////////////////

'''
# Plot histogram
total_shots = sum(qpu_counts.values())
probabilities = {k: v / total_shots for k, v in qpu_counts.items()}
#fig = plot_histogram(probs, bar_labels=False, figsize=(10, 6))  # Perbesar figure agar tidak terlalu padat
#fig = plot_histogram(qpu_counts, bar_labels=False, figsize=(10, 6))  # Perbesar figure agar tidak terlalu padat

fig = plot_histogram(probabilities, bar_labels=False, figsize=(10, 6))  # Perbesar figure agar tidak terlalu padat

# Ambil axis dan tambah label manual
ax = plt.gca()

# Tambahkan nilai probabilitas di atas setiap batang
for rect in ax.patches:
    height = rect.get_height()
    if height > 0.005:  # Jangan tampilkan label jika terlalu kecil
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            height,
            f'{height:.2f}',  # Format probabilitas
            ha='center',
            va='bottom',
            fontsize=8,
            rotation=0
        )

# Atur label sumbu X agar tidak bertumpuk
plt.xticks(rotation=45, fontsize=8)
plt.tight_layout()
plt.xlabel('Measurement Outcome (Bitstring)')
plt.ylabel('Probability')
plt.show()
'''