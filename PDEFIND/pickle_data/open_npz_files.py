from numpy import load

data = load('Greedy_01_NCV_bisplines_bins_fisher_prune_deg_2.npz')
lst = data.files
for item in lst:
    print(item)
    print(data[item])

