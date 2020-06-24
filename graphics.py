import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from Bio import Phylo

def graph_domains(t_domains, keys):
    """
    Represent the most common domains for the set of proteins 
    that align with the query sequence.
    """
    if len(t_domains) % 2 == 0:
        nrows = len(t_domains)//2
        ncols = len(t_domains)//nrows 
    else:
        nrows = len(t_domains)//2
        ncols = len(t_domains)//nrows + 1

    plt.figure(figsize=(9, 7))

    for i in range(len(t_domains)):
        number = nrows * 100 + ncols * 10 + 1 + i
        uniq_p = set(t_domains[i])
        uniq_p = list(uniq_p)

        d_p = dict()
        for item in uniq_p:
            d_p[item] = t_domains[i].count(item)

        l_p = list()
        for key,value in d_p.items():
            l_p.append((value,key))
        l_p.sort(reverse=True)

        domains = [item[1] for item in l_p]
        counts  = [item[0] for item in l_p]

        plt.subplot(number)
        plt.bar(domains[:5], counts[:5])
        plt.xticks(rotation=35, ha = "right", fontsize=6)
        plt.tight_layout()
        plt.ylabel("Number of proteins")
        plt.xlabel("Type of domain")
        plt.title(keys[i])

    plt.suptitle('Most common domains', fontsize=12, style = "italic", 
                                        weight = "bold", in_layout = True)
    plt.savefig("Results/Figures/Domains.png")
    plt.show()

    return "Results/Figures/Domains.png"
    
def graph_blastp(input_file, querys):
    """
    Plot the value of coverage and identity for each hit.
    The e value would be represented by the size of each point.
    """
    fig, ax = plt.subplots()

    for query in querys:
        input_handle = open(input_file, "r")
        data_list = list()
        for line in input_handle:
            columns = line.split("\t")

            if columns[0] == query:
                    data_list.append((float(columns[3]), float(columns[2]), float(columns[4])))

        input_handle.close()

        data_list.sort()

        e = [i[2] for i in data_list]
        for i in range(len(e)): #smaller evalue means the size of the point is smaller too
            if e[i] == 0.0 or e[i] < 1e-150:
                e[i] = 25
            elif e[i] < 1e-100:
                e[i] = 50
            elif e[i] < 1e-50:
                e[i] = 75
            elif e[i] < 1e-30:
                e[i] = 100
            else:
                e[i] = 150

        data = {'x': [i[0] for i in data_list],
        'y': [i[1] for i in data_list],
        'e': e
        }
        
        ax.scatter('x', 'y', s='e', data=data, label=query, alpha=0.6, edgecolors='none')
    
    #create a legend with all evalues
    e = [25, 50, 75, 100, 150]
    x = np.linspace(70,75,5)
    y = np.linspace(50,55,5)
    scatter = ax.scatter(x, y, s=e, alpha=0) #this is not shown in the final graphic

    #add a legend with evalues
    handles, labels = scatter.legend_elements(prop = "sizes", alpha=0.3)
    labels = ['< 1e-150', '< 1e-100', '< 1e-50', '< 1e-30', '> 1e-30']
    legend2 = ax.legend(handles, labels, loc="upper right", title="E-value")
    ax.add_artist(legend2)
    
    ax.legend(loc="upper left")
    plt.xlim(right = 100)
    plt.ylim(top = 100)
    plt.xlabel("Coverage")
    plt.ylabel("Identity")
    plt.title("Blastp results", fontsize=12, style = "italic",
                                weight = "bold", in_layout = True)
    plt.savefig("Results/Figures/Blastp.png")
    plt.show()

    return "Results/Figures/Blastp.png"

def graph_tree(querys):
    """
    Drow phylogenetic trees for each query.
    """
    
    for i in range(len(querys)):
        tree = Phylo.read("Results/" + querys[i] + "_tree.nw", "newick")
        tree.ladderize() 
        matplotlib.rc('font', size = 8)
        fig = plt.figure(figsize = (12,12))
        axes = fig.add_subplot(1,1,1)
        
        Phylo.draw(tree, axes = axes, do_show=False)
        plt.savefig("Results/Figures/" + querys[i] + "_tree.png")
        plt.close()
    return "Results/Figures/[query]_tree.png"