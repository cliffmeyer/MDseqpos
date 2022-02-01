"""
Script for plotting motif logos to png file
"""
import matplotlib.pyplot as plt
import pandas as pd
import logomaker 

def motif_list_to_dict(motif):
    """
    args:
        - motif: position-wise list of A,C,G,T prob lists eg [[0.5,0.5,0,0],[0,0,0.5,0.5]]
    returns:
        - dictionary eg {'A':[0.5,0.0],'C':[0.5,0.0],'G':[0.0,0.5],[0.0,0.5]}
    """
    INDEX_A,INDEX_C,INDEX_G,INDEX_T=0,1,2,3
    d = {'A':[],'C':[],'G':[],'T':[]}
    for pos in motif:
        d['A'].append(pos[INDEX_A])
        d['C'].append(pos[INDEX_C])
        d['G'].append(pos[INDEX_G])
        d['T'].append(pos[INDEX_T])
    return d


def pts_to_inches(pts):
    """
    Converts pts to inches
    """
    INCHES_PER_PT = 0.0138889
    return pts*INCHES_PER_PT


def plot_motif_info(*motifs,filename='test.png'):
    """
    Plot motif PWM and reverse complement in usual information format.
    args:
        - motifs: position-wise lists of A,C,G,T probability lists eg [[0.5,0.5,0,0],[0,0,0.5,0.5]]
        - filename: name of output png file
    """
    FIG_WIDTH_PTS = 660
    FIG_HEIGHT_PTS = 300
    num_rows = len(motifs)
    num_cols = 1

    SMALL_SIZE = 12
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 20

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig,axs = plt.subplots(num_rows,num_cols,tight_layout=True)
    if num_rows == 1:
        axs = [axs]

    for i,ax in enumerate(axs):
        motif_prob_dict = motif_list_to_dict(motifs[i])
        motif_prob_df = pd.DataFrame.from_dict(motif_prob_dict)
        motif_info_df = logomaker.transform_matrix(df=motif_prob_df,from_type='probability',to_type='information')

        motif_len = motif_info_df.shape[0]
        logomaker.Logo(motif_info_df,
            ax=ax,
            color_scheme= {'A':'lawngreen','C':'royalblue','G':'orange','T':'red'},
            font_name='Helvetica')

        # style using Axes methods
        ax.set_ylim([0,2])
        ax.set_ylabel('Information Content', labelpad=1)
        ax.set_xlabel('Position', labelpad=1)
        ax.set_yticks([0,0.5,1,1.5,2])
        ax.set_xticks(list(range(motif_len)))
        ax.set_xticklabels('%d'%x for x in list(range(1,motif_len+1)))
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fig.set_figwidth(pts_to_inches(FIG_WIDTH_PTS))
    fig.set_figheight(pts_to_inches(FIG_HEIGHT_PTS))
    fig.savefig(filename)
    plt.close(fig) 

def test():
    motif = [[1,0,0,0],[0.5,0.5,0,0],[0,0.5,0.5,0],[0,0,0,1]]
    motif_rc = [[0,0,0,1],[0,0.5,0.5,0],[0.5,0.5,0,0],[1,0,0,0]]
    plot_motif_info(5*motif,5*motif_rc,filename='test.png')


if __name__ == '__main__':
    test()


