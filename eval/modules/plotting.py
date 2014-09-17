import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy as sp


def plot(data, plot_info, fname=None, ax=None, labels=None):
    """This is a wrapper function that handles the data plotting"""

    if plot_info[2] == 'log10':
        data = sp.log10(data + 1)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 16), dpi=300)

    if len(data.shape) > 1:

        cmap = plt.get_cmap('jet')
        norm = plt.Normalize(0, data.shape[0])  
        label_handles = []
        if plot_info[1] == 'plot':
            for f in range(data.shape[0]):
                if labels is not None:
                    tmp, = ax.plot(sp.arange(data.shape[1]), data[f, :], color=cmap(norm(f)), label=labels[f])
                    label_handles.append(tmp)
                else:
                    ax.plot(sp.arange(data.shape[1]), data[f, :], color=cmap(norm(f)))
        elif plot_info[1] == 'bar':
            w = 0.8 / data.shape[0]
            for f in range(data.shape[0]):
                if labels is not None:
                    tmp = ax.bar(sp.arange(data.shape[1]) + 0.1 + f*w, data[f, :], w, color=cmap(norm(f)), label=labels[f])
                    label_handles.append(tmp)
                else:
                    ax.bar(sp.arange(data.shape[1]) + 0.1 + f*w, data[f, :], w, color=cmap(norm(f)))
            ax.set_xticks(sp.arange(data.shape[1]) + 0.5)
            ax.set_xticklabels(sp.arange(data.shape[1]))
            #ax.set_xlim([0.5, ax.get_xlim()[1]])
        if labels is not None:
            #ax.legend(label_handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            ax.legend(label_handles, labels, loc=0, fontsize='small')
    else:
        if plot_info[1] == 'plot':
            ax.plot(sp.arange(data.shape[0]), data)
        elif plot_info[1] == 'bar':
            ax.bar(sp.arange(data.shape[0]) + 0.1, data)
            ax.set_xticks(sp.arange(data.shape[0]) + 0.5)
            ax.set_xticklabels(sp.arange(data.shape[0]))
            #ax.set_xlim([0.5, ax.get_xlim()[1]])
    ax.set_xlabel(plot_info[3]) 
    if len(plot_info[2]) > 0:
        ax.set_ylabel('%s %s' % (plot_info[2], plot_info[4]))
    else:
        ax.set_ylabel(plot_info[4]) 
    ax.set_title(plot_info[5])


