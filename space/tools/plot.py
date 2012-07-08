import sys

class Plot(object):
    """Represents a Plot object, containing the data to plot and describing attributes. """

    def __init__(self, title, data, xlog=False, ylog=False, xlabel='', ylabel='', plottype='lines', yrange=None, xrange=None):
        self.title = title
        self.keytitle = title
        self.data = data
        self.xlog = xlog
        self.ylog = ylog
        self.plottype = plottype
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xrange = xrange
        self.yrange = yrange
        self.plotfile = ''
        ### sanitize input
        if self.xlog and self.xrange != None:
            self.xrange = (max(self.xrange[0], 1), self.xrange[1])
        if self.ylog and self.yrange != None:
            self.yrange = (max(self.yrange[0], 1), self.yrange[1])
            
        
    def _plot(self, plot, data, spread=False, number=1, plot_count=1):
        """Writes the raw plotting data into the plotting stream."""
       
        if type(self.data) == dict:
            if len(self.data.keys()) > 0:
                plot += '\"-\" with %s title \"%s\",' % (self.plottype, self.keytitle)
                index = self.data.keys()
                index.sort()
                for k in index:
                    if spread:
                        data  += (str((k * (number + 2)) + plot_count) + ' ' + str(self.data[k]) + '\n')
                    else:
                        data  += (str(k) + ' ' + str(self.data[k]) + '\n')
                data += 'e\n'
            else:
                plot += '\"-\" with %s title \"%s (no data)\",' % (self.plottype, self.keytitle)
                data += '0 0\ne\n'
        elif type(self.data) == list:
            if len(self.data) > 0:
                plot += '\"-\" with %s title \"%s\",' % (self.plottype, self.keytitle)
                value_offset = 0
                #xtics = ''
                for value in self.data:
                    if spread:
                        data  += (str((value_offset * (number + 2)) + plot_count) + ' ' + str(value) + '\n')
                        value_offset += 1
                    else:
                        data  += (str(value) + '\n')
                data += 'e\n'
            else:
                plot += '\"-\" with %s title \"%s (no data)\",' % (self.plottype, self.keytitle)
                data += '0 0\ne\n'
        
        return (plot, data)

    def set_keytitle(self, title):
        """Sets the string presented in the plot's key to >>title<<. """
        self.keytitle = title

    def set_plotfile(self, fname):
        """Associates a filename to the plotlist. """
        self.plotfile = fname

    def set_plottype(self, ptype='lines'):
        """Associates a filename to the plotlist. """
        self.plottype = ptype

    def len(self):
        """Return length of data list."""
        return len(self.data)

    def plot(self, plotpath, terminal='png'):
        """Plots the Plot object."""
        import os
        import subprocess
        import tempfile
        
        plot_tmp_name = tempfile.mkstemp()[1]
        try:
            plot_tmp = open(plot_tmp_name, 'w')
        except:
            print >> sys.stderr, 'ERROR: Could not open %s for writing' % plot_tmp_name
            sys.exit(1)

        print >> plot_tmp, 'set terminal %s size 1024,1024 enhanced' % terminal
        print >> plot_tmp, 'set output \"' + plotpath + '\"'
        print >> plot_tmp, 'reset'
        print >> plot_tmp, 'set title \"%s\"' % self.title
        print >> plot_tmp, 'set xlabel \"%s\"' % self.xlabel
        print >> plot_tmp, 'set ylabel \"%s\"' % self.ylabel
        print >> plot_tmp, 'set key'

        if self.xlog:
            print >> plot_tmp, 'set logscale x'
        if self.ylog:
            print >> plot_tmp, 'set logscale y'
        if self.xrange != None:
            print >> plot_tmp, 'set xrange [%i:%i]' % (self.xrange[0], self.xrange[1])
        if self.yrange != None:
            print >> plot_tmp, 'set yrange [%i:%i]' % (self.yrange[0], self.yrange[1])


        plot = 'plot '
        data = '\n'

        (plot, data) = self._plot(plot, data)

        print >> plot_tmp, plot[:-1], data

        plot_tmp.close()
        subprocess.call(['gnuplot', plot_tmp_name])
        os.remove(plot_tmp_name)
 
class Multiplot(object):
    """Object storing several data objects, under one title."""

    def __init__(self, plot=None, spread=False):
        self.list = []
        if plot != None:
            self.add(plot)
        else:
            self.title = ''
            self.xlabel = ''
            self.ylabel = ''
            self.xlog = False
            self.ylog = False
            self.xrange = None
            self.yrange = None
            self.plottype = 'lines'
        self.spread = spread

    def add(self, plot):
        """Add new Plot object to Multiplot."""
        
        if len(self.list) == 0:
            self.title = plot.title
            self.xlabel = plot.xlabel
            self.ylabel = plot.ylabel
            self.xlog = plot.xlog
            self.ylog = plot.ylog
            self.xrange = plot.xrange
            self.yrange = plot.yrange
            self.plottype = plot.plottype
        self.list.append(plot)

    def remove(self, plot):
        """Remove given Plot from the Multiplot object"""

        if plot in self.list:
            self.list.remove(plot)

    def _plot(self, plot_tmp, column=0, colsize=1.0, row=0, rowsize=1.0):
        """Writes the raw data for plotting the whole multiplot to the plot stream."""

        row_offset = 1.0 - rowsize

        plot = 'plot '
        data = '\n'

        plot_count = 1
        print >> plot_tmp, 'reset'
        print >> plot_tmp, 'set title \"%s\"' % self.title
        print >> plot_tmp, 'set size %f, %f' % (colsize, rowsize)
        print >> plot_tmp, 'set origin %f, %f' % (float(column) * float(colsize), row_offset - (float(row) * float(rowsize)))
        print >> plot_tmp, 'set xlabel \"%s\"'  % self.xlabel
        print >> plot_tmp, 'set ylabel \"%s\"' % self.ylabel
        print >> plot_tmp, 'set key'

        if self.xlog:
            print >> plot_tmp, 'set logscale x'
        if self.ylog:
            print >> plot_tmp, 'set logscale y'

        if self.xrange != None:
            if self.spread:
                print >> plot_tmp, 'set xrange [%i:%i]' % (self.xrange[0] - 1, self.xrange[1] + 1)
            else:
                print >> plot_tmp, 'set xrange [%i:%i]' % (self.xrange[0], self.xrange[1])
        if self.yrange != None:
            print >> plot_tmp, 'set yrange [%i:%i]' % (self.yrange[0], self.yrange[1])

        if self.spread:
            print >> plot_tmp, 'set boxwidth 0.8 absolute'
            
            ### modify xrange
            if self.xrange != None:
                self.xrange = (self.xrange[0] - 1, self.xrange[1] + 1)

            ### set tics
            max_len = 0
            for pl in self.list:
                max_len = max(max_len, pl.len())
            xtics = '('
            for i in range(max_len):
                xtics += ('\"' + str(i + 1) + '\" ' + str(((len(self.list) + 2) * (i + 2)) - ((len(self.list) + 3) / 2)) + ',')
            xtics = (xtics[:-1] + ')')
            print >> plot_tmp, 'set xtics %s' % xtics
            print >> plot_tmp, 'set xrange [%i:%i]' % (len(self.list) - 1, (len(self.list) + 2) * (max_len + 1) + 1)

        for pl in self.list:
            if pl.plottype == 'boxes':
                pl.plottype = 'boxes fs solid 0.7'
            (plot, data) = pl._plot(plot, data, self.spread, len(self.list), plot_count)
            plot_count += 1

        print >> plot_tmp, plot[:-1], data

    def get_list(self):
        """Return list of plot objects. """
        return self.list

    def set_spread(self, spread=True):
        """Set flag: spread histograms."""
        self.spread = spread

    def set_xlog(self, xlog=True):
        """Set flag: log for x axis."""
        self.xlog = xlog
        if self.xlog and self.xrange != None:
            self.xrange = (max(self.xrange[0], 1), self.xrange[1])

    def set_ylog(self, ylog=True):
        """Set flag: log for y axis."""
        self.ylog = ylog
        if self.ylog and self.yrange != None:
            self.yrange = (max(self.yrange[0], 1), self.yrange[1])

    def set_xrange(self, xrange=(0, 1)):
        """Redefine xrange."""
        self.xrange = xrange

    def set_yrange(self, yrange=(0, 1)):
        """Redefine yrange."""
        self.yrange = yrange

    def set_xlabel(self, label=''):
        """Redefine xlabel."""
        self.xlabel = label

    def set_ylabel(self, label=''):
        """Redefine ylabel."""
        self.ylabel = label

    def set_plottype(self, ptype='lines'):
        """Redefine plottype."""
        self.plottype = ptype
        for pl in self.list:
            pl.set_plottype = ptype

    def plot(self, plotpath, terminal='png'):
        """Plots the Multiplot object. """        
        import os
        import subprocess
        import tempfile
        
        plot_tmp_name = tempfile.mkstemp()[1]
        try:
            plot_tmp = open(plot_tmp_name, 'w')
        except:
            print >> sys.stderr, 'ERROR: Could not open %s for writing' % plot_tmp_name
            sys.exit(1)

        print >> plot_tmp, 'set terminal %s size 1024,1024 enhanced' % terminal
        print >> plot_tmp, 'set output \"' + plotpath + '\"'
        print >> plot_tmp, 'set multiplot'
    
        self._plot(plot_tmp)        

        print >> plot_tmp, 'unset multiplot'

        plot_tmp.close()
        subprocess.call(['gnuplot', plot_tmp_name])
        #subprocess.call(['cat', plot_tmp_name])
        os.remove(plot_tmp_name)
 

class Plotlist(object):
    """Object storing a list of different Plot objects"""
    
    def __init__(self, plot=None):
        self.list = []
        if plot != None:
            self.list.append(plot)

    def set_plotfile(self, fname):
        """Sets plotfile-attribute of all plots in the list to fname."""
        for pl in self.list:
            pl.set_plotfile(fname)

    def append(self, plot):
        """Append new Plot object to list. """
        self.list.append(plot)
    
    def iter(self):
        """Return list's iterator. """
        return self.list.__iter__()

    def remove(self, plot):
        """Remove given Plot object from list, if list contains it."""
        if plot in self.list:
            self.list.remove(plot)

    def _plot(self, plot_tmp, spread=False, column=0, colsize=1.0, row=0, rowsize=1.0, rows=2):
        """Writes the raw data for plotting the whole list to the plot stream."""

        from math import ceil

        cell = 0
        column = 0
        row = 0

        colsize = float(1) / ceil(float(len(self.list)) / rows) 
        rowsize = float(1) / float(rows)
        row_offset = 1.0 - rowsize

        self.list.sort(key=lambda x: x.title)

        for pl in self.list:
            plot = 'plot '
            data = '\n'

            print >> plot_tmp, 'reset'
            print >> plot_tmp, 'set title \"%s\"' % pl.title
            print >> plot_tmp, 'set size %f, %f' % (colsize, rowsize)
            print >> plot_tmp, 'set origin %f, %f' % (float(column) * float(colsize), row_offset - (float(row) * float(rowsize)))
            print >> plot_tmp, 'set xlabel \"%s\"' % pl.xlabel
            print >> plot_tmp, 'set ylabel \"%s\"' % pl.ylabel
            print >> plot_tmp, 'set key'

            if pl.xlog:
                print >> plot_tmp, 'set logscale x'
            if pl.ylog:
                print >> plot_tmp, 'set logscale y'

            (plot, data) = pl._plot(plot, data, spread, len(self.list))

            print >> plot_tmp, plot[:-1], data
        
            cell += 1
            column = (cell / rows )
            row = (cell % rows)

    def plot(self, plotpath, spread=False, rows=2, terminal='png'):
        """Plots all objects from the list. """        
        import os
        import subprocess
        import tempfile
        from math import ceil
        
        plot_tmp_name = tempfile.mkstemp()[1]
        try:
            plot_tmp = open(plot_tmp_name, 'w')
        except:
            print >> sys.stderr, 'ERROR: Could not open %s for writing' % plot_tmp_name
            sys.exit(1)

        _unit = ''
        _color = ''
        if terminal == 'png':
            _factor = 512
        elif terminal == 'postscript':
            _factor = 20
            _color = 'enhanced color'
            _unit = 'cm'
        else:
            _factor = 100
        print >> plot_tmp, 'set terminal %s size %i%s,%i%s %s' % (terminal, ceil(float(len(self.list)) / float(rows)) * 512, _unit, 2 * _factor, _unit, _color) #1024'
        print >> plot_tmp, 'set output \"' + plotpath + '\"'
        print >> plot_tmp, 'set multiplot'
    
        self._plot(plot_tmp, spread, rows)        

        print >> plot_tmp, 'unset multiplot'

        plot_tmp.close()
        subprocess.call(['gnuplot', plot_tmp_name])
        os.remove(plot_tmp_name)

class Plotdict(object):
    """Object storing differnt Multiplots in a dictionary. Keys are plot-titles and values are Multiplots. """    

    def __init__(self, plot=None):
        self.dict = dict()
        if plot != None:
            self.dict[plot.title] = Multiplot(plot)

    def __getitem__(self, key):
        return self.dict[key]

    def add(self, plot, label=None):
        """Add new plot to the dictionary or merge to existing plotlist."""
        if label != None:
            plot.set_keytitle(label)

        try:
            self.dict[plot.title].add(plot)
        except KeyError:
            self.dict[plot.title] = Multiplot(plot)

    def add_plotlist(self, plotlist, label=None):
        """Add new plotlist to the dictionary or merge elements, if plotlist with same title already exists. The plotlist is under a key does NOT have to be unique."""

        for plot in plotlist.iter(): 
            self.add(plot, label)

    def remove(self, plot):
        """Remove plotlist from dictionary, if title exists in the key list."""
        if self.dict.has_key(plot.title):
            del self.dict[plot.title]

    def keys(self):
        """Return keys of Plotdict."""
        return self.dict.keys()

    def plot(self, plotpath, rows=2, terminal='png'):
        """Plot all plotlists in the dictionary."""

        import os
        import subprocess
        from math import ceil
        import tempfile
        
        plot_tmp_name = tempfile.mkstemp()[1]
        try:
            plot_tmp = open(plot_tmp_name, 'w')
        except:
            print >> sys.stderr, 'ERROR: Could not open %s for writing' % plot_tmp_name
            sys.exit(1)

        _unit = ''
        _color = ''
        if terminal == 'png':
            _factor = 512
        elif terminal == 'postscript':
            _factor = 20
            _unit = 'cm'
            _color = 'enhanced color'
        else:
            _factor = 100
        print >> plot_tmp, 'set terminal %s size %i%s,%i%s %s' % (terminal, ceil(float(len(self.dict.keys())) / float(rows)) * 512, _unit, 2 * _factor, _unit, _color) #1024'
        print >> plot_tmp, 'set output \"' + plotpath + '\"'
        print >> plot_tmp, 'set multiplot'

        cell = 0
        column = 0
        row = 0

        colsize = float(1) / ceil(float(len(self.dict.keys())) / rows) 
        rowsize = float(1) / float(rows)

        for key in self.dict.keys():
            self.dict[key]._plot(plot_tmp, column, colsize, row, rowsize)
            cell += 1
            column = (cell / rows )
            row = (cell % rows)

        print >> plot_tmp, 'unset multiplot'

        plot_tmp.close()
        subprocess.call(['gnuplot', plot_tmp_name])
        os.remove(plot_tmp_name)

