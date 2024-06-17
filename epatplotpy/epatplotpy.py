

# epatplot_graph.py

import matplotlib.pyplot as plt
import os
import glob
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import sys
import warnings
from astropy.io import fits
from pysas.pyutils import pyutils
from datetime import datetime



def load_files(path_to_file):
    '''
    Loads the temporary text files into numpy 1D arrays.
    
    Args:
        path_to_file: the path to one file.
        
    Output:
        data: the data as an array.
    '''

    try:
        data = np.loadtxt(path_to_file)
    except FileNotFoundError:
        raise FileNotFoundError('Could not open file {}. Quitting...'.format(path_to_file))
    
    return data


def return_text():
    '''
    Returns the basic text obtained by reading some temporary files.

    Args:
        None
    Output:
        txt: the text ready and formatted.
    '''

    
    try:
        with open('info_plot.txt') as ip:
            txt = ip.read()
        os.remove('info_plot.txt')
    except FileNotFoundError:
        warnings.warn('Could not open info_plot.txt')
        txt = ''
    
    try:
        with open('err_info.txt') as ei:
            txt2 = ei.read()
        txt = txt + txt2
        os.remove('err_info.txt')
    except FileNotFoundError:
        warnings.warn('Could not open err_info.txt.')
        txt = txt + ''

    return txt


def gather_info(in_file):
    '''
    Gathers info related to the input file to be used in the text.

    Args:
        in_file: the input file.

    Output:
        instr, obsid, expid, src_n, ra, dec, timedel, object_name
    '''

    try:
        with fits.open(in_file) as ff:
            telescope = pyutils.get_key_word(ff, 'TELESCOP')
            instr = pyutils.get_key_word(ff, 'INSTRUME')
            submode = pyutils.get_key_word(ff, 'SUBMODE')
            filter_id = pyutils.get_key_word(ff, 'FILTER')
            livetime = pyutils.get_key_word(ff, 'LIVETIME')
            if livetime == 'unknwon':
                livetime = -999
            rev = pyutils.get_key_word(ff, 'REVOLUT')
            obj_name = pyutils.get_key_word(ff, 'OBJECT')
    except FileNotFoundError:
        warnings.warn('Could not open file {}.'.format(in_file))
        return('unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown')

    
    return(telescope, rev, instr, submode, filter_id, livetime, obj_name)


def get_new_frame():
    '''
    Returns the coordinates for the new axis in order to plot the models.

    Args:
        none

    Output:
        tuple containing (min_x, max_x, min_y, max_y)
    '''

    with open('tmp_wcs') as tw:
        lines = tw.readlines()

    for line in lines:
        if 'xmin' in line:
            min_x = float(line.split(':')[1])
        if 'xmax' in line:
            max_x = float(line.split(':')[1])
        if 'ymin' in line:
            min_y = float(line.split(':')[1])
        if 'ymax' in line:
            max_y = float(line.split(':')[1])

    return (min_x, max_x, min_y, max_y)


def main(filename, device, infile, xlabel_tag):
    '''
    Main function.

    Args:
        filename: the name of the output file.
        device: the original display input.
        infile: the original input file.
        xlabel_tag: the tag that will define the label in the X axis.

    Output:
        none.
    '''
   
    ax2_txt = return_text()
    #text = return_text(infile, GTI_file)
    x_data_files = glob.glob('xdatahis_*')
    y_data_files = glob.glob('ydatahis_*')
    x_data_files.sort()
    y_data_files.sort()
    x_hist_files = glob.glob('xhis*')
    y_hist_files = glob.glob('yhis*')
    x_hist_files.sort()
    y_hist_files.sort()
    x_line_files = glob.glob('x*line*')
    y_line_files = glob.glob('y*line*')
    x_line_files.sort()
    y_line_files.sort()

    if xlabel_tag == 'ADU1':
      xlabel = 'ADU Channel'
    elif xlabel_tag == 'ADU5':
        xlabel = 'ADU Channel [5 eV]'
    elif xlabel_tag == 'PIEV':
        xlabel = 'PI Channel [eV]'


    telescope, rev, instr, submode, filter_id, livetime, obj_name = gather_info(infile)

    # Get the range parameters if they have been set
    if os.path.isfile('tmp_wcs'):
        min_x, max_x, min_y, max_y = get_new_frame()

    if len(x_hist_files) != 0:
        n_pages = 2
    else:
        n_pages = 1

    with PdfPages(filename) as pdf:
        for pages in range(0, n_pages):
            plt.style.use('bmh')
            if pages == 0:
                # first plot
                fig, ax = plt.subplots(1)
                fig.text(0.9, 0.90, os.path.basename(infile), ha='right', 
                                                        fontsize = 6)
                #ax.set_title('{}/{} {} {}'.format(telescope, instr, submode, filter_id), loc = 'left', fontsize = 10)
                fig.text(0.05, 0.93, '{} / {} {} {}'.format(telescope, instr, submode, filter_id), fontsize = 10, weight='bold')
                fig.text(0.05, 0.97, '{} {}, pattern statistics:'.format(telescope, instr), fontsize = 6)
                fig.text(0.9, 0.97, pdf, fontsize = 5, ha='right')
                if livetime > 0:
                    fig.text(0.7, 0.84, 'Livetime: {} (s).'.format(round(livetime, 1)), fontsize = 7)
                if rev != 'unknown':
                    fig.text(0.9, 0.94, 'Rev.{}'.format(rev),  ha='right', fontsize = 6)
                if obj_name != 'unknown':
                    fig.text(0.9, 0.92, obj_name, fontsize = 6,  ha='right')
                try:
                    with open('energy_f_tmp') as en_f:
                        e_text = en_f.read()
                    e_text = e_text.strip().rstrip()
                    fig.text(0.18, 0.13, e_text.strip().rstrip(), fontsize = 6, color = 'grey')
                    os.remove('energy_f_tmp')
                except FileNotFoundError:
                    warnings.warn('Could not open "energy_f_tmp.')
                    pass

                axlabs=[] # RDS: make a list of the data labels for plot 1
                axlabels_dict = dict()
                if len(x_data_files) > 0:
                    for i in range(0, len(x_data_files)):
                        x = load_files(x_data_files[i])
                        y = load_files(y_data_files[i])
                        y = np.power([10] * len(y), y)
                        x = np.power([10] * len(x), x)
                        n_label = x_data_files[i][-1]
                        ax.step(x, y, label = n_label, linewidth = 0.7)
                        axlabs.append(n_label)
                        axlabels_dict.update({n_label: 'C{}'.format(i)})
                #fig.legend(loc = 'right')
                ax.tick_params(axis = 'x', bottom=True, top = True)
                ax.tick_params(axis = 'y', left=True, right = True)
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_xlabel(xlabel)
                # RDS: set the xrange as needed 
                if os.path.isfile('tmp_wcs'):
                    ax.set_xlim(np.power([10], min_x), np.power([10], max_x))
                    ax.set_ylim(np.power([10], min_y), 
                                        np.power([10], max_y))
                #ax.set_xlim(80.0,21000.0)
# RDS: Put pattern legend on axis rather than in a legend box - axis 1
                # Put the legend of each coloured line on the Y axis
                yboxes=[]
                # Loop over each pattern label and construct the full text
                for i in range(0, len(axlabs)):
                    yboxes.append(TextArea(axlabs[i]+"/adu   ", 
                              textprops=dict(color=axlabels_dict[axlabs[i]], 
                              size=9,rotation=90,ha='left',va='bottom')))
#
                # If data files existed and hance labels exist then 
                # add them into the legend
                if (len(yboxes) > 0):
                    ybox = VPacker(children=yboxes,align="bottom", 
                                                        pad=0, sep=5)
#                  # Display the legend
                    anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0., 
                       frameon=False, bbox_to_anchor=(-0.08, 0.0), 
                       bbox_transform=ax.transAxes, borderpad=0.)
                    ax.add_artist(anchored_ybox)
# End RDS
            else: # The second page
                fig, ax = plt.subplots(1)
                colour_id = 0
                hist_flag = False
                labels_dict = dict()
                
                # Loop over each hist file
                for i in range(0, len(x_hist_files)):

                    n_label = x_hist_files[i].split('_')[2]
                    if n_label == 'sd':
                        n_label = 's+d'
                    labels_dict.update({n_label: 'C{}'.format(colour_id)})
                    colour_id = colour_id + 1

                    # Ignore if the file is empty
                    if os.path.getsize(x_hist_files[i]):
                       x = load_files(x_hist_files[i])
                       y = load_files(y_hist_files[i])
                       x = np.power([10] * len(x), x)
                       if len(x) > 0:
                          hist_flag = True
                       ax.step(x, y, label = n_label, linewidth = 0.7, color = labels_dict[n_label])
                    
                # RDS: set the xrange as needed 
                if os.path.isfile('tmp_wcs'):
                    ax.set_xlim(np.power([10], min_x), np.power([10], max_x))
                    ax.set_ylim(-0.05, 1.05) # Fixed vals
                    os.remove('tmp_wcs')
                #.set_xlim(80.0,21000.0)
                ax.tick_params(axis = 'y', left=True, right = True)
                ax.tick_params(axis = 'x', bottom=True, top = True)

                ax2 = ax.twiny()
                for i in range(0, len(x_line_files)):
                    x = load_files(x_line_files[i])
                    y = load_files(y_line_files[i])
                    n_label = x_line_files[i].split('_')[2]
                    if n_label == 'sd':
                        n_label = 's+d'
                    x = np.power([10] * len(x), x)
                    ax2.plot(x, y, color = labels_dict[n_label], linewidth = 0.7)
                
                #fig.legend(loc = 'right')    ! Remove the legend plot - now on axis
                ax.set_xscale('log')
                ax.grid(None)
                
                hist_flag=True # Get rid

                if hist_flag:
                    ax2.axis('off')
                    ax.set_xlabel(xlabel)
                else:
                    ax.set_xticks([], minor = True)
#                    ax.tick_params(axis = 'x', bottom=True, top = True)
                    ax2.tick_params(axis = 'x', top = True)
                    #ax2.xaxis.tick_bottom()
                    ax2.set_xlabel(xlabel)
                    ax2.xaxis.set_label_position('bottom')
#                    ax.tick_params(axis = 'x', which = 'both', bottom=False, top = False, labelbottom = False)

                
                # Show the SRC coordinates
                try:
                    with open('src_info_tmp') as src_f:
                        e_text = src_f.read()
                        e_text = e_text.replace('\x00', '')
                    fig.text(0.10, 0.13, e_text, fontsize = 6, color = 'grey')
                    os.remove('src_info_tmp')
                except FileNotFoundError:
                    pass

# RS: Write the username, date, time info
                try:
                    uname = os.getlogin( ) 
                    now = datetime.now()
                    dt_string = now.strftime("%d-%b-%Y %H:%M")
                    fig.text(0.98,0.02,'{}    {}'.format(
                     uname,dt_string), fontsize = 5, ha='right')
                except:
                    pass
# RDS: Make the Y axis values smaller
                #ax2.tick_params(axis = 'y', which='major', labelsize=2)
# End RDS

# RDS: Put pattern legend on axis rather than in a legend box
                # Get the labels
                mylabs = list(labels_dict.keys())
                # Put the legend of each coloured line on the Y axis
                yboxes=[]
                # Loop over each pattern label and construct the full text
                for i in range(0, len(mylabs)):
                    yboxes.append(TextArea(mylabs[i]+"       ", 
                              textprops=dict(color=labels_dict[mylabs[i]], 
                              size=12,rotation=90,ha='left',va='bottom')))

                ybox = VPacker(children=yboxes,align="bottom", 
                                                        pad=0, sep=5)
                # Display the legend
                anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0., 
                       frameon=False, bbox_to_anchor=(-0.08, 0.15), 
                       bbox_transform=ax.transAxes, borderpad=0.)
                ax2.add_artist(anchored_ybox)
# End RDS

                ax2.set_xscale('log')
                ax.set_xscale('log')
                fig.text(0.16, 0.375, ax2_txt, fontsize = 7)

                plt.grid(False)

            plt.grid(None)

            if 'XW' in device.upper():
                plt.show()
            else:
                pdf.savefig()

            plt.close()
            

    for f in x_data_files:
        os.remove(f)
    for f in y_data_files:
        os.remove(f)
    for f in x_hist_files:
        os.remove(f)
    for f in y_hist_files:
        os.remove(f)
    for f in x_line_files:
        os.remove(f)
    for f in y_line_files:
        os.remove(f)

    return None