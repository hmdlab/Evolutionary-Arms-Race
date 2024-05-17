import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import baltic as bt

def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):

    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.01
        # *** is p < 0.001
        # etc.
        text = 'n.s.'
        p_list = [0.001, 0.01, 0.05]
        asterix = ['***', '**', '*']

        for p, a in zip(p_list, asterix):

            if data < p:

                text = a
                break
        

        #while data < p:
        #    text += '*'
        #    p /= 10.

        #    if maxasterix and len(text) == maxasterix:
        #        break

        #if len(text) == 0:
        #    text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)


def return_yaxis(self,x_attr=None,y_attr=None,target=None,size=None,colour=None,circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,
            zorder=None,outline=None,outline_size=None,outline_colour=None,**kwargs):

    if target==None: target=lambda k: k.is_leaf()
    if x_attr==None: x_attr=lambda k:k.x
    if y_attr==None: y_attr=lambda k:k.y
    if size==None: size=40
    if colour==None: colour='k'
    if zorder==None: zorder=3

    if outline==None: outline=True
    if outline_size==None: outline_size=lambda k: size(k)*2 if callable(size) else size*2
    if outline_colour==None: outline_colour=lambda k: 'k'

    if inwardSpace<0: inwardSpace-=self.treeHeight

    circ_s=circStart*math.pi*2
    circ=circFrac*math.pi*2

    allXs=list(map(x_attr,self.Objects))
    if normaliseHeight==None: normaliseHeight=lambda value: (value-min(allXs))/(max(allXs)-min(allXs))
    linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

    df=[]
    names=[]
    xs=[]
    ys=[]
    colours=[]
    sizes=[]

    outline_xs=[]
    outline_ys=[]
    outline_colours=[]
    outline_sizes=[]
    
    for k in filter(target,self.Objects):
        xs.append(x_attr(k))
        ys.append(y_attr(k))
        colours.append(colour(k)) if callable(colour) else colours.append(colour)
        sizes.append(size(k)) if callable(size) else sizes.append(size)

        if outline:
            outline_xs.append(xs[-1])
            outline_ys.append(ys[-1])
            outline_colours.append(outline_colour(k)) if callable(outline_colour) else outline_colours.append(outline_colour)
            outline_sizes.append(outline_size(k)) if callable(outline_size) else outline_sizes.append(outline_size)

        df.append([k.name, x_attr(k), y_attr(k), colours[-1]])

    df = pd.DataFrame(df, columns=['repeat adjusted name', 'x', 'y', 'color'])

    return df


# return peak summit position in repeat alignment
def return_position_in_repeat_alignment(name, position, map_dict):

    if name not in map_dict.keys() or position==0:
        return np.nan
    
    else:
        return_posi = map_dict[name][position]
        
        if return_posi != '-':
            return return_posi
        
        else:
            posi_plus = position
            posi_minus = position

            return_posi = map_dict[name][posi_plus]
            while return_posi == '-':

                posi_plus += 1

                if posi_plus in map_dict[name].keys():

                    return_posi = map_dict[name][posi_plus]

                else:
                    break
            
            while return_posi == '-':

                posi_minus -= 1

                if posi_minus in map_dict[name].keys():

                    return_posi = map_dict[name][posi_minus]

                else:
                    print(name)
                    break
            
            return float(return_posi)