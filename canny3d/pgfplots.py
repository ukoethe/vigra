import os,sys
import numpy as npy
import pdb


def writeTable(table,filename,names=None, asPairs=False):
    f = file(filename,"w")
    if names != None:
        for name in names:
            f.write("%-10s\t" % name)
        f.write("\n")
    if asPairs:
        for row in table:
            f.write("(")
            f.write("%-2s" % str(row[0]))
            for val in row[1:]:
                f.write(",%-10s" % str(val).replace(',','.')
)
            f.write(") ")
    else:
        for row in table:
            for val in row:
                f.write("%.10f\t" % val
)
            f.write("\n")

    f.close()

def pairs(table):
    ret = ""
    for row in table:
        ret += "("
        ret +="%-2s" % str(row[0])
        for val in row[1:]:
            ret +=",%-10s" % str(val)
        ret +=") "

    return ret


def plot(data,namebase,
         lnames=None,width=None,height=None,
         plotoptions=None,extraXTick=None,
         crange=None, name=None,title=None,yticksLabels=True):
    
    names = ["data%i" % i for i in range(len(data))]
    
    writeTable(npy.array(data).T,
               namebase + ".dat",
               names=names)

    f = open(namebase + ".tikz", "w") 
        

    if not yticksLabels:
        f.write(r"""\pgfplotsset{every y tick label/.append style={%
        inner sep=0pt,outer sep=0pt,minimum size=0pt}}%""")
    else:
        f.write(r"""\pgfplotsset{every y tick label/.append style={%
        inner sep=0pt,outer sep=1pt,minimum size=0pt}}%""")        
    f.write(r"""
    \pgfplotsset{every axis/.append style={
    font=\footnotesize,
    thin,
    tick style={ultra thin}}}
    \pgfplotsset{every axis label/.append style={
    draw,
    fill=green!20}}    
    %\pgfplotsset{
    %every axis plot post/.append style=
    %{opacity=0.5}}
\begin{tikzpicture}[remember picture]
    \begin{axis} [
    title style={text height=1.5ex,text depth=.25ex},""")
    if title!=None:
        f.write("\ntitle=%s,%%" % title)
    if name!=None:
        f.write("\nname=%s," % name)
    if width!=None:
        f.write("\nwidth=%s," % width)
    if height!=None:
        f.write("\nheight=%s," % height)
    if crange!=None:
        f.write("\nymin=%0.2f,%%" % (float(crange[0])))
        f.write("\nymax=%0.2f,%%" % (float(crange[1])))
    if not yticksLabels:
        f.write("\nyticklabel=\yticklabeltextnone,%%")
##    if ptype=="ybar":
##        f.write("\nybar,")
##        f.write("\nbar width=%s," % barwidth)
##    if ptype=="const":
##        f.write("\nconst plot,")
    if extraXTick != None:
        f.write("\nextra x ticks={%f}," % extraXTick['val'])
        f.write(r"""extra x tick style={grid=major,
        /pgfplots/tick label style={
        outer sep = 0.3cm
        }},""")
        f.write("extra x tick labels={%s}," % extraXTick['lab'])

    f.write(r"""
    legend style={nodes={anchor=west}}
    %%extra x ticks={-0.1}
    ]""")

    for i in range(len(data)-1):
        f.write(r"""
        \addplot%s
        table[x=data0,y=data%i] {%s};""" % \
                ("["+plotoptions+"]" if plotoptions != None else "",\
                 i+1,namebase + ".dat"))


    if lnames != None:
        lstr = lnames[0]
        for name in lnames[1:]:
            lstr += "," + name
        f.write(r"\legend{%s}" % lstr)

    f.write(r"""
    \end{axis}%%
    \end{tikzpicture}%%""")



def hist(data,namebase,
         lnames=None,nbins=100,brange=None,crange=None,normed=0,width=None,
         ptype="ybar",barwidth=None,height=None,drawStats=False,title=None,
         plotOptions=None,yticksLabels=True):

    
    if brange == None:
        brange = data[0].min(),data[0].max()

    hdatas = []
    hbins = None
    for i,d in enumerate(data):
        hdata,bins = npy.histogram(d,
                                   bins=nbins,
                                   normed=normed,
                                   range=brange)
        if hbins == None:
            hbins = bins
        else:
            assert (hbins == bins).all()

        print hdata.sum()
        hdatas.append(hdata)

    names = ["data%i" % i for i in range(len(hdatas))]
    names.append("bins")
    hdatas.append(hbins)    
    writeTable(npy.array(hdatas).T,
               namebase + ".dat",
               names=names)

    f = open(namebase + ".tikz", "w") 
        
    
    f.write(r"""
    \begin{tikzpicture}%
    \def\pgfplots@ytickalignnum{0}
    \pgfplotsset{every axis/.append style={%
    font=\footnotesize,%
    thin,%
    tick style={ultra thin}}}%
    \pgfplotsset{every axis label/.append style={%
    draw,%
    fill=green!20}}%
    \pgfplotsset{every x tick label/.append style={%
    /pgf/number format/.cd,fixed,precision=1,zerofill=false}}%
    \pgfplotsset{every tick/.append style={thin,black!70}}
    \pgfplotsset{every tick label/.append style={%""")
    #f.write("\nymajorticks=%s,"% ("true" if yticksLabels else "false"))
    f.write(r"""
    font=\tiny}}%""")
##     \pgfplotsset{%
##     every axis plot post/.append style={opacity=0.5}}%""")
    if not yticksLabels:
        f.write(r"""
        \pgfplotsset{every y tick label/.append style={%
        inner sep=0pt,outer sep=0pt,minimum size=0pt,/pgf/number format/.cd,fixed,precision=0,zerofill=false}}%""")
    else:
        f.write(r"""
        \pgfplotsset{every y tick label/.append style={%
        inner sep=0pt,outer sep=1pt,minimum size=0pt,/pgf/number format/.cd,int detect}}%""")        
    f.write(r"""
    \begin{axis}[%
    outer sep=0pt,
    tickwidth=2pt,
    title style={text height=1.5ex,text depth=.25ex},""")
    if title!=None:
        f.write("\ntitle=%s,%%" % title)
    if not yticksLabels:
        f.write("\nyticklabel=\yticklabeltextnone,%%")
    if width!=None:
        f.write("\nwidth=%s,%%" % width)
    if height!=None:
        f.write("\nheight=%s,%%" % height)
    if ptype=="ybar":
        f.write("\nybar,%%")
    if barwidth!=None:
        f.write("\nbar width=%s," % barwidth)
    else:
        if width != None:
            f.write("\nbar width=%s/%d," % (width,2*nbins))
    if ptype=="const":
        f.write("\nconst plot,%%")
    if crange!=None:
        f.write("\nymin=%0.2f,%%" % (float(crange[0])))
        f.write("\nymax=%0.2f,%%" % (float(crange[1])))
    f.write("\nenlargelimits=true,")
    f.write("\ntick align=inside,")

    
    f.write(r"""
    area cycle list,%
    legend style={nodes={anchor=west}},%
    name=mplot%
    %%extra x ticks={-0.1}
    ]%""")

    for i in range(len(hdatas)-1):
        f.write(r"""
        \addplot%s %%
        table[x=bins,y=data%i] {%s};%%""" % \
                ("["+plotOptions+"]%%" if plotOptions != None else "",\
                 i,namebase + ".dat"))

    if lnames != None:
        lstr = lnames[0]
        for name in lnames[1:]:
            lstr += "," + name
        f.write(r"\legend{%s}" % lstr)

    f.write(r"""
    \end{axis}%%""")
    if drawStats:
        assert len(data) == 1
        data = data[0]
        f.write(r"""
        \ifthenelse{\equal{\withTables}{1}}{%%
        \node[anchor=north,outer sep=1cm,
        draw=gray,
        drop shadow={opacity=0.25,shadow scale=0.95},
        fill=blue!4,
        rounded corners=1pt] at (mplot.south) {%%
        {\tiny
        \begin {tabular}{ll}
        mean& \pgfmathprintnumber{%s}\\
        \rowcolor[rgb]{1,1,1}
        var & \pgfmathprintnumber{%s}\\
        min & \pgfmathprintnumber{%s}\\
        \rowcolor[rgb]{1,1,1}
        max & \pgfmathprintnumber{%s}\\
        \end {tabular}}
        };}{}%%
        """ % (data.mean(),data.var(),data.min(),data.max()))

    f.write("\n\end{tikzpicture}%%")


def plots():
    ret = []
    for vals in ep.eindices:
        ret.append(npy.array(vals[1][0][1].values()))
    pgfplots.hist(ret,"../thesis/figs/hist2",bins=30,lnames=["\\tns{averageGradient}","\\tns{faceMean}","\\tns{emd}"],width="\\textwidth")

def bars(data,namebase,
         lnames=None,brange=None,crange=None,width=None,
         ptype="ybar",barwidth=None,height=None,title=None,
         plotOptions=None,yticksLabels=True):

    if brange == None:
        brange = data[0].min(),data[0].max()

    names = ["data%i" % i for i in range(len(data))]
    vals = pairs(npy.array(data).T)

    f = open(namebase + ".tikz", "w") 
        
    f.write(r"""
    \begin{tikzpicture}[remember picture]%
    \def\pgfplots@ytickalignnum{0}
    \pgfplotsset{every axis/.append style={%
    font=\footnotesize,%
    thin,%
    tick style={ultra thin}}}%
    \pgfplotsset{every axis label/.append style={%
    draw,%
    fill=green!20}}%""")
##     \pgfplotsset{every tick/.append style={thin,black!70}}
##     \pgfplotsset{every tick label/.append style={%""")
    #f.write("\nymajorticks=%s,"% ("true" if yticksLabels else "false"))
##     f.write(r"""
##     font=\tiny}}%""")
##     \pgfplotsset{%
##     every axis plot post/.append style={opacity=0.5}}%""")
##     if not yticksLabels:
##         f.write(r"""
##         \pgfplotsset{every y tick label/.append style={%
##         inner sep=0pt,outer sep=0pt,minimum size=0pt}}%""")
##     else:
##         f.write(r"""
##         \pgfplotsset{every y tick label/.append style={%
##         inner sep=0pt,outer sep=1pt,minimum size=0pt}}%""")        
    f.write(r"""
    \begin{axis}[%
    x tick label as interval,
    outer sep=0pt,
    tickwidth=2pt,
    title style={text height=1.5ex,text depth=.25ex},""")
    if title!=None:
        f.write("\ntitle=%s,%%" % title)
    if not yticksLabels:
        f.write("\nyticklabel=\yticklabeltextnone,%%")
    if width!=None:
        f.write("\nwidth=%s,%%" % width)
    if height!=None:
        f.write("\nheight=%s,%%" % height)
##     if ptype=="ybar":
##         f.write("\nybar,%%")
##     if barwidth!=None:
##         f.write("\nbar width=%s," % barwidth)
##     else:
##         if width != None:
##             f.write("\nbar width=%s/%d," % (width,2*nbins))
##     if ptype=="const":
##         f.write("\nconst plot,%%")
    if crange!=None:
        f.write("\nymin=%0.2f,%%" % (float(crange[0])))
        f.write("\nymax=%0.2f,%%" % (float(crange[1])))
    f.write("\nenlargelimits=0.03,")
    f.write("\ntick align=inside,")

    
    f.write(r"""
%    xtick=data,
    ybar interval=0.8, 
    legend style={nodes={anchor=west}},%
    name=mplot%,
    x tick label style={
    rotate=90,anchor=east,
    /pgf/number format/.cd,sci,sci subscript
    }
    ]%""")

    f.write(r"""
    \addplot[draw=black,fill=green!80!white] coordinates
    {%s};%%""" % (vals))
##            ("["+plotOptions+"]%%" if plotOptions != None else "",\


    if lnames != None:
        lstr = lnames[0]
        for name in lnames[1:]:
            lstr += "," + name
        f.write(r"\legend{%s}" % lstr)

    f.write(r"""
    \end{axis}%%""")
    f.write("\n\end{tikzpicture}%%")

