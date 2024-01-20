#!/usr/bin/env python3


# https://matplotlib.org/gallery/user_interfaces/embedding_in_tk_sgskip.html
# https://matplotlib.org/gallery/animation/simple_anim.html
# https://qiita.com/nv-h/items/92feeb34338c09c6d2a2

import serial

import tkinter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.animation as animation

import numpy as np


def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    readSerial.close()


def init():  # only required for blitting to give a clean slate.
    line.set_ydata(np.sin(x))
    return line,

data_dict={}
ys_cur=None
def animate(i):
    global ys_cur
    lines = readSerial.read_all() # 1 line (upto '\n')
    if len(lines)>0:
        for l in lines.decode('utf-8').split("\n")[:-1]:
            arr=l.split(",")
            print(l)
            if len(arr)>=2:
                x=int(arr[0].strip())
                y=float(arr[1].strip())
                data_dict[x]=y
    items=list(sorted(data_dict.items()))
    #print(items)
    if len(items)>0:
        xs=[x/180*np.pi for x,y in items]+[items[0][0]]
        ys=[y for x,y in items]+[items[0][1]]
        if ys_cur is None or len(ys_cur)<len(ys):
            ys_cur=ys
        else:
            ys_cur=[ys[i]+0.5*(ys_cur[i]-ys[i]) for i in range(len(ys))]
        line.set_xdata(xs)  # update the data.
        line.set_ydata(ys_cur)  # update the data.
    return line,


readSerial = serial.Serial('/dev/ttyUSB0',115200, timeout=3)
root = tkinter.Tk()
root.wm_title("Embedding in Tk anim")

fig = Figure()
# FuncAnimationより前に呼ぶ必要がある
canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.

x = np.arange(0, 3, 0.01)  # x軸(固定の値)
l = np.arange(0, 8, 0.01)  # 表示期間(FuncAnimationで指定する関数の引数になる)
plt = fig.add_subplot(111, polar = True)
plt.set_ylim([0.0, 5.0])
line, = plt.plot(x, np.sin(x))

ani = animation.FuncAnimation(fig, animate, l,
    init_func=init, interval=10, blit=True,
    )

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

button = tkinter.Button(master=root, text="Quit", command=_quit)
button.pack()

tkinter.mainloop()

