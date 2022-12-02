
from telnetlib import theNULL
from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as figure
from tkinter import *
import csv
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = './' + sys.argv[1]
    else:
        filename = '../render/md_render.txt'

N = 100
radius = 0.0
T = 0.0
delta_t = 1.0
frames = T/delta_t
bound_x = 0
bound_y = 0
cor_points = 0





i=0
j=0
head_readed = False
read_corl = False
with open(filename, newline='\n' ) as f:
    reader = csv.reader(f,delimiter = ';')
    for row in reader:
        if (head_readed == False):
            bound_x = float(row[0])
            bound_y = float(row[1])
            N = int(row[2])
            radius = float(row[3])
            steps = int(row[4])
            delta_t = float(row[5])
            save_every_frame = int(row[6])
            cor_points = int(row[7])
            frames = int(steps/save_every_frame)
            frame_x = np.zeros(int(N*(frames+1))).reshape(N,int(frames+1))
            frame_y = np.zeros(int(N*(frames+1))).reshape(N,int(frames+1))
            
            current_frame = np.zeros(int(frames+1))
            current_t = np.zeros(int(frames+1))
            potential_e = np.zeros(int(frames+1))
            cinetic_e = np.zeros(int(frames+1))
            sqared_displacment = np.zeros(int(frames))
            corl_x = np.zeros(int(row[7]))
            corl = np.zeros(int(row[7]))
            head_readed = True
            continue

        if ((i == N) and (not read_corl)):
          current_frame[j] = float(row[1])
          cinetic_e[j] = float(row[2])
          potential_e[j] = float(row[3])
          sqared_displacment[j]= float(row[4])
          current_t[j] = float(row[5])
          i=0
          j=j+1
          continue

        if (float(row[0]) >= 1487.0):
            read_corl = True
            i=0
            continue

        if (read_corl):
            corl_x[i]=float(row[0])
            corl[i] = float(row[1])
            i=i+1
            continue

        frame_x[i,j] = float(row[0])
        frame_y[i,j] = float(row[1])
        i=i+1

full_e = cinetic_e+potential_e

mashtab = 800/bound_x

traked = 1
tragectory_x = frame_x[traked,:]*mashtab
tragectory_y = frame_y[traked,:]*mashtab

plt.figure(figsize=(13,5), dpi=120)
plt.subplot(1,2,1)
plt.title('Correlation function')
plt.plot(corl_x,corl,'b.')
#plt.plot(np.ones(3)*radius,np.linspace(0,1e-1,3),'r--')
plt.plot(np.ones(3)*bound_x/2,np.linspace(0,np.amax(corl),3),'g--')

plt.subplot(1,2,2)
plt.title('Squared displacment')
plt.plot(np.arange(0,sqared_displacment.__len__())*delta_t,sqared_displacment[:] , "b-")

plt.show()

root = Tk()
c = Canvas(root,width=int(bound_x*mashtab+400),heigh=int(bound_y*mashtab+20))
c.pack()

for k in range(int(frames)):
    fr = k
    c.create_oval(tragectory_x[fr]-2,tragectory_y[fr]-2, tragectory_x[fr]+2,tragectory_y[fr]+2, width=0, fill="red")

c.create_line(bound_x*mashtab,0,bound_x*mashtab,bound_y*mashtab+20,width=5)
c.create_line(bound_x*mashtab,bound_y*mashtab,0,bound_y*mashtab,width=5)
c.create_text(bound_x*mashtab+5,155, text=("Particles: " + str(N)), fill="black", anchor="nw", font=('Helvetica 15 bold'))

for i in range(frames):
    c.delete("del")
    c.create_text(bound_x*mashtab+5,5, text=("Time: " + str(current_frame[i]) + '*' +  str(delta_t)),
     fill="black", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+5,35, text=("Temperature: " + str(current_t[i])), fill="black", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+5,65, text=("Cinetic Energy: " + str(cinetic_e[i])), fill="black", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+5,95, text=("Potencial Energy: " + str(potential_e[i])), fill="black", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+5,125, text=("Suadred displacment: " + str(sqared_displacment[i])), fill="black", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_line(0,bound_y*mashtab+15,(int(i/(frames/1000))/1000)*bound_x*mashtab,bound_y*mashtab+15,width=20, fill="green", tags="del")



    for j in range(N):
        if j==traked:
            c.create_oval((frame_x[j,i]-radius)*mashtab,(frame_y[j,i]-radius)*mashtab,(frame_x[j,i]+radius)*mashtab,(frame_y[j,i]+radius)*mashtab, fill='white', tags="del")
        else:
            c.create_oval((frame_x[j,i]-radius)*mashtab,(frame_y[j,i]-radius)*mashtab,(frame_x[j,i]+radius)*mashtab,(frame_y[j,i]+radius)*mashtab, fill='black', tags="del")
       # c.create_oval((frame_x[j,i]-4*radius)*mashtab,(frame_y[j,i]-4*radius)*mashtab,(frame_x[j,i]+4*radius)*mashtab,(frame_y[j,i]+4*radius)*mashtab)

    root.update()
    #sleep(delta_t*mashtab)
    sleep(0.005)
    
    if (i%(frames/100) == 0):
        print(int(i/(frames/100)), '%')
        if (int(i/(frames/100)) == 99):
            print("Done!")
            root.quit()
           # root.destroy()

root.mainloop()


