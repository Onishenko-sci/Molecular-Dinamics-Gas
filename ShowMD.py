from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as figure
from tkinter import *
import csv
import sys
from scipy.optimize import curve_fit

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = './' + sys.argv[1]
    else:
        filename = '../render/md_render.txt'

i=0
j=0
head_readed = False
read_corelation_data = False
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
            frame_rate = int(row[6])
            cor_points = int(row[7])
            frames = int(steps/frame_rate)
            frame_x = np.zeros(int(N*(frames+1))).reshape(N,int(frames+1))
            frame_y = np.zeros(int(N*(frames+1))).reshape(N,int(frames+1))
            
            current_frame = np.zeros(int(frames+1))
            current_t = np.zeros(int(frames))
            potential_e = np.zeros(int(frames+1))
            cinetic_e = np.zeros(int(frames+1))
            sqared_displacment = np.zeros(int(frames))
            corl_x = np.zeros(int(row[7]))
            corl = np.zeros(int(row[7]))
            head_readed = True
            continue

        if ((i == N) and (not read_corelation_data)):
          current_frame[j] = float(row[1])
          cinetic_e[j] = float(row[2])
          potential_e[j] = float(row[3])
          sqared_displacment[j]= float(row[4])
          current_t[j] = float(row[5])
          i=0
          j=j+1
          continue

        if (float(row[0]) == 101):
            read_corelation_data = True
            i=0
            continue

        if (read_corelation_data):
            corl_x[i]=float(row[0])
            corl[i] = float(row[1])
            i=i+1
            continue

        frame_x[i,j] = float(row[0])
        frame_y[i,j] = float(row[1])
        i=i+1

mashtab = 800/bound_x
radius = radius/2

traked = N-1 #Number of tracked atom
tragectory_x = frame_x[traked,:]*mashtab
tragectory_y = frame_y[traked,:]*mashtab

def sqr(x,a):
    return a*x**2
X = np.arange(0,len(sqared_displacment))*frame_rate
X = X*delta_t*(10**9)
Y = sqared_displacment/(radius*2)**2
#[0:int(steps/(frame_rate*2))]

args, covar = curve_fit(sqr,X,Y)
Y_predicted = sqr(X,args[0])

plt.figure(figsize=(11,5))
plt.subplot(1,2,1)
plt.title('Correlation function')
plt.ylabel('Correlation function')
plt.xlabel('Distance between particles, Angstrom')
plt.plot(corl_x*(10**10),corl,'b.')
plt.plot(np.array([2*radius,2*radius])*(10**10),np.linspace(0,corl.max(),2),'g--')

plt.subplot(1,2,2)
plt.title('Mean squared displacement')
plt.ylabel('Mean squared displacement, meters')
plt.xlabel('time, ns')
plt.plot(X,Y , "b.")
plt.plot(X,Y_predicted , "r--")

plt.show()

total_e = cinetic_e+potential_e

plt.figure(figsize=(10,7), dpi=120)
plt.subplot(2,1,1)
plt.plot(np.arange(0,len(current_t)), current_t, 'b.')
plt.subplot(2,1,2)
plt.plot(np.arange(0,len(potential_e)), potential_e, 'b-')
plt.plot(np.arange(0,len(cinetic_e)), cinetic_e, 'r-')
plt.plot(np.arange(0,len(total_e)), total_e, 'g-')
plt.show()

scene_x = 5
scene_y = 5
density = (N)/(bound_x*bound_y)
density = density*(0.34*10**(-9))**2

root = Tk()
c = Canvas(root,width=int(bound_x*mashtab+500),heigh=int(bound_y*mashtab+30), bg="#303030", highlightthickness=0) 
c.pack()

for k in range(int(frames)):
    fr = k
    c.create_oval(tragectory_x[fr]-2+5,tragectory_y[fr]-2+5, tragectory_x[fr]+2+5,tragectory_y[fr]+2+5, width=0, fill="white")

c.create_line(bound_x*mashtab+5,5,bound_x*mashtab+5,bound_y*mashtab+5,width=5,fill="#BDD4F1")
c.create_line(bound_x*mashtab+5,bound_y*mashtab+5,5,bound_y*mashtab+5,width=5,fill="#BDD4F1")
c.create_line(5,5,5,bound_y*mashtab+5,width=5,fill="#BDD4F1")
c.create_line(5,5,bound_x*mashtab+5,5,width=5,fill="#BDD4F1")
c.create_text(bound_x*mashtab+10,bound_y*mashtab-15, text=("Particles: " + str(N)), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'))
c.create_text(bound_x*mashtab+10,bound_y*mashtab-45, text=("Density: " + str(density)[0:5]), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'))

for i in range(frames):
    c.delete("del")
    c.create_text(bound_x*mashtab+10,5,   text=("Time: " + str(current_frame[i]) + '*' +  str(delta_t)), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10,35,  text=("Temperature: " + str(current_t[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10,65,  text=("Kinetic Energy: " + str(cinetic_e[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10,95,  text=("Potential Energy: " + str(potential_e[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10,125,  text=("Total Energy: " + str(potential_e[i]+cinetic_e[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")
    c.create_text(bound_x*mashtab+10,155, text=("Mean squared displacement: " + str(sqared_displacment[i])), fill="#BDD4F1", anchor="nw", font=('Helvetica 15 bold'), tags="del")

    c.create_line(0,bound_y*mashtab+25,(int(i/(frames/1000))/1000)*(bound_x*mashtab+500),bound_y*mashtab+25,width=10, fill="#6BAABF", tags="del")



    for j in range(N):
        if j==traked:
            c.create_oval(scene_x + (frame_x[j,i]-radius)*mashtab,
                          scene_y + (frame_y[j,i]-radius)*mashtab,
                          scene_x + (frame_x[j,i]+radius)*mashtab,
                          scene_y + (frame_y[j,i]+radius)*mashtab, 
                          fill='#FB0000', tags="del", activefill="blue",outline="#FB0000")
        else:
            c.create_oval(scene_x + (frame_x[j,i]-radius)*mashtab,
                          scene_y + (frame_y[j,i]-radius)*mashtab,
                          scene_x + (frame_x[j,i]+radius)*mashtab,
                          scene_y + (frame_y[j,i]+radius)*mashtab,
                          fill="#6BAABF", tags="del", outline="#6BAABF")

    root.update()
    #sleep(0.001)
    
    if (i%(frames/100) == 0):
        print(int(i/(frames/100)), '%')
        if (int(i/(frames/100)) == 99):
            print("Done!")
            root.quit()


root.mainloop()


