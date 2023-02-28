'''When using this script in a published work, please cite the following paper.
Iliev, S.; Tsibranska, S.; Kichev, I.; Tcholakova, S.; Denkov, N.; Ivanova, A. Computational procedure for analysis of crystallites in polycrystalline solids of quasilinear molecules. Molecules 2023
'''
import numpy as np
import pandas as pd
from tkinter import filedialog
from tkinter import *
import sys

class StdRed(object):
    def __init__(self, textwid):
        self.text_space = textwid

    def write(self, texts):
        self.text_space.config(state=NORMAL)
        self.text_space.insert(END, texts)
        self.text_space.see(END)
        self.text_space.update_idletasks()

root = Tk()
root.configure(bg="black")
root.title("Slopes")

text = Text(root, width=150, font=("Helvetica", 10), highlightcolor="#33FFC9", state=NORMAL, bg="#BAB8B8", fg="black")
text.grid(row=1, columnspan=3)
sys.stdout = StdRed(text)

class main_program:
    @staticmethod
    def open_directory():
        file_path = filedialog.askopenfilenames(filetypes=[("Text files", ".dat"),
                                                           ("Text files", ".txt"),
                                                           ("Text files", ".xvg"),
                                                           ("All files", "*.*")])
        for filenames in file_path:
            df = pd.read_csv(filenames, sep='\t', header=None)
            df = df.drop(df.columns[0], axis=1)
            df = np.array(df)

            arglen = len(filenames)
            output1 = filenames[:arglen - 4]
            output1 = output1 + "-slopes.csv"

            fin_diff = []
            slopes = []    

            num_rows = len(df)
            num_cols = len(df[0])
            num_mols = int(num_cols/6)
            print("rows = " + str(num_rows) + " ; " "columns = " + str(num_cols) + " ; " "number of molecules = " + str(num_mols))

            for k in range(num_rows):
                for i in range(0, num_cols, 6):
                    fin_diff.append(df[k][i+3] - df[k][i])
                    fin_diff.append(df[k][i+4] - df[k][i+1])
                    fin_diff.append(df[k][i+5] - df[k][i+2])

            for j in range(0, len(fin_diff), 3):
                slopes.append(abs(fin_diff[j]/np.sqrt(fin_diff[j]**2 + fin_diff[j+1]**2 + fin_diff[j+2]**2)))
                
                slopes.append(abs(fin_diff[j+1]/np.sqrt(fin_diff[j]**2 + fin_diff[j+1]**2 + fin_diff[j+2]**2)))
                
                slopes.append(abs(fin_diff[j+2]/np.sqrt(fin_diff[j]**2 + fin_diff[j+1]**2 + fin_diff[j+2]**2)))

            angles = np.rad2deg(np.arccos(slopes))
            row_ang = len(angles)
            x = []
            y = []
            z = []
            av_angles = []
            
            for i in range(0, int(num_mols*3), 3):
                for j in range(i, row_ang, int(3*num_mols)):
                    x.append(angles[j])
                    y.append(angles[j+1])
                    z.append(angles[j+2])

            x = np.array(x)
            x_av = np.mean(x.reshape(-1, num_rows), axis=1)
            x_std = np.std(x.reshape(-1, num_rows), axis=1, ddof=1)

            y = np.array(y)
            y_av = np.mean(y.reshape(-1, num_rows), axis=1)
            y_std = np.std(y.reshape(-1, num_rows), axis=1, ddof=1)

            z = np.array(z)
            z_av = np.mean(z.reshape(-1, num_rows), axis=1)
            z_std = np.std(z.reshape(-1, num_rows), axis=1, ddof=1)

            output = pd.DataFrame({"X": x_av, "X-St.Dev": x_std, "Y": y_av, "Y-St.Dev": y_std, "Z": z_av, "Z-St.Dev": z_std})
            output.to_csv(output1, index=False)
        
            print(output1 + " has been created!")


# Define the browse button
browse_txt = StringVar()
browse_btn = Button(root, textvariable=browse_txt, command=lambda: main_program.open_directory(),
                    font=("Helvetica", 14), bg="#33FFC9", fg="black", width=15, height=1)
browse_txt.set("Browse")
browse_btn.grid(row=0, column=1)

root.mainloop()