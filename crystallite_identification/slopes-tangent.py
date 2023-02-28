'''If you are using this script, please cite the following paper.
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
        # self.text_space.config(state=DISABLED)

#    def flush(self):
#        pass

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
                if fin_diff[j] < 1e-6 and fin_diff[j] >= 0: 
                    fin_diff[j] = fin_diff[j] + 1e-6
                if fin_diff[j] > -1e-6 and fin_diff[j] < 0: 
                    fin_diff[j] = fin_diff[j] - 1e-6
                slopes.append(fin_diff[j+1]/fin_diff[j])
                
                if fin_diff[j+1] < 1e-6 and fin_diff[j+1] >= 0:
                     fin_diff[j+1] = fin_diff[j+1] + 1e-6
                if fin_diff[j+1] > -1e-6 and fin_diff[j+1] < 0: 
                    fin_diff[j+1] = fin_diff[j+1] - 1e-6
                slopes.append(fin_diff[j+2]/fin_diff[j+1])
                
                if fin_diff[j+2] < 1e-6 and fin_diff[j+2] >= 0: 
                    fin_diff[j+2] = fin_diff[j+2] + 1e-6
                if fin_diff[j+2] > -1e-6 and fin_diff[j+2] < 0: 
                    fin_diff[j+2] = fin_diff[j+2] - 1e-6
                slopes.append(fin_diff[j]/fin_diff[j+2])

            angles = np.arctan(slopes)*57.2958
            row_ang = len(angles)
            yx = []
            zy = []
            xz = []
            av_angles = []
            
            for i in range(0, int(num_mols*3), 3):
                for j in range(i, row_ang, int(3*num_mols)):
                    yx.append(angles[j])
                    zy.append(angles[j+1])
                    xz.append(angles[j+2])

            yx = np.array(yx)
            yx_av = np.mean(yx.reshape(-1, num_rows), axis=1)
            yx_std = np.std(yx.reshape(-1, num_rows), axis=1, ddof=1)

            zy = np.array(zy)
            zy_av = np.mean(zy.reshape(-1, num_rows), axis=1)
            zy_std = np.std(zy.reshape(-1, num_rows), axis=1, ddof=1)

            xz = np.array(xz)
            xz_av = np.mean(xz.reshape(-1, num_rows), axis=1)
            xz_std = np.std(xz.reshape(-1, num_rows), axis=1, ddof=1)

            output = pd.DataFrame({"YX": yx_av, "YX-St.Dev": yx_std, "ZY": zy_av, "ZY-St.Dev": zy_std, "XZ": xz_av, "XZ-St.Dev": xz_std})
            output.to_csv(output1, index=False)
        
            print(output1 + " has been created!")


# Define the browse button
browse_txt = StringVar()
browse_btn = Button(root, textvariable=browse_txt, command=lambda: main_program.open_directory(),
                    font=("Helvetica", 14), bg="#33FFC9", fg="black", width=15, height=1)
browse_txt.set("Browse")
browse_btn.grid(row=0, column=1)

root.mainloop()