
import os
import numpy as np

cwd = os.getcwd()
basepath = 'outputs/'
dirs = os.listdir(basepath)

for dir in dirs:
    if dir != '.DS_Store':
        subdirs = os.listdir(basepath + dir)
        for files in subdirs:
            if files == 'output_passive.txt':
                path = basepath + dir + '/' + filesfile
                checker = 1
                while checker == 1:
                    with open(path, "r+") as file:
                        datalen = np.loadtxt(path,usecols=0)
                        dat_si = len(datalen)
                        data = np.loadtxt(path,skiprows=dat_si-1)
                        if data.shape[0] == 5:
                            if len(data)>1:
                                if ( data[1] != 40.0):
                                    lines = file.readlines()
                                    lines = lines[:-1]  # Remove the last line
                                    file.seek(0)  # Rewind to the beginning
                                    file.truncate()
                                    file.writelines(lines)  # Overwrite with the modified lines     
                                else:
                                    checker = 0
                            else:
                                checker = 0
                                print('Failed to load data')
                        else:
                            lines = file.readlines()
                            lines = lines[:-1]  # Remove the last line
                            file.seek(0)  # Rewind to the beginning
                            file.truncate()
                            file.writelines(lines)  # Overwrite with the modified lines
                            #print('removed line')