# -*- coding: utf-8 -*-
filename1 = 'C:/Users/kolesar/Documents/cv/51_Schönaich bei Böblingen/04_depth_profiel_analysis/RBS measurements for depth profile/141211/1_S.TXT'
def get_table2(filename1):
    dataFile = open(filename1, 'r')
    channel1 = []
    counts1 = []
    energy1 = []
    header1 = dataFile.readline()
    header2 = dataFile.readline()
    header3 = dataFile.readline()
    header4 = dataFile.readline()
    header5 = dataFile.readline()
    header6 = dataFile.readline()
    header7 = dataFile.readline()
    header8 = dataFile.readline()
    header9 = dataFile.readline()
    header10 = dataFile.readline()
    header11 = dataFile.readline()
    header12 = dataFile.readline()
    header13 = dataFile.readline()
    header14 = dataFile.readline()
    header15 = dataFile.readline()
    header16 = dataFile.readline()
    header17 = dataFile.readline()
    header18 = dataFile.readline()
    header19 = dataFile.readline()
    header20 = dataFile.readline()
    header21 = dataFile.readline()
    header22 = dataFile.readline()
    header23 = dataFile.readline()
    header24 = dataFile.readline()
    header25 = dataFile.readline()
    header26 = dataFile.readline()
    header27 = dataFile.readline()
    header28 = dataFile.readline()
    header29 = dataFile.readline()
    header30 = dataFile.readline()
    header31 = dataFile.readline()
    header32 = dataFile.readline()
    header33 = dataFile.readline()
    header34 = dataFile.readline()
    header35 = dataFile.readline()
    header36 = dataFile.readline()
    header37 = dataFile.readline()
    header38 = dataFile.readline()
    for line in dataFile:
        x1, y1,y2,y3, y4, y5, y6, y7  = line.split()
        channel1.append(float(x1))
        counts1.append(float(y1))
    dataFile.close()
    return (channel1, counts1)