import pandas as pd
import pickle
from analysis_tools import summarize_bpm as sbpm
# from analysis_tools import summarize_wpm as swpm

def collectresults(resultfile,fdrcut,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,validationfile):

    pklin = open(resultfile,"rb")
    rfd = pickle.load(pklin)
    pklin.close()

    fdrBPM, fdrWPM, fdrPATH = rfd.fdrbpm2, rfd.fdrwpm2, rfd.fdrpath2

    vfdrBPM = fdrBPM[fdrBPM<=fdrcut].dropna()
    vfdrWPM = fdrWPM[fdrWPM<=fdrcut].dropna()
    vfdrPATH = fdrPATH[fdrPATH<=fdrcut].dropna()

    # print(vfdrBPM)
    # print(vfdrWPM)
    # print(vfdrPATH)

    pcutoff, netcut = 1, 1

    sbpm.summarize_bpm(ssmfile,bpmindfile,resultfile,validationfile,snppathwayfile,fdrcut,pcutoff,netcut)

    return
