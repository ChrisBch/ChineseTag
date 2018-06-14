# coding: utf-8
import numpy as np
import pandas as pd

from tkinter import *
import tkinter.messagebox as messagebox
import tkinter.filedialog as filedialog
import time

GetInfo = ""

class Application(Frame):
#This is a GUI class, we can build a Interface by this
#recommed that the save file is the directory which is the path that you save your every matrix
#e.g. D:\SavePOS\...
#the selectTrain is the file that you need to train the module
#e.g. D:\Train\train.txt
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        self.SaveLabel = Label(self, text = "SelectSave")
        self.SaveLabel.pack(fill=X)
        self.fileinput = Entry(self)
        self.fileinput.pack(fill=X)
        self.CoruLabel = Label(self, text = "SelectTrain")
        self.CoruLabel.pack(fill=X)
        self.Coruinput = Entry(self)
        self.Coruinput.pack(fill=X)
        self.txtinput = Label(self, text = "Test")
        self.txtinput.pack(fill=X)        
        self.txtinput = Entry(self)
        self.txtinput.pack(fill=X)        
        self.alertButton = Button(self, text='Start', command=self.start)
        self.alertButton.pack(fill=X)
        self.alertButton = Button(self, text='Help', command=self.help)
        self.alertButton.pack(fill=X)
        self.text = Text(self)
        self.text.pack()


    def start(self):
        global GetInfo
        GetInfo = self.txtinput.get()        
        result = Viterbi(self.fileinput.get(),self.Coruinput.get(),self.fileinput.get(),self.fileinput.get(),self.fileinput.get())
        self.text.insert(END,str(result))


    def help(self):
        HelpText = "n   普通名词 \nnt  时间名词 \nnd  方位名词\nnl  处所名词\nnh  人名\nnhf 姓\nnhs 名\nns  地名\nnn  族名\nni  机构名\nnz  其他专名\nv   动词" \
            "\nvu  能愿动词\na   形容词\nf   区别词\nm   数词　　\nq   量词\nd   副词\nr   代词\np   介词\nc   连词\nu   助词\ne   叹词\no   拟声词\ni   习用语\nj   缩略语" \
            "\nh   前接成分\nk   后接成分\ng   语素字\nx   非语素字\nw   标点符号\nws  非汉字字符串\nwu  其他未知的符号\nvd  趋向动词\nvl  联系动词"
        messagebox.showinfo("List of POS",message=HelpText)
        
        
def GetCoreDictionary(filepath="YourCorpus",csvpath="your save path"):
#This function can return the emission-frequency-matrix in dataframe format
#This function will be called by the GenerEmissionProb function, you may check it for more details
#csvpath is transfer-frequency-matrix storage path in your computer
#filepath is the path to your computer for storing text that can be trained
    EmissionMatrix = pd.DataFrame(columns=["n","nt","nd","nl","nh","nhf","nhs","ns","nn","ni","nz","v","vd","vl","vu",
                                     "a","f","m","q","d","r","p","c","u","e","o","i","j","h","k","g","x","w","ws","wu","mq",""])
    with open(filepath,"r",encoding="utf-8-sig") as reader:
        CoreDict={}
        for index,line in enumerate(reader):
            TempWordList = []
            TempPOSList = []
            flag = 0
            if index%1000==0:
                app.text.insert(END,"we have done "+str(index)+" lines to generate Emission Matrix"+"\n")
                app.update_idletasks()
            for index2,word in enumerate(line):
                if (index2 not in [0,1,2,3]) and word!="\n":
                    if word==" ":
                        TempWord = "".join(TempWordList)
                        TempPOS = "".join(TempPOSList[1:])
                        if TempWord not in CoreDict.keys():
                            CoreDict[TempWord] = {}
                        if TempPOS not in CoreDict[TempWord].keys():
                            CoreDict[TempWord][TempPOS] = 1
                        else:
                            CoreDict[TempWord][TempPOS] = CoreDict[TempWord][TempPOS] + 1
                        TempPOSList = []
                        TempWordList = []
                        flag = 0
                    elif word=="/":
                        TempPOSList.append(word)
                        flag = 1
                    else:
                        if flag==1:
                            TempPOSList.append(word)
                        else:
                            TempWordList.append(word)    
    for key in CoreDict.keys():
        EmissionMatrix.loc[key] = np.zeros(37)
        for sub_key in CoreDict[key].keys():
            EmissionMatrix[sub_key][key] = CoreDict[key][sub_key] 
    EmissionMatrix.drop(columns=[""],inplace=True)
    EmissionMatrix.loc["NTK"] = np.zeros(36) 
    EmissionMatrix.to_csv(csvpath+r"\EmissonFreq.csv")
    app.text.insert(END,"we have finished to generate Emisson matrix \n")
    app.update_idletasks()    
    return EmissionMatrix


def GetTranMatrix(filepath="YourCorpus",csvpath="your save path"):
#This function can return the transfer-frequency-matrix in dataframe format
#This function will be called by the GenerTranProb function, you may check it for more details
#csvpath is transfer-frequency-matrix storage path in your computer
#filepath is the path to your computer for storing text that can be trained
#It is recommended to add "r" before the path if you use windows. e.g., r"D:\example"
    TranMatrix = pd.DataFrame(data=np.zeros([37,37]),
                              index=["n","nt","nd","nl","nh","nhf","nhs","ns","nn","ni","nz","v","vd","vl","vu",
                                     "a","f","m","q","d","r","p","c","u","e","o","i","j","h","k","g","x","w","ws","wu","mq",""],
                              columns=["n","nt","nd","nl","nh","nhf","nhs","ns","nn","ni","nz","v","vd","vl","vu",
                                     "a","f","m","q","d","r","p","c","u","e","o","i","j","h","k","g","x","w","ws","wu","mq",""])
    with open(filepath,"r",encoding="utf-8-sig") as Reader:
        for index,line in enumerate(Reader):
            if index%1000==0:
                app.text.insert(END,"we have processed "+str(index)+" lines to generate Transfer Matrix"+"\n")
                app.update_idletasks()
            flag = 0
            TempPOSList = []
            OldPOS = None
            for index2,word in enumerate(line):
                if (index2 not in [0,1,2,3]) and word!="\n":
                    if word==" ":
                        if TempPOSList!=[]:
                            TempPOS = "".join(TempPOSList[1:])
                            TempPOSList = []
                            if OldPOS!=None:
                                TranMatrix[TempPOS][OldPOS] = TranMatrix[TempPOS][OldPOS] + 1
                            OldPOS = TempPOS
                            flag = 0
                    elif word=="/":
                        TempPOSList.append(word)
                        flag = 1
                    else:
                        if flag==1:
                            TempPOSList.append(word)
    TranMatrix.drop(index=[""],columns=[""],inplace=True)
    TranMatrix.to_csv(csvpath+r"\TranFreq.csv")
    return TranMatrix


def GenerTranProb(csvpath1="your save path",filepath="YourCorpus",csvpath2="your tran save path"):
#This function can return the transfer-probability-matrix in dataframe format
#csvpath1 is transfer-probability-matrix storage path in your computer
#csvpath2 is transfer-frequency-matrix storage path in your computer
#filepath is the path to your computer for storing text that can be trained
    TranMatrix = GetTranMatrix(filepath,csvpath2)
    TranMatrix['Col_sum'] = TranMatrix.apply(lambda x: x.sum(), axis=1)
    for r in range(36):
        for i in range(36):
            TranMatrix.iloc[r,i] = TranMatrix.iloc[r,i]/TranMatrix.iloc[r,36]
    TranMatrix.drop(columns=["Col_sum"],inplace=True)
    GenerTranProb = TranMatrix
    GenerTranProb.to_csv(csvpath1+r"\TranProb.csv")
    return GenerTranProb    


def GenerEmissionProb(csvpath1="your save path",filepath="YourCorpus",csvpath2="your emission save path"):
#This function can return the emissionr-probability-matrix in dataframe format
#csvpath1 is emission-probability-matrix storage path in your computer
#csvpath2 is emission-frequency-matrix // path in your computer
#filepath is the path to your computer for storing text that can be trained
    EmissionMatrix = GetCoreDictionary(filepath,csvpath2)
    EmissionMatrix['Col_sum'] = EmissionMatrix.apply(lambda x: x.sum(), axis=1)
    for r in range(len(EmissionMatrix)):
        if r%1000 == 0:
            app.text.insert(END,"we have calculated "+str(r)+" words' Emisson Probility Matrix"+"\n")
            app.update_idletasks()
        for i in range(36):
            EmissionMatrix.iloc[r,i] = EmissionMatrix.iloc[r,i]/EmissionMatrix.iloc[r,36]
    EmissionMatrix.drop(columns=["Col_sum"],inplace=True)
    GenerEmissionProb = EmissionMatrix
    GenerEmissionProb.to_csv(csvpath1+r"\EmissonProb.csv")
    app.text.insert(END,"we have finished to generate Emisson Probability Metrix\n")
    return GenerEmissionProb


def GetWord(sentence):
    WordList = []
    TempWordList1 = []
    for word in sentence:
        if word==" " :
            TempWord = "".join(TempWordList1)
            if TempWord in WordList:
                pass
            else: 
                WordList.append(TempWord)
            TempWordList1 = []
        else:
            TempWordList1.append(word)
    return WordList


def Viterbi(csvpathTran1="your tranmatrx save path",
            filepath="file",csvpathTran2="prob",
            csvpathEmission1="save",csvpathEmission2="prob"):
#This function can return the POS result in list format
#csvpathTran1(or 2) is transfer-frequency(probability)-matrix storage path in your computer
#csvpathEmisson1(or 2) is emission-frequency(probability)-matrix storage path in your computer
#filepath is the path to your computer for storing text that can be trained
    POPList = ["n","nt","nd","nl","nh","nhf","nhs","ns","nn","ni","nz","v","vd","vl","vu",
     "a","f","m","q","d","r","p","c","u","e","o","i","j","h","k","g","x","w","ws","wu","mq"] 
    DictOldPOS = {}
    EmissionProb = GenerEmissionProb(csvpathEmission1,filepath,csvpathEmission2)
    TranProb = GenerTranProb(csvpathTran1,filepath,csvpathTran2)
    WordList = GetWord(GetInfo)
    app.text.insert(END,"received the test sentence, we can start test now"+"\n")
    app.update_idletasks()
    NodeMatrix = [["" for i in range(36)] for i in range(2*len(WordList))]
    start = WordList.pop(0)
    j = 0
    try:
        for i in POPList:
            if EmissionProb[i][start] != 0:
                j = j + 1
                DictOldPOS[i] = EmissionProb[i][start]
                NodeMatrix[0][j-1] = i
        for i in range(1,len(WordList)+1):
            j = 0
            k = 0
            DictTempPOS = {}
            DictNewPOS = {}
            for m in POPList:
                if EmissionProb[m][WordList[i-1]] != 0:
                    j = j + 1
                    DictNewPOS[m] = EmissionProb[m][WordList[i-1]]
                    NodeMatrix[2*i][j-1] = m
            for key_1 in DictNewPOS.keys():
                TranFlag = 0
                TempBestList = []
                for key_2 in DictOldPOS.keys():
                    TempTranProb = TranProb[key_2][key_1]*DictOldPOS[key_2]*DictNewPOS[key_1]
                    if TempTranProb>=TranFlag:
                        TranFlag = TempTranProb
                        TempBestList.append(key_2)
                k = k + 1
                NodeMatrix[2*i-1][k-1] = TempBestList[-1]
                DictTempPOS[key_1] = TranFlag
            DictOldPOS = {}
            DictOldPOS = DictTempPOS
        KeyMax = max(DictOldPOS.keys())
        CountKey = 0
        for key in DictOldPOS.keys():
            CountKey = CountKey + 1
            if key==KeyMax:
                break
        BestList = [KeyMax,NodeMatrix[-3][CountKey-1]]
        count = CountKey
        for i in reversed(range(1,len(WordList))):
            TempKey = NodeMatrix[2*i-1][count-1]
            BestList.append(TempKey)
            count = 0
            for j in NodeMatrix[2*i-2]:
                count = count + 1
                if j == TempKey:
                    break
        BestList.reverse()
        return BestList
    except KeyError:
        return "there is no such word in dictionary"


if __name__ == "__main__":
    app = Application()  
    app.master.geometry('500x300+500+200')
    app.master.title('Chinese Tagging System')
    app.mainloop()