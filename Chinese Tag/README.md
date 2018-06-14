# ChineseTag

Module.py can run in the IDE which has a GUI can make the training become easy.
This Chinese tagging system is based on the HMM model and is decoded using the Viterbi algorithm.

P.S. I have already packaged the module with pyinstaller as a module.exe, but the file is too large to be uploaded to github.
Anyone who wants to share this module can click below.

https://drive.google.com/drive/folders/1np74U3aCEjE6ioCOl7rO0sHTyXbuQ2Mb?usp=sharing

After you execute that exe, in the SelectSave you can enter the path to save the transmit and transfer matrices generated during the training (according to HMM model requirements) which are csv format and can be opened with excel. 

In the Slectfile you should enter your “path + file” to allow the system to select your data source. 

Finally, enter the sentence you want to test. Each word in this sentence should be separated by a space. 

Then click the start button to start training the model.


The data I use in training comes from the following address,(thanks for the author's sharing) 
And maybe you need to change some code if you want to use other format of data.

https://github.com/liwenzhu/corpusZh

If you are interested in that system or anything of my research, 
please contact with me :)

email: chris.yuan.ece@gmail.com
