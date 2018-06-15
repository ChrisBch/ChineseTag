#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define BUCKETCOUNT 16 
/**************************************************
Author: ChenHan Yuan

Note:
First, the program defines a hash table, 
Second,the GetXXTable function is defined to obtain the mnemonic symbol table and the register table. 
Then, the program defines two modes
(1 corresponds to a program with block (2.11) and 2 corresponds to no block (2.9)). 
After that,this is two rounds of traversal of assembly language.
And generate the corresponding Table and destination code
***************************************************/ 
struct hashEntry
{
    const char* key;
    char* value;
    struct hashEntry* next;
};
typedef struct hashEntry entry;

struct hashTable
{
    entry bucket[BUCKETCOUNT]; 
};
typedef struct hashTable table;

void initHashTable(table* t)
{
    int i;
    if (t == NULL)return;

    for (i = 0; i < BUCKETCOUNT; ++i) {
        t->bucket[i].key = NULL;
        t->bucket[i].value = NULL;
        t->bucket[i].next = NULL;
    }
}

int keyToIndex(const char* key)
{
    int index , len , i;
    if (key == NULL)return -1;

    len = strlen(key);
    index = (int)key[0];
    for (i = 1; i<len; ++i) {
        index *= 1103515245 + (int)key[i];
    }
    index >>= 27;
    index &= (BUCKETCOUNT - 1);
    return index;
}

char* strDup(const char* str)
{
    int len;
    char* ret;
    if (str == NULL)return NULL;

    len = strlen(str);
    ret = (char*)malloc(len + 1);
    if (ret != NULL) {
        memcpy(ret , str , len);
        ret[len] = '\0';
    }
    return ret;
}

int insertEntry(table* t , const char* key , const char* value)
{
    int index , vlen1 , vlen2;
    entry* e , *ep;

    if (t == NULL || key == NULL || value == NULL) {
        return -1;
    }

    index = keyToIndex(key);
    if (t->bucket[index].key == NULL) {
        t->bucket[index].key = strDup(key);
        t->bucket[index].value = strDup(value);
    }
    else {
        e = ep = &(t->bucket[index]);
        while (e != NULL) { 
            if (strcmp(e->key , key) == 0) {
                vlen1 = strlen(value);
                vlen2 = strlen(e->value);
                if (vlen1 > vlen2) {
                    free(e->value);
                    e->value = (char*)malloc(vlen1 + 1);
                }
                memcpy(e->value , value , vlen1 + 1);
                return index;  
            }
            ep = e;
            e = e->next;
        } 
        e = (entry*)malloc(sizeof (entry));
        e->key = strDup(key);
        e->value = strDup(value);
        e->next = NULL;
        ep->next = e;
    }
    return index;
}

const char* findValueByKey(const table* t , const char* key)
{
    int index;
    const entry* e;
    if (t == NULL || key == NULL) {
        return NULL;
    }
    index = keyToIndex(key);
    e = &(t->bucket[index]);
    if (e->key == NULL) return NULL;
    while (e != NULL) {
        if (0 == strcmp(key , e->key)) {
            return e->value;   
        }
        e = e->next;
    }
    return NULL;
}

//You can change the following file path
void SaveTable(table* t,int mode)
{
	FILE *fp = NULL;
	if(mode==0)
	{
		fp = fopen("D:\\OpTable.txt","w+");
		fputs("Index MemCode Format|OpCode \n", fp);
	}
	if (mode==1){
		fclose(fp);
		fp = fopen("D:\\SymTable.txt","w+");
		fputs("Index SymName Address \n", fp);
	}
	if (mode==2){
		fclose(fp);
		fp = fopen("D:\\LitTable.txt","w+");
		fputs("Index LitName Address \n", fp);
	}
	if (mode==3){
		fclose(fp);
		fp = fopen("D:\\RegTable.txt","w+");
		fputs("Index RegName RegNum \n", fp);//in %x format
	}
    int i;
    int k=0;
    entry* e;
    if (t == NULL)return;
    for (i = 0; i<BUCKETCOUNT; ++i) {
        e = &(t->bucket[i]);
        while (e->key != NULL) {
        	k = k + 1;
            fprintf(fp,"%d %s %s \n" ,k,e->key,e->value);
            if (e->next == NULL)break;
            e = e->next;
        }
    }
    fclose(fp);
}

table GetOpTable()
{
	table Optable;
	initHashTable(&Optable);
	insertEntry(&Optable,"ADD","3/4|18");
	insertEntry(&Optable,"ADDF","3/4|58");
	insertEntry(&Optable,"ADDR","2|90");
	insertEntry(&Optable,"AND","3/4|40");
	insertEntry(&Optable,"CLEAR","2|B4");
	insertEntry(&Optable,"COMP","3/4|28");
	insertEntry(&Optable,"COMPF","3/4|88");
	insertEntry(&Optable,"COMPR","2|A0");
	insertEntry(&Optable,"DIV","1|24");
	insertEntry(&Optable,"DIVF","1|64");
	insertEntry(&Optable,"DIVR","1|9C");
	insertEntry(&Optable,"FIX","3/4|C4");
	insertEntry(&Optable,"FLOAT","3/4|C0");
	insertEntry(&Optable,"HIO","3/4|F4");
	insertEntry(&Optable,"J","3/4|3C");
	insertEntry(&Optable,"JEQ","3/4|30");
	insertEntry(&Optable,"JGT","3/4|34");
	insertEntry(&Optable,"JLT","3/4|38");
	insertEntry(&Optable,"JSUB","3/4|48");
	insertEntry(&Optable,"LDA","3/4|00");
	insertEntry(&Optable,"LDB","3/4|68");
	insertEntry(&Optable,"LDCH","3/4|50");
	insertEntry(&Optable,"LDF","3/4|70");
	insertEntry(&Optable,"LDL","3/4|08");
	insertEntry(&Optable,"LDS","3/4|6C");
	insertEntry(&Optable,"LDT","3/4|74");
	insertEntry(&Optable,"LDX","3/4|04");
	insertEntry(&Optable,"LPS","3/4|E0");
	insertEntry(&Optable,"UML","3/4|20");
	insertEntry(&Optable,"MULF","3/4|60");
	insertEntry(&Optable,"MULR","2|98");
	insertEntry(&Optable,"NORM","1|C8");
	insertEntry(&Optable,"OR","3/4|44");
	insertEntry(&Optable,"RD","3/4|D8");
	insertEntry(&Optable,"RMO","2|AC");
	insertEntry(&Optable,"RSUB","3/4|4C");
	insertEntry(&Optable,"SHIFTL","2|A4");
	insertEntry(&Optable,"SHIFTR","2|A8");
	insertEntry(&Optable,"SIO","1|F0");
	insertEntry(&Optable,"SSK","3/4|EC");
	insertEntry(&Optable,"STA","3/4|0C");
	insertEntry(&Optable,"STB","3/4|78");
	insertEntry(&Optable,"STCH","3/4|54");
	insertEntry(&Optable,"STF","3/4|80");
	insertEntry(&Optable,"STI","3/4|D4");
	insertEntry(&Optable,"STL","3/4|14");
	insertEntry(&Optable,"STS","3/4|7C");
	insertEntry(&Optable,"STSW","3/4|E8");
	insertEntry(&Optable,"STT","3/4|84");
	insertEntry(&Optable,"STX","3/4|10");
	insertEntry(&Optable,"SUB","3/4|1C");
	insertEntry(&Optable,"SUBF","3/4|5C");
	insertEntry(&Optable,"SUBR","2|94");
	insertEntry(&Optable,"SVC","2|B0");
	insertEntry(&Optable,"TD","3/4|E0");
	insertEntry(&Optable,"TIO","1|F8");
	insertEntry(&Optable,"TIX","3/4|2C");
	insertEntry(&Optable,"TIXR","2|B8");
	insertEntry(&Optable,"WD","3/4|DC");
	SaveTable(&Optable,0);
	return Optable;
}

table GetRegTable()
{
	table Regtable;
	initHashTable(&Regtable);
	insertEntry(&Regtable,"A","00");
	insertEntry(&Regtable,"X","10");
	insertEntry(&Regtable,"L","20");
	insertEntry(&Regtable,"PC","80");
	insertEntry(&Regtable,"SW","90");
	insertEntry(&Regtable,"S","40");
	insertEntry(&Regtable,"T","50");
	insertEntry(&Regtable,"F","60");
	SaveTable(&Regtable,3);
	return Regtable;
}

//You can change the following default start address of program block 
table GetBloTable()
{
	table Blotable;
	initHashTable(&Blotable);
	insertEntry(&Blotable,"CDATA","66");
	insertEntry(&Blotable,"CBLKS","71");
	insertEntry(&Blotable,"defau","00");
	return Blotable;
 } 
 
int main()
{
	table OP = GetOpTable();
	table REG= GetRegTable();
	table BloTable = GetBloTable();
	table SymTable;
	initHashTable(&SymTable);
	table LitTable;
	initHashTable(&LitTable);
	table BloNumTable;
	initHashTable(&BloNumTable);
	table BloTableNum;
	initHashTable(&BloTableNum);
	char mode; 
	char buffer[50];
	char TempLoc[6];
	char TempSym[6];
	char TempOp[6];
	char TempOpSym[8];
	char TempOpSym2[7];
	char TempExpSym[2];
	char TempBase[8];
	char Temp[]="        ";
	char Temp2[]="        ";
	int TempInt;
	int TempInt2;
	int TempInt3;
	int LtorgFlag;
	int LocCounter;
	int BaseAddr;
	int BloCount[3]={0x00,0x00,0x00};
	int BloCounter;
	char TempBloCounter[]="\0\0";
	char TempUse[]="\0\0\0\0\0";
	char EquLit[7]="\0\0";
	char EquLit2[6]="\0\0";
	int FlagA=0;
	int FlagB=0;
	printf("which mode would you like to program:\n1:with blocks\n0:no block\n");
	mode = getchar();
	FILE *fpR;
	//You can change the path of source program file  
	if(mode=='0')
	fpR = fopen("D:\\srcpro.txt","r");
	else
	fpR = fopen("D:\\srcpro2.txt","r");
	while (!feof(fpR)){
		fgets(buffer, sizeof(buffer), fpR);
		strncpy(TempOp,buffer+8,6);
		strncpy(TempSym,buffer,6);
		strncpy(TempOpSym,buffer+16,8);
		strncpy(TempOpSym2,buffer+25,7);
		strncpy(TempExpSym,buffer+7,1);
		for(int j=0;TempSym[j]!='\0';j++){
		if((TempSym[j]==' ')||(TempSym[j]=='l'))
			TempSym[j]='\0';
		}
		for(int j=0;TempOp[j]!='\0';j++){
		if(TempOp[j]==' ')
			TempOp[j]='\0';
		}
		for(int j=0;TempExpSym[j]!='\0';j++){
		if(TempExpSym[j]==' ')
			TempExpSym[j]='\0';
		}
		for(int j=0;TempOpSym[j]!='\0';j++){
		if(TempOpSym[j]==' ')
			TempOpSym[j]='\0';
		}
		for(int j=0;TempOpSym2[j]!='\0';j++){
		if(TempOpSym2[j]==' ')
			TempOpSym2[j]='\0';
		}
		if (strcmp(TempSym,"COPY")==0){
			LocCounter = 0x00;
			itoa(BloCounter,TempBloCounter,10);
			insertEntry(&BloNumTable,"defau",TempBloCounter);
			insertEntry(&BloTableNum,TempBloCounter,"defau");
			strncpy(TempUse,"defau",5);
		}
		if (strcmp(TempOp,"USE")==0){
			if(TempOpSym[0]=='\0'){
			 BloCount[atoi(findValueByKey(&BloNumTable,TempUse))]= LocCounter;
			 LocCounter = BloCount[atoi(findValueByKey(&BloNumTable,"defau"))];
			 strncpy(TempUse,"defau",5);
			}
			else if((FlagA==1)&&(FlagB==1)){
			  BloCount[atoi(findValueByKey(&BloNumTable,TempUse))]= LocCounter;
			  LocCounter= BloCount[atoi(findValueByKey(&BloNumTable,TempOpSym))];
			  strncpy(TempUse,TempOpSym,5);	
			}
			else {	
			  BloCount[atoi(findValueByKey(&BloNumTable,TempUse))]= LocCounter;
			  LocCounter = 0x00;
			  BloCounter = BloCounter+1;
			  itoa(BloCounter,TempBloCounter,10);
			  insertEntry(&BloNumTable,TempOpSym,TempBloCounter);
			  insertEntry(&BloTableNum,TempBloCounter,TempOpSym);
			  strncpy(TempUse,TempOpSym,5);
			  if(strcmp(TempOpSym,"CDATA")==0)
			  FlagA=1;
			  if(strcmp(TempOpSym,"CBLKS")==0)
			  FlagB=1;
			}
		}
		else if(findValueByKey(&OP,TempOp)!=NULL){
			strncpy(Temp,findValueByKey(&OP,TempOp),1);
			LtorgFlag = 0;
			if(buffer[15]=='='){
				if(EquLit[0]=='\0')
				strncpy(EquLit,TempOpSym,6);
				else
				strncpy(EquLit2,TempOpSym,5);
			}
			if((TempSym[0]!=' ')&&(LtorgFlag==0)){
        	    /*a bug about why add a space*/ 
        	    if (mode=='0'){
        	    	itoa(LocCounter,TempLoc,16);
        	    	insertEntry(&SymTable,TempSym,TempLoc);
				}
        	    else{
        	    	itoa(LocCounter,TempLoc,16);
        	    	strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	    	insertEntry(&SymTable,TempSym,TempLoc);
				}
		    }
			if(TempExpSym[0]=='+')
			LocCounter = LocCounter + 0x04;
			else if(Temp[0]=='2')
			LocCounter = LocCounter + 0x02;
			else
			LocCounter = LocCounter + 0x03;
		}
		else if(strcmp(TempOp,"START")==0)
		LocCounter = LocCounter;
		else if(strcmp(TempOp,"END")==0)
		printf("we have finished the first circle\n");
		else if(strcmp(TempOp,"BASE")==0)
		{
			strncpy(TempBase,TempOpSym,8);
			for(int j=0;TempBase[j]!='\0';j++){
		    if(TempBase[j]==' ')
			  TempBase[j]='\0';
		    }
		}
		else if(strcmp(TempOp,"LTORG")==0){
			LtorgFlag = 1;
			if(EquLit[0]!='\0'){
			  if(mode=='0'){
			    itoa(LocCounter,TempLoc,16);
			    insertEntry(&LitTable,EquLit,TempLoc);	
		      }
		      else{
		        itoa(LocCounter,TempLoc,16);
		        strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	    insertEntry(&LitTable,EquLit,TempLoc);
			  }
			LocCounter = LocCounter + strlen(EquLit)-3;
		    }
		    if(EquLit2[0]!='\0'){
			  if(mode=='0'){
			    itoa(LocCounter,TempLoc,16);
			    insertEntry(&LitTable,EquLit2,TempLoc);	
		      }
		      else{
		        itoa(LocCounter,TempLoc,16);
		        strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	    insertEntry(&LitTable,EquLit2,TempLoc);
			  }
			LocCounter = LocCounter + strlen(EquLit2)-3;
		    }
		}
		else if(strcmp(TempOp,"WORD")==0)
		{
			if(strcmp(TempBase,TempSym)==0)
			BaseAddr =LocCounter;
			if((LtorgFlag)||(mode=='1'))
			{
			 if(mode=='0'){
			   itoa(LocCounter,TempLoc,16);
			   insertEntry(&LitTable,TempSym,TempLoc);	
		     }
		     else{
		       itoa(LocCounter,TempLoc,16);
		       strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	   insertEntry(&LitTable,TempSym,TempLoc);
			 }
			}
			LocCounter = LocCounter + 0x03;	
		}
		else if(strcmp(TempOp,"RESW")==0)
		{
			if(strcmp(TempBase,TempSym)==0)
			BaseAddr =LocCounter;
			if((LtorgFlag)||(mode=='1'))
			{
			  if(mode=='0'){
			   itoa(LocCounter,TempLoc,16);
			   insertEntry(&LitTable,TempSym,TempLoc);	
		      }
		      else{
		       itoa(LocCounter,TempLoc,16);
		       strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	   insertEntry(&LitTable,TempSym,TempLoc);
			  }	
			}
			LocCounter = LocCounter + 3*atoi(TempOpSym);
		}
	    else if(strcmp(TempOp,"RESB")==0)
		{
			if(strcmp(TempBase,TempSym)==0)
			BaseAddr =LocCounter;
			if((LtorgFlag)||(mode=='1'))
			{
			  if(mode=='0'){
			   itoa(LocCounter,TempLoc,16);
			   insertEntry(&LitTable,TempSym,TempLoc);	
		      }
		      else{
		       itoa(LocCounter,TempLoc,16);
		       strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	   insertEntry(&LitTable,TempSym,TempLoc);
			  }	
			}
			LocCounter = LocCounter + 1*atoi(TempOpSym);
		    
		}
		else if(strcmp(TempOp,"BYTE")==0){
			if(strcmp(TempBase,TempSym)==0)
			BaseAddr =LocCounter;
			if(mode=='0'){
			 if(LtorgFlag){
			  itoa(LocCounter,TempLoc,16);
			  insertEntry(&LitTable,TempSym,TempLoc);	
			 }
			 else{
			  itoa(LocCounter,TempLoc,16);
			  insertEntry(&SymTable,TempSym,TempLoc);
			 }
		    }
		    else{
		    	if (TempSym[0]!='\0'){
		    	  itoa(LocCounter,TempLoc,16);
		    	  strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
			      insertEntry(&SymTable,TempSym,TempLoc);
				}
				else{
				  itoa(LocCounter,TempLoc,16);
				  strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
			      insertEntry(&LitTable,TempSym,TempLoc);	
				}
			}
			LocCounter = LocCounter + (strlen(TempOpSym)-3)/2;
		}
        else if(strcmp(TempOp,"EQU")==0)
        {
        	if(strcmp(TempBase,TempSym)==0)
			BaseAddr =LocCounter;
        	if((LtorgFlag)||(mode=='1'))
        	{
        	 strncpy(Temp,TempOpSym,1);
        	 if(Temp[0]=='*'){
        	 	if(mode=='0'){
			     itoa(LocCounter,TempLoc,16);
			     insertEntry(&LitTable,TempSym,TempLoc);	
		        }
		        else{
		         itoa(LocCounter,TempLoc,16);
		         strcat(TempLoc,findValueByKey(&BloNumTable,TempUse));
        	     insertEntry(&LitTable,TempSym,TempLoc);
			    }	
			 }
			 else{
			 	strncpy(Temp,buffer+24,1);
			 	if(Temp[0]=='-')
			 	{
			 	 char AddOpSym[8];
			 	 strncpy(AddOpSym,TempOpSym2,6);
			 	  if(mode=='0'){
			 	 	itoa(atoi(findValueByKey(&LitTable,TempOpSym))-atoi(findValueByKey(&LitTable,AddOpSym)),AddOpSym,10);
			 	    insertEntry(&LitTable,TempSym,AddOpSym);
			      }
			      else
			      {
			     	char a[5];
			     	char b[3];
			     	strncpy(a,findValueByKey(&LitTable,TempOpSym),4);
			     	strncpy(b,findValueByKey(&LitTable,AddOpSym),1);
			     	itoa(atoi(a)-atoi(b),AddOpSym,10);
					strcat(AddOpSym,findValueByKey(&BloNumTable,TempUse));
			 	    insertEntry(&LitTable,TempSym,AddOpSym);
				  }
				 }
			    }	
			}
		}
		else
		  printf("error:No such assembler synax\n");
	}
	SaveTable(&SymTable,1);
	SaveTable(&LitTable,2);
	int line=0;
	int ObjCode;
	char Change[8];
	char Chan[6];
	char Dex1[11]="0x";
	char BloFlag[]="\0";
	fclose(fpR);
	FILE *fpR2;
	FlagA=0;
	FlagB=0;
	int Modi=0;
	int ModiCounter[4]={0x00,0x00,0x00};
	char start[3];
	int ProLenth=LocCounter;
	//You can change the path of storaging Table file  
	if(mode=='0')
	fpR2 = fopen("D:\\srcpro.txt","r");
	else
	fpR2 = fopen("D:\\srcpro2.txt","r");
	FILE *fpWOR;
	fpWOR = fopen("D:\\original.txt","w");
	FILE *fpWRC;
	fpWRC = fopen("D:\\record.txt","w+");
	while (!feof(fpR2)){
		line = line + 1;
		fgets(buffer, sizeof(buffer), fpR2);
		strncpy(TempOp,buffer+8,6);
		strncpy(TempSym,buffer,6);
		strncpy(TempOpSym,buffer+16,8);
		strncpy(TempOpSym2,buffer+25,7);
		for(int j=0;TempSym[j]!='\0';j++){
		if(TempSym[j]==' ')
			TempSym[j]='\0';
		}
		for(int j=0;TempOp[j]!='\0';j++){
		if(TempOp[j]==' ')
			TempOp[j]='\0';
		}
		for(int j=0;TempOpSym[j]!='\0';j++){
		if(TempOpSym[j]==' ')
			TempOpSym[j]='\0';
		}
		for(int j=0;TempOpSym2[j]!='\0';j++){
		if(TempOpSym2[j]==' ')
			TempOpSym2[j]='\0';
		}
		if(strcmp(TempOp,"START")==0){
			fprintf(fpWRC,"H%6s%.6x%.6x\n",TempSym,atoi(TempOpSym),LocCounter);
			strncpy(start,TempOpSym,2);
			LocCounter = 0x00;
			strncpy(TempUse,"defau",5);
			fprintf(fpWOR,"%d   %s\n",line,buffer );	
		}
		else if (strcmp(TempOp,"USE")==0){
			fprintf(fpWOR,"%d   %s\n",line,buffer );
			if(TempOpSym[0]=='\0'){
			 BloCount[atoi(findValueByKey(&BloNumTable,TempUse))]= LocCounter;
			 LocCounter = BloCount[atoi(findValueByKey(&BloNumTable,"defau"))];
			 strncpy(TempUse,"defau",5);
			}
			else if((FlagA==1)&&(FlagB==1)){
			  BloCount[atoi(findValueByKey(&BloNumTable,TempUse))]= LocCounter;
			  LocCounter= BloCount[atoi(findValueByKey(&BloNumTable,TempOpSym))];
			  strncpy(TempUse,TempOpSym,5);	
			}
			else {	
			  BloCount[atoi(findValueByKey(&BloNumTable,TempUse))]= LocCounter;
			  LocCounter = 0x00;
			  BloCounter = BloCounter+1;
			  strncpy(TempUse,TempOpSym,5);
			  if(strcmp(TempOpSym,"CDATA")==0)
			  FlagA=1;
			  if(strcmp(TempOpSym,"CBLKS")==0)
			  FlagB=1;
			}
		}
		else if(findValueByKey(&OP,TempOp)!=NULL){
			if(buffer[7]=='+'){
				LocCounter = LocCounter + 0x04;
				strncpy(Change,findValueByKey(&OP,TempOp)+4,2);
				strcat(Dex1,Change);
				TempInt3 = strtol(Dex1, NULL, 16);
				for(int j=2;Dex1[j]!='\0';j++)
			    Dex1[j]='\0';
		        for(int j=0;Change[j]!='\0';j++)
			    Change[j]='\0';
				if(buffer[15]==' '){
					strncpy(Change,findValueByKey(&SymTable,TempOpSym),4);
				    strcat(Dex1,Change);
				    TempInt = strtol(Dex1, NULL, 16);
				    for(int j=2;Dex1[j]!='\0';j++)
			         Dex1[j]='\0';
		            for(int j=0;Change[j]!='\0';j++)
			         Change[j]='\0';
					ObjCode = ((TempInt3+0x03)<<24)+(0x01<<20)+TempInt;
					if(mode=='0'){
						ModiCounter[Modi] = LocCounter;
						Modi = Modi+1;
					}
				}
				if(buffer[15]=='#'){
					if(findValueByKey(&SymTable,TempOpSym)!=NULL){
						strncpy(Change,findValueByKey(&SymTable,TempOpSym),4);
				        strcat(Dex1,Change);
					    TempInt = strtol(Dex1, NULL, 16);
					    for(int j=2;Dex1[j]!='\0';j++)
			             Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			             Change[j]='\0';
			            ObjCode = ((TempInt3+0x01)<<24)+(0x01<<20)+TempInt;
					}	
					else if(findValueByKey(&LitTable,TempOpSym)!=NULL){ 
						strncpy(Change,findValueByKey(&LitTable,TempOpSym),4);
				        strcat(Dex1,Change);
					    TempInt = strtol(Dex1, NULL, 16);
				 	    for(int j=2;Dex1[j]!='\0';j++)
			             Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			             Change[j]='\0';
					    ObjCode = ((TempInt3+0x01)<<24)+(0x01<<20)+TempInt;
					}                                                      
					else{
						strcat(Dex1,TempSym);
						TempInt = strtol(Dex1, NULL, 16);
					    for(int j=2;Dex1[j]!='\0';j++)
			             Dex1[j]='\0';
						ObjCode = ((TempInt3+0x01)<<24)+(0x01<<20)+TempInt; 
				    }
				}
				if(buffer[15]=='@'){
					strncpy(Change,findValueByKey(&SymTable,TempSym),4);
				    strcat(Dex1,Change);
					TempInt = strtol(Dex1, NULL, 16);
					for(int j=2;Dex1[j]!='\0';j++)
			         Dex1[j]='\0';
		            for(int j=0;Change[j]!='\0';j++)
			         Change[j]='\0';
					if(((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800))
					ObjCode = ((TempInt3+0x02)<<24)+(0x03<<20)+TempInt-LocCounter;
					else
					ObjCode = ((TempInt3+0x02)<<24)+(0x03<<20)+TempInt-BaseAddr;
				}
				fprintf(fpWRC,"T%.6x04%x\n",LocCounter,ObjCode);
				fprintf(fpWOR,"%d   %s   %x\n",line,buffer,ObjCode );
			}
			else{
				strncpy(Change,findValueByKey(&OP,TempOp)+4,2);
				strcat(Dex1,Change);
				TempInt3 = strtol(Dex1, NULL, 16);
				for(int j=2;Dex1[j]!='\0';j++)
			    Dex1[j]='\0';
		        for(int j=0;Change[j]!='\0';j++)
			    Change[j]='\0';
				if(findValueByKey(&OP,TempOp)[0]=='2'){
					LocCounter = LocCounter + 0x02;
					strncpy(Change,findValueByKey(&OP,TempOp)+2,2);
				    strcat(Dex1,Change);
				    TempInt3 = strtol(Dex1, NULL, 16);
				    for(int j=2;Dex1[j]!='\0';j++)
			        Dex1[j]='\0';
		            for(int j=0;Change[j]!='\0';j++)
			        Change[j]='\0';
					if(TempOpSym2[0]=='\0'){
						strncpy(Change,findValueByKey(&REG,TempOpSym),2);
				        strcat(Dex1,Change);
				        TempInt2 = strtol(Dex1, NULL, 16);
				        for(int j=2;Dex1[j]!='\0';j++)
			            Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			            Change[j]='\0';
						ObjCode = (TempInt3<<8)+TempInt2;
					}
					else
					{
						strncpy(Change,findValueByKey(&REG,TempOpSym),1);
				        strcat(Dex1,Change);
				        TempInt = strtol(Dex1, NULL, 16);
				        for(int j=2;Dex1[j]!='\0';j++)
			            Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			            Change[j]='\0';
			            strncpy(Change,findValueByKey(&REG,TempOpSym2),1);
				        strcat(Dex1,Change);
				        TempInt2 = strtol(Dex1, NULL, 16);
				        for(int j=2;Dex1[j]!='\0';j++)
			            Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			            Change[j]='\0';
						ObjCode = (TempInt3<<8)+(TempInt<<4)+TempInt2;
					}
					fprintf(fpWRC,"T%.6x02%.4x\n",LocCounter,ObjCode);
				    fprintf(fpWOR,"%d   %s   %x\n",line,buffer,ObjCode );
				}
				else{
					LocCounter = LocCounter + 0x03;
					if(buffer[15]==' '){
						if(findValueByKey(&SymTable,TempOpSym)!=NULL){
							strncpy(Change,findValueByKey(&SymTable,TempOpSym),4);
				            strcat(Dex1,Change);
				            TempInt = strtol(Dex1, NULL, 16);
				            for(int j=2;Dex1[j]!='\0';j++)
			                    Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                    Change[j]='\0';
							if(mode=='1'){
								strncpy(Chan,findValueByKey(&SymTable,TempOpSym),5);
								for(int i=strlen(Chan)-1;Chan[i]!='\0';i++)
								BloFlag[0]=Chan[i];
								strncpy(Change,findValueByKey(&SymTable,TempOpSym),strlen(Chan)-1);
				                strcat(Dex1,Change);
				                TempInt = strtol(Dex1, NULL, 16);
				                for(int j=2;Dex1[j]!='\0';j++)
			                     Dex1[j]='\0';
		                        for(int j=0;Change[j]!='\0';j++)
			                     Change[j]='\0';
			                    strncpy(Change,findValueByKey(&BloTable,findValueByKey(&BloTableNum,BloFlag)),4);
			                    strcat(Dex1,Change);
			                    TempInt = TempInt+strtol(Dex1, NULL, 16);
			                    for(int j=2;Dex1[j]!='\0';j++)
			                     Dex1[j]='\0';
		                        for(int j=0;Change[j]!='\0';j++)
			                     Change[j]='\0';
			                    for(int j=0;Chan[j]!='\0';j++)
			                     Chan[j]='\0';
							}
						}
						else{
							strncpy(Change,findValueByKey(&LitTable,TempOpSym),4);
				            strcat(Dex1,Change);
				            TempInt = strtol(Dex1, NULL, 16);
				 	        for(int j=2;Dex1[j]!='\0';j++)
			                Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                Change[j]='\0';
			                if(mode=='1'){
								strncpy(Chan,findValueByKey(&LitTable,TempOpSym),5);
								for(int i=strlen(Chan)-1;Chan[i]!='\0';i++)
								BloFlag[0]=Chan[i];
								strncpy(Change,findValueByKey(&LitTable,TempOpSym),strlen(Chan)-1);
				                strcat(Dex1,Change);
				                TempInt = strtol(Dex1, NULL, 16);
				                for(int j=2;Dex1[j]!='\0';j++)
			                     Dex1[j]='\0';
		                        for(int j=0;Change[j]!='\0';j++)
			                     Change[j]='\0';
			                    strncpy(Change,findValueByKey(&BloTable,findValueByKey(&BloTableNum,BloFlag)),4);
			                    strcat(Dex1,Change);
			                    TempInt = TempInt+strtol(Dex1, NULL, 16);
			                    for(int j=2;Dex1[j]!='\0';j++)
			                     Dex1[j]='\0';
		                        for(int j=0;Change[j]!='\0';j++)
			                     Change[j]='\0';
			                    for(int j=0;Chan[j]!='\0';j++)
			                     Chan[j]='\0';
			                
							}
						}
						if(buffer[16]==' ')
						ObjCode = ((TempInt3+0x03)<<16)+(0x00<<12);
						else{
						  if(TempOpSym2[0]!='\0')
						  {
							if((((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800))||(mode==1))
								ObjCode = ((TempInt3+0x03)<<16)+(0x0a<<12)+TempInt-LocCounter;
							else
								ObjCode = ((TempInt3+0x03)<<16)+(0x0c<<12)+TempInt-BaseAddr;
						  } 
						  else{
						  	if(mode=='1')
						  	{
						  		ObjCode = ((TempInt3+0x03)<<16)+(0x02<<12)+TempInt-LocCounter;
						  		if((TempInt-LocCounter)<0x00)
									  ObjCode = (0x01<<12)+ObjCode;
						    }
						  	else{
							  if(((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800)){
								  ObjCode = ((TempInt3+0x03)<<16)+(0x02<<12)+TempInt-LocCounter;
								  if((TempInt-LocCounter)<0x00)
									  ObjCode = (0x01<<12)+ObjCode;
							  }
							  else{
								  ObjCode = ((TempInt3+0x03)<<16)+(0x04<<12)+TempInt-BaseAddr;
								  if((TempInt-BaseAddr)<0x00)
					        	  ObjCode = (0x01<<12)+ObjCode;
							  }
						    }
						  }
					    }
					}
				    if(buffer[15]=='#'){
					    if(findValueByKey(&LitTable,TempOpSym)!=NULL){
				            strncpy(Change,findValueByKey(&LitTable,TempOpSym),4);
				            strcat(Dex1,Change);
				            TempInt = strtol(Dex1, NULL, 16);
				 	        for(int j=2;Dex1[j]!='\0';j++)
			                Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                Change[j]='\0';
			                if(((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800))
							ObjCode = ((TempInt3+0x01)<<16)+(0x02<<12)+TempInt-LocCounter;
						    else
							ObjCode = ((TempInt3+0x01)<<16)+(0x04<<12)+TempInt-BaseAddr;	
						}	
					    else if(findValueByKey(&SymTable,TempOpSym)!=NULL){
					    	strncpy(Change,findValueByKey(&SymTable,TempOpSym),4);
				            strcat(Dex1,Change);
				            TempInt = strtol(Dex1, NULL, 16);
				 	        for(int j=2;Dex1[j]!='\0';j++)
			                Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                Change[j]='\0';
					        if(((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800))
					        ObjCode = ((TempInt3+0x01)<<16)+(0x02<<12)+TempInt-LocCounter;
						    else
							ObjCode = ((TempInt3+0x01)<<16)+(0x04<<12)+TempInt-BaseAddr;
						}
						else{
							strcat(Dex1,TempOpSym);
				            TempInt = strtol(Dex1, NULL, 16);
				 	        for(int j=2;Dex1[j]!='\0';j++)
			                Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                Change[j]='\0';
							ObjCode = ((TempInt3+0x01)<<16)+(0x00<<12)+TempInt;
						}    
				    }
				    if(buffer[15]=='@'){
				    	strncpy(Change,findValueByKey(&LitTable,TempOpSym),6);
				        strcat(Dex1,Change);
				        TempInt = strtol(Dex1, NULL, 16);
				 	    for(int j=2;Dex1[j]!='\0';j++)
			            Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			            Change[j]='\0';
			            if(mode=='1'){
							strncpy(Chan,findValueByKey(&LitTable,TempOpSym),5);
							for(int i=strlen(Chan)-1;Chan[i]!='\0';i++)
							BloFlag[0]=Chan[i];
							strncpy(Change,findValueByKey(&LitTable,TempOpSym),strlen(Chan)-1);
				            strcat(Dex1,Change);
				            TempInt = strtol(Dex1, NULL, 16);
				            for(int j=2;Dex1[j]!='\0';j++)
			                    Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                    Change[j]='\0';
			                strncpy(Change,findValueByKey(&BloTable,findValueByKey(&BloTableNum,BloFlag)),4);
			                strcat(Dex1,Change);
			                TempInt = TempInt+strtol(Dex1, NULL, 16);
			                for(int j=2;Dex1[j]!='\0';j++)
			                    Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                    Change[j]='\0';
			                for(int j=0;Chan[j]!='\0';j++)
			                    Chan[j]='\0';
						}
					    if((((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800))||(mode==1))
							ObjCode = ((TempInt3+0x02)<<16)+(0x02<<12)+TempInt-LocCounter;
						else
							ObjCode = ((TempInt3+0x02)<<16)+(0x04<<12)+TempInt-BaseAddr;
				    }
				    if(buffer[15]=='='){
				    	if(findValueByKey(&LitTable,TempOpSym)!=NULL)
				    	strncpy(Change,findValueByKey(&LitTable,TempOpSym),6);
				    	else
				    	{
				    	 char tem[5];
				    	 itoa(ProLenth,tem,16);
				    	 strncpy(Change,tem,6);
						 }
				        strcat(Dex1,Change);
				        TempInt = strtol(Dex1, NULL, 16);
				 	    for(int j=2;Dex1[j]!='\0';j++)
			            Dex1[j]='\0';
		                for(int j=0;Change[j]!='\0';j++)
			            Change[j]='\0';
			            if(mode=='1'){
							strncpy(Chan,findValueByKey(&LitTable,TempOpSym),5);
							for(int i=strlen(Chan)-1;Chan[i]!='\0';i++)
							BloFlag[0]=Chan[i];
							strncpy(Change,findValueByKey(&LitTable,TempOpSym),strlen(Chan)-1);
							strcat(Dex1,Change);
				            TempInt = strtol(Dex1, NULL, 16);
				            for(int j=2;Dex1[j]!='\0';j++)
			                    Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                    Change[j]='\0';
			                strncpy(Change,findValueByKey(&BloTable,findValueByKey(&BloTableNum,BloFlag)),4);
			                strcat(Dex1,Change);
			                TempInt = TempInt+strtol(Dex1, NULL, 16);
			                for(int j=2;Dex1[j]!='\0';j++)
			                    Dex1[j]='\0';
		                    for(int j=0;Change[j]!='\0';j++)
			                    Change[j]='\0';
			                for(int j=0;Chan[j]!='\0';j++)
			                    Chan[j]='\0';
						}
					    if((((TempInt-LocCounter)<=0x7FF)&&((TempInt-LocCounter)>=-0x800))||(mode==1)){
					    	ObjCode = ((TempInt3+0x03)<<16)+(0x02<<12)+TempInt-LocCounter; 
						}
						else
						ObjCode = ((TempInt3+0x03)<<16)+(0x04<<12)+TempInt-BaseAddr;
					}
					fprintf(fpWRC,"T%.6x03%x\n",LocCounter,ObjCode);
				    fprintf(fpWOR,"%d   %s   %.6x\n",line,buffer,ObjCode );
				}
			}
		}
		else if(strcmp(TempOp,"END")==0){
			if(mode=='1'){
		 	 if(EquLit2[0]!='\0'){
			 LocCounter = LocCounter + strlen(EquLit2)-3;
			 char tem2[4];
			 for(int i=2;i<=strlen(EquLit2)-2;i++)
			 tem2[i-2]=EquLit2[i];
			 fprintf(fpWRC,"T%.6x%.2d%s\n",LocCounter,strlen(EquLit2)-3,tem2);
		     fprintf(fpWOR,"%d   *     =%s   \n%s\n",line,EquLit2,tem2);
		     }
		    }
		    line = line+1;
		    fprintf(fpWOR,"%d   %s   \n",line,buffer);
		}
		else if(strcmp(TempOp,"LTORG")==0)
		{
		 fprintf(fpWOR,"%d   %s   \n",line,buffer);
		 if(EquLit[0]!='\0'){
		 	
			line=line+1;
			LocCounter = LocCounter + strlen(EquLit)-3;
			char tem1[4];
			for(int i=2;i<=strlen(EquLit);i++)
			tem1[i-2]=EquLit[i];
			fprintf(fpWRC,"T%.6x%.2d%s\n",LocCounter,strlen(EquLit)-3,"454F46");
		    fprintf(fpWOR,"%d   *     =%s   %.6\n%s\n",line,EquLit,"454F46" );
		    }
		}
		else if(findValueByKey(&LitTable,TempSym)!=NULL)
		{
			if(strcmp(TempOp,"WORD")==0)
			LocCounter = LocCounter + 0x03;
		    if(strcmp(TempOp,"RESW")==0)
			LocCounter = LocCounter + 3*atoi(TempOpSym);
	        if(strcmp(TempOp,"RESB")==0)
	        LocCounter = LocCounter + 1*atoi(TempOpSym);
		    if(strcmp(TempOp,"BYTE")==0)
		    LocCounter = LocCounter + (strlen(TempOpSym)-3)/2;
			fprintf(fpWOR,"%d   %s  \n",line,buffer);
		}
		else if(strcmp(TempOp,"BASE")==0)
		fprintf(fpWOR,"%d   %s   \n",line,buffer);
		else if(strcmp(TempOp,"EQU")==0)
		fprintf(fpWOR,"%d   %s   \n",line,buffer);
		else if(findValueByKey(&SymTable,TempSym)!=NULL){
			char char1[3];
			for(int i=2;;i++){
			 if(ispunct(TempOpSym[i])!=0)
			 break;
			 char1[i-2] =TempOpSym[i];
			 }
			LocCounter = LocCounter + (strlen(TempOpSym)-3)/2;
			fprintf(fpWOR,"%d   %s  %s\n",line,buffer,char1);
			fprintf(fpWRC,"T%.6x%.2d%s\n",LocCounter,(strlen(TempOpSym)-3)/2,char1);
		}
		else if(buffer[7]=='='){
			int ascii;
			char Combine[3];
			char Combine2[8];
			ascii=TempOp[2];
			itoa(ascii,Combine2,16);
			for(int i=3;;i++){
				if(ispunct(TempOp[i])!=0)
				break;
				ascii=TempOp[i];
				itoa(ascii,Combine,16);
				strcat(Combine2,Combine);
				}
			LocCounter = LocCounter + strlen(TempOp)-3;
			fprintf(fpWRC,"T%.6x03%s\n",LocCounter,Combine2);
			fprintf(fpWOR,"%d   %s  %s\n",line,buffer,Combine2);
		}
		else
		printf("error\n");
	}
	//you may change the threshold of i in case of much more command need to modify.
	if (mode=='0'){
	 for(int i=0;i<=2;i++)
	 fprintf(fpWRC,"M%.6x05\n",ModiCounter[i]-0x03);	
	}
	fprintf(fpWRC,"E%.6d",atoi(start));
	printf("we have finished second circle");
	return 0;
}

  
