
import sys, re


class WSString:

    def __init__(self, val):

        # in reality this is going to have to be a very complex method
        # that creates the str based on a number of different possible
        # values

        if (isinstance(val,dict)):

            self.str = self.hashMapToString(val)

        elif (isinstance(val,list)):
        
            self.str = self.listToString(val)
            
        #elif (isinstance(val,str)):

        else :
            
            self.str = val

    # public methods

    def getHashMap(self):

        return self.stringToHashMap(self.str)

    def getArrayList(self):

        return self.stringToList(self.str)


    # For legacy. remove
        
    def getStructHM(self):

        return self.stringToHashMap(self.str)


    def getStruct(self):

        return self.stringToList(self.str)



    # attempt to write the method properly
    # what does the index do?

    def listToString(self, lst, deep=0, index=0):

        strg = ''
        lstSize = len(lst)
        subIndex = 0
        fixDeep = deep

        strg += '<Arraylist'+`fixDeep`+'.'+`index`+'>'
        deep += 1

        for i in range(0, lstSize):

            if isinstance(lst[i],dict):
                submap = lst[i]
                if len(submap.keys() > 0):
                    subIndex += 1
                    subString = self.hashMapToString(submap, deep, subIndex)
                    strg += subString

            elif isinstance(lst[i],list):
                sublst = lst[i]
                if len(sublst > 0):
                    subIndex += 1
                    subString = self.listToString(sublst, deep, subIndex)
                    strg += subString
                    
            # not sure why we don't use the types here
            else:
                strg += lst[i]
                        

            if i < lstSize -1:
                strg += '<#>'
                    
        strg += '<Arraylist'+`fixDeep`+'.'+`index`+'/>'

        return strg
                    

    # attempt to write the method properly


    def hashMapToString(self, hashmap, deep=0, index=0):

        strg = ''
        subIndex = 0
        fixDeep = deep

        # is this really correct??

        if len(hashmap.keys()) > 0:

            set = hashmap.keys()
            nbKeys = len(set)

            #print set, nbKeys

            strg += '<Hashmap'+`fixDeep`+'.'+`index`+'>'
            deep += 1

            for i in range(nbKeys):

                #print 'looping ', fixDeep, i
            
                strg += set[i]
                strg += '<=>'

                if hashmap[set[i]] == None:
                    strg += 'valIsNullPointer'

                else:

                    #print 'appending ', hashmap[set[i]] , ' (from ',  set[i], ')' 
                    if isinstance(hashmap[set[i]],str):
                        strg += hashmap[set[i]]
                        
                    elif isinstance(hashmap[set[i]],dict):
                        submap = hashmap[set[i]]
                        if len(submap.keys()) > 0:
                            subIndex += 1
                            subString = self.hashMapToString(submap, deep, subIndex)
                            strg += subString
                        
                    elif isinstance(hashmap[set[i]],list):
                        sublist = hashmap[set[i]]
                        if len(sublist) > 0:
                            subIndex += 1
                            subString = self.listToString(sublist, deep, subIndex)
                            strg += subString

                    elif isinstance(hashmap[set[i]],bool):
                        if hashmap[set[i]] == True:
                            strg += 'booleanTrue'
                        else:
                            strg += 'booleanFalse'

                    elif isinstance(hashmap[set[i]],int):
                        strg += 'valIsInteger'
                        strg += `hashmap[set[i]]`

                    elif isinstance(hashmap[set[i]],float):
                        strg += 'valIsFloat'
                        strg += hashmap[set[i]]

                if i < nbKeys -1:
                    strg += '<;>'
                    
            strg += '<Hashmap'+`fixDeep`+'.'+`index`+'/>'

        return strg
                    

                


    # attempt to write the method properly


    def stringToList(self, strg):

        lst = []
        regex = ''

        if len(strg) > 0:            
            regex = '<Arraylist[0-9]+.[0-9]+>'
            pattern = re.compile(regex)
            matcher = pattern.search(strg, 0)
            #print ' first match ', matcher.start()

            if matcher.start() == 0:
                startTag = matcher.group()
                #print 'group', startTag
                endTag = startTag[:-1] + '/>'
                #print 'group', endTag
                startSeq = matcher.end()
                endSeq= strg.find(endTag,startSeq)
                strg=strg[startSeq:endSeq]

                while len(strg) > 0:

                    cutAt = strg.find('<#>')

                    #print ' value from ', strg

                    if strg[:8] == '<Hashmap':
                        value = self.stringToHashMap(strg)
                        #print ' getting hashmap from ', strg
                        lst.append(value)
                
                        regex = '<Hashmap[0-9]+.[0-9]+>'
                        pattern = re.compile(regex)
                        matcher = pattern.search(strg, 0)
                        endTag = matcher.group()[:-1] + '/>'
                        endSeq= strg.find(endTag,startSeq)
                        strg=strg[endSeq + len(endTag):]
                        # print ' trimming string to: ', strg, startSeq, endSeq

                        cutAt = strg.find('<#>')
                        if cutAt == 0:
                            strg = strg[3:]

                        #print ' returning ', strg
                        
                    elif strg[:10] == '<Arraylist':
                        value = self.stringToList(strg)
                        lst.append(value)
                            
                        regex = '<Arraylist[0-9]+.[0-9]+>'
                        pattern = re.compile(regex)
                        matcher = pattern.search(strg, 0)
                        endTag = matcher.group()[:-1] + '/>'
                        endSeq= strg.find(endTag,startSeq)
                        strg=strg[endSeq + len(endTag):]
                            
                        cutAt = strg.find('<#>')
                        if cutAt == 0:
                            strg = strg[3:]

                    elif strg[:13] == 'valIsInteger':
                        cutAt = strg.find('<#>')
                        if cutAt < 0:
                            value = strg[13:]
                            strg = ''
                        else:
                            value = strg[13:cutAt]
                            strg = strg[cutAt+3:]
                        lst.append(value)

                        # String
                    else:
                        cutAt = strg.find('<#>')
                        if cutAt < 0:
                            value = strg
                            strg = ''
                        else:
                            value = strg[:cutAt]
                            strg = strg[cutAt+3:]

                        lst.append(value)
                    
        return lst
                

    # attempt to write the method properly

    def stringToHashMap(self, strg):

        map = {}
        regex = ''

        if len(strg) > 0:            
            regex = '<Hashmap[0-9]+.[0-9]+>'
            pattern = re.compile(regex)
            matcher = pattern.search(strg, 0)
            #print ' first match ', matcher.start()

            if matcher.start() == 0:
                startTag = matcher.group()
                #print 'group', startTag
                endTag = startTag[:-1] + '/>'
                #print 'group', endTag
                startSeq = matcher.end()
                endSeq= strg.find(endTag,startSeq)
                strg=strg[startSeq:endSeq]

                # print strg

                while len(strg) > 0:

                    cutAt = strg.find('<=>')
                    key = strg[0:cutAt]
                    # print 'KEY VALUE ', key
                    strg=strg[cutAt+3:]

                    #print ' value from ', strg

                    # print strg[:8]
                    if strg[:8] == '<Hashmap':
                        value = self.stringToHashMap(strg)
                        #print ' getting hashmap from ', strg
                        map[key]=value
                
                        regex = '<Hashmap[0-9]+.[0-9]+>'
                        pattern = re.compile(regex)
                        matcher = pattern.search(strg, 0)
                        endTag = matcher.group()[:-1] + '/>'
                        endSeq= strg.find(endTag,startSeq)
                        strg=strg[endSeq + len(endTag):]
                        #print ' trimming string to: ', strg, startSeq, endSeq

                        cutAt = strg.find('<;>')
                        if cutAt == 0:
                            strg = strg[3:]

                        #print ' returning ', strg
                        
                    elif strg[:10] == '<Arraylist':
                        value = self.stringToList(strg)
                        map[key]=value
                            
                        regex = '<Arraylist[0-9]+.[0-9]+>'
                        #print 'ArrayList group', matcher.group()
                        pattern = re.compile(regex)
                        matcher = pattern.search(strg, 0)
                        endTag = matcher.group()[:-1] + '/>'
                        endSeq= strg.find(endTag,startSeq)
                        strg=strg[endSeq + len(endTag):]
                            
                        cutAt = strg.find('<;>')
                        if cutAt == 0:
                            strg = strg[3:]

                    elif strg[:12] == 'valIsInteger':
                        cutAt = strg.find('<;>')
                        if cutAt < 0:
                            value = strg[12:]
                            strg = ''
                        else:
                            value = strg[12:cutAt]
                            strg = strg[cutAt+3:]
                        map[key] = value

                        # String
                    else:
                        cutAt = strg.find('<;>')
                        if cutAt < 0:
                            value = strg
                            strg = ''
                        else:
                            value = strg[:cutAt]
                            strg = strg[cutAt+3:]
                        map[key] = value
                    
        return map
        
                

