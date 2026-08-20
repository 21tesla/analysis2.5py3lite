import os, stat

def getDayStamp():

  import time
  
  dayStamp = time.strftime("%Y-%b-%d", time.gmtime(time.time()))
  
  return dayStamp

def getFileFromFtp(ftpSite,ftpDir,ftpFileName,localFileName):

  import ftplib

  FTP = ftplib.FTP

  ftp = FTP(ftpSite)
  ftp.login()
  ftp.cwd(ftpDir)

  fout = open(localFileName,'w')
  ftp.retrbinary("RETR %s" % ftpFileName, fout.write)
  fout.close()

  ftp.quit()

def getDataFromHttp(urlLocation):
  
  import urllib
  
  r1 = urllib.urlopen(urlLocation)
  code = r1.getcode()
  if (code < 400):
    data = r1.read()
    r1.close()
  else:
    data = None
    print "URL %s returned code %d" % (urlLocation, code)

  return data

def getTextFromHttp(urlLocation):
  
  data = getDataFromHttp(urlLocation)
  text = data.split("\n")
  
  return text

def getReferenceFileFromFtp(ftpSite,ftpDir,ftpFileName,localFtpFileName):

  # Might need to create this, based on function below... probably best to plug in code for file download, keep one def for workflow.
  pass
  
def getReferenceTextFileFromHttp(urlLocation,localFilePath, refText = "", isGzipped = False):
  
  """
  Function tries to get info online if required (updates every day), if this does not
  work it will use the latest available data from a locally saved file...
  """
  
  useLocalFile = False
  stampedLocalFilePath = "%s_%s" % (localFilePath,getDayStamp())
  
  if os.path.exists(stampedLocalFilePath):

    print "Using up-to-date local %s file..." % refText
  
  else:
    
    #
    # Get latest file - either to be used as information source or be deleted...
    #
    
    (path,baseName) = os.path.split(localFilePath)
    baseName = "%s_" % baseName
    
    lastStampedLocalFilePath = None
    
    localFiles = os.listdir(path)
    
    for localFile in localFiles:
      if localFile[:len(baseName)] == baseName:
        lastStampedLocalFilePath = os.path.join(path,localFile)
        break
  
    #
    # Now try to download, and save latest file
    #
    
    
    try:
      data = getDataFromHttp(urlLocation)
    except:
      data = None

    if data:
      
      if isGzipped:
        saveLocalFilePath = "%s.gz" % stampedLocalFilePath
      else:
        saveLocalFilePath = stampedLocalFilePath
      
      fout = open(saveLocalFilePath,'w')
      fout.write(data)
      fout.close()
      
      if lastStampedLocalFilePath:
        os.remove(lastStampedLocalFilePath)
        addText = ' (removed old file)'
      else:
        addText = ''
        
      print "Downloaded up-to-date file for %s, saved locally%s..." % (refText,addText)
      
      # Unpack
      if isGzipped:
        os.spawnlp(os.P_WAIT, 'gunzip', 'gunzip', saveLocalFilePath)
    
    else:

      if lastStampedLocalFilePath:
        stampedLocalFilePath = lastStampedLocalFilePath
        lastDayStamp = localFile[len(baseName):]
        print "Download failed for %s, using file stamped on day %s" % (refText,lastDayStamp)
      
      else:
        print "Error: could not get reference data for %s!" % refText
        return None

  #
  # Now get the data from the file
  #
  
  fin = open(stampedLocalFilePath)
  dataLines = fin.readlines()
  fin.close()
  
  return dataLines


"""

Code below from http://code.activestate.com/recipes/146306/

Allows uploading of file to web server, waits for result. See Coco.py for example for usage.

"""

def posturl(url, fields, files):

  import urlparse

  urlparts = urlparse.urlsplit(url)
  return post_multipart(urlparts[1],urlparts[2],fields,files)

def post_multipart(host, selector, fields, files):

  """
  Post fields and files to an http host as multipart/form-data.
  fields is a sequence of (name, value) elements for regular form fields.
  files is a sequence of (name, filename, value) elements for data to be uploaded as files
  Return the server's response page.
  """
  content_type, body = encode_multipart_formdata(fields, files)
  
  """
  # TODO replace by urllib2!?!? See http://docs.python.org/library/urllib2.html
  import urllib2
  print url
  request = urllib2.Request(url)
  request.add_header('content-type', content_type)
  request.add_header('content-length', str(len(body)))
  request.add_data(body)
  
  urlHandle = urllib2.urlopen(request)
  
  return urlHandle.geturl()
  
  """
  import httplib
  h = httplib.HTTP(host)
  h.putrequest('POST', selector)
  h.putheader('content-type', content_type)
  h.putheader('content-length', str(len(body)))
  h.endheaders()
  h.send(body)
  errcode, errmsg, headers = h.getreply()
  return h.file.read()

def encode_multipart_formdata(fields, files):
  """
  fields is a sequence of (name, value) elements for regular form fields.
  files is a sequence of (name, filename, value) elements for data to be uploaded as files
  Return (content_type, body) ready for httplib.HTTP instance
  """
  BOUNDARY = '----------ThIs_Is_tHe_bouNdaRY_$'
  CRLF = '\r\n'
  L = []
  for (key, value) in fields:
      L.append('--' + BOUNDARY)
      L.append('Content-Disposition: form-data; name="%s"' % key)
      L.append('')
      L.append(value)
  for (key, filename, value) in files:
      L.append('--' + BOUNDARY)
      L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename))
      L.append('Content-Type: %s' % get_content_type(filename))
      L.append('')
      L.append(value)
  L.append('--' + BOUNDARY + '--')
  L.append('')
  body = CRLF.join(L)
  content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
  return content_type, body

def get_content_type(filename):

  import mimetypes
  
  return mimetypes.guess_type(filename)[0] or 'application/octet-stream'


"""

Code below from http://peerit.blogspot.com/2007/07/multipartposthandler-doesnt-work-for.html

"""
####
# 02/2006 Will Holcomb <wholcomb@gmail.com>
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# 7/26/07 Slightly modified by Brian Schneider  

import urllib
import urllib2
import mimetools, mimetypes
from cStringIO import StringIO

class Callable:
    def __init__(self, anycallable):
        self.__call__ = anycallable

# Controls how sequences are uncoded. If true, elements may be given multiple values by
#  assigning a sequence.
doseq = 1

class MultipartPostHandler(urllib2.BaseHandler):
    handler_order = urllib2.HTTPHandler.handler_order - 10 # needs to run first

    def http_request(self, request):
        data = request.get_data()
        if data is not None and type(data) != str:
            v_files = []
            v_vars = []
            try:
                 for(key, value) in data.items():
                     if type(value) == file:
                         v_files.append((key, value))
                     else:
                         v_vars.append((key, value))
            except TypeError:
                systype, value, traceback = sys.exc_info()
                raise TypeError, "not a valid non-string sequence or mapping object", traceback

            if len(v_files) == 0:
                data = urllib.urlencode(v_vars, doseq)
            else:
                boundary, data = self.multipart_encode(v_vars, v_files)

                contenttype = 'multipart/form-data; boundary=%s' % boundary
                if(request.has_header('Content-Type')
                   and request.get_header('Content-Type').find('multipart/form-data') != 0):
                    print "Replacing %s with %s" % (request.get_header('content-type'), 'multipart/form-data')
                request.add_unredirected_header('Content-Type', contenttype)

            request.add_data(data)
        
        return request

    def multipart_encode(vars, files, boundary = None, buf = None):
        if boundary is None:
            boundary = mimetools.choose_boundary()
        if buf is None:
            buf = StringIO()
        for(key, value) in vars:
            buf.write('--%s\r\n' % boundary)
            buf.write('Content-Disposition: form-data; name="%s"' % key)
            buf.write('\r\n\r\n' + value + '\r\n')
        for(key, fd) in files:
            file_size = os.fstat(fd.fileno())[stat.ST_SIZE]
            filename = fd.name.split('/')[-1]
            contenttype = mimetypes.guess_type(filename)[0] or 'application/octet-stream'
            buf.write('--%s\r\n' % boundary)
            buf.write('Content-Disposition: form-data; name="%s"; filename="%s"\r\n' % (key, filename))
            buf.write('Content-Type: %s\r\n' % contenttype)
            # buffer += 'Content-Length: %s\r\n' % file_size
            fd.seek(0)
            buf.write('\r\n' + fd.read() + '\r\n')
        buf.write('--' + boundary + '--\r\n\r\n')
        buf = buf.getvalue()
        return boundary, buf
    multipart_encode = Callable(multipart_encode)

    https_request = http_request
