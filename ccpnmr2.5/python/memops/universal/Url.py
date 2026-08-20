"""
======================COPYRIGHT/LICENSE START==========================

Url.py: Url code for CCPN

Copyright (C) 2003-2010  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""


def fetchUrl(url, values=None, headers=None, timeout=None):

    return _fetchUrl(url, data=values, headers=headers, timeout=timeout)

    # import socket
    # import urllib
    # import urllib2
    #
    # if not headers:
    #     headers = {}
    #
    # try:
    #     # from Python 2.6 there is a timeout option in urlopen()
    #     # but for now assume Python 2.5 compatibility
    #     oldTimeout = socket.getdefaulttimeout()
    #     socket.setdefaulttimeout(timeout)
    #     if values:
    #         data = {}
    #         for key in values:
    #             value = values[key]
    #             if isinstance(value, unicode):
    #                 value = value.encode('utf-8')
    #             data[key] = value
    #         data = urllib.urlencode(data)
    #     else:
    #         data = None
    #     request = urllib2.Request(url, data, headers)
    #     response = urllib2.urlopen(request)
    #     result = response.read()
    # finally:
    #     socket.setdefaulttimeout(oldTimeout)
    #
    # return result


# uploadFiles() is slightly modified version of code from Michael Foord
# and that said:

# Copyright Michael Foord, 2004 & 2005.
# Released subject to the BSD License
# Please see http://www.voidspace.org.uk/documents/BSD-LICENSE.txt
# Based on http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/146306
# With inspiration from urllib2_file.py by Fabien Seisen, http://fabien.seisen.org/python/ (without his twiddly bits)
# It actually uses my upload script located at http://www.voidspace.xennos.com

def uploadFiles(url, fileFields, fields=None, boundary=None):
    """Uploads regular fields and files to specified url
    url is the Url to send data to.
    fileFields is a sequence of (name, fileName) elements, or the dictionary
    equivalent, for file form fields.
    fields is a sequence of (name, value) elements, or the dictionary equivalent,
    for regular form fields.
    Returns response."""

    import mimetypes
    import mimetools
    import os
    import urllib2

    if not fields:
        fields = ()

    if not boundary:
        boundary = '-----' + mimetools.choose_boundary() + '-----'

    CRLF = '\r\n'
    xx = []
    if isinstance(fields, dict):
        fields = fields.items()
    for (key, value) in fields:
        xx.append('--' + boundary)
        xx.append('Content-Disposition: form-data; name="%s"' % key)
        xx.append('')
        xx.append(str(value))

    if isinstance(fileFields, dict):
        fileFields = fileFields.items()
    for (key, fileName) in fileFields:
        fp = open(fileName, 'rb')
        value = fp.read()
        fp.close()

        fileName = os.path.basename(fileName)
        fileType = mimetypes.guess_type(fileName)[0] or 'application/octet-stream'
        xx.append('--' + boundary)
        xx.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, fileName))
        xx.append('Content-Type: %s' % fileType)
        xx.append('')
        xx.append(value)

    xx.append('--' + boundary + '--')
    xx.append('')
    body = CRLF.join(xx)

    contentType = 'multipart/form-data; boundary=%s' % boundary
    headers = {
        'Content-type'  : contentType,
        'Content-length': str(len(body)),
        }

    request = urllib2.Request(url, body, headers)
    handle = urllib2.urlopen(request)
    try:
        result = handle.read()
    finally:
        handle.close()

    return result


def uploadFile(url, fileKey, fileName, fields=None, boundary=None):
    """Uploads regular fields and file to specified url
    url is the Url to send data to.
    fileKey is the form key for the file.
    fileName is the full path for the file.
    fields is a sequence of (name, value) elements, or the dictionary equivalent,
    for regular form fields.
    Returns response."""

    return uploadFiles(url, ((fileKey, fileName),), fields, boundary)


def fetchHttpResponse(url, method='GET', data=None, headers=None):
    """Generate an http, and return the response
    """
    import ssl
    import os
    import certifi
    import urllib
    import urllib3.contrib.pyopenssl
    from urllib import urlencode

    if not headers:
        headers = {'Content-type': 'application/x-www-form-urlencoded;charset=UTF-8',
                    'User-Agent': 'ccpn-v2.5.0'
        }
    body = urlencode(data).encode('utf-8') if data else None

    # urllib3.contrib.pyopenssl.inject_into_urllib3()

    from urllib3.util import ssl_
    if ssl_.IS_PYOPENSSL:
        import urllib3.contrib.pyopenssl
        urllib3.contrib.pyopenssl.extract_from_urllib3()

    context = ssl._create_unverified_context()
    # context = ssl.create_default_context(purpose=ssl.Purpose.SERVER_AUTH,
    #                                      cafile=None,
    #                                      capath=None)

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    # context = ssl.SSLContext(ssl.PROTOCOL_TLS)
    # context.verify_mode = ssl.CERT_NONE
    # context.options |= ssl.OP_NO_SSLv2 | ssl.OP_NO_SSLv3 | ssl.OP_NO_COMPRESSION
    # context.load_default_certs()

    # _http_pid = os.getpid()
    # _http = urllib3.PoolManager(ssl_context=context)

    # create the options list for creating an http connection
    options = {#'cert_reqs': 'NONE',
               #'ca_certs' : certifi.where(),
               'timeout'  : urllib3.Timeout(connect=3.0, read=3.0),
               'retries'  : urllib3.Retry(1, redirect=False),
               'ssl_context'  : context
               }

    def _getProxyIn(proxyDict):
        """Get the first occurrence of a proxy type in the supplied dict
        """
        # define a list of proxy identifiers
        proxyCheckList = ['HTTPS_PROXY', 'https', 'HTTP_PROXY', 'http']
        for pCheck in proxyCheckList:
            proxyUrl = proxyDict.get(pCheck, None)
            if proxyUrl:
                return proxyUrl

    proxyUrl = _getProxyIn(os.environ) or _getProxyIn(urllib.getproxies())
    # proxy may still not be defined
    if proxyUrl:
        http = urllib3.ProxyManager(proxyUrl, **options)
    else:
        http = urllib3.PoolManager(**options)

    # generate an http request
    response = http.request(method, url,
                            headers=headers,
                            body=body,
                            preload_content=False)
    response.release_conn()

    # return the http response
    return response


def _fetchUrl(url, data=None, headers=None, timeout=2.0, proxySettings=None, decodeResponse=True):
    """Fetch url request from the server
    """
    import logging

    urllib3_logger = logging.getLogger('urllib3')
    urllib3_logger.setLevel(logging.CRITICAL)

    # if not proxySettings:
    #
    #     # read the proxySettings from the preferences
    #     from ccpn.util.UserPreferences import UserPreferences
    #
    #     _userPreferences = UserPreferences(readPreferences=True)
    #     if _userPreferences.proxyDefined:
    #         proxyNames = ['useProxy', 'proxyAddress', 'proxyPort', 'useProxyPassword',
    #                       'proxyUsername', 'proxyPassword']
    #         proxySettings = {}
    #         for name in proxyNames:
    #             proxySettings[name] = _userPreferences._getPreferencesParameter(name)

    response = fetchHttpResponse(url, method='POST', data=data, headers=headers)

    # if response:
    #     ll = len(response.data)
    #     print('>>>>>>responseUrl', proxySettings, response.data[0:min(ll, 20)])

    return response.data.decode('utf-8') if decodeResponse else response


def uploadFile(url, fileName, data=None):
    import os

    if not data:
        data = {}

    with open(fileName, 'rb') as fp:
        fileData = fp.read()

    data['fileName'] = os.path.basename(fileName)
    data['fileData'] = fileData

    try:
        return fetchUrl(url, data)
    except:
        return None
