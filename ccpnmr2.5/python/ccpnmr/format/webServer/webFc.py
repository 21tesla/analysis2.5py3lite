# Import the cherrypy and mako modules. If everything
# went ok in the install these should be
# available on your system

import cherrypy
from mako.template import Template
from mako.lookup import TemplateLookup

#
# Import os.path for relative imports
#

import os.path

#
# Import authentication stuff from Util.py
#

from ccpnmr.format.webServer.Util import AuthController, require, is_current_session_user

#
# Secure convert/ layer, session specific
#

"""

All Python objects are hung off SessionInfo class (see Util.py)

Current security: one browser connects, user has to log in. This session will keep on working until user logs
out (or 1 hour of inactivity).
"""

class Convert(object):
 
  #
  # Local configuration information - is used for all methods of Convert class. This basically
  # defines a 'secure' section based on a session.
  #
  
  _cp_config = {
  'tools.auth.on': True,
  'auth.require': [is_current_session_user]  # Conditions for user to access these pages
  }
    
  # The @cherrypy.expose tells the server that this
  # method can be accessed from a webpage. We need to
  # put this before any method that we want to access
  # from a webpage

  @cherrypy.expose
  def index(self,*args,**keywords):
      # mylookup.get_template is mako functionality that
      # specifies the file you want to render. In this case
      # it pulls the file 'fcTest.html' from '/fchtml'

    if 'doRealSubmit' in args:
      formatConvert(args)
      self.nextFcPage()
    elif 'update' in args:
      regenerateThisPage()
    else:
      for key in keywords.keys():\
        print key, keywords[key]
      for arg in args:\
        print arg
       #createThisPageFirstTime()

    # now use the keywords to render the updated page
    # 
    mytemplate = mylookup.get_template("fcTest.html")
      # 'mytemplate.render()' is for rendering the file i.e.
      # doing any substitutions in it, and in this case,
      # creating a straight html file
    return mytemplate.render()

  @cherrypy.expose
  def fcOutputPage(self):
      # The 2 lines below get the page 'fcTestOutput.html'
      # and render it in the browser while passing the two
      # variables fcName and fcWorld
      # to the html file

      mytemplate=mylookup.get_template("fcTestOutput.html")
      return mytemplate.render(fcName="FormatConverter",
                               fcWorld="planet of cherries")


  # in case submit is pressed, same page is reloaded, so in case fcInputPage 
  #  submits data we can use the following construction:
  # 
  def fcInputPage(self,*args):
    if 'doRealSubmit' in args:
      formatConvert(args)
      self.nextFcPage()
    elif 'update' in args:
      regenerateThisPage()
    else:
      createThisPageFirstTime()


  # see if this will work in a general way:
  @cherrypy.expose
  def goToPage(self,pageName):
      # The 2 lines below get the page 'fcTestOutput.html'
      # and render it in the browser while passing the two
      # variables fcName and fcWorld
      # to the html file

      mytemplate=mylookup.get_template(pageName)
      return mytemplate.render()


class ConvertSession(Convert):
 
  # TODO TODO
  # Session based - no user authentication. Still need to do on_login for new session
  # and on_logout when finishing session (ability to delete data). Also mail link
  # to email address... always or if wanted?
  #
  
  _cp_config = {}
  
  #
  # TODO probably need different index page here, set sessionInfo if not done yet
  # No checking of any user stuff, just set session. On 'logout', set sessionInfo to 
  # None and delete directory, so all info is gone...
  #


class Root(object):
  
  #
  # Authorisation controller, set up in auth
  #
  
  auth = AuthController()

  @cherrypy.expose
  def index(self):
  
    return """
<html>
<head>
    <title>FormatConverter web server</title>
    <link rel="stylesheet" type="text/css" href="css/style.css" type="text/css"></link>
    <script type="application/javascript" src="js/some.js"></script>
</head>
<html>
<body>
This is the public section
</body>
</html>
"""

  @cherrypy.expose
  def fcLinks(self):
      mytemplate=mylookup.get_template("fcLinks.html")
      return mytemplate.render()

  @cherrypy.expose
  def fcReference(self):
      mytemplate=mylookup.get_template("fcReference.html")
      return mytemplate.render()
        
if __name__ == '__main__':
  
  current_dir = os.path.dirname(os.path.abspath(__file__))

  # TemplateLookup is a Mako object that specifies where mako
  # is to find the html files you want to render. I recommend
  # going through the Basic usage of Mako docs to get started:
  # http://www.makotemplates.org/docs/usage.html. Mako does not
  # need a full physical path and is able to understand a
  # relative path from the location of the python file it has
  # been imported into. Therefore specifying 'arjunhtml' below
  # is enough and 'arjunhtml' is the subdirectory. So the
  # python file sits in ..../fcpyfiles then the 'fchtml'
  # directory sits in ..../fcpyfiles/fchtml for the
  # below to work. You need to change this to suit your own
  # directory structure. Or for a quick start, just create a
  # subdirectory below where you put this python file and put
  # the sub directory name instead of 'fchtml'

  mylookup = TemplateLookup(directories=['data/html'])  # Note that this is now global to script

  #
  # Database setup
  #
  
  from ccpnmr.format.webServer.Util import passwordFile
  
  if not os.path.exists(passwordFile):
    import anydbm
    userDb = anydbm.open(passwordFile,'c')
    userDb[''] = ''
    userDb.close()  
  
  #
  # More configuration, installation specific.
  #
  
  conf = {'/css': {'tools.staticdir.on': True,
                   'tools.staticdir.dir': os.path.join(current_dir, 'data', 'css')},
          '/images': {'tools.staticdir.on': True,
                      'tools.staticdir.dir': os.path.join(current_dir, 'data', 'images')},                           
          '/convert/css': {'tools.staticdir.on': True,
                           'tools.staticdir.dir': os.path.join(current_dir, 'data', 'css')},
          '/convert/images': {'tools.staticdir.on': True,
                              'tools.staticdir.dir': os.path.join(current_dir, 'data', 'images')},              
          '/convertSession/css': {'tools.staticdir.on': True,
                                  'tools.staticdir.dir': os.path.join(current_dir, 'data', 'css')},
          '/convertSession/images': {'tools.staticdir.on': True,
                                     'tools.staticdir.dir': os.path.join(current_dir, 'data', 'images')},              
          '/auth/css': {'tools.staticdir.on': True,
                           'tools.staticdir.dir': os.path.join(current_dir, 'data', 'css')},
          '/auth/images': {'tools.staticdir.on': True,
                              'tools.staticdir.dir': os.path.join(current_dir, 'data', 'images')},              
         }

  #
  # Connect different classes together via root object
  #
  
  root = Root()
  root.convert = Convert()
  root.convertSession = ConvertSession()
  cherrypy.tree.mount(root, config=conf)

  #
  # Set up site-wide config (means all errors go to data/site.log in this directory)
  #
  
  cherrypy.config.update({'environment': 'production',
                          'log.error_file': 'data/site.log',
                          'log.screen': True,
                          'server.socket_host': "localhost",
                          'server.socket_port': 8080,
                          'server.thread_pool': 10,
                          'tools.sessions.on': True,
                          })

  
  #
  # Start the server
  #
  
  cherrypy.quickstart(root, '/', config = conf)
