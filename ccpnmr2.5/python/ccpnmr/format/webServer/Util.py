#
# Controller to provide login and logout actions
#
# Modified from http://tools.cherrypy.org/wiki/AuthenticationAndAccessRestrictions
#

import os, anydbm, glob

# TODO this should go in Constants.py
passwordFile = 'data/users.db'
SESSION_KEY = 'formatConversion'
localUserData = os.path.join('data','user')

userPresentMessage = u"You cannot login because another user is currently logged in from this browser. Log out the other user if you want to continue."

class SessionInfo:

  def __init__(self,username,sessionId):
    
    self.username = username
    self.id = sessionId
    
    #
    # TODO: might have to create username alias with no @ in it if using email address!
    #
    
    self.dataDir = os.path.join(localUserData,"%s.%s" % (username,sessionId))
  
  def startFormatConversion(self):

    #
    # Import formatConversion code
    #

    from ccpnmr.format.general.Conversion import FormatConversion
  
    self.formatConversion = FormatConversion()
    
  # Add other methods/attributes when required!

def fetch_users():
  # Password encrypted in here, so cannot give them 'back' to the user!
  # 
  # TODO TODO Open this db at highest server level?!?! Is that possible with 10 threads?
  #
  usersDict = {}
  userDb = anydbm.open(passwordFile,'r')
  for user in userDb.keys():
    usersDict[user] = userDb[user]
  userDb.close()
  
  return usersDict

#
# Security stuff - could in principle implement this in a way so that password is always 
# transferred in its encrypted form (except when making new user). Need to encrypt via javaScript
# in this case on client side - see http://www.movable-type.co.uk/scripts/sha1.html
#
#import sha
# Encrypts password
#def sha_password_encrypter(password):
# 	  return sha.new(password).hexdigest()  # Creates an sha object, then hexdigest() it

# Gets password for user name
#def fetch_password(username):
#  users = fetch_users()
#  if users.has_key(username):
#    password = sha_password_encrypter(users[username])
#  else:
#    password = None
    
#  return password
import cherrypy

def check_auth(*args, **kwargs):
  """A tool that looks in config for 'auth.require'. If found and it
  is not None, a login is required and the entry is evaluated as a list of
  conditions that the user must fulfill"""
  conditions = cherrypy.request.config.get('auth.require', None)
  if conditions is not None:
  
    sessionInfo = cherrypy.session.get(SESSION_KEY, None)

    if sessionInfo:
      cherrypy.request.login = sessionInfo.username
      for condition in conditions:
        # A condition is just a callable that returns true or false
        if not condition():
          print "NO GOOD"
          raise cherrypy.HTTPRedirect("/auth/login")
    else:
      raise cherrypy.HTTPRedirect("/auth/login")
    
cherrypy.tools.auth = cherrypy.Tool('before_handler', check_auth)

def require(*conditions):
    """A decorator that appends conditions to the auth.require config
    variable."""
    def decorate(f):
        if not hasattr(f, '_cp_config'):
            f._cp_config = dict()
        if 'auth.require' not in f._cp_config:
            f._cp_config['auth.require'] = []
        f._cp_config['auth.require'].extend(conditions)
        return f
    return decorate

def is_current_session_user():
  
  # Redundant, but keep anyway
  
  validUser = False
  
  if cherrypy.session.has_key(SESSION_KEY) and cherrypy.session[SESSION_KEY].username == cherrypy.request.login:
    validUser = True

  return validUser

def check_credentials(username,password):

  usersDict = fetch_users()
  
  error_msg = u'Wrong username/password - please try again.'
  
  if usersDict.has_key(username):
    if password == usersDict[username]:
      if (not cherrypy.session.has_key(SESSION_KEY)) or (not cherrypy.session[SESSION_KEY]):
        error_msg = None
      elif cherrypy.session[SESSION_KEY].username == username:
        error_msg = None
      else:
        error_msg = userPresentMessage
  
  return error_msg
  
def check_newuser(username,password):

  usersDict = fetch_users()
  
  error_msg = None
  
  if usersDict.has_key(username):
    if password == usersDict[username]:
      error_msg = None
    else:
      error_msg = u"Username already exists with different password. Please login or use different username."
  else:
    userDb = anydbm.open(passwordFile,'w')
    if userDb.has_key(''):
      del(userDb[''])
    userDb[username] = password
    userDb.close()
      
  return error_msg

class AuthController(object):
    
    #
    #
    # TODO: the on_login and on_logout have to be made session-only as well (so if no username should still work)
    #       probably need to go elsewhere. Separate username related and session-only code.
    #
    #
    
    
    def on_login(self, username):
    
      """
      Called on successful login
      
      Here make a directory for this particular username/session combination
      
      Used for file upload, CCPN project, and generated files (to be downloaded).
      
      WARNING/TODO: if browser crashes, and relogin, have to make sure find existing directories first, reset everything here!!
      
      """
      
      session = cherrypy.session[SESSION_KEY]
      
      localSessionDirs = glob.glob(os.path.join(localUserData,"%s.*" % session.username))
      
      if localSessionDirs:
      
        #
        # TODO here have to ask if want to continue existing session?
        #
        
        if len(localSessionDirs) == 1: # This should always be true!
          
          os.rename(localSessionDirs[0],session.dataDir)
      
      #
      # Now make sure session directory exists
      #
      
      if not os.path.exists(session.dataDir):
        os.mkdir(session.dataDir)
        
      #
      # TODO here probably have to do some initialisation, like load existing data, ...
      #
      
      
    def on_logout(self, username):
      """
      Called on logout
      
      Here remove all user data from directory, if required!
      
      """
    
    def get_loginform(self, username, msg="Enter login information", from_page="/", checkType = 'login'):
        
        if checkType == 'login':
          title='FormatConversion login'
        else:
          title="FormatConversion new user"
          
        return """
<html>
<head>
  <title>%(title)s</title>  
</head>
<body>
    <h1>%(title)s</h1>  
    <form method="post" action="/auth/login">
    <input type="hidden" name="checkType" value="%(checkType)s" />
    <input type="hidden" name="from_page" value="%(from_page)s" />
    %(msg)s<br />
    Email address (user name): <input type="text" name="username" value="%(username)s" /><br />
    Password: <input type="password" name="password" /><br />
    <input type="submit" value="Log in" />
</body>
</html>""" % locals()

    @cherrypy.expose
    def new(self,username=None,password=None):
      
      return self.login(username = username, password = password, checkType = 'newUser', from_page = "/auth/new")
    
    @cherrypy.expose
    def login(self, username=None, password=None, from_page="/", checkType = 'login'):

      if username is None or password is None:
        return self.get_loginform("", from_page=from_page,checkType=checkType)

      if checkType == 'login':
        error_msg = check_credentials(username, password)
      else:
        error_msg = check_newuser(username,password)

      if error_msg:
        return self.get_loginform(username, error_msg, from_page = from_page,checkType=checkType)
      else:
        if cherrypy.session.has_key(SESSION_KEY) and cherrypy.session[SESSION_KEY]:
          # Doublecheck
          if cherrypy.session[SESSION_KEY].username != username:
            return self.get_loginform(username, userPresentMessage, from_page = from_page,checkType=checkType)
        else:
          sessionKey = cherrypy.request.headers['Cookie']
          (void,sessionId) = sessionKey.split('=')
          cherrypy.session[SESSION_KEY] = SessionInfo(username,sessionId)
          
        cherrypy.request.login = username
        self.on_login(username)
        raise cherrypy.HTTPRedirect("/convert")
    
    @cherrypy.expose
    def logout(self, from_page="/"):
      sess = cherrypy.session
      username = sess[SESSION_KEY].username
      sess[SESSION_KEY] = None
      if username:
        cherrypy.request.login = None
        self.on_logout(username)
      raise cherrypy.HTTPRedirect("/")
